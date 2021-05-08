#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/traffic-control-module.h"
#include "ns3/flow-monitor-module.h"
#include "ns3/applications-module.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

using namespace ns3;

class Lcg {
public:
	using seed_t = uint32_t;
	Lcg(seed_t seed, seed_t m = 100, seed_t a = 13, seed_t c = 1) 
		: seed(seed), m(m), a(a), c(c) {}

	seed_t operator()() {
		seed = ((a * seed) + c) % m;
		return seed;
	}

	seed_t seed, m, a, c;
};

template<typename T>
double normalize(T x, T min, T max) {
	return static_cast<double>(x - min) / static_cast<double>(max - min);
}

double expDist(double value, double lambda) {
	return -log(1 - value) / lambda;
}

void testPrng() {
	Lcg lcg(1);

	// ref: https://en.wikipedia.org/wiki/Linear_congruential_generator
	lcg.m = std::numeric_limits<Lcg::seed_t>::max();
	lcg.c = 0u;

	// our class: time ./waf --run scratch/project  2.54s user 0.20s system 107% cpu 2.547 total
	// their class: time ./waf --run scratch/project  2.59s user 0.17s system 107% cpu 2.556 total

	//Ptr<UniformRandomVariable> urv = CreateObject<UniformRandomVariable> ();
	Ptr<ExponentialRandomVariable> erv = CreateObject<ExponentialRandomVariable> ();

	constexpr size_t n = 1000;
	constexpr double lambda = 300.;

	std::vector<double> ourResults(n);

	std::generate(ourResults.begin(), ourResults.end(), [&]() {
		return expDist(normalize(lcg(), 0u, lcg.m), lambda);
		//return expDist(lcg(), lambda);
	});

	std::vector<double> theirResults(n);

	std::generate(theirResults.begin(), theirResults.end(), [&]() {
		return erv->GetValue(1.0/lambda, 1.);
	});

	{
		std::ofstream ourFile("our_results.txt");
		for(auto f : ourResults) {
			ourFile << f << ' ';
		}
	}

	{
		std::ofstream theirFile("their_results.txt");
		for(auto f : theirResults) {
			theirFile << f << ' ';
		}
	}

	/*
	std::cout << "Our results: {";
	for(auto f : ourResults) {
		std::cout << f << ", ";
	}
	std::cout << "\b\b}\n";

	std::cout << "Their results: {";
	for(auto f : theirResults) {
		std::cout << f << ", ";
	}
	std::cout << "\b\b}\n";
	*/
}

class NetworkTopology {
public:
	struct Binding {
		uint32_t first, second;
		std::string datarate;
	};

	NetworkTopology(std::initializer_list<Binding> list) {
		uint32_t max = 0;
		for(auto &l : list) {
			if(l.first > max) {
				max = l.first;
			}
			if(l.second > max) {
				max = l.second;
			}
		}
		nodes.Create(max + 1);
		for(const auto &b : list) {
			devices.Add(connect(b));
		}

		stack.Install(nodes);

		std::string queueSize = "10000";

		tch.SetRootQueueDisc ("ns3::FifoQueueDisc", "MaxSize", StringValue (queueSize+"p"));

		/*
		for(auto & device : devices) {
			qdiscs.Add(tch.Install(device));
		}*/
		qdiscs.Add(tch.Install(devices));

		const std::string base = "10.0.";
		const std::string end = ".0";
		const std::string mask = "255.255.255.0";

		for(uint32_t i = 0; i < devices.GetN()/2; i++) {
			ipv4.SetBase(Ipv4Address((base + std::to_string(i) + end).c_str()), mask.c_str());
			ipv4container.Add(ipv4.Assign({devices.Get(i*2), devices.Get(i*2+1)}));
		}

		Ipv4GlobalRoutingHelper::PopulateRoutingTables ();
	}

	NetDeviceContainer connect(const Binding &binding) {
		PointToPointHelper ptp;
		ptp.SetDeviceAttribute("DataRate", StringValue(binding.datarate));
		ptp.SetQueue ("ns3::DropTailQueue", "MaxSize", StringValue ("1p"));
		return ptp.Install(nodes.Get(binding.first), nodes.Get(binding.second));
	}


	Ipv4AddressHelper ipv4;
	Ipv4InterfaceContainer ipv4container;
	QueueDiscContainer qdiscs;
	TrafficControlHelper tch;
	InternetStackHelper stack;
	NetDeviceContainer devices;
	NodeContainer nodes;
};

int main(int argc, char** argv) {
	NetworkTopology topology = {
		{0, 4, "5Mbps"},	// A <-> E
		{1, 5, "5Mbps"},	// B <-> F
		{2, 5, "5Mbps"},	// C <-> F
		{3, 6, "5Mbps"},	// D <-> G
		{4, 6, "5Mbps"},	// E <-> G
		{5, 6, "8Mbps"},	// F <-> G
		{6, 7, "10Mbps"},	// G <-> Server
		{6, 8, "8Mbps"},	// G <-> Router
	};

	/*
	 *  N Queues
	 *
	 *  1 -- 2  1
	 *        \ |
	 *          5 -- 1
	 *        / |
	 *  1 -- 3  1
	 *       |
	 *       1
	 *
	 *  1 + 2 + 1 + 5 + 1 + 1 + 3 + 1 + 1 = 16
	 */


	/*
	 *  Interface index to node name (10.0.*.*)
	 *
	 *    0.1      0.2
	 *  A ------------ E    S
	 *             4.1  \   | 6.2
	 *                   \  |
	 *                    \ |
	 *                 4.2 \| 6.1   7.1            7.2
	 *                      G ------------------------ R
	 *                 5.2 / \ 3.2
	 *                    /   \
	 *                   /     \
	 *             5.1  /       \ 3.1
	 *  B ------------ F         D
	 *    1.1      1.2 | 2.2
	 *                 |
	 *                 |
	 *                 | 2.1
	 *                 C
	 */

	/*
	for(uint32_t i = 0; i < topology.ipv4container.GetN(); i++) {
		std::cout << i << '\n';
		auto nint = topology.ipv4container.Get(i).first->GetNInterfaces();
		for(uint32_t j = 0; j < nint; j++) {
			std::cout << '\t' << topology.ipv4container.Get(i).first->GetAddress(j, 0) << '\n';
		}
	}*/
	
    uint16_t port_number = 9;  
    ApplicationContainer server_apps;
	UdpServerHelper serverS (port_number);
	server_apps.Add(serverS.Install(topology.nodes.Get(7)));

	Ptr<UdpServer> S1 = serverS.GetServer();

	TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");

	Ptr<Socket> source1 = Socket::CreateSocket(topology.nodes.Get(7),tid);
	InetSocketAddress remote1 = InetSocketAddress (topology.ipv4container.GetAddress(15), port_number);
	source1->Connect (remote1);

	for(uint32_t i = 0; i < topology.ipv4container.GetN();i++) {
		std::cout << topology.ipv4container.GetAddress(i) << '\n';
	}
}
