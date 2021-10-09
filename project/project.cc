#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/csma-module.h"
#include "ns3/traffic-control-module.h"
#include "ns3/flow-monitor-module.h"
#include "ns3/applications-module.h"
#include "ns3/netanim-module.h"
#include "ns3/mobility-helper.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <vector>

bool doCustomPrng = false;

using LcgFunc = std::function<double()>;

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
	Ptr<ExponentialRandomVariable> erv = CreateObject<ExponentialRandomVariable>();

	constexpr size_t n = 1000;
	constexpr double lambda = 300.;

	std::vector<double> ourResults(n);

	std::generate(ourResults.begin(), ourResults.end(), [&]() {
			//return normalize(lcg(), 0u, lcg.m);
		return expDist(normalize(lcg(), 0u, lcg.m), lambda);
	});

	std::vector<double> theirResults(n);

	erv->SetAttribute("Mean", DoubleValue(1./lambda));
	std::generate(theirResults.begin(), theirResults.end(), [&]() {
			//return urv->GetValue(0., 1.);
			//return erv->GetValue(1.0/lambda, 1.);
		return erv->GetValue();
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

std::map<Address, int> packetDistribution;

static void receivedMsg (Ptr<Socket> router, Ptr<Socket> client, Ptr<const Packet> p, const Address &srcAddress , const Address &dstAddress)
{
	//std::cout << "::::: A packet received at the Server! Time:   " << Simulator::Now ().GetSeconds () << std::endl;
	packetDistribution[srcAddress]++;

	Ptr<UniformRandomVariable> rand=CreateObject<UniformRandomVariable>();

	if(rand->GetValue(0.0,1.0) <= 0.3) {
		//std::cout << "::::: Transmitting from Server to Router   "  << std::endl;
		router->Send (Create<Packet> (p->GetSize()));
	}
	else {	// P(0.7)
		//std::cout << "::::: Transmitting from GW to Controller   "  << std::endl;
		client->SendTo(Create<Packet> (p->GetSize()),0,srcAddress);
	}
}

static void durationCallback(Ptr<OutputStreamWrapper> wrapper, Ptr< const Packet > packet, Ptr< NetDevice > txDevice, Ptr< NetDevice > rxDevice, Time duration, Time lastBitTime)
{
	*wrapper->GetStream() << duration.As(Time::MS) << '\n';
}

static void sojournCallback(Ptr<OutputStreamWrapper> wrapper, ns3::Time t)
{

    *wrapper->GetStream()<< t.As(Time::MS) << std::endl;

}

uint32_t globalN = 0; // hehe
uint32_t globalSum = 0;

void tcPacketsInQueue (QueueDiscContainer qdiscs, Ptr<OutputStreamWrapper> stream)
{
	auto disc = qdiscs.Get(qdiscs.GetN() - 4);
	//auto disc = qdiscs.Get(0);
	uint32_t currentSum = disc->GetNPackets();
	//std::cout << currentSum << '\n';
	globalSum += currentSum;
	globalN++;
	*stream->GetStream() << Simulator::Now().GetSeconds () << "\t" << currentSum << std::endl;
}

static void generateTraffic (Ptr<Socket> socket, Ptr<ExponentialRandomVariable> randomSize,	Ptr<ExponentialRandomVariable> randomTime, LcgFunc customRandomSize, LcgFunc customRandomTime)
{
	uint32_t pktSize = doCustomPrng ? customRandomSize() : randomSize->GetInteger(); //Get random value for packet size
	//std::cout << "::::: A packet is generate at Node "<< socket->GetNode ()->GetId () << " with size " << pktSize <<" bytes ! Time:   " << Simulator::Now ().GetSeconds () << std::endl;

	
	// We make sure that the message is at least 12 bytes. The minimum length of the UDP header. We would get error otherwise.
	if(pktSize<12){
		pktSize=12;
	}
	
	socket->Send(Create<Packet>(pktSize));

	Time pktInterval = Seconds(doCustomPrng ? customRandomTime() : randomTime->GetValue()); //Get random value for next packet generation time 
	Simulator::Schedule(pktInterval, &generateTraffic, socket, randomSize, randomTime, customRandomSize, customRandomTime); //Schedule next packet generation
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
		bindings.assign(list);
	}
	
	void installAsPtp() {
		for(const auto &b : bindings) {
			devices.Add(connect(b));
		}

		stack.Install(nodes);

		std::string queueSize = "10000";

		tch.SetRootQueueDisc ("ns3::FifoQueueDisc", "MaxSize", StringValue (queueSize+"p"));

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

	void installAsCsma() {
		CsmaHelper csma;
		for(const auto& b : bindings) {
			csma.SetChannelAttribute("DataRate", StringValue(b.datarate));
			csma.SetQueue ("ns3::DropTailQueue", "MaxSize", StringValue ("1p"));
			NodeContainer newNodes;
			newNodes.Add(nodes.Get(b.first));
			newNodes.Add(nodes.Get(b.second));
			auto installation = csma.Install(newNodes);
			devices.Add(installation);
		}

		stack.Install(nodes);

		std::string queueSize = "10000";

		tch.SetRootQueueDisc ("ns3::FifoQueueDisc", "MaxSize", StringValue (queueSize+"p"));

		qdiscs.Add(tch.Install(devices));

		const std::string base = "10.0.";
		const std::string end = ".0";
		const std::string mask = "255.255.255.0";

		for(uint32_t i = 0; i < devices.GetN()/2; i++) {
			ipv4.SetBase(Ipv4Address((base + std::to_string(i) + end).c_str()), mask.c_str());
			ipv4container.Add(ipv4.Assign({devices.Get(i*2), devices.Get(i*2+1)}));
		}
		std::cout << "Contains: " << ipv4container.GetN() << '\n';

		Ipv4GlobalRoutingHelper::PopulateRoutingTables ();
	}

	NetDeviceContainer connect(const Binding &binding) {
		
		PointToPointHelper ptp;
		ptp.SetDeviceAttribute("DataRate", StringValue(binding.datarate));
		ptp.SetQueue ("ns3::DropTailQueue", "MaxSize", StringValue ("1p"));
		return ptp.Install(nodes.Get(binding.first), nodes.Get(binding.second));
		
		/*
		CsmaHelper csma;
		csma.SetChannelAttribute("DataRate", StringValue(binding.datarate));
		csma.SetQueue("ns3::DropTailQueue", "MaxSize", StringValue("1p"));
		NodeContainer csmaNodes;
		csmaNodes.Add(nodes.Get(binding.first));
		csmaNodes.Add(nodes.Get(binding.second));
		return csma.Install(nodes);
		*/
	}

	uint32_t setSource(uint32_t nodeId, uint32_t addressId, uint32_t port) {
		TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
		Ptr<Socket> source = Socket::CreateSocket(nodes.Get(nodeId), tid);
		auto target = ipv4container.GetAddress(addressId);
		InetSocketAddress remote = InetSocketAddress(target, port);
		source->Connect(remote);
		sockets.push_back(source);
		return sockets.size() -1;
	}

	uint32_t setSource(uint32_t nodeId, uint32_t port) {

		TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
		Ptr<Socket> source = Socket::CreateSocket(nodes.Get(nodeId), tid);
		sockets.push_back(source);
		return sockets.size() -1;
	}

	void installServer(uint32_t nodeId, uint32_t port) {
		UdpServerHelper server(port);
		servers.Add(server.Install(nodes.Get(nodeId)));
	}

	void startSendingSomePackets(uint32_t socketId, double time, double size, LcgFunc customRandomSize, LcgFunc customRandomTime) {
		//Mean inter-transmission time
		Ptr<ExponentialRandomVariable> randomTime = CreateObject<ExponentialRandomVariable> ();
		randomTime->SetAttribute ("Mean", DoubleValue (time));

		//Mean packet time q8^)
		//mean = 100; // 100 Bytes
		Ptr<ExponentialRandomVariable> randomSize = CreateObject<ExponentialRandomVariable> ();

		randomSize->SetAttribute ("Mean", DoubleValue (size));

		auto source = sockets[socketId];


		Simulator::ScheduleWithContext (source->GetNode ()->GetId (), Seconds (1.0), &generateTraffic, source, randomSize, randomTime, customRandomSize, customRandomTime);
	}

	std::vector<Binding> bindings;
    ApplicationContainer servers;
	std::vector<Ptr<Socket>> sockets;
	Ipv4AddressHelper ipv4;
	Ipv4InterfaceContainer ipv4container;
	QueueDiscContainer qdiscs;
	TrafficControlHelper tch;
	InternetStackHelper stack;
	NetDeviceContainer devices;
	NodeContainer nodes;
};

int main(int argc, char** argv) {
	//testPrng();
	//return 0;
	
	uint32_t seed = 12;
	bool doCsma = false;
	CommandLine cmd;
	cmd.AddValue("csma", "Use CSMA instead of Point-to-Point for connections", doCsma);
	cmd.AddValue("seed", "Seed to use", seed);
	cmd.AddValue("useCustomPrng", "Use custom prng instead of the ns3 version", doCustomPrng);
	cmd.Parse(argc, argv);

	RngSeedManager::SetSeed(seed);
	Time::SetResolution(Time::Unit::US);
	std::cout << "Using seed: " << seed << '\n';

	Lcg customRandomTime(seed, std::numeric_limits<Lcg::seed_t>::max());
	customRandomTime.c = 0u;
	Lcg customRandomSize = customRandomTime;

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

	if(doCsma) {
		std::cout << "Using CSMA...\n";
		topology.installAsCsma();
	} else {
		std::cout << "Using PtP...\n";
		topology.installAsPtp();
	}

	if(doCustomPrng) {
		std::cout << "Using custom prng...\n";
	} else {
		std::cout << "Using ns3 prng...\n";
	}

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
	
    uint16_t port_number = 9;  

	// Set up sources

	// From server to router
	auto routerSocketIndex = topology.setSource(7, 15, port_number);

	// From server to clients
	auto clientSocketIndex = topology.setSource(7, port_number);

	// From client to server
	auto aSocketIndex = topology.setSource(0, 13, port_number);
	auto bSocketIndex = topology.setSource(1, 13, port_number);
	auto cSocketIndex = topology.setSource(2, 13, port_number);
	auto dSocketIndex = topology.setSource(3, 13, port_number);

	// Set up servers
	topology.installServer(0, port_number);
	topology.installServer(1, port_number);
	topology.installServer(2, port_number);
	topology.installServer(3, port_number);
	topology.installServer(7, port_number);
	topology.installServer(8, port_number);

	auto routerSocket = topology.sockets[routerSocketIndex];
	auto clientSocket = topology.sockets[clientSocketIndex];

	// When server recieves a message
	auto server = topology.servers.Get(4);
	server->TraceConnectWithoutContext ("RxWithAddresses", MakeBoundCallback (&receivedMsg, routerSocket, clientSocket));

	auto funcCustomRandomTimeA = [&]() {
		constexpr double time = 1./0.002;
		return expDist(normalize(customRandomTime(), 0u, std::numeric_limits<Lcg::seed_t>::max()), time);
		//return expDist(customRandomTime(), time);
	};

	auto funcCustomRandomTimeB = [&]() {
		constexpr double time = 1./0.002;
		return expDist(normalize(customRandomTime(), 0u, std::numeric_limits<Lcg::seed_t>::max()), time);
		//return expDist(customRandomTime(), time);
	};

	auto funcCustomRandomTimeC = [&]() {
		constexpr double time = 1./0.0005;
		return expDist(normalize(customRandomTime(), 0u, std::numeric_limits<Lcg::seed_t>::max()), time);
		//return expDist(customRandomTime(), time);
	};

	auto funcCustomRandomTimeD = [&]() {
		constexpr double time = 1./0.001;
		return expDist(normalize(customRandomTime(), 0u, std::numeric_limits<Lcg::seed_t>::max()), time);
		//return expDist(customRandomTime(), time);
	};

	auto funcCustomRandomSize = [&]() {
		constexpr double size = 1./100.;
		return expDist(normalize(customRandomSize(), 0u, std::numeric_limits<Lcg::seed_t>::max()), size);
		//return expDist(customRandomSize(), size);
	};

	topology.startSendingSomePackets(aSocketIndex, 0.002, 100., funcCustomRandomSize, funcCustomRandomTimeA);
	topology.startSendingSomePackets(bSocketIndex, 0.002, 100., funcCustomRandomSize, funcCustomRandomTimeB);
	topology.startSendingSomePackets(cSocketIndex, 0.0005, 100., funcCustomRandomSize, funcCustomRandomTimeC);
	topology.startSendingSomePackets(dSocketIndex, 0.001, 100., funcCustomRandomSize, funcCustomRandomTimeD);

	MobilityHelper help;
	help.SetMobilityModel("ns3::ConstantPositionMobilityModel");
	help.InstallAll();

	double simulationTime = 5.0+1.0; // +1 because of warmup time

	AnimationInterface anim("out.xml");
	anim.EnablePacketMetadata(true);
	anim.SetConstantPosition(topology.nodes.Get(0), -100, -100);
	anim.SetConstantPosition(topology.nodes.Get(1), -100, -0);
	anim.SetConstantPosition(topology.nodes.Get(2), -50, +100);
	anim.SetConstantPosition(topology.nodes.Get(3),   0, +100);
	anim.SetConstantPosition(topology.nodes.Get(4), -50, -100);
	anim.SetConstantPosition(topology.nodes.Get(5), -50, -0);
	anim.SetConstantPosition(topology.nodes.Get(6), 0, -0);
	anim.SetConstantPosition(topology.nodes.Get(7), 0, -100);
	anim.SetConstantPosition(topology.nodes.Get(8), 100, -0);

	anim.UpdateNodeDescription(topology.nodes.Get(0), "A");
	anim.UpdateNodeDescription(topology.nodes.Get(1), "B");
	anim.UpdateNodeDescription(topology.nodes.Get(2), "C");
	anim.UpdateNodeDescription(topology.nodes.Get(3), "D");
	anim.UpdateNodeDescription(topology.nodes.Get(4), "E");
	anim.UpdateNodeDescription(topology.nodes.Get(5), "F");
	anim.UpdateNodeDescription(topology.nodes.Get(6), "G");
	anim.UpdateNodeDescription(topology.nodes.Get(7), "Server");
	anim.UpdateNodeDescription(topology.nodes.Get(8), "Router");
	anim.EnableQueueCounters(Seconds(0), Seconds(simulationTime), Seconds(0.001));


	AsciiTraceHelper asciiTraceHelper;
	Ptr<OutputStreamWrapper> stream = asciiTraceHelper.CreateFileStream ("project_queue.tr");

	auto transmissionStreamFG = asciiTraceHelper.CreateFileStream("transmission_time_f_g.tr");
	auto sojournStreamFG = asciiTraceHelper.CreateFileStream("sojourn_time_f_g.tr");
	auto deviceF = topology.devices.Get(5);
	auto queueF = topology.qdiscs.Get(10);
	deviceF->GetChannel()->TraceConnectWithoutContext("TxRxPointToPoint", MakeBoundCallback(&durationCallback, transmissionStreamFG));
	queueF->TraceConnectWithoutContext("SojournTime", MakeBoundCallback(&sojournCallback, sojournStreamFG));

	auto transmissionStreamAS = asciiTraceHelper.CreateFileStream("transmission_time_a_s.tr");
	auto sojournStreamAS = asciiTraceHelper.CreateFileStream("sojourn_time_a_s.tr");
	auto deviceA = topology.devices.Get(0);
	auto queueA = topology.qdiscs.Get(0);
	deviceA->GetChannel()->TraceConnectWithoutContext("TxRxPointToPoint", MakeBoundCallback(&durationCallback, transmissionStreamAS));
	queueA->TraceConnectWithoutContext("SojournTime", MakeBoundCallback(&sojournCallback, sojournStreamAS));

	auto transmissionStreamSA = asciiTraceHelper.CreateFileStream("transmission_time_s_a.tr");
	auto sojournStreamSA = asciiTraceHelper.CreateFileStream("sojourn_time_s_a.tr");
	auto deviceS = topology.devices.Get(6);
	auto queueS = topology.qdiscs.Get(12);
	deviceS->GetChannel()->TraceConnectWithoutContext("TxRxPointToPoint", MakeBoundCallback(&durationCallback, transmissionStreamSA));
	queueS->TraceConnectWithoutContext("SojournTime", MakeBoundCallback(&sojournCallback, sojournStreamSA));


	//for (float t=1.0; t < simulationTime; t+=0.001) {
	for (float t=1.0; t < simulationTime; t+=0.00001) {
		Simulator::Schedule(Seconds(t), &tcPacketsInQueue, topology.qdiscs, stream);
	}

	Simulator::Stop(Seconds(simulationTime));
	FlowMonitorHelper flowmon;
	Ptr<FlowMonitor> monitor = flowmon.InstallAll();

	Simulator::Run();
	auto stats = monitor->GetFlowStats();
	auto probe = monitor->GetAllProbes();

	std::cout << "Collected stats for " << stats.size() << " occurences\n";

	auto nB = ++packetDistribution.begin();
	auto nC = ++++packetDistribution.begin();

	auto sumBandC = nB->second + nC->second;

	const auto avgDelay = stats[6].delaySum / sumBandC / simulationTime; // Delay
	std::cout << "F avgDelay: " << avgDelay << '\n';

	for(const auto& i : packetDistribution) {
		std::cout << i.first << ' ' << i.second << '\n';
	}

	Simulator::Destroy();
	flowmon.SerializeToXmlFile("flowmon_stats.xml", false, false);

	std::cout << globalSum << ' ' << globalN << '\n';
	std::cout << static_cast<double>(globalSum) / static_cast<double>(globalN) << '\n';

	/*
	std::cout << "Sum of all enqueued:\n";
	for(auto it = topology.qdiscs.Begin(); it != topology.qdiscs.End(); ++it) {
		std::cout <<  std::distance(topology.qdiscs.Begin(), it) << ": " << (*it)->GetStats().nTotalEnqueuedPackets << '\n';
		std::cout <<  std::distance(topology.qdiscs.Begin(), it) << ": " << (*it)->GetStats().nTotalDequeuedPackets << '\n';
		//sum += (*it)->GetStats().nTotalEnqueuedPackets;
	}
	//std::cout << sum / globalN << '\n';
	*/
}
