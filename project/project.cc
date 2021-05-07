#include "ns3/core-module.h"

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

int main(int argc, char** argv) {
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
