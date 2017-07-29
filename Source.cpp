#include <iostream>
#include <array>
#include "rode.h"
#include <functional> // für std::function
#include <vector>
#include <ctime> //um performance zu messen
#include <cmath>//um Funktionen zu testen in testalgorithms()

template<class T, class M>
void printTrajectory(std::pair<M,T> trajectory) {
	std::cout << "t:\tu:\n";
	for (size_t i = 0; i < (trajectory.second).size(); i++) {
		std::cout << (trajectory.first)[i] << "\t";
		for (size_t j = 0; j < (trajectory.second)[i].size(); j++) {
			std::cout << (trajectory.second)[i][j] << "\t";
		}
		std::cout << "\n";
	}
}

template<class T, class M>
void saveTrajectory(std::pair<M, T> trajectory) {

	FILE *f;

	errno_t errorCode = fopen_s(&f, "C:/Users/Tristan_local/Documents/Visual Studio 2015/Projects/ODE/rode/trajectories/test.txt", "w");

	if (f == NULL) {

		printf("Error opening file!\n");
		std::cin.get();
		exit(1);
	}

	for (size_t i = 0; i < (trajectory.second).size(); i++) {
		fprintf(f, "%lf\t", (trajectory.first)[i]);
		for (size_t j = 0; j < (trajectory.second)[i].size(); j++) {
			fprintf(f, "%lf\t", trajectory.second[i][j]);
		}
		fprintf(f, "\n");
	}

	printf("Everything worked!\n");
	fclose(f);
}


double f1(std::vector<double> u, double t)  {
	return 0.02*u[0]*u[1];
}


double f2(std::vector<double> u, double t)  {
	return u[1];
}

void testlotkavolterra() {
	double gamma = 1.5;
	double s = 0.01;
	double c = 0.8;
	double delta = 2;
	double r = s*c;
	double sq = 0.02;
	double cq = 0.4;
	double inhibition = 0.001;

	auto u1 = [&](std::vector<double>&& u, double t) 
	{return gamma*u[0] * (1 - (u[0] * inhibition)) - s*u[0] * u[1]; };

	auto u2 = [&](std::vector<double>&& u, double t)
	{return -delta*u[1] + r*u[0] * u[1]; };

	std::vector<std::function<double(std::vector<double>, double)>> functions;
	functions.push_back(u1);
	functions.push_back(u2);
	//std::array<[]double(*)(std::vector<double>, double), 2> a = { u1,u2 };
	std::array<double, 2> tspan = { 0,10 };
	std::array<double, 2> u0 = { 200,100 };

	
	clock_t begin = std::clock();
	r::ode4(functions, tspan, u0, 0.001);

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Zeit ODE1:" << elapsed_secs << "\n";

	begin = std::clock();
	r::ode1_old(functions, tspan, u0, 0.001);
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "Zeit ODE2:" << elapsed_secs << "\n";

	

//	printTrajectory(r::ode2(functions, tspan, u0, 0.001));
//	printTrajectory(r::ode4(functions, tspan, u0, 0.001));
	//printTrajectory(r::ode1(functions, { 0,30 }, { 200,100 }, 0.001));
	//printTrajectory(r::ode1(functions, tspan, u0, 0.001));

//	printTrajectory(r::ode1(functions, { 0,30 }, { 200,100 }, 0.001));
//	saveTrajectory(r::ode2(functions, tspan, u0, 0.001));
	//printTrajectory(r::ode1({ f1,f2 }, { 1,5 }, { -1,2 }, 1));
}

void testalgorithms() {
	double h = 0.001;
	double t_start = 0;
	double t = t_start;

	std::array<double, 2> tspan = { t_start,1 };
	std::array<double, 2> u0 = { 1,1 };
	auto u1 = [&](std::vector<double>&& u, double t) {
		return u[0];
	};
	auto u2 = [&](std::vector<double>&& u, double t) {
		return -u[1];
	};
	std::vector<std::function<double(std::vector<double>, double)>> functions;
	functions.push_back(u1);
	functions.push_back(u2);

	auto r1 = [&](double t) {
		return exp(t);
	};
	auto r2 = [&](double t) {
		return exp(-t);
	};
	std::vector<std::function<double(double)>> real_functions;
	real_functions.push_back(r1);
	real_functions.push_back(r2);
	std::cout << "Globaler Fehler (u-y)/h:\n";
	auto result = r::ode4(functions, tspan, u0, h);
	std::cout << "t:\tu:\n";
	for (size_t i = 0; i < (result.second).size(); i++) {
		std::cout << (result.first)[i] << "\t";
		for (size_t j = 0; j < (result.second)[i].size(); j++) {
			std::cout << (result.second)[i][j] - (real_functions[j](t)) << "\t";
		}
		std::cout << "\n";
		t = t + h;
	}
//	printTrajectory(result);
}

void sample_ode() {
	double h = 0.2;
	double t_start = 0;
	double t = t_start;
	double t_end = 1;

	std::array<double, 2> tspan = { t_start,t_end };
	std::array<double, 2> u0 = { 1,1 };

	auto u1 = [&](std::vector<double>&& u, double t) {
		return t;
	};
	auto u2 = [&](std::vector<double>&& u, double t) {
		return -u[1];
	};
	std::vector<std::function<double(std::vector<double>, double)>> functions;
	functions.push_back(u1);
	functions.push_back(u2);

	auto result = r::ode4(functions, tspan, u0, h);

		printTrajectory(result);
}

int main() {

	testlotkavolterra();
//	testalgorithms();
//	sample_ode();

	std::cin.get();
}