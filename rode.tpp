#include <vector>
#include <limits> //für Maschinengenauigkeit
#include<cmath>
#include <numeric>


//	template<class fContainer, class numberContainer, class number>
template<class fContainer, class numberContainer, class initContainer, class T>
auto r::ode1_old(const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h) {
	//typedef fContainer::value_type function;


	typename numberContainer::const_iterator it_tspan = tspan.begin();
	typename fContainer::const_iterator it_fContainer = u.begin();
	typename initContainer::const_iterator it_initContainer = u0.begin();

	T t = *it_tspan + h;
	T t_end = *(++it_tspan);
	std::vector<initContainer::value_type> temp(u0.size());
	std::vector<std::vector<initContainer::value_type>> trajectory_u(1, std::vector<initContainer::value_type>(u0.begin(), u0.end()));
	std::vector<initContainer::value_type> trajectory_t(1, t);

	for (size_t i = 1; t <= t_end; i++) {
		size_t j = 0;
		for (it_fContainer = u.begin(); it_fContainer != u.end(); ++it_fContainer) {
			temp[j] = trajectory_u[i - 1][j] + h*(*it_fContainer)(trajectory_u[i - 1], t);
			j++;
		}
		trajectory_u.push_back(std::vector<initContainer::value_type>(u0.size()));
		std::copy(temp.begin(), temp.end(), (trajectory_u[i]).begin());
		std::fill(temp.begin(), temp.end(), 0);
		trajectory_t.push_back(t);
		t = t + h;
	}

	auto trajectory = std::make_pair(trajectory_t, trajectory_u);
	it_fContainer = u.begin();
	return trajectory;
}
/*
	template<class fContainer, class numberContainer, class initContainer, class T>
	auto r::ode1(const std::initializer_list<fContainer>& u, const std::initializer_list<numberContainer>& tspan, const std::initializer_list<initContainer>& u0, T h){
	return r::ode1((u),(tspan),u0,h);
	}
*/

template<class fContainer, class numberContainer, class initContainer, class T, class M>
M r::ode1_old(fContainer& u, const std::initializer_list<numberContainer>& tspan, const std::initializer_list<initContainer>& u0, T h) {
	return r::ode1_old((u), (tspan), u0, h);
}

/*Template für explizite Runge Kutta Verfahren ohne Schrittweitensteuerung.*/
template<class fContainer, class numberContainer, class initContainer, class T, class alg, class rk_Container1D, class rk_Container2D >
auto r::Explicit_RKTemplate(r::RKTableInterface<alg, rk_Container1D, rk_Container2D>& butcher_tableau, const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h) {
	//fContainer& u contains the differential equation functions
	//numberContainer tspan contains the endpoints of the interval t_start, t_end
	//initContainer& u0 contains the initial values
	//t contains the step size

	typename rk_Container2D::const_iterator it_rk_Table_a = butcher_tableau.getA().begin();

	const size_t u_size = u.size();
	const auto t_start = *tspan.begin();
	const auto t_end = *(++tspan.begin());

	const auto t_difference = abs((t_start - t_end) / h);

	std::vector<std::vector<initContainer::value_type>> trajectory_u(1, std::vector<initContainer::value_type>(u0.begin(), u0.end()));
	std::vector<initContainer::value_type> trajectory_t(1, t_start);

	//Check if intervall is multiple of step size to properly prereserve vector

	if (floor(t_difference) == t_difference) {
		trajectory_u.reserve(t_difference);
		trajectory_t.reserve(t_difference);

	}
	else {
		trajectory_u.reserve(static_cast<size_t>(t_difference) + 1);
		trajectory_t.reserve(static_cast<size_t>(t_difference) + 1);
	}

	//Array with temporary solution of next step
	std::vector<initContainer::value_type> u_t(trajectory_u[0].begin(), trajectory_u[0].end());

	//Number of Runge-Kutta stages
	const size_t L = (*it_rk_Table_a).size();

	//Solutions of intermediate locations between step sizes
	std::vector<std::vector<initContainer::value_type>> y(L, std::vector<initContainer::value_type>(u_size, 0));

	//Sum of the intermediate runge kutta stages
	std::vector<initContainer::value_type> sum(u0.size(), 0);

	auto& tk_L = butcher_tableau.getC();
	auto& A = butcher_tableau.getA();
	auto& B = butcher_tableau.getB_P();

	bool intervall_overshoot = false;
	bool one_last_repetition = true;
	auto t = t_start;
	size_t iteration = 0;

	while (t < t_end && one_last_repetition) {

		//Entscheide, ob loop nocheinmal wiederholt werden muss
		one_last_repetition = !intervall_overshoot;

		if (t + h > t_end) {
			intervall_overshoot = true;
			h = t_end - t;
		}
		//Berechne Stützstellen
		for (size_t i = 0; i < L; i++) {
			sum = u_t;
			for (size_t k = 0; k < u_size; k++) {
				for (size_t j = 0; j < i; j++) {
					//	if(k==0)
					//		std::cout << "s: " << i << "j: " << j << "Funktion:" << k <<"prior_sum[k]:"<<sum[k]<<"updated_sum[k]: "<< y[j][k] * A[i][j] <<"y[i][k]: "<< y[j][k]<< "t: "<<t<<"\n";
					sum[k] += h*y[j][k] * A[i][j];
				}
			}
			for (size_t k = 0; k < u_size; k++) {
				y[i][k] = u[k](sum, t + h*tk_L[i]);
			}
		}
		iteration++;
		t = t + h;
		//Berechne u_{t+1}
		for (size_t i = 0; i < u_size; i++) {
			for (size_t j = 0; j < L; j++) {
				u_t[i] += h*y[j][i] * B[j];
			}
		}
		trajectory_u.push_back(u_t);
		trajectory_t.push_back(t);
	}
	return std::make_pair(trajectory_t, trajectory_u);
}

/*Template für explizite Runge Kutta Verfahren ohne Schrittweitensteuerung.*/
template<class fContainer, class numberContainer, class initContainer, class alg, class rk_Container1D, class rk_Container2D, class measure>
auto r::Explicit_variable_RKTemplate(
	r::RKTableInterface<alg, rk_Container1D, rk_Container2D>& butcher_tableau,
	const fContainer& u,
	const numberContainer& tspan,
	const initContainer& u0,
	const numberContainer& h_values,
	measure rule) {
	//fContainer& u contains the differential equation functions
	//numberContainer tspan contains the endpoints of the interval t_start, t_end
	//initContainer& u0 contains the initial values
	//t contains the step size

	typename rk_Container2D::const_iterator it_rk_Table_a = butcher_tableau.getA().begin();

	const size_t u_size = u.size();
	const auto t_start = *tspan.begin();
	const auto t_end = *(++tspan.begin());
	numberContainer::value_type h_min=std::min(h_values[0],h_values[1]);
	numberContainer::value_type h_max = std::max(h_values[0], h_values[1]);
	numberContainer::value_type h = h_max;
	numberContainer::value_type h_current = h;
	const auto t_difference = abs((t_start - t_end) / h_min);//Es wird eventuell sehr viel reserviert!

	std::vector<std::vector<initContainer::value_type>> trajectory_u(1, std::vector<initContainer::value_type>(u0.begin(), u0.end()));
	std::vector<initContainer::value_type> trajectory_t(1, t_start);

	//Check if intervall is multiple of step size to properly prereserve vector

	if (floor(t_difference) == t_difference) {
		trajectory_u.reserve(t_difference);
		trajectory_t.reserve(t_difference);

	}
	else {
		trajectory_u.reserve(static_cast<size_t>(t_difference) + 1);
		trajectory_t.reserve(static_cast<size_t>(t_difference) + 1);
	}

	//Array with temporary solution of next step
	std::vector<initContainer::value_type> u_t;
	std::vector<initContainer::value_type> error(u_size, 0);
	//Number of Runge-Kutta stages
	const size_t L = (*it_rk_Table_a).size();

	//Solutions of intermediate locations between step sizes
	std::vector<std::vector<initContainer::value_type>> y(L, std::vector<initContainer::value_type>(u_size, 0));

	//Sum of the intermediate runge kutta stages
	std::vector<initContainer::value_type> sum(u0.size(), 0);

	auto& tk_L = butcher_tableau.getC();
	auto& A = butcher_tableau.getA();
	auto& B_P = butcher_tableau.getB_P();
	auto& B_Plus = butcher_tableau.getB_Plus();
	const auto p = butcher_tableau.getOrder();
	bool intervall_overshoot = false;
	bool one_last_repetition = true;
	auto t = t_start;
	size_t iteration = 0;

	//std::array<>

	while (t < t_end && one_last_repetition) {
		u_t = trajectory_u[iteration];
		//Entscheide, ob loop nocheinmal wiederholt werden muss
		one_last_repetition = !intervall_overshoot;

		if (t + h > t_end) {
			intervall_overshoot = true;
			h = t_end - t;
		}
		//Berechne Stützstellen
		for (size_t i = 0; i < L; i++) {
			sum = u_t;
			for (size_t k = 0; k < u_size; k++) {
				for (size_t j = 0; j < i; j++) {
					//	if(k==0)
					//		std::cout << "s: " << i << "j: " << j << "Funktion:" << k <<"prior_sum[k]:"<<sum[k]<<"updated_sum[k]: "<< y[j][k] * A[i][j] <<"y[i][k]: "<< y[j][k]<< "t: "<<t<<"\n";
					sum[k] += h*y[j][k] * A[i][j];
				}
			}
			for (size_t k = 0; k < u_size; k++) {
				y[i][k] = u[k](sum, t + h*tk_L[i]);
			}
		}

		//Berechne u_{t+1}
		for (size_t i = 0; i < u_size; i++) {
			for (size_t j = 0; j < L; j++) {
				u_t[i] += h*y[j][i] * B_Plus[j];
				error[i] += (B_Plus[j] - B_P[j])*y[j][i];
				//std::cout << "error[" << i << "]: " << error[i] << "\n";
			}
		}
		h_current = h;
		bool stepSizeOK = rule(error, h);
		fill(error.begin(), error.end(), 0);
		if (stepSizeOK) {
			iteration++;
			t = t + h_current;
			trajectory_u.push_back(u_t);
			trajectory_t.push_back(t);
		}
	}
	return std::make_pair(trajectory_t, trajectory_u);
}


/*Funktionsaufrufe*/
template<class alg, class rk_Container1D, class rk_Container2D, class T>
auto f2(r::RKTableInterface<alg, rk_Container1D, rk_Container2D> rk, T h) {
	std::vector<double> trajectory_t;
	std::vector<std::vector<double>> trajectory_u;
	std::cout << "Erfolg!\n";
	return std::make_pair(trajectory_t, trajectory_u);
}

//Euler
template<class fContainer, class numberContainer, class initContainer, class T>
auto r::ode1(const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h) {
	r::euler_table<std::array<double, 1>, std::array<std::array<double, 1>, 1>> euler;

	//	std::vector<double> trajectory_t;
	//	std::vector<std::vector<double>> trajectory_u;
	//	return std::make_pair(trajectory_t,trajectory_u);
	//	return	f2(euler,0.5);
	return r::Explicit_RKTemplate(euler, u, tspan, u0, h);
	//	return Explicit_RKTemplate2(euler,u,u0,h);
}

//Heun
template<class fContainer, class numberContainer, class initContainer, class T>
auto r::ode2(const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h) {
	r::heun_table<std::array<double, 2>, std::array<std::array<double, 2>, 2>> heun;
	return r::Explicit_RKTemplate(heun, u, tspan, u0, h);
}

template<class fContainer, class numberContainer, class initContainer, class T>
auto r::ode4(const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h) {
	r::dopri<std::array<double, 7>, std::array<std::array<double, 7>, 7>> dopri;
	return r::Explicit_RKTemplate(dopri, u, tspan, u0, h);

}

//Dopri
template<class fContainer, class numberContainer, class initContainer>
auto r::ode45(const fContainer& u, const numberContainer& tspan, const initContainer& u0, const numberContainer& h_values) {
	r::dopri<std::array<double, 7>, std::array<std::array<double, 7>, 7>> dopri;
	auto t_difference = abs(tspan[1] - tspan[0]);
//	initContainer::value_type tol = t_difference*0.1; //ändern
	auto tol = 0.00001;
	const auto p = dopri.do_getOrder();

	numberContainer::value_type h_min = std::min(h_values[0], h_values[1]);
	numberContainer::value_type h_max = std::max(h_values[0], h_values[1]);

	auto rule = [&](const auto local_d_error, auto& h) { //C++ 2014 erlaubt keyword "lambda" in Lambdas
		auto c = 1;
		auto k = 2;
		auto error_l2 = std::pow(std::inner_product(local_d_error.begin(), local_d_error.end(), local_d_error.begin(), 0.0),0.5);
		auto h_opt = pow(tol / (t_difference*error_l2), static_cast<double>(1) / p)*h;
		if (h_opt<c*h) {
			if (h_opt < h_min) {
				std::cout << "Abbruch! h_opt: " << h_opt<<" h_min: "<<h_min;
				throw std::invalid_argument("Current step size below h_min.");
			}
			h = std::min(h_opt*0.5, h_max);
			return false;
		}
		else {
			//td::cout<<"tau: " << error_l2 << "\n;";
			//std::cout << "opt Berechnung:"<<(error_l2)/(pow(h,p))*pow(h,p)<<"hopt"<<h_opt<<"\n";
			h = std::min(k*h,h_max);
			return true;
		}

	};

	return r::Explicit_variable_RKTemplate(dopri, u, tspan, u0, h_values,rule);
}


/*
template<class fContainer, class numberContainer, class initContainer, class T, class M>
M r::ode2(const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h){
euler<std::array<double,1>,std::array<std::array<double,1>,1>> euler_table;
return r::Explicit_RKTemplate(euler_table,u,tspan,u0,h);
}

*/