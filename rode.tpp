#include <vector>
#include <limits> //für Maschinengenauigkeit
#include<cmath>

	
//	template<class fContainer, class numberContainer, class number>
template<class fContainer, class numberContainer, class initContainer, class T>
auto r::ode1_old(const fContainer& u, const numberContainer& tspan,const initContainer& u0, T h){
	//typedef fContainer::value_type function;


	typename numberContainer::const_iterator it_tspan = tspan.begin();
	typename fContainer::const_iterator it_fContainer = u.begin();
	typename initContainer::const_iterator it_initContainer = u0.begin();

	T t = *it_tspan;
	T t_end = *(++it_tspan);
	std::vector<initContainer::value_type> temp(u0.size());
	std::vector<std::vector<initContainer::value_type>> trajectory_u(1,std::vector<initContainer::value_type>(u0.begin(),u0.end()));
	std::vector<initContainer::value_type> trajectory_t(1,t);
	
	for (size_t i=1; t <= t_end; i++) {
	size_t j=0;
		for (it_fContainer = u.begin(); it_fContainer != u.end(); ++it_fContainer) {
			temp[j]=trajectory_u[i-1][j]+h*(*it_fContainer)(trajectory_u[i-1],t);
			j++;
		}
		trajectory_u.push_back(std::vector<initContainer::value_type>(u0.size()));
		std::copy(temp.begin(),temp.end(),(trajectory_u[i]).begin());
		std::fill(temp.begin(),temp.end(),0);
		trajectory_t.push_back(t);
		t = t + h;

	}
	
	auto trajectory=std::make_pair(trajectory_t,trajectory_u);
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
	M r::ode1_old(fContainer& u, const std::initializer_list<numberContainer>& tspan, const std::initializer_list<initContainer>& u0, T h){
	return r::ode1_old((u),(tspan),u0,h);
	}

	/*Template für explizite Runge Kutta Verfahren ohne Schrittweitensteuerung.*/
	template<class fContainer, class numberContainer, class initContainer, class T, class alg, class rk_Container1D, class rk_Container2D >
	auto r::Explicit_RKTemplate(r::RKTableInterface<alg,rk_Container1D, rk_Container2D>& butcher_tableau, const fContainer& u, const numberContainer& tspan,const initContainer& u0, T h){
	
	typename numberContainer::const_iterator it_tspan = tspan.begin();
	typename fContainer::const_iterator it_fContainer = u.begin();
	typename initContainer::const_iterator it_initContainer = u0.begin();
	typename rk_Container2D::const_iterator it_rk_Table_a = butcher_tableau.getA().begin();
	typename rk_Container1D::const_iterator it_B_P=butcher_tableau.getB_P().begin();
	const initContainer::value_type t_start = *it_tspan;
	initContainer::value_type t=t_start+h;
	const initContainer::value_type t_end = *(++it_tspan);
	if(t_start>t_end){
	h=h*-1;
	}
	std::vector<std::vector<initContainer::value_type>> trajectory_u(1,std::vector<initContainer::value_type>(u0.begin(),u0.end()));
	std::vector<initContainer::value_type> u1(trajectory_u[0].begin(),trajectory_u[0].end());
	std::vector<initContainer::value_type> trajectory_t(1,t_start);
	const size_t L=(*it_rk_Table_a).size();

	std::vector<std::vector<initContainer::value_type>> y_K_L(L, std::vector<initContainer::value_type>(u0.size(),0));

	bool end_not_reached=true;
	double precision=std::numeric_limits<double>::denorm_min();
	std::vector<initContainer::value_type> temp(u0.size(),0);

	rk_Container1D tk_L;

	for (size_t j=0;j<L;j++){
	tk_L[j]=t+(*(butcher_tableau.getC().begin()+j))*h;
	}

	for (size_t i=1; end_not_reached; i++){
	if(abs(t_end-t)<2*precision){
	end_not_reached=false;
	}
		for (size_t s_1=0;s_1<L;s_1++){
			for (size_t s_2=0;s_2<=s_1;s_2++){
				rk_Container2D::const_iterator currentAs_1_it=butcher_tableau.getA().begin();
				for (size_t func_i=0;func_i<u.size();func_i++){
				temp[func_i]=temp[func_i]+(*((*(currentAs_1_it+s_2)).begin()))*(*std::next(it_fContainer,func_i))(y_K_L[s_1],tk_L[s_2]);			
				}						
			}
			for(size_t func_i=0;func_i<u.size();func_i++){
				y_K_L[s_1][func_i]=trajectory_u[i-1][func_i]+h*temp[func_i];		
			}
			for(size_t func_i2=0;func_i2<u.size();func_i2++){
//				u1[func_i2]=u1[func_i2]+h*(*std::next(it_B_P,s_1))*(*std::next(it_fContainer,func_i2))(y_K_L[s_1],tk_L[s_1]);
				u1[func_i2]=u1[func_i2]+h*(*std::next(butcher_tableau.getB_P().begin(),s_1))*(*std::next(it_fContainer,func_i2))(y_K_L[s_1],tk_L[s_1]);				
			}
			std::fill(temp.begin(),temp.end(),0);
		}
		trajectory_u.push_back(std::vector<initContainer::value_type>(u0.size()));
		std::copy(u1.begin(),u1.end(),(trajectory_u[i]).begin());
		trajectory_t.push_back(t);
		t=t+h;
		if(t>t_end && t-h<t_end){
		t=t-(t-t_end);
		}
}
/*
*/
//	std::cout<<"Success!";
	return std::make_pair(trajectory_t,trajectory_u);		
}

/*Funktionsaufrufe*/
template<class alg, class rk_Container1D, class rk_Container2D, class T>
auto f2(r::RKTableInterface<alg, rk_Container1D, rk_Container2D> rk, T h) {
	std::vector<double> trajectory_t;
	std::vector<std::vector<double>> trajectory_u;
	std::cout<<"Erfolg!\n";
	return std::make_pair(trajectory_t,trajectory_u);
}

//Euler
template<class fContainer, class numberContainer, class initContainer, class T>
auto r::ode1(const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h){
	r::euler_table<std::array<double,1>,std::array<std::array<double,1>,1>> euler;

//	std::vector<double> trajectory_t;
//	std::vector<std::vector<double>> trajectory_u;
//	return std::make_pair(trajectory_t,trajectory_u);
//	return	f2(euler,0.5);
	return r::Explicit_RKTemplate(euler,u,tspan,u0,h);
//	return Explicit_RKTemplate2(euler,u,u0,h);
}

//Heun
template<class fContainer, class numberContainer, class initContainer, class T>
auto r::ode2(const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h){
	r::heun_table<std::array<double,2>,std::array<std::array<double,2>,2>> heun;
	return r::Explicit_RKTemplate(heun,u,tspan,u0,h);
}

template<class fContainer, class numberContainer, class initContainer, class T>
auto r::ode4(const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h){
	r::dopri<std::array<double,7>,std::array<std::array<double,7>,7>> dopri;
	return r::Explicit_RKTemplate(dopri,u,tspan,u0,h);

}

/*
template<class fContainer, class numberContainer, class initContainer, class T, class M>
M r::ode2(const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h){
euler<std::array<double,1>,std::array<std::array<double,1>,1>> euler_table;
return r::Explicit_RKTemplate(euler_table,u,tspan,u0,h);
}

*/