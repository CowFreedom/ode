#pragma once

#ifndef rode

namespace r {

	/*Interfaces*/
	/*Runge Kutta Interface*/
	template<class Derived>
	class enable_down_cast {
	private:
		typedef enable_down_cast Base;
	public:
		Derived const* self() const {
			return static_cast<Derived const*>(this);
		}
		Derived* self() {
			return static_cast<Derived*>(this);
		}
	protected:
		~enable_down_cast() = default;
	};

	template<class Impl, class rk_Container1D, class rk_Container2D>
	class RKTableInterface :public enable_down_cast<Impl> {

	private:
		using enable_down_cast<Impl>::self;

	public:
		rk_Container2D getA() { return self()->do_getA(); }
		rk_Container1D getC() { return self()->do_getC(); }
		rk_Container1D getB_P() { return self()->do_getB_P(); }
		rk_Container1D getB_Plus() { return self()->do_getB_Plus(); }
		unsigned int getOrder() { return self()->do_getOrder(); }

		void setB_Plus(rk_Container1D&) { self()->do_setB_Plus(); };

		void setB_P(rk_Container1D&) { self()->do_setB_P(); };

		template<class T>
		void setC(const std::initializer_list<T>&&) { self()->do_setC(const std::initializer_list<T>&&); };

		template<class T>
		void setA(const std::initializer_list<std::initializer_list<T>>) { self()->do_setA(const std::initializer_list<std::initializer_list<T>>); };

		template<class T>
		void setB_P(const std::initializer_list<T>&&) { self()->do_setB_P(const std::initializer_list<T>&&); };

		template<class T>
		void setB_Plus(const std::initializer_list<T>&&) { self()->do_setB_Plus(const std::initializer_list<T>&&); };


	protected:
		~RKTableInterface() = default;
		rk_Container2D a;
		rk_Container1D b_p;
		rk_Container1D b_plus;
		rk_Container1D c;

	};

	/*Template für explizite Runge Kutta Verfahren ohne Schrittweitensteuerung.*/
	template<class fContainer, class numberContainer, class initContainer, class T, class alg, class rk_Container1D, class rk_Container2D >
	auto Explicit_RKTemplate(RKTableInterface<alg,rk_Container1D, rk_Container2D>& butcher_tableau, const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h);

	/*
	Algorithmus: Euler (explizit)
	Ordnung des Verfahrens: 1
	Adaptive Schrittweite: Nein

	Diese Funktion nutzt den expliziten Euler, ein Verfahren erster Ordnung.
	Aufgrund der vergleichsweise geringen Genauigkeit sollte dieser Algorithmus nur für einfachste ODE's
	genutzt werden oder wenn Performance wichtig ist.*/
	template<class fContainer, class numberContainer, class initContainer, class T>
	auto ode1_old(const fContainer& u, const numberContainer& tspan,const initContainer& u0, T h);

//	template<class fContainer, class numberContainer, class initContainer, class T>
//	auto ode1(const std::initializer_list<fContainer>& u, const std::initializer_list<numberContainer>& tspan, const std::initializer_list<initContainer>& u0, T h);

	template<class fContainer, class numberContainer, class initContainer, class T, class M>
	M ode1_old(fContainer& u, const std::initializer_list<numberContainer>& tspan, const std::initializer_list<initContainer>& u0, T h);

	template<class fContainer, class numberContainer, class initContainer, class T>
	auto ode1(const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h);

	//Heun
	template<class fContainer, class numberContainer, class initContainer, class T>
	auto ode2(const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h);


	template<class fContainer, class numberContainer, class initContainer, class T>
	auto ode4(const fContainer& u, const numberContainer& tspan, const initContainer& u0, T h);
}




#include "rk_tables.tpp"

#include "rode.tpp"

#endif