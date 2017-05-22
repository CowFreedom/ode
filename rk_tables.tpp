/*Dies sind die Butcher-Tableaus für verschiedene Runge Kutta Verfahren.
Wichtig ist, dass jede neu definierte Klasse von RKTableInterface erbt, welches ein
Interface darstellt.*/

#ifndef rk_tables





/*explizite normale Tabellen*/

//***Euler start***//
namespace r {
template<class rk_Container1D,class rk_Container2D>
class euler_table : public r::RKTableInterface<euler_table<rk_Container1D, rk_Container2D>,rk_Container1D, rk_Container2D> {
	friend class r::RKTableInterface<euler_table, rk_Container1D, rk_Container2D>;
public:
	euler_table();
	rk_Container2D do_getA() { return a; }
	rk_Container1D do_getB_P() { return b_p; };
	rk_Container1D do_getB_Plus() { return b_plus; };
	rk_Container1D do_getC() { return c; };
	unsigned int do_getOrder(){ return 1; }

	template<class T>
	void do_setC(const std::initializer_list<T>&& ct) {/* c = std::move(ct);*/ std::copy(ct.begin(), ct.end(), c.begin()); };

	template<class T>
	void do_setA(const std::initializer_list<std::initializer_list<T>>);

	template<class T>
	void do_setB_P(const std::initializer_list<T>&& bt) { /*b_p = std::move(bt); */std::copy(bt.begin(), bt.end(), b_p.begin());};

	template<class T>
	void do_setB_Plus(const std::initializer_list<T>&& bpt)=delete;
};

template<class rk_Container1D, class rk_Container2D>
template<class T>
void euler_table<rk_Container1D, rk_Container2D>::do_setA(const std::initializer_list<std::initializer_list<T>> at) {
	rk_Container2D::iterator temp = a.begin();
	for (std::initializer_list<std::initializer_list<T>>::iterator iter = at.begin(); iter != at.end(); iter++) {
		std::copy((*iter).begin(), (*iter).end(), (*temp).begin());
		temp++;
	}
}



template<class rk_Container1D, class rk_Container2D>
euler_table<rk_Container1D, rk_Container2D>::euler_table() {
	do_setC<rk_Container1D::value_type>({ 0 });
	do_setB_P<rk_Container1D::value_type>({1});
	do_setA<rk_Container1D::value_type>({{0}});
}
//***Euler end***//

//***Heun Start***//

template<class rk_Container1D,class rk_Container2D>
class heun_table : public r::RKTableInterface<heun_table<rk_Container1D, rk_Container2D>,rk_Container1D, rk_Container2D> {
	friend class r::RKTableInterface<heun_table, rk_Container1D, rk_Container2D>;
public:
	heun_table();
	rk_Container2D do_getA() { return a; }
	rk_Container1D do_getB_P() { return b_p; };
	rk_Container1D do_getB_Plus() { return b_plus; };
	rk_Container1D do_getC() { return c; };
	unsigned int do_getOrder(){ return 1; }

	template<class T>
	void do_setC(const std::initializer_list<T>&& ct) {/* c = std::move(ct);*/ std::copy(ct.begin(), ct.end(), c.begin()); };

	template<class T>
	void do_setA(const std::initializer_list<std::initializer_list<T>>);

	template<class T>
	void do_setB_P(const std::initializer_list<T>&& bt) { /*b_p = std::move(bt); */std::copy(bt.begin(), bt.end(), b_p.begin());};

	template<class T>
	void do_setB_Plus(const std::initializer_list<T>&& bpt)=delete;
};

template<class rk_Container1D, class rk_Container2D>
template<class T>
void heun_table<rk_Container1D, rk_Container2D>::do_setA(const std::initializer_list<std::initializer_list<T>> at) {
	rk_Container2D::iterator temp = a.begin();
	for (std::initializer_list<std::initializer_list<T>>::iterator iter = at.begin(); iter != at.end(); iter++) {
		std::copy((*iter).begin(), (*iter).end(), (*temp).begin());
		temp++;
	}
}

template<class rk_Container1D, class rk_Container2D>
heun_table<rk_Container1D, rk_Container2D>::heun_table() {
	do_setC<rk_Container1D::value_type>({ 0,1 });
	do_setB_P<rk_Container1D::value_type>({0.5,0.5});
	do_setA<rk_Container1D::value_type>({{0,0},{1,0}});
//	ralstons method do_setC<rk_Container1D::value_type>({ 0,2/3 });
//	do_setB_P<rk_Container1D::value_type>({static_cast<double>(1/4),static_cast<double>(3)/4});
//	do_setA<rk_Container1D::value_type>({{0,0},{static_cast<double>(2)/3,0}});
}


/*Eingebettete Tabellen*/

template<class rk_Container1D,class rk_Container2D>
class dopri : public r::RKTableInterface<dopri<rk_Container1D, rk_Container2D>,rk_Container1D, rk_Container2D> {
	friend class r::RKTableInterface<dopri, rk_Container1D, rk_Container2D>;
public:
	dopri();
	rk_Container2D do_getA() { return a; }
	rk_Container1D do_getB_P() { return b_p; };
	rk_Container1D do_getB_Plus() { return b_plus; };
	rk_Container1D do_getC() { return c; };
	unsigned int do_getOrder(){ return 5; }

	template<class T>
	void do_setC(const std::initializer_list<T>&& ct) {/* c = std::move(ct);*/ std::copy(ct.begin(), ct.end(), c.begin()); };

	template<class T>
	void do_setA(const std::initializer_list<std::initializer_list<T>>);

	template<class T>
	void do_setB_P(const std::initializer_list<T>&& bt) { /*b_p = std::move(bt); */std::copy(bt.begin(), bt.end(), b_p.begin());};

	template<class T>
	void do_setB_Plus(const std::initializer_list<T>&& bpt) { /*b_plus = std::move(bpt);*/std::copy(bpt.begin(), bpt.end(), b_plus.begin());};
};

/*Reihenfolge der Templates Deklerationen ist wichtig, siehe*/
template<class rk_Container1D, class rk_Container2D>
template<class T>
void dopri<rk_Container1D, rk_Container2D>::do_setA(const std::initializer_list<std::initializer_list<T>> at) {
	std::array<double, 7> test = { 0.075,0.225,0,0,0,0,0 };
	rk_Container2D::iterator temp = a.begin();
	for (std::initializer_list<std::initializer_list<T>>::iterator iter = at.begin(); iter != at.end(); iter++) {
		std::copy((*iter).begin(), (*iter).end(), (*temp).begin());
		temp++;
	}
}



template<class rk_Container1D, class rk_Container2D>
dopri<rk_Container1D, rk_Container2D>::dopri() {
	do_setC<rk_Container1D::value_type>({ 0,0.2,0.3,0.8,(8 / static_cast<double>(9)),1,1 });
	do_setB_Plus<rk_Container1D::value_type>({(35/static_cast<double>(384)),0,(500/static_cast<double>(1113)),(125/static_cast<double>(192)),(-2187/static_cast<double>(6784)),(11/84),0});
	do_setB_P<rk_Container1D::value_type>({(5179/static_cast<double>(57600)),0,(7571/static_cast<double>(16695)),0.6140625,(-92097/static_cast<double>(339200)),(187/static_cast<double>(2100)),0.025 });
	do_setA<rk_Container1D::value_type>({  {0,0,0,0,0,0,0},
	{0.2,0,0,0,0,0,0},
	{0.075,0.225,0,0,0,0,0},
	{44 / static_cast<double>(45),static_cast<double>(-56) / 15,static_cast<double>(32 / 9),0,0,0,0},
	{(static_cast<double>(19372) / 6561),(static_cast<double>(-25360) / 2187),(static_cast<double>(64448) / 6561),(static_cast<double>(-212) / 729),0,0,0},
	{ (static_cast<double>(9017) / 3168),(static_cast<double>(-355) / 33),(static_cast<double>(46732) / 5247),(static_cast<double>(49) / 176),(static_cast<double>(-5103) / 18656),0,0 },
	{ (static_cast<double>(35) / 384),0,(500 / static_cast<double>(1113)),(static_cast<double>(125) / 192),(static_cast<double>(-2187) / 6784),(static_cast<double>(11) / 84),0 }
	});
}
}
#endif