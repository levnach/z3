#include "math/dd/dd_pdd.h"
#include "math/dd/pdd_eval.h"
#include "math/dd/pdd_interval.h"

namespace dd {
static void test1() {
    pdd_manager m(3);
    pdd v0 = m.mk_var(0);
    pdd v1 = m.mk_var(1);
    pdd v2 = m.mk_var(2);
    std::cout << v0 << "\n";
    std::cout << v1 << "\n";
    std::cout << v2 << "\n";
    pdd c1 = v0 * v1 * v2;
    pdd c2 = v2 * v0 * v1;
    std::cout << c1 << "\n";
    SASSERT(c1 == c2);

    c1 = v0 + v1 + v2;
    c2 = v2 + v1 + v0;
    std::cout << c1 << "\n";
    SASSERT(c1 == c2);

    c1 = (v0+v1) * v2;
    c2 = (v0*v2) + (v1*v2);
    std::cout << c1 << "\n";
    SASSERT(c1 == c2);
    c1 = (c1 + 3) + 1;
    c2 = (c2 + 1) + 3;
    std::cout << c1 << "\n";
    SASSERT(c1 == c2);
    c1 = v0 - v1;
    c2 = v1 - v0;
    std::cout << c1 << " " << c2 << "\n";

    c1 = v1*v2;
    c2 = (v0*v2) + (v2*v2);
    pdd c3 = m.zero();
    VERIFY(m.try_spoly(c1, c2, c3));
    std::cout << c1 << " " << c2 << " spoly: " << c3 << "\n";

    c1 = v1*v2;
    c2 = (v0*v2) + (v1*v1);
    VERIFY(m.try_spoly(c1, c2, c3));
    std::cout << c1 << " " << c2 << " spoly: " << c3 << "\n";

    c1 = (v0*v1) - (v0*v0);
    c2 = (v0*v1*(v2 + v0)) + v2;
    c3 = c2.reduce(c1);
    std::cout << c1 << " " << c2 << " reduce: " << c3 << "\n";
}

static void test2() {
    std::cout << "\ntest2\n";
    // a(b^2)cd + abc + bcd + bc + cd + 3 reduce by  bc
    pdd_manager m(4);
    pdd a = m.mk_var(0);
    pdd b = m.mk_var(1);
    pdd c = m.mk_var(2);
    pdd d = m.mk_var(3);
    pdd e = (a * b * b * c * d) + (2*a*b*c) + (b*c*d) + (b*c) + (c*d) + 3;
    std::cout << e << "\n";
    pdd f = b * c;
    pdd r_ef = m.reduce(e, f);
    m.display(std::cout);
    std::cout << "result of reduce " << e << " by " << f << " is " << r_ef <<  "\n";
    pdd r_fe = m.reduce(f, e);
    std::cout << "result of reduce " << f << " by " << e << " is " << r_fe << "\n" ;
    VERIFY(r_fe == f);
}

static void test3() {
    std::cout << "\ntest3\n";
    pdd_manager m(4);
    pdd a = m.mk_var(0);
    pdd b = m.mk_var(1);
    pdd c = m.mk_var(2);
    pdd d = m.mk_var(3);
        
    pdd e = a + c;
    for (unsigned i = 0; i < 5; i++) {
        e = e * e;
    }
    e = e * b;
    std::cout << e << "\n";
}

static void test_reset() {
    std::cout << "\ntest reset\n";
    pdd_manager m(4);
        pdd a = m.mk_var(0);
        pdd b = m.mk_var(1);
        pdd c = m.mk_var(2);
        pdd d = m.mk_var(3);
        std::cout << (a + b)*(c + d) << "\n";

        unsigned_vector l2v;
        for (unsigned i = 0; i < 4; ++i) 
            l2v.push_back(3 - i);
        m.reset(l2v);
        a = m.mk_var(0);
        b = m.mk_var(1);
        c = m.mk_var(2);
        d = m.mk_var(3);
        std::cout << (a + b)*(c + d) << "\n";
    }

static void  test7() {
    std::cout << "\ntest7\n";
    
    pdd_manager m(5);
    pdd a_ = m.mk_var(0);
  
    pdd a = m.mk_var(1);
    pdd b = m.mk_var(2);
    pdd c = m.mk_var(3);
    pdd d = m.mk_var(4);
    pdd_eval ev(m);
    std::function<rational (unsigned)> f = [](unsigned j){ return rational(j); };
    ev.var2val() = f;
    pdd p = a * d * c + a + d;
    VERIFY(ev(p) == rational(1*4*3 + 1 + 4));
    p = a * d * d * c + a + d * c;
    VERIFY(ev(p) == rational(1*4*4*3 + 1 + 4*3));
}
static void  test8() {
    std::cout << "\ntest8\n";
    
    pdd_manager m(5);
    pdd a_ = m.mk_var(0);
  
    pdd a = m.mk_var(1);
    pdd b = m.mk_var(2);
    pdd c = m.mk_var(3);
    pdd d = m.mk_var(4);
    reslimit lim;
    pdd_interval ev(m, lim);

    dep_intervals intervals(lim);
    typedef dep_intervals::interval interval;
    std::function<interval (unsigned, bool )> f = [&intervals](unsigned j, bool){
                                                      interval a;
                                                      intervals.set_interval_for_scalar(a, rational(j));
                                                      intervals.set_upper(a, rational(j + 1));                                                                                return a;
                                                  };
    ev.var2interval() = f;
    pdd p = a * d + d;
    interval i = ev.get_interval<dep_intervals::with_deps_t::without_deps>(p);
    VERIFY(intervals.lower(i) == rational(1*4 + 4));
    VERIFY(intervals.upper(i) == rational((1+1)*(4+1) + (4+1)));
}
>>>>>>> add pdd_interval to evaluate intervals of pdd expressions
}

void tst_pdd() {
    dd::test1();
    dd::test2();
    dd::test3();
    dd::test_reset();
    dd::test7();
    dd::test8();
}
