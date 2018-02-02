/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Lev Nachmanson
*/
#pragma once
#include <map>
#include "util/trace.h"
#include "util/lp/stacked_value.h"
#include "util/lp/stacked_map.h"
namespace lp {
enum class endpoint_kind { START, END, STEND }; 
// represents the set of disjoint intervals of integer number
template <typename T>
class integer_domain {
#ifdef Z3DEBUG
    //    std::set<int> m_domain;
#endif
    struct endpoint { int m_start_expl; int m_end_expl;
        endpoint() : m_start_expl(-1), m_end_expl(-1) {}
        endpoint(int s, int e) : m_start_expl(s), m_end_expl(e) {}
        endpoint_kind kind() const {
            lp_assert(m_start_expl != -1 || m_end_expl != -1);
            if (m_end_expl == -1) {
                return endpoint_kind::START;
            }
            if (m_start_expl == -1) {
                return endpoint_kind::END;
            }
            
            return endpoint_kind::STEND;
        }
        bool operator==(const endpoint& e) const {
            return m_start_expl == e.m_start_expl && m_end_expl == e.m_end_expl;
        }
        bool operator!=(const endpoint & e) const { return !(*this == e); }

        void print(std::ostream & out) const {
            if (m_start_expl != -1 && m_end_expl != -1)
                out << "(" <<  m_start_expl << ", " << m_end_expl << ")";
            else {
                if (m_start_expl != -1) {
                    out << "(" << m_start_expl << ")";
                }
                else if (m_end_expl != -1) {
                    out << "(" << m_end_expl << ")";
                }
            }
        }
    };
    stacked_map<T, endpoint>                 m_endpoints;
    stacked_value<bool>                      m_empty;
    typedef typename std::map<T, endpoint>::iterator iter;
    typedef typename std::map<T, endpoint>::const_iterator const_iter;
    typedef typename std::map<T, endpoint>::reverse_iterator riter;
    
public:
    // the default constructor creates a set containing all integer numbers
    integer_domain() : integer_domain(false) {}


    // if is_empty = false then the constructor creates a set containing all integer numbers,
    // otherwise it creates an empty set
    integer_domain(bool is_empty) : m_empty(is_empty) {
#if Z3DEBUG
        // if (!is_empty) {
        //     for (int i = 0; i <= 100; i++)
        //         m_domain.insert(i);
        // }
#endif
    }
    
    void init_to_contain_all() {
        m_empty = false;
        m_endpoints.clear();
#if Z3DEBUG
        // for (int i = 0; i <= 100; i++)
        //     m_domain.insert(i);
#endif
    }

    // copy constructor 
    integer_domain(const integer_domain<T> & t) :
#if Z3DEBUG
        //        m_domain(t.m_domain),
#endif
        m_endpoints(t.m_endpoints),
        m_empty(t.m_empty)
    {
    }

    // needed for debug only
    void restore_domain() {
#if Z3DEBUG
        // for (int i = 0; i <= 100; i++)
        //     if (contains(i))
        //         m_domain.insert(i);
        //     else
        //         m_domain.erase(i);
#endif
    }
    
    bool operator==(const integer_domain<T> & t) {
        return m_empty == t.m_empty && m_endpoints() == t.m_endpoints();
    }
    bool contains_all() const {
        return !m_empty  && m_endpoints.empty();
    }

    bool contains(const T & x) const {
        if (contains_all())
            return true;
        bool neg_inf;
        const_iter l;
        bool found_left_point = get_left_point(x, neg_inf, l);
        if (!found_left_point)
            return has_neg_inf();
        if (neg_inf)
            return true;
        if (pos(l) == x)
            return true;
        return is_proper_start(l);
    }

    void handle_right_point_in_union(const_iter &r, const T &y) {
        if (pos(r) == y) {
            if (is_proper_start(r))
                erase(r);
            else
                set_end(y);
        }
        else if (pos(r) == y + 1) {
            if (is_proper_start(r)) {
                erase(r);
            }
            else {
                set_end(r->first);
            }
        }
        else if (!is_proper_end(r))
            set_end(y);

        lp_assert(is_correct());
    }
    void handle_left_point_in_union(const_iter& l, const T &x, const T & y) {
        if (pos(l) == x || pos(l) + 1 == x) {
            if (is_proper_end(l)) {
                l++;
                erase(std::prev(l));
            }
            else {
                if (is_one_point_interval(l)) {
                    set_start(pos(l));
                }
            }
        }
        else {
            if (!is_proper_start(l)) {
                set_start(x);
            }
        }
            
        while (l!= m_endpoints.end() && pos(l) <= x)
            l++;
        remove_from_the_left(y, l);
    }

    void unite_with_interval(const T& x, const T& y) {
        lp_assert(false);
        //         TRACE("disj_intervals", tout << "unite_with_interval(" << x << ", " << y << ")\n";);
        // #if Z3DEBUG
        //         // for (int i = std::max(x, 0); i <= std::min(100, y); i++)
        //         //     m_domain.insert(i);
        // #endif

        //         lp_assert(x <= y);
        //         if (x == y) {
        //             unite_with_one_point_interval(x);
        //             return;
        //         }

        //         const_iter l, r;
        //         bool neg_inf, pos_inf;
        //         bool found_left_point = get_left_point(x, neg_inf, l);
        //         bool found_right_point = get_right_point(y, pos_inf, r);
        //         m_empty = false;

        //         if (!found_left_point) {
        //             if (!found_right_point) {
        //                 m_endpoints.clear();
        //                 set_start(x);
        //                 set_end(y);
        //                 return;
        //             }
        //             // found_right_point is true
        //             if (pos_inf) {
        //                 m_endpoints.clear();
        //                 set_start(x);
        //                 return;
        //             }
        //             remove_from_the_left(y);
                        
        //             if (pos(m_endpoints.begin()) == y || pos(m_endpoints.begin()) == y + 1) {
        //                 if (is_proper_start(m_endpoints.begin()))
        //                     m_endpoints.erase(m_endpoints.begin());
        //                 else 
        //                     set_end(pos(m_endpoints.begin()));
        //                 set_start(x);
        //             }
        //             else {
        //                 lp_assert(pos(m_endpoints.begin()) > y + 1);
        //                 if (is_start(m_endpoints.begin()))
        //                     set_end(y);
        //                 set_start(x);
        //             }
        //             return;
        //         }

        //         lp_assert(found_left_point);
        //         if (!found_right_point) {
        //             bool neg_inf = has_neg_inf();
        //             remove_from_the_right(x);
        //             if (m_endpoints.empty()) {
        //                 if (!neg_inf)
        //                     set_start_end(x, y);
        //                 else
        //                     set_end(y);
        //                 return;
        //             }
        //             if (pos(m_endpoints.rbegin()) == x || pos(m_endpoints.rbegin()) == x - 1) {
        //                 if (is_proper_end(m_endpoints.rbegin())) {
        //                     m_endpoints.erase(m_endpoints.rbegin());
        //                 }
        //                 else if (is_one_point_interval(m_endpoints.rbegin())) {
        //                     set_start(pos(m_endpoints.rbegin()));
        //                 }
        //                 set_end(y);
        //             }
        //             else {
        //                 if (is_end(m_endpoints.rbegin())) {
        //                     set_start(x);
        //                     set_end(y);
        //                 }
        //                 else {
        //                     set_end(y);
        //                 }
        //             }
        //             return;
        //         }

        //         // found_right_point and found_left_point
        //         if (!neg_inf)
        //             handle_left_point_in_union(l, x, y);
        //         else {
        //             remove_from_the_left(y);
        //         }
        //         if (!pos_inf)
        //             handle_right_point_in_union(r, y);
        //         else
        //             remove_from_the_right(x);
    }

    bool has_pos_inf() const {
        if (m_empty)
            return false;
        
        if (m_endpoints.empty())
            return true;
        
        lp_assert(m_endpoints.rbegin() != m_endpoints.rend());
        return m_endpoints.rbegin()->second.kind() == endpoint_kind::START;
    }

    bool has_neg_inf() const {
        if (m_empty)
            return false;

        if (m_endpoints.empty())
            return true;
        auto it = m_endpoints.begin();
        return is_proper_end(it->second.kind());//m_endpoints.begin());
    }
    
    bool is_correct() const {
        if (m_empty) {
            if (m_endpoints.size() > 0) {
                TRACE("disj_intervals", tout << "is empty is true but m_endpoints.size() = " << m_endpoints.size() << std::endl;);
                return false;
            }
            return true;
        }
        bool expect_end;
        bool prev = false;
        T prev_x;
        for (auto t : m_endpoints()) {
            if (prev && ((expect_end && !is_end(t.second)) || (!expect_end && !is_start(t.second)))) {
                TRACE("disj_intervals", tout << "x = " << t.first << "\n";);
                if (expect_end) {
                    TRACE("disj_intervals", tout << "expecting an interval end\n";);
                } else {
                    TRACE("disj_intervals", tout << "expecting an interval start\n";);
                }
                return false;
            }

            if (prev) {
                if (t.first - prev_x <= 1 && !expect_end) {
                    TRACE("disj_intervals", tout << "the sequence is not increasing or the gap is too small: " << prev_x << ", " << t.first << std::endl;);
                    return false;
                }
            } 
            if (t.second.kind() == endpoint_kind::STEND) {
                expect_end = false; // swallow a point interval
            }
            else {
                if (prev)
                    expect_end = !expect_end;
                else
                    expect_end = is_start(t.second);
            }
            prev = true;
            prev_x = t.first;
        }
#if Z3DEBUG
        // for (int i = 0; i <= 100; i++ ) {
        //     if ( (m_domain.find(i) != m_domain.end()) != contains(i)) {
        //         TRACE("disj_intervals", tout << "incorrect value of contains(" << i << ") is = " << contains(i) << std::endl;);
        //         return false;
        //     }
        // }
#endif
        return true;
    }
public:
    void print(std::ostream & out) const {
        if (m_empty) {
            out << "empty\n";
            return;
        }
        if (m_endpoints.empty()){
            out << "[-oo,oo]\n";
            return;
        }
        bool first = true;
        for (auto t : m_endpoints()) {
            if (first) {
                if (t.second.kind() == endpoint_kind::END) {
                    out << "[-oo," << t.first; t.second.print(out); tout << "]";
                }
                else if (t.second.kind() == endpoint_kind::START) {
                    out << "[" << t.first; t.second.print(out); tout << ",";
                } else if (t.second.kind() == endpoint_kind::STEND) {
                    out << "[" << t.first; t.second.print(out); tout << "]";
                }
                first = false;
            } else {
                if (t.second.kind() == endpoint_kind::START) {
                    out << "[" << t.first; t.second.print(out); tout << ",";
                }
                else if (t.second.kind() == endpoint_kind::END) {
                    out << t.first; t.second.print(out); tout << "]";
                }
                else if (t.second.kind() == endpoint_kind::STEND) {
                    out << "[" << t.first; t.second.print(out); tout << "]";;
                }
            }
        }
        if (has_pos_inf())
            out << "oo]";
        
        out << "\n";
    }
    
    void push() { m_endpoints.push(); m_empty.push(); }
    void pop() { m_endpoints.pop(); m_empty.pop(); }
    void pop(unsigned k) { while(k--) pop(); }

    bool intersect_with_bound(const T & x, bool is_lower, unsigned explanation) {
        return is_lower? intersect_with_lower_bound(x, explanation) : intersect_with_upper_bound(x, explanation);
    }
    // we intersect the existing set with the half open to the right interval
    // returns true if the domain changes
    bool intersect_with_lower_bound(const T& x, unsigned explanation) {
#ifdef Z3DEBUG
        // for (int i = 0; i < x; i++)
        //     m_domain.erase(i);
#endif
        TRACE("disj_intervals", tout << "intersect_with_lower_bound(" << x << ")\n";);

        if (m_empty)
            return false;
        if (m_endpoints.empty()) {
            set_start(x, explanation);
            return true;
        }
        bool pos_inf = has_pos_inf();
        auto it = m_endpoints.begin();
        while (it != m_endpoints.end() && pos(it) < x) {
            m_endpoints.erase(it);
            it = m_endpoints.begin();
        }
        if (m_endpoints.empty()) {
            if (!pos_inf) {
                m_empty = true;
                return true;
            } 
            set_start(x, explanation);
            return true;
        }
        lp_assert(pos(it) >= x);
        if (pos(it) == x) {
            if (is_proper_end(it)) {
                set_start(x, explanation);
                return true;
            }
        }
        else { // x(it) > x
            if (is_proper_end(it)) {
                set_start(x, explanation);
                return true;
            }
        }
        return false;
    }
public:
    bool intersection_with_upper_bound_is_empty(const T& x) const {
        if (has_neg_inf())
            return false;
        if (m_empty)
            return true;
        T b;
        lp_assert(get_lower_bound(b));
        get_lower_bound(b);
        return x < b;
    }

    bool intersection_with_lower_bound_is_empty(const T& x) const {
        if (has_pos_inf())
            return false;
        if (m_empty)
            return true;
        T b;
        lp_assert(get_upper_bound(b));
        get_upper_bound(b);
        return x > b;
    }

    // we intersect the existing set with the half open interval
    // returns true if there is a change
    bool intersect_with_upper_bound(const T& x, unsigned explanation) {
#ifdef Z3DEBUG
        // for (int i = 100; i > x; i--)
        //     m_domain.erase(i);
#endif
        TRACE("disj_intervals", tout << "intersect_with_upper_bound(" << x << ")\n";);
        if (m_empty)
            return false;
        if (m_endpoints.empty()) {
            set_end(x, explanation);
            return true;
        }
        bool neg_inf = has_neg_inf();
        auto it = m_endpoints.rbegin();

        bool change = false;
        while (!m_endpoints.empty() && pos(it) > x) {
            m_endpoints.erase(std::prev(m_endpoints.end()));
            it = m_endpoints.rbegin();
            change = true;
        }
        if (m_endpoints.empty()) {
            if (!neg_inf) {
                m_empty = true;
                return true;
            }
            set_end(x, explanation);
            change = true;
        }
        lp_assert(pos(it) <= x);
        if (pos(it) == x) {
            if (is_one_point_interval(it)) {} 
            else if (is_proper_end(it)) {}
            else {// is_proper_start(it->second)
                set_end(x, explanation);
                change = true;
            }
        }
        else { // pos(it) < x} 
            if (is_proper_start(it)) {
                set_end(x, explanation);
                change = true;
            }
        }
        lp_assert(is_correct());
        return change;
    }
public:
    void intersect_with_interval(const T& x, const T & y) {
#ifdef Z3DEBUG
        // for (int i = 0; i <= 100; i++)
        //     if (i < x || i > y)
        //         m_domain.erase(i);
#endif

        TRACE("disj_intervals", tout << "intersect_with_interval(" << x << ", " << y <<")\n";);
        if (m_empty)
            return;
        lp_assert(x <= y);
        intersect_with_lower_bound(x);
        intersect_with_upper_bound(y);
    }

    // add an intervar [x, inf]
    void unite_with_interval_x_pos_inf(const T& x) {
        lp_assert(false);
        //         if (contains_all())
        //             return;
        // #if Z3DEBUG
        //         // for (int i = x; i <= 100; i++)
        //         //     m_domain.insert(i);
        // #endif
        //         TRACE("disj_intervals", tout << "unite_with_interval_x_pos_inf(" << x << ")\n";);
        //         if (m_empty) {
        //             set_start(x);
        //             m_empty = false;
        //             return;
        //         }
        //         bool neg_inf = has_neg_inf();
        //         remove_from_the_right(x);
        //         if (m_endpoints.empty()) {
        //             if (!neg_inf)
        //                 set_start(x);
        //             return;
        //         }
        //         auto it = m_endpoints.rbegin();
        //         lp_assert(pos(it) <= x);
        //         if (pos(it) == x) {
        //             if (is_proper_end(it)) {
        //                 m_endpoints.erase(x);
        //             } else {
        //                 set_start(x);
        //             }
        //         } else if (pos(it) == x - 1 && is_end(it)) {
        //             if (is_proper_start(it)) {
        //                 // do nothing
        //             }
        //             else if (is_proper_end(it)) {
        //                 m_endpoints.erase(it);
        //             }
        //             else {
        //                 lp_assert(is_one_point_interval(it));
        //                 set_start(it);
        //             }
        //         } else {
        //             if (!has_pos_inf())
        //                 set_start(x);
        //         }
    }

    // add an interval [-inf, x]
    void unite_with_interval_neg_inf_x(const T& x) {
        lp_assert(false); // not implemented
        // #if Z3DEBUG
        //         // for (int i = 0; i <= x; i++)
        //         //     m_domain.insert(i);
        // #endif
        //         TRACE("disj_intervals", tout << "unite_with_interval_neg_inf_x(" << x << ")\n";);
        //         if (m_empty) {
        //             set_end(x);
        //             m_empty = false;
        //             return;
        //         }
        //         bool pos_inf;
        //         const_iter r;
        //         bool found_right_point = get_right_point(x, pos_inf, r);
        //         if (!found_right_point) {
        //             m_endpoints.clear();
        //             set_end(x);
        //             return;
        //         }
        //         if (pos_inf) {
        //             m_endpoints.clear();
        //             return;
        //         }
        //         lp_assert(pos(r) >= x);
        //         if (pos(r) == x || pos(r) == x + 1) {
        //             if (is_proper_start(r))
        //                 erase(r);
        //             else if (is_one_point_interval(r)) {
        //                 set_end(pos(r));
        //             } // do nothing for the proper end
        //         } else {
        //             if (!is_proper_end(r))
        //                 set_end(x);
        //         }
        
        //         while (!m_endpoints.empty() && m_endpoints.begin()->first < x) {
        //             m_endpoints.erase(m_endpoints.begin());
        //         }
        //         lp_assert(is_correct());
    }

private:
    bool is_start(endpoint_kind x) const { return x == endpoint_kind::START || x == endpoint_kind::STEND; }
    bool is_start(const iter & it) const { return is_start(it->second);  }
    bool is_start(const const_iter & it) const { return is_start(it->second);  }
    bool is_start(const riter & it) const { return is_start(it->second);  }
    bool is_start(const endpoint & e) const { return is_start(e.kind());  }
   
    bool is_end(endpoint_kind x) const { return x == endpoint_kind::END || x == endpoint_kind::STEND; }
    bool is_end(const iter & it) const { return is_end(it->second);  }
    bool is_end(const const_iter & it) const { return is_end(it->second); }
    bool is_end(const riter & it) const {  return is_end(it->second); }
    bool is_end(const endpoint& e) const {  return is_end(e.kind()); }

    
    T pos(const iter & it) const {
        return it->first;
    }
    T pos(const const_iter & it) const {
        return it->first;
    }
    T pos(const riter & it) const {
        return it->first;
    }

    T bound_kind(iter & it) const {
        return it->second;
    }

    T bound_kind(riter & it) const {
        return it->second;
    }

    bool is_proper_start(endpoint_kind x) const { return x == endpoint_kind::START; }
    bool is_proper_start(const riter &x) const { return is_proper_start(x->second);}
    bool is_proper_start(const iter &x) const { return is_proper_start(x->second);}
    bool is_proper_start(const const_iter &x) const { return is_proper_start(x->second);}
    bool is_proper_start(const endpoint &x) const { return is_proper_start(x.kind());}
        
    bool is_proper_end(endpoint_kind x) const { return x == endpoint_kind::END; }
    bool is_proper_end(const iter & it) const { return is_proper_end(it->second.kind()); }
    bool is_proper_end(const const_iter & it) const { return is_proper_end(it->second); }
    bool is_proper_end(const riter & it) const { return is_proper_end(it->second); }
    bool is_proper_end(const endpoint & x) const { return is_proper_end(x.kind()); }

    bool is_one_point_interval(const endpoint & x) const { return is_one_point_interval(x.kind()); }
    bool is_one_point_interval(endpoint_kind x) const { return x == endpoint_kind::STEND; }
    bool is_one_point_interval(const iter & it) const { return is_one_point_interval(it->second); }
    bool is_one_point_interval(const const_iter & it) const {
        return is_one_point_interval(it->second);
    }
    bool is_one_point_interval(const riter & it) const {
        return is_one_point_interval(it->second);
    }
    
    void erase(const iter& it) {
        m_endpoints.erase(it);
    }
    void erase(const const_iter& it) {
        m_endpoints.erase(it);
    }

    void erase(const riter& it) {
        m_endpoints.erase(it);
    }

    void erase(const T& x) {
        m_endpoints.erase(x);
    }
    
    /*    void set_one_point_interval(const T& x, unsigned explanation) {
          auto it = m_endpoints().find(x);
          set_one_point_interval(it, explanation);
          }*/

    void set_start(const_iter &t, unsigned explanation) {
        lp_assert(t != m_endpoints.end());
        endpoint e = t->second;
        e.m_start_expl = explanation;
        m_endpoints[t->first] = e;
    }

    void set_start(const T& x, unsigned explanation) {
        endpoint e = get_endpoint(x);
        e.m_start_expl = explanation;
        m_endpoints[x] = e;
    }

    endpoint get_endpoint(const T& x) const {
        auto it = m_endpoints().find(x);
        if (it == m_endpoints().end())
            return endpoint();
        return it->second;
    }
    
    void set_end(const T& x, unsigned expl) {
        endpoint e = get_endpoint(x);
        e.m_end_expl = expl;
        m_endpoints[x] = e;
    }

    void set_end(const_iter& t, unsigned expl) {
        endpoint e = t->second;
        e.m_end_expl = expl;
        m_endpoints[t->first] =  e;
    }

    void set_end(riter& t, unsigned explanation) {
        endpoint e = t->second;
        e.m_end_expl = expl;
        m_endpoints[t->first] =  e;
    }

private:
    /*    void set_start_end(const T& x, const T & y, unsigned expl) {
          set_start(x, expl);
          set_end(y, expl);
          }*/

    void unite_with_one_point_interval(const T &x) {
        TRACE("disj_intervals", tout << "unite_with_one_point_interval(" << x << ")\n";);
        if (m_empty) {
            m_empty = false;
            set_one_point_interval(x);
            return;
        }
        bool has_left_neig, has_right_neig;
        bool neg_inf, pos_inf;
        const_iter l, r;
        has_left_neig = get_left_point(x, neg_inf, l);
        has_right_neig = get_right_point(x, pos_inf, r);
        if (!has_left_neig) {
            if (!has_right_neig) {
                set_one_point_interval(x);
                return;
            }
            lp_assert(!neg_inf);
            // if (pos(r) == x ) nothing happens
            if (pos(r) == x + 1) {
                if (is_one_point_interval(r)) {
                    set_start(x);
                    set_end(x + 1);
                } else {
                    lp_assert(is_proper_start(r));
                    erase(r);
                    set_start(x);
                }
            } else {
                lp_assert(pos(r) > x);
                set_one_point_interval(x);
            }
            return;
        }

        // has_left_neig
        if (!has_right_neig) {
            if (pos_inf)
                return;
            // if (pos(l) == x nothing happens
            if (pos(l) == x - 1) {
                if (is_proper_end(l)) {
                    erase(l);
                    set_end(x);
                } else {
                    lp_assert(is_one_point_interval(l));
                    set_start(l->first);
                    set_end(x);
                }
            } else {
                lp_assert(pos(l) < x - 1);
                set_one_point_interval(x);
            }
            return;
        }
        // has both neighbors
        if (neg_inf || pos_inf)
            return;
        if (pos(l) == x || pos(r) == x)
            return;

        // now the cases pos(l) == pos(r) or  pos(l) + 1 == pos(r) are impossible
        if (is_proper_start(l)) {
            lp_assert(is_proper_end(r));
            return;
        }
            
        if (pos(l) + 2 < pos(r)) { // we can glue only to one neighbor
            if (pos(l) + 1 == x) {
                if (is_proper_end(l)) {
                    erase(l);
                    set_end(x);
                }
                else {
                    lp_assert(!is_proper_start(l));
                    set_start(pos(l));
                    set_end(x);
                }
            }
            else if (x + 1 == pos(r)) {
                if (is_proper_start(r)) {
                    erase(r);
                }
                else if (is_one_point_interval(r)) {
                    set_end(pos(r));
                } 
                set_start(x);
            }
            else {
                set_one_point_interval(x);
            }
        } else {
            lp_assert(pos(l) + 2 == pos(r));
            lp_assert(pos(l) + 1 == x); // x is just in between l and r
            if (is_proper_end(l)) {
                erase(l);
            } else {
                lp_assert(is_one_point_interval(l));
                set_start(l->first);
            }
            if (is_proper_start(r)) {
                erase(r);
            } else {
                lp_assert(is_one_point_interval(r));
                set_end(r->first);
            }
        }

    }
    // return false if there are no points y in m_endpoints with pos(y) <= x
    // and negative infiniti is not true
    bool get_left_point(const T& x, bool &neg_inf, const_iter & l) const {
        if (m_empty) {
            neg_inf = false;
            return false;
        }

        if (m_endpoints.empty()) {
            neg_inf = true;
            return true;
        }

        l = m_endpoints.lower_bound(x);
        if (l == m_endpoints.end()) {
            l--;
            neg_inf = false;
            return true;
        }
        if (pos(l) == x) {
            neg_inf = false;
            return true;
        }
        if (l == m_endpoints.begin())
            return neg_inf = has_neg_inf();
        l--;
        neg_inf = false;
        return true;
    }
    
    // return false iff there are no points y in m_endpoints with pos(y) >= x, and positive infinity is not true
    bool get_right_point(const T& x, bool &pos_inf, const_iter & r) const {
        if (m_empty) {
            pos_inf = false;
            return false;
        }
            
        if (m_endpoints.empty()) {
            pos_inf = true;
            return true;
        }
        
        r = m_endpoints.lower_bound(x);
        if (r == m_endpoints.end()) {
            return pos_inf = has_pos_inf();;
        }
        
        pos_inf = false;
        return true;
    }
    


    // inserts x, if needed and return a pointer to the leftmost point that is a candidate to removal
    /*
      iter insert_start_for_unite_with_interval(const T& x) {
      auto lbx = m_endpoints.lower_bound(x);
      if (lbx == m_endpoints.end()) {
      if (!has_pos_inf()) {
      T lastx = pos(m_endpoints.rbegin());
      if (lastx + 1 < x) 
      set_start_end(x, y);
      else {
      m_endpoints.erase(lastx);
      set_end(y);
      }
      }
      return;
      }
      if (pos(lbx) == x) {
      if(is_proper_end(lbx)) {
      m_endpoints.erase(x);
      } else { set_start(x);}
      }
      if (pos(lbx) > x) {
      if (pos(lbx) > y + 1) {
      if (is_end(lbx))
      return;
      set_start_end(x, y);
      return;
      }
      if (pos(lpx) == y + 1) {
      if (is_end(lbx))
      }
                
      }
        
      lp_assert(false); // not implemented
      }
    */

    void remove_from_the_left(const T & y ) {
        while (!m_endpoints.empty() && pos(m_endpoints.begin()) < y) {
            m_endpoints.erase(m_endpoints.begin());
        }
    }

    void remove_from_the_left(const T & y, const_iter& l ) {
        while (l!= m_endpoints.end() && pos(l) < y) {
            l++;
            erase(std::prev(l));
        }
    }

    
    void remove_from_the_right(const T & x) {
        while (!m_endpoints.empty() && pos(m_endpoints.rbegin()) > x) {
            m_endpoints.erase(m_endpoints.rbegin());
        }
    }

    void remove_from_the_right(const T & x, riter & r) {
        while (!m_endpoints.empty() && pos(r) > x) {
            r--;
            m_endpoints.erase(std::next(r));
        }
    }
public:
    bool get_lower_bound_with_expl(T& b, unsigned & expl) const {
        if (m_empty)
            return false;
        if (has_neg_inf())
            return false;
        expl = m_endpoints.begin()->second.m_start_expl;
        if (expl == static_cast<unsigned>(-1))
            return false;
        b = pos(m_endpoints.begin());
        return true;
    }

    bool get_lower_bound(T& b) const {
        if (m_empty)
            return false;
        if (has_neg_inf())
            return false;
        b = pos(m_endpoints.begin());
        return true;
    }

    int get_lower_bound_expl() const {
        if (m_empty)
            return -1;
        if (has_neg_inf())
            return -1;
        return m_endpoints.begin()->second.m_start_expl;
    }

    int get_upper_bound_expl() const {
        if (m_empty)
            return -1;
        if (has_pos_inf())
            return -1;
        return m_endpoints.rbegin()->second.m_end_expl;
    }
    
    bool get_upper_bound_with_expl(T& b, unsigned & expl) const {
        if (m_empty)
            return false;
        if (has_pos_inf())
            return false;
        expl = m_endpoints.rbegin()->second.m_end_expl;
        if (expl == static_cast<unsigned>(-1))
            return false;
        b = m_endpoints.rbegin()->first;
        return true;
    }

    bool get_upper_bound_and_kind_with_expl(T& b, endpoint_kind & kind, unsigned & expl) const {
        if (m_empty)
            return false;
        if (has_pos_inf())
            return false;
        b = m_endpoints.rbegin()->first;
        kind = m_endpoints.rbegin()->second.kind();
        expl = m_endpoints.rbegin()->second.m_explanation;
        return true;
    }

    bool get_upper_bound(T& b) const {
        if (m_empty)
            return false;
        if (has_pos_inf())
            return false;
        b = m_endpoints.rbegin()->first;
        return true;
    }

    bool is_empty() const { return m_empty; }

    bool is_fixed() const {
        if (has_pos_inf() || has_neg_inf())
            return false;
        T l;
        get_lower_bound(l);
        T u;
        get_upper_bound(u);
        return l==u;
    }


    bool improves_with_lower_bound(const T & v) const {
        T b;
        bool lower_bound_exists = get_lower_bound(b);
        return (!lower_bound_exists || v > b) &&
            !intersection_with_lower_bound_is_empty(v);
    }

    bool improves_with_upper_bound(const T & v) const {
        T b;
        bool upper_bound_exists = get_upper_bound(b);
        return (!upper_bound_exists || v < b) &&
            !intersection_with_upper_bound_is_empty(v);
    }
    
    // returns true if adding the bound b narrows the domain, but does not make it empty
    bool improves(const T & v, bool is_lower_bound) const {
        if (is_lower_bound)
            return improves_with_lower_bound(v);
        return improves_with_upper_bound(v);
    }
};
}
