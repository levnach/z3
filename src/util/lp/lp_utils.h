/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
  This file should be present in z3 and in Lean.
  The idea is that it is only one different file in Lean and z3 source inside of LP
*/
#pragma once
#define lp_for_z3
#include <string>
#ifdef lp_for_z3
#include "../rational.h"
#include "../sstream.h"
#include "../z3_exception.h"
namespace lp {
    inline void throw_exception(const std::string & str) {
         throw default_exception(str);
    }
    template <typename T> class numeric_traits {};
    typedef rational mpq;
    typedef z3_exception exception;
    #ifdef LEAN_DEBUG
    inline void lean_assert(bool b) {}
    #else
#define lean_assert(_x_) {}
    #endif
    inline void lean_unreachable() { lean_assert(false); }
    template<>
    class numeric_traits<double> {
    public:
        static bool precise() { return false; }
        static double g_zero;
        static double const & zero() { return g_zero; }
        static double g_one;
        static double const & one() { return g_zero; }
        static bool is_zero(double v) { return v == 0.0; }
        static double const & get_double(double const & d) { return d;}
        static double log(double const & d) { NOT_IMPLEMENTED_YET(); return d;}
    };
    
    template<>
    class numeric_traits<rational> {
    public:
        static bool precise() { return true; }
        static rational const & zero() { return rational::zero(); }
        static rational const & one() { return rational::one(); }
        static bool is_zero(const rational & v) { return v.is_zero(); }
        static double const  get_double(const rational  & d) { return d.get_double();}
        static rational log(rational const& r) { UNREACHABLE(); return r; }
 };
}
namespace std {
template<>
struct hash<rational> {
    inline size_t operator()(const rational & v) const {
        return v.hash();
    }
};
}

template <class T>
inline void hash_combine(std::size_t & seed, const T & v) {
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
template<typename S, typename T> struct hash<pair<S, T>> {
    inline size_t operator()(const pair<S, T> & v) const {
        size_t seed = 0;
        hash_combine(seed, v.first);
        hash_combine(seed, v.second);
        return seed;
    }
};
}



#else // end of lp_for_z3
#include <utility>
#include <functional>
#include "util/numerics/mpq.h"
#ifdef __CLANG__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmismatched-tags"
#endif
namespace std {
template<>
struct hash<lean::mpq> {
    inline size_t operator()(const lean::mpq & v) const {
        return v.hash();
    }
};
}

template <class T>
inline void hash_combine(std::size_t & seed, const T & v) {
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
template<typename S, typename T> struct hash<pair<S, T>> {
    inline size_t operator()(const pair<S, T> & v) const {
        size_t seed = 0;
        hash_combine(seed, v.first);
        hash_combine(seed, v.second);
        return seed;
    }
};
}
#ifdef __CLANG__
#pragma clang diagnostic pop
#endif
#endif
