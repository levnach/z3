/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
*/

#pragma once
// this class implements a vector with some stack functionality
#include <vector>
#include <set>
#include <stack>
#include <unordered_map>
#include "util/lp/stacked_map.h"
namespace lean {


template <typename B> class stacked_vector {
    struct delta {
        std::unordered_map<unsigned,B> m_original_changed;
        unsigned m_prev_size;
    };
    std::vector<B> m_vec;
    std::stack<delta> m_stack;
public:
    class ref {
        stacked_vector<B> & m_vec;
        unsigned m_key;
    public:
        ref(stacked_vector<B> & m, unsigned  key) :m_vec(m), m_key(key) {}
        ref & operator=(const B & b) {
            m_vec.replace(m_key, b);
            return *this;
        }
        ref & operator=(const ref & b) { lean_assert(false); return *this; }
        operator const B&() const {
            return m_vec.m_vec[m_key];
        }
    };
private:
    void replace(unsigned a, const B & b)  {
        if (!m_stack.empty()) {
            if (a < m_stack.top().m_prev_size) {
                if (m_vec[a] != b) {
                    auto & orig_changed= m_stack.top().m_original_changed;
                    auto itt = orig_changed.find(a);
                    if (itt == orig_changed.end()) {
                        orig_changed.emplace(a, m_vec[a]);
                    } else if (itt->second == b) {
                        orig_changed.erase(itt);
                    }
                }
            }
        }
        m_vec[a] = b;
    }
public:
    
    ref operator[] (unsigned  a) {
        lean_assert(a < size());
        return ref(*this, a);
    }

    const B & operator[](unsigned  a) const {
        lean_assert(a < size());
        auto it = m_vec.find(a);
        if (it == m_vec.end()) {
            lean_assert(false);
        }

        return it->second;
    }

    
    
    unsigned size() const {
        return m_vec.size();
    }

    void push() {
        delta d;
        d.m_prev_size = size();
        m_stack.push(d);
    }
    
    void pop() {
        pop(1);
    }
    void pop(unsigned k) {
        while (k-- > 0) {
            if (m_stack.empty())
                return;
            delta & d = m_stack.top();
            m_vec.resize(d.m_prev_size);

            for (auto & t: d.m_original_changed) {
                m_vec[t.first] = t.second;
            }
            m_stack.pop();
        }
    }

    void push_back(const B& b) {
        if (m_stack.empty()) {
            m_vec.push_back(b);
            return;
        }
        delta & d = m_stack.top();
        int i = m_vec.size();
        if (i < d.m_prev_size) {
            auto it = d.m_original_changed.find(i);
            lean_assert(it != d.m_original_changed.end());
            if (it->second == b) {
                d.m_original_changed.erase(it);
            }
        }
        m_vec.push_back(b);
    }

    void clear() {
        if (m_stack.empty()) {
            m_vec.clear();
            return;
        }

        delta & d = m_stack.top();
        int i = 0;
        for (auto & l: m_vec) {
            if (i >= d.m_prev_size) continue;
            auto it = d.m_original_changed.find(i);
            if (it == d.m_original_changed.end())
                d.m_original_changed.emplace(i, l);
            i++;
        }
        m_vec.clear();
    }
    
};
}
