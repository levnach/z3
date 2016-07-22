/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
*/

#pragma once
// this class implements a map with some stack functionality
#include <unordered_map>
#include <set>
#include <stack>
namespace lean {


template <typename A, typename B,
          typename Hash = std::hash<A>,
          typename KeyEqual = std::equal_to<A>,
          typename Allocator = std::allocator< std::pair<const A, B> >
          > class stacked_map {
    struct delta {
        std::unordered_set<A, Hash, KeyEqual> m_new;
        std::unordered_map<A,B, Hash, KeyEqual, Allocator> m_original_changed;
    };
    std::unordered_map<A,B,Hash, KeyEqual, Allocator> m_map;
    std::stack<delta> m_stack;
public:
    class ref {
        stacked_map<A,B,Hash, KeyEqual, Allocator> & m_map;
        const A & m_key;
    public:
        ref(stacked_map<A,B,Hash, KeyEqual, Allocator> & m, const A & key) :m_map(m), m_key(key) {}
        ref & operator=(const B & b) {
            m_map.emplace_replace(m_key, b);
            return *this;
        }
        ref & operator=(const ref & b) { lean_assert(false); return *this; }
        operator const B&() const {
            auto it = m_map.m_map.find(m_key);
            lean_assert(it != m_map.m_map.end());
            return it->second;
        }
    };
private:
    void emplace_replace(const A & a, const B & b)  {
        if (!m_stack.empty()) {
            auto it = m_map.find(a);
            if (it == m_map.end()) {
                m_stack.top().m_new.insert(a);
                m_map.emplace(a, b);
            } else if (it->second != b) {
                auto & orig_changed= m_stack.top().m_original_changed;
                auto itt = orig_changed.find(a);
                if (itt == orig_changed.end()) {
                    orig_changed.emplace(a, it->second);
                } else if (itt->second == b) {
                    orig_changed.erase(itt);
                }
                it->second = b;
            }            
        } else { // there is no stack: just emplace or replace
            m_map[a] = b;
        }
    }
public:    
    ref operator[] (const A & a) {
        return ref(*this, a);
    }

    const B & operator[]( const A & a) const {
        auto it = m_map.find(a);
        if (it == m_map.end()) {
            lean_assert(false);
        }

        return it->second;
    }

    bool try_get_value(A& key, B& val) const  {
        auto it = m_map.find(key);
        if (it == m_map.end())
            return false;
        
        val = it->second;
        return true;
    }
    bool try_get_value(A&& key, B& val) const  {
        auto it = m_map.find(std::move(key));
        if (it == m_map.end())
            return false;
        
        val = it->second;
        return true;
    }
    
    unsigned size() const {
        return m_map.size();
    }

    bool contains(A & key) const {
        return m_map.find(key) != m_map.end();
    }

    bool contains(A && key) const {
        return m_map.find(std::move(key)) != m_map.end();
    }
    
    void push() {
        delta d;
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
            for (auto & t : d.m_new) {
                m_map.erase(t);
            }
            for (auto & t: d.m_original_changed) {
                m_map[t.first] = t.second;
            }
            m_stack.pop();
        }
    }

    const std::unordered_map<A, B,Hash, KeyEqual, Allocator>& operator()() const { return m_map;}
};
}
