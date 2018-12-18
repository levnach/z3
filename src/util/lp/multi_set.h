/*++
  Copyright (c) 2017 Microsoft Corporation

  Module Name:

  <name>

  Abstract:

  <abstract>

  Author:
  Nikolaj Bjorner (nbjorner)
  Lev Nachmanson (levnach)

  Revision History:


  --*/
#include <unordered_map>

namespace nla {
template <typename T>
struct multi_set {
    std::unordered_map<T, int> m;

    template <typename K>
    multi_set(const K& c) {
        for (const auto & v : c) {
            add(v);
        }
    }
    void add(const T& e) {
        auto i = m.find(e);
        if (i == m.end()) {
            m[e] = 1;
        } else {
            i->second++;
        }
    }

    void remove(const T& e) {
        auto i = m.find(e);
        if (i == m.end()) {
            return; // noop
        } else {
            i->second--;
            if (i->second == 0)
                m.erase(i);
        }
    }
    bool contains(const T& e) const {
        return m.find(e) != m.end();
    } 
};
}
