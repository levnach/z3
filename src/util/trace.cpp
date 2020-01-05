/*++
Copyright (c) 2006 Microsoft Corporation

Module Name:

    trace.cpp

Abstract:

    Trace message support.

Author:

    Leonardo de Moura (leonardo) 2006-09-13.

Revision History:

--*/
#include "util/trace.h"
#include "util/str_hashtable.h"

#ifdef _TRACE

std::ofstream tout(".z3-trace"); 

static int g_enable_all_trace_tags = false;
static str_hashtable* g_enabled_trace_tags = nullptr;

static str_hashtable& get_enabled_trace_tags() {
    if (!g_enabled_trace_tags) {
        g_enabled_trace_tags = alloc(str_hashtable);
    }
    return *g_enabled_trace_tags;
}

void finalize_trace() {
    dealloc(g_enabled_trace_tags);
    g_enabled_trace_tags = nullptr;
}

void enable_trace(const char * tag) {
    int v = atoi(tag);
    if (v) {
        g_enable_all_trace_tags = v;
    } else {
        get_enabled_trace_tags().insert(tag);
    }
}

void enable_all_trace(bool flag) {
    g_enable_all_trace_tags = flag;
}

void disable_trace(const char * tag) {
    get_enabled_trace_tags().erase(tag);
}

bool is_trace_enabled(const char * tag) {
    if (g_enable_all_trace_tags) {
        if (get_enabled_trace_tags().contains(tag))
            return false;
        std::string stag(tag);
        std::string bug("bug");
        if (stag.find(bug) != std::string::npos)
            return false;
        g_enable_all_trace_tags--;
        if (g_enable_all_trace_tags == 0)
            exit(3);
        if (g_enable_all_trace_tags % 5)
            return false;
            
    }
    return g_enable_all_trace_tags || 
        (g_enabled_trace_tags && get_enabled_trace_tags().contains(tag));
}

void close_trace() {
    tout.close();
}

void open_trace() {
    tout.open(".z3-trace");
}

#endif
