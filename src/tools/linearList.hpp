#ifndef LINEARLIST
#define LINEARLIST

#include "type.hpp"
#include <vector>
#include <cstdio>

// struct linearListNode {
//     v_size v;
//     int nxt;

//     linearListNode():v(0), nxt(-1) {}
// };

// struct listHeader {
//     int bg, ed, sz;
// };

// class linearList {
//     linearListNode * pool;
//     v_size n;
//     std::vector<listHeader> heads;

// public:
//     linearList(v_size n_) {
//         n = n_;
//         pool = new linearListNode[n];
//     }

//     insert()
// };

// class memoryPool {
//     v_size * bg;
//     v_size m, d;


//     memoryPool(v_size m_, v_size d_) {
//         //Graph edge number
//         m = m_; d = d_;
//         bg = new v_size[m];
//     }

//     v_size * alloc() {
        
//     }
// };

class clique {
    v_size *bg = nullptr;
    v_size minD;
    v_size d;
    v_size *pMinD = nullptr;
    v_size l;

public:
    clique(v_size d_) {
        bg = new v_size[d_];
        d = minD = d_;
        l = 0;
    }

    ~clique() { delete [] bg; }

    void push_back(v_size v, v_size degree) {
        bg[l++] = v;
        if(degree < minD) {
            minD = degree; pMinD = bg + l - 1;
        }
    }

    v_size VectexWithMinD() { return *pMinD; }

    v_size * begin() { return bg; }
    v_size * end() {return bg + l;}

    bool empty() {return l == 0;}

    void clear() {l = 0; minD = d;}

    void print() {
        for(v_size i = 0; i < l; i++) {
            printf("%u ", bg[i]);
        }
        printf("\n");
    }

    v_size size() {return l;}

    v_size operator [] (v_size v) {
        return bg[v];
    }
};

#endif