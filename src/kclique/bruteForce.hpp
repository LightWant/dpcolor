#ifndef BRUTEFORCE_HPP
#define BRUTEFORCE_HPP

#include "../graph/graph.hpp"
#include "../tools/type.hpp"
#include "../tools/hopstotchHash.hpp"
#include "../tools/linearSet.hpp"

#include <cassert>
#include <tuple>
#include <random>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>

using Pair = std::pair<v_size, v_size>;
using std::vector;
#define pP first
#define pR second

class bruteForce {
private:
    Graph * g;
    v_size ** S;
    v_size k;
    hopstotchHash * hashTable;

public:
    bruteForce(Graph * g_) {g = g_;}
    ~bruteForce() {
        delete [] hashTable;
        for(v_size i = 0; i < k; i++) delete [] S[i];
        delete [] S;
    }
    bool cc(v_size u, v_size v) {
        return hashTable[u].contain(v);
    }

    void previousWork() {
        hashTable = new hopstotchHash[g->vCnt];
        for(v_size u = 0; u < g->vCnt; u++) {
            if(g->pIdx[u + 1] == g->pIdx[u]) continue;
            hashTable[u].build(g->pEdge + g->pIdx[u], 
                g->pIdx[u + 1] - g->pIdx[u]);
        }

        S = new v_size*[k];
        for(v_size i = 0; i < k; i++) {
            S[i] = new v_size[g->degeneracy];
        }
    }

    void run(v_size k_) {
        k = k_;

        double t = clock();

        previousWork();
        double sum = 0.0;

        for(v_size u = 0; u < g->vCnt; u++) {
            v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];
            if(outDegree + 1 < k) continue;

            memcpy(S[0], g->pEdge + g->pIdx2[u], sizeof(v_size) * outDegree);
            
            double num = search(outDegree, 0);
            sum += num;
            
           // printf("%u %.0f\n", u, num);
        }

        printf("time:%.2fs\n", (clock() - t)/CLOCKS_PER_SEC);
        printf("%u-Clique:%.0f\n", k, sum);

        fflush(stdout);
    }

    double search(v_size n, v_size d) {
        if(n == 0) return 0.0;

        if(d + 2 == k) {
            return n;
        }

        double ans = 0.0;
        for(v_size i = 0; i < n; i++) {
            v_size u = S[d][i];
            
            v_size l = 0;
            for(v_size j = i + 1; j < n; j++) {
                if(cc(u, S[d][j])) {
                    S[d + 1][l++] = S[d][j];
                }
            }

            ans += search(l, d + 1);
        }

        return ans;
    }
};

#endif