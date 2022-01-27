#ifndef SHADOWBASEDONCOLORS_H
#define SHADOWBASEDONCOLORS_H


#include <cassert>
#include <tuple>
#include <random>
#include <vector>
#include <algorithm>
#include <random>

using Pair = std::pair<v_size, v_size>;
using std::vector;

struct sampleBasedOnColors {
//染色，dp，采样
//输入 按照颜色分类，颜色总数
    Graph * g = nullptr;
    v_size * pEdge = nullptr;
    v_size * pIdx = nullptr;
    double ** dp = nullptr;
    double * memoryPool = nullptr;
    v_size k;
    v_size * clique = nullptr;
    v_size * a = nullptr;
    hopstotchHash * hashTable;
    std::default_random_engine e;
    // std::function<bool(v_size, v_size)> connect;

    void init(Graph * g_, v_size k_, hopstotchHash * hashTable_) {
        g = g_; k = k_;
        clique = new v_size[k];
        hashTable = hashTable_;

        dp = new double*[g->cc + 1];
        memoryPool = new double[(g->cc + 1) * (k+1)]();
        v_size p = 0;
        for(v_size i = 0; i < g->cc + 1; i++) {
            dp[i] = memoryPool + p;
            p += k + 1;
        }

        a = new v_size[g->cc]();

        pEdge = new v_size[g->degeneracy*g->degeneracy];
        pIdx = new v_size[g->degeneracy];
    }

    ~sampleBasedOnColors() {
        if(clique != nullptr) delete [] clique;
        if(memoryPool != nullptr) delete [] memoryPool;
        if(dp != nullptr) delete [] dp;
        if(a != nullptr) delete [] a;
        if(pEdge != nullptr) delete [] pEdge;
        if(pIdx != nullptr) delete [] pIdx;
    }

    void computeDP(v_size u) {
        for(v_size i = 0; i < g->cc; i++) a[i] = 0;
        for(v_size i = g->pIdx2[u]; i < g->pIdx[u+1]; i++) {
            v_size v = g->pEdge[i];
            a[g->color[v]]++;
        }

        dp[0][0] = 1.0;
        for(v_size i = 1; i <= g->cc; i++) {
            dp[i][0] = 1.0;
            for(v_size j = 1; j <= k && j <= i; j++) {
                dp[i][j] = dp[i - 1][j - 1]*a[i - 1] + dp[i - 1][j];
            }
        }
    }

    void sortVbyColor(v_size u) {
        v_size deg = g->pIdx[u+1] - g->pIdx2[u];
        memcpy(pEdge, g->pEdge + g->pIdx2[u], 
            sizeof(v_size)*deg);
        
        auto cmp = [&](v_size a, v_size b) {
            return g->color[a] < g->color[b];
        };
        std::sort(pEdge, pEdge + deg, cmp);

        pIdx[g->color[pEdge[0]]] = 0;
// v_size tmp = 1;
        for(v_size i = 1; i < deg; i++) {
            if(g->color[pEdge[i]] != g->color[pEdge[i-1]]) {
                pIdx[g->color[pEdge[i]]] = i;
// assert(a[g->color[pEdge[i-1]]] == tmp);
// tmp = 1;
            }
            // else tmp++;
        }
    }

    bool connect(v_size u, v_size v) {
        return hashTable[u].contain(v);
    }
    v_size randChoose(v_size c) {
        v_size deg = a[c];
        return pEdge[pIdx[c] + std::rand() % deg];
    }

    int sampleOneTime() {
        v_size i = g->cc;
        v_size j = k;
        v_size l = 0;
        
        std::uniform_real_distribution<double> u(0, 1);

        while(l < k) {
            while(i > 0 && j > 0) {
                if(u(e) + 1e-8 < dp[i-1][j-1]*a[i-1]/dp[i][j]) {
                    clique[l++] = randChoose(i-1);
                    j--;
                }
                i--;
            }

            if(i == 0 && j > 0) {
                i = g->cc; j = k; l = 0;
            }
        }

        // if(l < k) return -1;

        for(v_size i = 0; i < l; i++) {
            for(v_size j = i + 1; j < l; j++) {
                if(!connect(clique[i], clique[j])) return 0;
            }
        }
        
        return 1;
    }

//     double sample(v_size u) {
//         // int sampleTimes = 1;
//         int sampleTimes = int(dp[g->cc][k]/1e9);
//         if(sampleTimes < 1) sampleTimes = 1;

//         int sampleTotalTimes = 0;
//         int t = 0;

//         computeDP(u);
// // printf("%.2f ", dp[g->cc][k]);
//         sortVbyColor(u);

//         while(sampleTotalTimes < sampleTimes) {
//             int ret = sampleOneTime();
            
//             // if(ret == -1) continue;
            
//             sampleTotalTimes++;
//             t += ret;
//         }

// // printf("%.2f ", 1.0 * t / sampleTimes);
//         return 1.0 * t / sampleTimes * dp[g->cc][k];
//     }
    double sample(v_size u) {
        // int sampleTimes = 1;
        int preT = 0, preTotals = 0;
        int t = 0, totals = 0;
        int cnt = 0;
        // constexpr int batchSize = 2;

        computeDP(u);
// printf("%.2f ", dp[g->cc][k]);
        sortVbyColor(u);

        int maxTime = 100;
        while(maxTime--) {

            t = 0; totals = 0;
            
            for(int i = 0; i < 2; i++) {
                int ret = sampleOneTime();
                // if(ret == -1) continue;
                totals++;
                t += ret;
            }
            

            if(preTotals == 0 || preT == 0) {
                preT += t; preTotals += totals;
            }
            else if(1.0*preT/preTotals-1.0*(preT+t)/(preTotals+totals) < 1e-7) {
                preT += t; preTotals += totals;
                cnt++;
                if(cnt > 2) {
                    cnt = 0;
                    break;
                }
            }
            else {
                preT += t; preTotals += totals;
                cnt = 0;
            }
        }

// printf("totals %d, %u\n", preTotals, (g->pIdx[u+1]-g->pIdx2[u])/10);
        return 1.0 * preT / preTotals * dp[g->cc][k];
    }

};

#endif