#ifndef CCCPATHPARALLEL2_H
#define CCCPAHTPARALLEL2_H

#include "../../graph/graph.hpp"
#include "../../tools/type.hpp"
#include "../../tools/hopstotchHash.hpp"
#include "../../tools/bitmap.hpp"
#include "../../tools/linearSet.hpp"

#include <cassert>
#include <tuple>
#include <random>
#include <vector>
#include <algorithm>
#include <random>

using Pair = std::pair<v_size, v_size>;
using std::vector;

// constexpr v_size batchSize = 50;

struct cccpathParallel {
    v_size sz;
    Graph * g;
    hopstotchHash * hashTable;
    v_size k;
    v_size threads;

    double * experiments;
    double sumW;
    double *** dp;

    v_size ** pEdge = nullptr;
    v_size ** pIdx = nullptr;
    v_size ** pColor = nullptr;
    v_size ** sortByColor = nullptr;
    v_size vCnt, eCnt;

    v_size ** clique = nullptr;
    v_size N = 5000000;

    v_size chunkSize = 32;
    v_size *st = nullptr, *ed = nullptr;// *p = nullptr;

    // std::default_random_engine e[200];
    // std::uniform_real_distribution<double>  ** uiDistribution;

    void init(v_size sz_, std::vector<v_size> & nodes, v_size N_=5000000) {
        sz = sz_;
        N = N_;
        experiments = new double[sz];
        
        st[0] = 0;
        e_size s = 0;
        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];
            v_size deg = g->pIdx[u + 1] - g->pIdx2[u];
            s += deg*deg + k*deg*deg;
        }
        s = s / threads;
        e_size ss = 0;
        v_size j = 1;
        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];
            v_size deg = g->pIdx[u + 1] - g->pIdx2[u];
            ss += deg*3+deg*deg+k*deg*deg;
            if(ss >= s * j) {
                st[j++] = i;
                if(j == threads) break;
            }
        }
        for(v_size i = 0; i+1 < threads; i++) {
            ed[i] = st[i + 1];
        }
        if(threads > 1)
            ed[threads - 1] = sz;

        #pragma omp parallel reduction(+:sumW)
        {
            int threadId = omp_get_thread_num();
            
            v_size i;

            while((i = __sync_fetch_and_add(&st[threadId], chunkSize)) < ed[threadId]) {
                v_size j = std::min(i + chunkSize, ed[threadId]);
                for(v_size x = i; x < j; x++) {
                    v_size u = nodes[x];
                    memset(pColor[threadId], 0, sizeof(v_size)*(g->cc+1));
                    v_size deg = g->pIdx[u+1] - g->pIdx2[u];
                    for(v_size i = 0; i < deg; i++) {
    // printf("%d %d\n", g->color[g->pEdge[g->pIdx2[u] + i]] + 1, g->cc);fflush(stdout);
                        pColor[threadId][ g->color[g->pEdge[g->pIdx2[u] + i]] + 1 ]++;
                    }

                    for(v_size i = 1; i < g->cc; i++) {
                        pColor[threadId][i] += pColor[threadId][i - 1];
                    }

                    for(v_size i = 0; i < deg; i++) {
                        v_size v = g->pEdge[g->pIdx2[u] + i];
                        sortByColor[threadId][ pColor[threadId][g->color[v]]++ ] = v;
                    }

                    memcpy(g->pEdge + g->pIdx2[u], sortByColor[threadId], sizeof(v_size)*deg);

                    computeDP(u, threadId);

                    double sumD = 0.0;
                    for(v_size j = 0; j < g->pIdx[u+1] - g->pIdx2[u]; j++) {
                        sumD += dp[threadId][j][k];
                    }

                    experiments[x] = sumD;
                    sumW += sumD;
                }
            }

            for(v_size tId = 0; tId < threads; tId++) {
                while((i = __sync_fetch_and_add(&st[tId], chunkSize)) < ed[tId]) {
                    v_size j = std::min(i + chunkSize, ed[tId]);
                    for(v_size x = i; x < j; x++) {
                        v_size u = nodes[x];
                        memset(pColor[threadId], 0, sizeof(v_size)*(g->cc+1));
                        v_size deg = g->pIdx[u+1] - g->pIdx2[u];
                        for(v_size i = 0; i < deg; i++) {
                            pColor[threadId][ g->color[g->pEdge[g->pIdx2[u] + i]] + 1 ]++;
                        }

                        for(v_size i = 1; i < g->cc; i++) {
                            pColor[threadId][i] += pColor[threadId][i - 1];
                        }

                        for(v_size i = 0; i < deg; i++) {
                            v_size v = g->pEdge[g->pIdx2[u] + i];
                            sortByColor[threadId][ pColor[threadId][g->color[v]]++ ] = v;
                        }

                        memcpy(g->pEdge + g->pIdx2[u], sortByColor[threadId], sizeof(v_size)*deg);

                        computeDP(u, threadId);

                        double sumD = 0.0;
                        for(v_size j = 0; j < g->pIdx[u+1] - g->pIdx2[u]; j++) {
                            sumD += dp[threadId][j][k];
                        }

                        experiments[x] = sumD;
                        sumW += sumD;
                    }
                }
            }

            delete [] sortByColor[threadId];
            delete [] pColor[threadId];
        }
    }

    void initForSingleNode(v_size k_, Graph * g_, 
        hopstotchHash * hashTable_, v_size threads_) {
        k = k_;
        g = g_;
        hashTable = hashTable_;
        threads = threads_;

        sumW = 0.0;
        clique = new v_size*[threads];
        for(v_size i = 0; i < threads; i++) {
            clique[i] = new v_size[k];
        }

        dp = new double**[threads];
        for(v_size t = 0; t < threads; t++) {
            // double * memoryPool = new double[g->degeneracy * (k+1)]();
            dp[t] = new double*[g->degeneracy];
            // v_size p = 0;

            for(v_size i = 0; i < g->degeneracy; i++) {
                dp[t][i] = new double[k+1];
                // p += k + 1;
            }
        }

        for(v_size t = 0; t < threads; t++) {
            for(v_size i = 0; i < g->degeneracy; i++) {
                dp[t][i][0] = 0;
                dp[t][i][1] = 1;
            }
        }

        pEdge = new v_size*[threads];
        pIdx = new v_size*[threads];
        for(v_size t = 0; t < threads; t++) {
            pEdge[t] = new v_size[g->degeneracy*g->degeneracy];
            pIdx[t] = new v_size[g->degeneracy + 1];
        }
        sortByColor = new v_size*[threads];
        for(v_size i = 0; i < threads; i++) {
            sortByColor[i] = new v_size[g->degeneracy + 1];
        }
        pColor = new v_size*[threads];
        for(v_size i = 0; i < threads; i++) {
            pColor[i] = new v_size[g->cc + 1];
        }

        st = new v_size[threads];
        ed = new v_size[threads];
    }

    ~cccpathParallel() {
        if(experiments != nullptr) delete [] experiments;
        for(v_size i = 0; i < threads; i++) {
            for(v_size j = 0; j < g->degeneracy; j++)
                delete [] dp[i][j];
            delete [] dp[i];
        }
        if(dp != nullptr) delete [] dp;

        for(v_size i = 0; i < threads; i++) {
            delete [] pEdge[i];
            delete [] pIdx[i];
        }
        if(pEdge != nullptr) delete [] pEdge;
        if(pIdx != nullptr) delete [] pIdx;

        for(v_size i = 0; i < threads; i++) {
            delete [] clique[i];
        }
        if(clique != nullptr) delete [] clique;
        if(sortByColor != nullptr) delete [] sortByColor;
        if(pColor != nullptr) delete [] pColor;

        if(st != nullptr) delete [] st;
        if(ed != nullptr) delete [] ed;
        // if(p != nullptr) delete [] p;
    }

    bool connect(v_size u, v_size v) {
        return hashTable[u].contain(v);
    }
    
    void computeDP(v_size u, v_size tId) {
        v_size outDegree = g->pIdx[u+1] - g->pIdx2[u];
        pIdx[tId][0] = 0;
        for(v_size i = 0; i < outDegree; i++) {
            v_size v = sortByColor[tId][i];
            pIdx[tId][i + 1] = pIdx[tId][i]; 

            for(v_size j = i + 1; j < outDegree; j++) {
                v_size w = sortByColor[tId][j];

                if(g->color[v] == g->color[w]) continue;
                if(hashTable[v].contain(w)) {
                    pEdge[tId][pIdx[tId][i + 1]++] = j;
                }
            }
        }

        for(v_size j = 2; j <= k; j++) {
            for(v_size i = 0; i < outDegree; i++) {
                dp[tId][i][j] = 0.0;
                for(v_size l = pIdx[tId][i]; l < pIdx[tId][i + 1]; l++) {
                    dp[tId][i][j] += dp[tId][pEdge[tId][l]][j - 1];
                }
            }
        }
    }

    int sampleOneTime(v_size id, v_size u, 
        std::uniform_real_distribution<double> & d,
        std::default_random_engine & e,
        v_size tId) {
        
        v_size preId = -1;

        double sumD = experiments[id];
        double x = d(e);
        
        double sumTmp = 0.0;
        for(v_size i = 0; i < g->pIdx[u+1] - g->pIdx2[u]; i++) {
            sumTmp += dp[tId][i][k];
            if(sumTmp + 1e-10 >= x * sumD) {
                clique[tId][0] = sortByColor[tId][ i ];
                preId = i;
                break;
            }
        }

        for(v_size i = 1; i < k; i++) {
            x = d(e);
            sumTmp = sumD = 0.0;
            
            for(v_size j = pIdx[tId][preId]; j < pIdx[tId][preId + 1]; j++) {
                sumD += dp[tId][pEdge[tId][j]][k - i];
            }

            for(v_size j = pIdx[tId][preId]; j < pIdx[tId][preId + 1]; j++) {
                sumTmp += dp[tId][pEdge[tId][j]][k - i];
                if(sumTmp + 1e-10 >= x * sumD) {
                    clique[tId][i] = sortByColor[tId][ pEdge[tId][j] ];;
                    preId = pEdge[tId][j];
                    break;
                }
            }

            for(v_size j = 0; j < i - 1; j++) {
                if(!connect(clique[tId][i], clique[tId][j])) {
                    return 0;
                }
            }
        }
        
        return 1;
    }

    double sample(std::vector<v_size> & nodes, e_size sampleTimes) {
        e_size t = 0;
        e_size sampleTotalTimes = 0;
        double ans = 0.0;

        std::random_device rd;
        std::default_random_engine e[100];
        for(v_size i = 0; i < threads; i++) {
            e[i].seed(rd());
        }

        st[0] = 0;
        e_size s = 0;
        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];
            v_size deg = g->pIdx[u + 1] - g->pIdx2[u];
            e_size expectedSampleTime 
                    = std::ceil(sampleTimes * experiments[i] / sumW);
            s += deg*deg + k*deg*deg+k*2*expectedSampleTime;
        }
        s = s / threads;
        e_size ss = 0;
        v_size j = 1;
        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];
            v_size deg = g->pIdx[u + 1] - g->pIdx2[u];
            e_size expectedSampleTime 
                    = std::ceil(sampleTimes * experiments[i] / sumW);
            ss += deg*deg + k*deg*deg+k*2*expectedSampleTime;
            if(ss >= s * j) {
                st[j++] = i;
                if(j == threads) break;
            }
        }
        for(v_size i = 0; i+1 < threads; i++) {
            ed[i] = st[i + 1];
        }
        if(threads > 1)
            ed[threads - 1] = sz;
        // work[i] = deg*deg+k*deg*deg+k*2*expectedSampleTime;

        #pragma omp parallel reduction(+:sampleTotalTimes, t, ans)
        {
            int threadId = omp_get_thread_num();
            // std::default_random_engine e(time(0)*(threadId+1));
            // std::default_random_engine e(seeds[threadId]);
            std::uniform_real_distribution<double> uiDistribution(0, 1);
            v_size i;

            while((i = __sync_fetch_and_add(&st[threadId], chunkSize)) < ed[threadId]) {
                v_size j = std::min(i + chunkSize, ed[threadId]);
                for(v_size x = i; x < j; x++) {
                    v_size u = nodes[x];

                    // e_size expectedSampleTime 
                    //     = std::round(sampleTimes * experiments[i] / sumW);
                    e_size expectedSampleTime 
                        = std::ceil(sampleTimes * experiments[x] / sumW);

                    if(expectedSampleTime == 0) continue; 
                    sortByColor[threadId] = g->pEdge + g->pIdx2[u];
                    computeDP(u, threadId);

                    v_size tt = 0;
                    for(v_size j = 0; j < expectedSampleTime; j++) {
                        tt += sampleOneTime(x, u, uiDistribution, e[threadId], threadId);
                    }

                    ans += 1.0*tt/expectedSampleTime*experiments[x];
                    sampleTotalTimes += expectedSampleTime;
                    t += tt;
                }
            }

            for(v_size tId = 0; tId < threads; tId++) {
                while((i = __sync_fetch_and_add(&st[tId], chunkSize)) < ed[tId]) {
                    v_size j = std::min(i + chunkSize, ed[tId]);
                    for(v_size x = i; x < j; x++) {
                        v_size u = nodes[x];

                    // e_size expectedSampleTime 
                    //     = std::round(sampleTimes * experiments[i] / sumW);
                        e_size expectedSampleTime 
                            = std::ceil(sampleTimes * experiments[x] / sumW);

                        if(expectedSampleTime == 0) continue; 
                        sortByColor[threadId] = g->pEdge + g->pIdx2[u];
                        computeDP(u, threadId);

                        v_size tt = 0;
                        for(v_size j = 0; j < expectedSampleTime; j++) {
                            tt += sampleOneTime(x, u, uiDistribution, e[threadId], threadId);
                        }

                        ans += 1.0*tt/expectedSampleTime*experiments[x];
                        sampleTotalTimes += expectedSampleTime;
                        t += tt;
                    }
                }
            }
        }
        // printf("sampleTimes %u\n", sampleTimes);
        // printf("sample rate %f\n", 1.0 * t / sampleTimes);
        printf("| %.6f %llu %.1f", 1.0 * t / sampleTotalTimes, sampleTotalTimes, sumW);
        // return 1.0 * t / sampleTotalTimes * sumW;
        return ans;
    }
};

#endif