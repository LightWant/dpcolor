#ifndef CCPARALLEL_H
#define CCPARALLEL_H

#include "../../graph/graph.hpp"
#include "../../tools/type.hpp"
#include "../../tools/hopstotchHash.hpp"
#include "../../tools/linearSet.hpp"


#include <cassert>
#include <tuple>
#include <random>
#include <vector>
#include <algorithm>
#include <random>

using Pair = std::pair<v_size, v_size>;
using std::vector;

struct shadowPlusColorParallel {
    v_size sz;
    Graph * g;
    hopstotchHash * hashTable;
    v_size k;
    double * experiments;
    double sumW, maxW;
    v_size maxD;
    double *** dp;
    
    v_size ** pEdge = nullptr;
    v_size ** pIdx = nullptr;

    v_size ** sortByColor = nullptr;
    v_size ** pColor = nullptr;

    v_size ** clique = nullptr;
    v_size ** a = nullptr;
    vector<std::uniform_int_distribution<int>> uForColor;

    v_size N = 5000000;
    int threads;

    void init(v_size sz_, std::vector<v_size> & nodes, v_size N_=5000000) {
        sz = sz_;
        N = N_;
        experiments = new double[sz];
        for(int i = 0; i < g->cc; i++) {
            uForColor.push_back(std::uniform_int_distribution<int>(0, i));
        }

        #pragma omp parallel reduction(+:sumW)
        {
            int threadId = omp_get_thread_num();
            
            #pragma omp for schedule(dynamic, 32) 
            for(v_size i = 0; i < sz; i++) {
                v_size u = nodes[i];
                
                // sortVbyColor(i, u, threadId);
                v_size deg = g->pIdx[u+1] - g->pIdx2[u];
                for(v_size c = 0; c < g->cc; c++) a[i][c] = 0;
                for(v_size j =  g->pIdx2[u]; j < g->pIdx[u + 1]; j++) {
                    v_size c = g->color[g->pEdge[j]];
                    a[i][c]++;
                }
                computeDP(i, threadId);
                // set(i, g->pIdx[u+1] - g->pIdx2[u], dp[threadId][g->cc][k]);
                experiments[i] = dp[threadId][g->cc][k];
                sumW += dp[threadId][g->cc][k];
            }
        }
printf("sumW %f\n", sumW);
    }

    void initForSingleNode(v_size k_, Graph * g_, 
        hopstotchHash * hashTable_, v_size threads_, v_size sz) {
        k = k_;
        g = g_;
        hashTable = hashTable_;
        threads = threads_;

        sumW = 0.0;
        maxD = 0;
        
        clique = new v_size*[threads + 1];
        for(v_size i = 0; i < threads; i++) {
            clique[i] = new v_size[k*1000];
        }

        dp = new double**[threads + 1];
        for(v_size t = 0; t < threads; t++) {
            dp[t] = new double*[g->cc + 1];

            for(v_size i = 0; i < g->cc + 1; i++) {
                dp[t][i] = new double[k+1];
            }
        }

        pEdge = new v_size*[threads];
        pIdx = new v_size*[threads];
        for(v_size t = 0; t < threads; t++) {
            pEdge[t] = new v_size[g->degeneracy*g->degeneracy];
            pIdx[t] = new v_size[g->cc + 5];
        }

        a = new v_size*[sz];
        for(int i = 0; i < sz; i++) {
            a[i] = new v_size[g->cc + 1]();
            // for(v_size j = 0; j < g->cc; j++) {
            //     a[i][j] = 0;
            // }
        }
        // sortByColor = new v_size*[threads];
        // for(v_size i = 0; i < threads; i++) {
        //     sortByColor[i] = new v_size[g->degeneracy + 1];
        // }
        // pColor = new v_size*[threads];
        // for(v_size i = 0; i < threads; i++) {
        //     pColor[i] = new v_size[g->cc + 1];
        // }
    }

    ~shadowPlusColorParallel() {
        if(experiments != nullptr) delete [] experiments;
        // if(memoryPool != nullptr) delete [] memoryPool;
        for(v_size i = 0; i < threads; i++) {
            for(v_size j = 0; j < g->cc + 1; j++)
                delete [] dp[i][j];
            delete [] dp[i];
        }
        if(dp != nullptr) delete [] dp;
        for(v_size i = 0; i < sz; i++) {
            delete [] a[i];
        }
        // delete [] a[0];
        if(a != nullptr) delete [] a;

        for(v_size i = 0; i < threads; i++) {
            delete [] pEdge[i];
            delete [] pIdx[i];
            delete [] clique[i];
        }
        if(pEdge != nullptr) delete [] pEdge;
        if(pIdx != nullptr) delete [] pIdx;
        if(clique != nullptr) delete [] clique;
    }

    bool connect(v_size u, v_size v) {
        return hashTable[u].contain(v);
    }

    void computeDP(v_size ii, v_size threadId) {
        dp[threadId][0][0] = 1.0;
        for(v_size i = 1; i <= g->cc; i++) {
            dp[threadId][i][0] = 1.0;
            for(v_size j = 1; j <= k && j <= i; j++) {
                dp[threadId][i][j] = 
                    dp[threadId][i - 1][j - 1]*a[ii][i - 1]
                    + dp[threadId][i - 1][j];
            }
        }
    }

    void sortVbyColor(v_size ii, v_size u, v_size threadId) {
        pIdx[threadId][0] = 0;
        for(v_size i = 1; i <= g->cc; i++) {
            pIdx[threadId][i] = pIdx[threadId][i - 1] + a[ii][i - 1];
        }
// int tmp = 0;
        for(v_size i = g->pIdx2[u]; i < g->pIdx[u + 1]; i++) {
            v_size v = g->pEdge[i];
            v_size c = g->color[v];
// if(c == 0) tmp++;
            pEdge[threadId][pIdx[threadId][c]] = v;
            pIdx[threadId][c]++;
        }
// if(tmp != a[ii][0]) {
//     printf("xxx %d %u\n", tmp, a[ii][0]);fflush(stdout);
// }
// if(pIdx[threadId][0] != a[ii][0])
//     printf("%u %u\n", pIdx[threadId][0],  a[ii][0]);
//         for(v_size c = 1; c < g->cc; c++) {
//             if(pIdx[threadId][c] != pIdx[threadId][c-1]+a[ii][c]) {
// printf("%u %u %u\n", pIdx[threadId][c], pIdx[threadId][c-1], a[ii][c]);
// fflush(stdout);
            // }
        // }

        pIdx[threadId][0] = 0;
        for(v_size i = 1; i <= g->cc; i++) {
            pIdx[threadId][i] = pIdx[threadId][i - 1] + a[ii][i - 1];
        }
    }

    int sampleOneTime(std::uniform_real_distribution<double> & udistri, 
        std::default_random_engine & e,
        v_size ii, v_size threadId, v_size degree) {
        v_size i = g->cc;
        v_size j = k;
        v_size l = 0;

        // while(l < k) {//unuseful loop; can eliminate
        while(i > 0 && j > 0) {
            if(udistri(e)+1e-8 < dp[threadId][i-1][j-1]*a[ii][i-1]/dp[threadId][i][j]) {
                clique[threadId][l++] = i - 1;
                j--;
            }
            i--;
        }
// printf("lk %u %u\n", l, k);fflush(stdout);
        assert(l == k);

        // auto rd = [&](int max) {   
        //     return uForColor[max](e);
        // };

        for(v_size i = 0; i < l; i++) {
            // clique[threadId][i] = randChoose(ii, clique[threadId][i], threadId);
            v_size c = clique[threadId][i];
            v_size deg = a[ii][c];
            v_size offset = deg * udistri(e)+1e-5;
            offset %= deg;
// if(pIdx[threadId][c] + offset >= degree) {
//     printf("error\n");fflush(stdout);
// }
            clique[threadId][i] = pEdge[threadId][pIdx[threadId][c] + offset];
        }

        for(v_size i = 0; i < l; i++) {
            for(v_size j = i + 1; j < l; j++) {
                if(!connect(clique[threadId][i], clique[threadId][j])) return 0;
            }
        }
        
        return 1;
    }



    double sample(std::vector<v_size> & nodes, e_size sampleTimes) {
        v_size t = 0;
        e_size sampleTotalTimes = 0;
        // double ans = 0.0;
        std::random_device rd;
        std::default_random_engine e[100];
        for(v_size i = 0; i < threads; i++) {
            e[i].seed(rd());
        }

        e_size * chunk = new e_size[sz]();
        v_size * pChunk = new v_size[threads + 1];
        pChunk[0] = 0;
        e_size sumCompute = 0;

        #pragma omp parallel for schedule(dynamic, 32) reduction(+:sumCompute, sampleTotalTimes)
        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];

            e_size expectedTime = std::round(sampleTimes * (experiments[i] / sumW));
            if(expectedTime == 0) {
                chunk[i] = 0;
                continue;
            }
// tmpp+=expectedTime;
            v_size outDegree = g->pIdx[u+1] - g->pIdx2[u];
            chunk[i] = expectedTime * k * k;
            chunk[i] += g->cc+outDegree+g->cc*k;
            chunk[i] += outDegree*(1+log(outDegree));
            
            sumCompute += chunk[i];
            sampleTotalTimes += expectedTime;
        }

        e_size chunkSz = (sumCompute+threads-1) / threads;
        v_size p = 1, tmp = 0;
        for(v_size i = 0; i < sz; i++) {
            tmp += chunk[i];
            if(tmp >= p * chunkSz) {
                pChunk[p] = i;
                p++;
                if(p == threads) break;
            }
        }
        while(p <= threads) {
            pChunk[p] = sz;
            p++;
        }

        v_size * pt = new v_size[threads];
        for(int i = 0; i < threads; i++) {
            pt[i] = pChunk[i];
        }
// double sumW2 = 0;
        #pragma omp parallel proc_bind(spread) reduction(+:t)
        {
            int threadId = omp_get_thread_num();
            std::uniform_real_distribution<double> uiDistribution(0, 1);
            v_size i;
            // v_size tt = 0;
            // #pragma omp for schedule(dynamic, 1) 
            // for(v_size i = 0; i < sz; i++) {
// double st = omp_get_wtime();
            while((i = __sync_fetch_and_add(&pt[threadId], 1)) < pChunk[threadId + 1]) {
            // for(v_size i = pChunk[threadId]; i < pChunk[threadId + 1]; i++) {
                v_size u = nodes[i];

                e_size expectedSampleTime 
                    = std::round(sampleTimes * experiments[i] / sumW);

                if(expectedSampleTime == 0) {
                    // i++;
                    continue; 
                }

                computeDP(i, threadId);
                sortVbyColor(i, u, threadId);
// sumW2 += dp[threadId][g->cc][k];
                // v_size t2 = 0;
                v_size tt = 0;
                for(v_size j = 0; j < expectedSampleTime; j++) {
                    tt += sampleOneTime(uiDistribution, e[threadId], i, threadId, g->pIdx[u+1]-g->pIdx2[u]);
                }
                t += tt;
                // __sync_fetch_and_add(&sampleTotalTimes, expectedSampleTime);
// tt += expectedSampleTime;
                // i++;
            }

            //stealing
            for(int j = 0; j < threads; j++) {
                if(j != threadId && pt[j] < pChunk[j + 1]) {
                    int i;
                    while((i = __sync_fetch_and_add(&pt[j], 1)) < pChunk[j + 1]) {
                        v_size u = nodes[i];
                        e_size expectedSampleTime 
                            = std::round(sampleTimes * experiments[i] / sumW);
                        if(expectedSampleTime == 0) {
                            continue; 
                        }
                        // sortVbyColor(i, u, threadId);
                        computeDP(i, threadId);
                        sortVbyColor(i, u, threadId);
                        v_size tt = 0;
                        for(v_size l = 0; l < expectedSampleTime; l++) {
                            tt += sampleOneTime(uiDistribution, e[threadId], i, threadId, g->pIdx[u+1]-g->pIdx2[u]);
                        }
                        t += tt;
                    }
                }
            }
// printf("%u %f\n", threadId, omp_get_wtime() - st);
            #pragma omp barrier

            if(sampleTotalTimes < sampleTimes) {
                e_size leftTimes = sampleTimes - sampleTotalTimes;
// printf("%llu\n", leftTimes);
                std::discrete_distribution<int> 
                    udistribution(experiments, experiments + sz);

                #pragma omp for
                for(e_size i = 0; i < leftTimes; i++) {
                    int id = udistribution(e[threadId]);
                    v_size u = nodes[id];
                    sortVbyColor(id, u, threadId);
                    computeDP(id, threadId);
                    
                    t += sampleOneTime(uiDistribution, e[threadId], id, threadId, g->pIdx[u+1]-g->pIdx2[u]);
                }
            }
        }

        delete [] chunk;
        delete [] pChunk;
        delete [] pt;

        // if(sampleTotalTimes < sampleTimes) sampleTotalTimes = sampleTimes;
        
        // printf("sampleTimes %u\n", sampleTimes);
        printf("| %f %u sum w%.0f\n", 1.0 * t / sampleTimes, t, sumW);
        // printf("diff %f\n", 1.0 * t / sampleTotalTimes * sumW-ans);
        return 1.0 * t / sampleTotalTimes * sumW;
    }
};

#endif

//     double sample(std::vector<v_size> & nodes, e_size sampleTimes) {
//         v_size t = 0;
//         e_size sampleTotalTimes = 0;
//         double ans = 0.0;
//         std::random_device rd;
//         std::default_random_engine e[100];
//         for(v_size i = 0; i < threads; i++) {
//             e[i].seed(rd());
//         }

//         #pragma omp parallel reduction(+:sampleTotalTimes, t, ans)
//         {
//             int threadId = omp_get_thread_num();
//             std::uniform_real_distribution<double> uiDistribution(0, 1);
//             // v_size tt = 0;
// // double st = omp_get_wtime();
//             #pragma omp for schedule(static, 32)
//             for(v_size i = 0; i < sz; i++) {
//             // while((i = __sync_fetch_and_add(&pt[threadId], 1)) < pChunk[threadId + 1]) {
//                 v_size u = nodes[i];

//                 e_size expectedSampleTime
//                     = std::round(sampleTimes * experiments[i] / sumW);

//                 if(expectedSampleTime == 0) {
//                     continue; 
//                 }

//                 computeDP(i, threadId);
//                 sortVbyColor(i, u, threadId);
//                 // v_size t2 = 0;
//                 v_size tt = 0;

//                 for(e_size j = 0; j < expectedSampleTime; j++) {
//                     tt += sampleOneTime(uiDistribution, e[threadId], i, threadId);
//                 }

//                 t += tt;
//                 sampleTotalTimes += expectedSampleTime;
//             }
// // printf("%f\n", omp_get_wtime() - st);
// //             #pragma omp barrier

// //             if(sampleTotalTimes < sampleTimes) {
// //                 v_size leftTimes = sampleTimes - sampleTotalTimes;
// // // printf("%llu\n", leftTimes);
// //                 std::discrete_distribution<int> 
// //                     udistribution(experiments, experiments + sz);

// //                 #pragma omp for
// //                 for(v_size i = 0; i < leftTimes; i++) {
// //                     int id = udistribution(e[threadId]);
// //                     v_size u = nodes[id];
// //                     computeDP(id, threadId);
// //                     sortVbyColor(id, u, threadId);
// //                     t += sampleOneTime(uiDistribution, e[threadId], id, threadId);
// //                 }
// //             }
//         }

//         // printf("sampleTimes %u\n", sampleTimes);
//         printf("| %f %llu", 1.0 * t / sampleTotalTimes, sampleTotalTimes);
//         // printf("diff %f\n", 1.0 * t / sampleTotalTimes * sumW-ans);
//         return 1.0 * t / sampleTotalTimes * sumW;
//     }
