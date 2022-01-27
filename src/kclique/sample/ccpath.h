#ifndef CCPATH_H
#define CCPAHT_H

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

// constexpr v_size batchSize = 50;

struct ccpath {
    v_size sz;
    Graph * g;
    hopstotchHash * hashTable;
    v_size k;
    double * experiments;
    double sumW, maxW;
    v_size maxD;
    double ** dp;
    double * memoryPool = nullptr;

    v_size * pEdge = nullptr;
    v_size * pIdx = nullptr;
    v_size vCnt, eCnt;

    v_size * clique = nullptr;
    std::default_random_engine e;
    e_size N = 5000000;

    void init(v_size sz_, std::vector<v_size> & nodes, e_size N_=5000000) {
        sz = sz_;
        N = N_;
        experiments = new double[sz];

        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];

            computeDP(u);

            double sumD = 0.0;
            for(v_size i = 0; i < g->pIdx[u+1] - g->pIdx2[u]; i++) {
                sumD += dp[i][k];
            }
// printf("%.0f\\", sumD);
            set(i, g->pIdx[u+1] - g->pIdx2[u], sumD);
        }
    }

    void initForSingleNode(v_size k_, Graph * g_, hopstotchHash * hashTable_) {
        k = k_;
        g = g_;
        hashTable = hashTable_;
        sumW = 0.0;
        maxD = 0;
        clique = new v_size[k];

        dp = new double*[g->degeneracy];
        memoryPool = new double[g->degeneracy * (k+1)]();
        v_size p = 0;
        for(v_size i = 0; i < g->degeneracy; i++) {
            dp[i] = memoryPool + p;
            p += k + 1;
        }

        for(v_size i = 0; i < g->degeneracy; i++) {
            dp[i][0] = 0;
            dp[i][1] = 1;
        }

        pEdge = new v_size[g->degeneracy*g->degeneracy];
        pIdx = new v_size[g->degeneracy];
    }

    ~ccpath() {
        if(experiments != nullptr) delete [] experiments;
        if(memoryPool != nullptr) delete [] memoryPool;
        if(dp != nullptr) delete [] dp;
        if(pEdge != nullptr) delete [] pEdge;
        if(pIdx != nullptr) delete [] pIdx;
        if(clique != nullptr) delete [] clique;
    }

    bool connect(v_size u, v_size v) {
        return hashTable[u].contain(v);
    }

    void set(v_size i, v_size vCnti, double experimentsi) {
        experiments[i] = experimentsi;
        sumW += experimentsi;
        maxD = std::max(maxD, vCnti);
        maxW = std::max(maxW, experimentsi);
    }
    
    void computeDP(v_size u) {
        v_size outDegree = g->pIdx[u+1] - g->pIdx2[u];

        pIdx[0] = 0;
        for(v_size i = 0; i < outDegree; i++) {
            v_size v = g->pEdge[g->pIdx2[u] + i];
            pIdx[i + 1] = pIdx[i];

            for(v_size j = i + 1; j < outDegree; j++) {
                v_size w = g->pEdge[g->pIdx2[u] + j];
                // if(g->color[w] < g->color[v]) {
                    if(hashTable[v].contain(w)) {
                        pEdge[pIdx[i + 1]++] = j;
                    }
                // }
            }     
        }

        for(v_size j = 2; j <= k; j++) {
            for(v_size i = 0; i < outDegree; i++) {
                dp[i][j] = 0.0;
                for(v_size l = pIdx[i]; l < pIdx[i + 1]; l++) {
                    dp[i][j] += dp[pEdge[l]][j - 1];
                }
            }
        }
    }

    v_size shrink(std::vector<v_size> & nodes) {
        v_size sz2 = sz;
        v_size l = 0;
        double newSumW = 0;

        while(l < sz2) {
            if(experiments[l] * 100 < maxW) {
                std::swap(nodes[l], nodes[--sz2]);
                std::swap(experiments[l], experiments[sz2]);
            }
            else {
                newSumW += experiments[l];
                l++;
            }
        }
        sumW = newSumW;

        printf("%f\n", maxW);
        printf("%u-%u-%u\n", l, sz2, sz);

        sz = l;

        return l + 1;
    }

    void print() {
        // int tmp[] = {0,0,0,0,0,0,0};
        printf("%u\n", sz);
        for(v_size i = 0; i < sz; i++) {
            printf("%.2f\n", experiments[i]);
        }
    }

    int sampleOneTime(v_size id, v_size u, 
        std::uniform_real_distribution<double> & d) {
        v_size preId = -1;

        double sumD = experiments[id];
        double x = d(e);
        
        double sumTmp = 0.0;
        for(v_size i = 0; i < g->pIdx[u+1] - g->pIdx2[u]; i++) {
            sumTmp += dp[i][k];
            if(sumTmp + 1e-10 >= x * sumD) {
                clique[0] = g->pEdge[g->pIdx2[u] + i];
                preId = i;
                break;
            }
        }

        for(v_size i = 1; i < k; i++) {
            sumTmp = sumD = 0.0;
            for(v_size j = pIdx[preId]; j < pIdx[preId + 1]; j++) {
                sumD += dp[pEdge[j]][k - i];
            }

            x = d(e);
            for(v_size j = pIdx[preId]; j < pIdx[preId + 1]; j++) {
                sumTmp += dp[pEdge[j]][k - i];
                if(sumTmp + 1e-10 >= x * sumD) {
                    clique[i] = g->pEdge[g->pIdx2[u] + pEdge[j]];
                    preId = pEdge[j];
                    break;
                }
            }

            for(v_size j = 0; j < i; j++) {
                if(!connect(clique[i], clique[j])) {
                    return 0;
                }
            }
        }
        
        return 1;
    }

    double sample(std::vector<v_size> & nodes, e_size sampleTimes) {
        e_size t = 0;
        e_size sampleTotalTimes = 0;
        std::default_random_engine generator;
        std::uniform_real_distribution<double> uiDistribution(0, 1);
        // double ans = 0.0;

        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];

            computeDP(u);

            e_size expectedSampleTime 
                = std::round(sampleTimes * experiments[i] / sumW);

            sampleTotalTimes += expectedSampleTime;

            // if(expectedSampleTime == 0) continue; 
            // t = 0;
            for(v_size j = 0; j < expectedSampleTime; j++) {
                t += sampleOneTime(i, u, uiDistribution);
            }
            // ans += 1.0*t/expectedSampleTime*sumW;
        }

        // printf("sampleTimes %u\n", sampleTimes);
        // printf("sample rate %f\n", 1.0 * t / sampleTimes);
     //           printf("| %.6f", 1.0 * t / sampleTimes);

	printf("| %.6f", 1.0 * t / sampleTotalTimes);

	return 1.0 * t / sampleTotalTimes * sumW;
        // return ans;
    }

    double sampleMultilayer(std::vector<v_size> & nodes, v_size sampleTimes) {
        v_size sampleTotalTimes = 0;
        std::default_random_engine generator;
        std::uniform_real_distribution<double> uiDistribution(0, 1);
        double ans = 0.0;

        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];

            computeDP(u);

            v_size expectedSampleTime 
                = std::round(sampleTimes * experiments[i] / sumW);

            if(expectedSampleTime == 0) continue;

            sampleTotalTimes += expectedSampleTime;

            v_size t = 0;
            for(v_size j = 0; j < expectedSampleTime; j++) {
                t += sampleOneTime(i, u, uiDistribution);
            }
            ans += 1.0 * t / expectedSampleTime * experiments[i];
        }

        return ans;
    }

     double sampleAuto(std::vector<v_size> & nodes) {
        v_size sampleTimes = 0;
        constexpr v_size batchSize = 50;
        int noChange = 0;
        
        v_size t = 0;
        std::default_random_engine generator(time(0));
		std::discrete_distribution<int> 
            distribution(experiments, experiments + sz);
        std::uniform_real_distribution<double> uiDistribution(0, 1);
        double eps = 1e-5;

        // if(8 <= k+1 && k+1 <= 9) eps /= 10;
        // else if(9 < k+1 && k+1 <= 13) eps /= 100;
        // else if(13 < k+1) eps /= 1000;

        v_size maxBatches = 50000000;
        while(maxBatches--) {
            v_size tmpT = 0;

            for(v_size i = 0; i < batchSize; i++) {
                int id = distribution(generator);
                v_size u = nodes[id];

                computeDP(u);

                int ret = sampleOneTime(id, u, uiDistribution);

                tmpT += ret;
            }

            if(sampleTimes == 0) {
                t += tmpT; sampleTimes += batchSize;
// if(tmpT == 0) printf("sample error\n");
                continue;
            }

            double preRate = 1.0*t/sampleTimes;
            double nowRate = 1.0*(t+tmpT)/(sampleTimes+batchSize);
            
            t += tmpT; sampleTimes += batchSize;

            if(preRate > 0 && abs(preRate - nowRate) < eps) {
                noChange++;
                if(noChange > 15) break;
            }
            else noChange = 0;
        }

        printf("sampleTimes %u\n", sampleTimes);
        printf("sample rate %f\n", 1.0 * t / sampleTimes);
        return 1.0 * t / sampleTimes * sumW;
    }
};

#endif

  // double sample(std::vector<v_size> & nodes) {
    //     v_size sampleTimes = sz*20;
    //     v_size t = 0;
    //     v_size sampleTotalTimes = 0;
    //     std::default_random_engine generator;
	// 	std::discrete_distribution<int> 
    //         distribution(experiments, experiments + sz);
    //     std::uniform_real_distribution<double> uiDistribution(0, 1);

    //     while(sampleTotalTimes < sampleTimes) {
    //         int id = distribution(generator);
    //         v_size u = nodes[id];

    //         computeDP(u);
    //         sortVbyColor(u);

    //         int ret = sampleOneTime(uiDistribution);

    //         sampleTotalTimes++;
    //         t += ret;
    //     }

    //     return 1.0 * t / sampleTimes * sumW;
    // }

//    
