#ifndef CC_H
#define CC_H

#include <cassert>
#include <tuple>
#include <random>
#include <vector>
#include <algorithm>
#include <random>

using Pair = std::pair<v_size, v_size>;
using std::vector;

constexpr v_size batchSize = 50;

struct shadowPlusColor {
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

    v_size * clique = nullptr;
    v_size * a = nullptr;
    std::default_random_engine e;
    v_size N = 5000000;

    void init(v_size sz_, std::vector<v_size> & nodes, v_size N_=5000000) {
        sz = sz_;
        N = N_;
        experiments = new double[sz];

        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];
            computeDP(u);
            set(i, g->pIdx[u+1] - g->pIdx2[u], dp[g->cc][k]);
        }
    }

    void initForSingleNode(v_size k_, Graph * g_, hopstotchHash * hashTable_) {
        k = k_;
        g = g_;
        hashTable = hashTable_;
        sumW = 0.0;
        maxD = 0;
        clique = new v_size[k];

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

    ~shadowPlusColor() {
        if(experiments != nullptr) delete [] experiments;
        if(memoryPool != nullptr) delete [] memoryPool;
        if(dp != nullptr) delete [] dp;
        if(a != nullptr) delete [] a;
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
    
    v_size randChoose(v_size c) {
        v_size deg = a[c];
        return pEdge[pIdx[c] + std::rand() % deg];
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

        for(v_size i = 1; i < deg; i++) {
            if(g->color[pEdge[i]] != g->color[pEdge[i-1]]) {
                pIdx[g->color[pEdge[i]]] = i;
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

    double sampleOneBatchForRate(v_size u) {
        int tmpT = 0;
        static std::uniform_real_distribution<double> uiDistribution(0, 1);
        for(v_size i = 0; i < 5; i++) {
            computeDP(u);
            sortVbyColor(u);

            int ret = sampleOneTime(uiDistribution);

            tmpT += ret;
        }

        return 1.0 * tmpT / 5;
    }

    int sampleOneTime(std::uniform_real_distribution<double> & u) {
        v_size i = g->cc;
        v_size j = k;
        v_size l = 0;

        // while(l < k) {//unuseful loop; can eliminate
            while(i > 0 && j > 0) {
                if(u(e) + 1e-8 < dp[i-1][j-1]*a[i-1]/dp[i][j]) {
                    clique[l++] = randChoose(i-1);
                    j--;
                }
                i--;
            }

        //     if(i == 0 && j > 0) {
        //         i = g->cc; j = k; l = 0;
        //     }
        // }

        // if(l < k) return -1;

        for(v_size i = 0; i < l; i++) {
            for(v_size j = i + 1; j < l; j++) {
                if(!connect(clique[i], clique[j])) return 0;
            }
        }
        
        return 1;
    }

    double sample(std::vector<v_size> & nodes) {
        v_size sampleTimes = sz*20;
        v_size t = 0;
        v_size sampleTotalTimes = 0;
        std::default_random_engine generator;
		std::discrete_distribution<int> 
            distribution(experiments, experiments + sz);
        std::uniform_real_distribution<double> uiDistribution(0, 1);

        while(sampleTotalTimes < sampleTimes) {
            int id = distribution(generator);
            v_size u = nodes[id];

            computeDP(u);
            sortVbyColor(u);

            int ret = sampleOneTime(uiDistribution);

            sampleTotalTimes++;
            t += ret;
        }

        return 1.0 * t / sampleTimes * sumW;
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

        if(8 <= k+1 && k+1 <= 9) eps /= 10;
        else if(9 < k+1 && k+1 <= 13) eps /= 100;
        else if(13 < k+1) eps /= 1000;

        v_size maxBatches = 50000000;
        while(maxBatches--) {
            v_size tmpT = 0;

            for(v_size i = 0; i < batchSize; i++) {
                int id = distribution(generator);
                v_size u = nodes[id];

                computeDP(u);
                sortVbyColor(u);

                int ret = sampleOneTime(uiDistribution);

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

    double sample(std::vector<v_size> & nodes, e_size sampleTimes) {
        v_size t = 0;
        e_size sampleTotalTimes = 0;
        std::default_random_engine generator;
        std::uniform_real_distribution<double> uiDistribution(0, 1);
        // double ans = 0.0;

        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];

            computeDP(u);
            sortVbyColor(u);

            e_size expectedSampleTime 
                = std::round(sampleTimes * experiments[i] / sumW);

            sampleTotalTimes += expectedSampleTime;

            // v_size t2 = 0;
            for(v_size j = 0; j < expectedSampleTime; j++) {
                int ret = sampleOneTime(uiDistribution);
                t += ret;
                // t2 += ret;
            }

            // ans += 1.0 * t2 / expectedSampleTime * dp[g->cc][k];
        }

        // printf("sampleTimes %u\n", sampleTimes);
        printf("| %f sumW%.0f\n", 1.0 * t / sampleTimes, sumW);
        // printf("diff %f\n", 1.0 * t / sampleTotalTimes * sumW-ans);
        return 1.0 * t / sampleTotalTimes * sumW;
    }
};

#endif