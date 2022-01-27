#ifndef SHADOWROW_H
#define SHADOWROW_H


#include <cassert>
#include <tuple>
#include <random>
#include <vector>
#include <algorithm>
#include <random>

using Pair = std::pair<v_size, v_size>;
using std::vector;

struct shadows {
    v_size sz;
    Bitmap ** graphs;
    Bitmap * memoryPool;
    v_size * l;
    e_size * eCnt;
    v_size * vCnt;
    double * experiments;
    double sumW;
    v_size maxD;
    v_size * rng;

    void init(v_size sz_, v_size dengeneracy, 
        vector<v_size> & nodes, Graph * g) {
        sz = sz_;
        graphs = new Bitmap*[sz];
        e_size total = 0;
        for(v_size i = 0; i < sz; i++) {

            v_size u = nodes[i];
            v_size outDegree = g->pIdx[u+1]- g->pIdx2[u];
            total += outDegree;
        }
        memoryPool = new Bitmap[total];
        v_size p = 0;

        for(v_size i = 0; i < sz; i++) {
            v_size u = nodes[i];
            v_size outDegree = g->pIdx[u+1]- g->pIdx2[u];
// printf("%u de%u ", i, outDegree);fflush(stdout);

            graphs[i] = memoryPool + p;
            p += outDegree;

            for(v_size j = 0; j < outDegree; j++) {
                graphs[i][j].resize(outDegree);
            }
        }
        l = new v_size[sz];
        eCnt = new e_size[sz];
        vCnt = new v_size[sz];
        experiments = new double[sz];
        sumW = 0.0;
        maxD = 0;
    }
    ~shadows() {
        // for(v_size i = 0; i < sz; i++) {
        //     delete [] graphs[i];
        // }
        delete [] memoryPool;
        delete [] graphs;
        delete [] l;
        delete [] eCnt;
        delete [] vCnt;
        delete [] experiments;
    }
    void set(v_size i, v_size li, 
        e_size eCnti, v_size vCnti, double experimentsi) {
        // graphs[i] = subGraph;
        l[i] = li;
        eCnt[i] = eCnti;
        vCnt[i] = vCnti;
        experiments[i] = experimentsi;
        sumW += experimentsi;
        maxD = std::max(maxD, vCnti);
    }
    double f(v_size x) {
        return pow(2.71828, x)/sqrt(acos(-1)*2*x*x*x*x*x);
    }
    double sample() {
        std::default_random_engine generator;
		std::discrete_distribution<int> distribution(experiments, experiments + sz);

        rng = new v_size[maxD];

        v_size Xr = 0;
        // v_size sampleTimes = std::max(sz, 50000u);
        v_size sampleTimes = f(l[0])*maxD*maxD;
        // e_size sampleTimes = 0;
        // for(v_size i = 0; i < sz; i++) {
        //     sampleTimes += eCnt[i];
        // }
printf("sptimes %u\n", sampleTimes);
        for(v_size i = 0; i < sampleTimes; i++) {
            int shadowId = distribution(generator);
            for(v_size i = 0; i < vCnt[shadowId]; i++) {
                rng[i] = i;
            }
            std::random_shuffle(rng, rng + vCnt[shadowId]);
            bool f = true;

            for(v_size j = 1; j < l[shadowId]; j++) {
                for(v_size m = 0; m < j; m++) {
                    if(!graphs[shadowId][rng[j]][rng[m]]) {
                        f = false; break;
                    }
                }
                if(!f) break;
            }

            if(f) Xr++;
        }

        delete [] rng;

        return Xr * sumW / sampleTimes;
    }
};

struct shadows2 {
    v_size sz;
    Graph * g;
    hopstotchHash * hashTable;
    v_size k;
    double * experiments;
    double sumW;
    v_size maxD;

    void init(v_size sz_, v_size k_, Graph * g_, hopstotchHash * hashTable_) {
        sz = sz_;
        k = k_;
        g = g_;
        experiments = new double[sz];
        sumW = 0.0;
        maxD = 0;
        hashTable = hashTable_;
    }
    ~shadows2() {
        delete [] experiments;
    }

    bool connect(v_size u, v_size v) {
        return hashTable[u].contain(v);
    }
    void set(v_size i, v_size vCnti, double experimentsi) {
        experiments[i] = experimentsi;
        sumW += experimentsi;
        maxD = std::max(maxD, vCnti);
    }
    double f(v_size x) {
        return pow(2.71828, x)/sqrt(acos(-1)*2*x*x*x*x*x);
    }
    double sample(std::vector<v_size> & nodes) {
        std::default_random_engine generator;
		std::discrete_distribution<int> distribution(experiments, experiments + sz);

        v_size *rng = new v_size[maxD];

        v_size Xr = 0;
        // v_size sampleTimes = std::max(sz, 50000u);
        v_size sampleTimes = f(k)*maxD*maxD;
        if(sampleTimes > 1000000) sampleTimes = 1000000;
        // e_size sampleTimes = 0;
        // for(v_size i = 0; i < sz; i++) {
        //     sampleTimes += eCnt[i];
        // }
printf("sptimes %u maxD %u\n", sampleTimes, maxD);
        for(v_size i = 0; i < sampleTimes; i++) {
            int id = distribution(generator);
            v_size u = nodes[id];

            v_size deg = g->pIdx[u+1]-g->pIdx2[u];
            for(v_size i = 0; i < deg; i++) {
                rng[i] = i;
            }
            std::random_shuffle(rng, rng + deg);
            bool f = true;

            for(v_size j = 1; j < k; j++) {

                v_size x = g->pEdge[g->pIdx2[u]+rng[j]];
                for(v_size m = 0; m < j; m++) {
                    v_size p1 = g->pIdx2[u]+rng[m];
                    v_size y = g->pEdge[p1];

                    if(!connect(x, y)) {
                        f = false; break;
                    }
                }
                if(!f) break;
            }

            if(f) Xr++;
        }

        delete [] rng;

        return Xr * sumW / sampleTimes;
    }
};

#endif