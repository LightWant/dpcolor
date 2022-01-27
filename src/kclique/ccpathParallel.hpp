#ifndef samplePlusExactParallelPARALLEL_HPP
#define samplePlusExactParallelPARALLEL_HPP

#include "../graph/graph.hpp"
#include "../tools/type.hpp"
#include "../tools/hopstotchHash.hpp"
#include "../tools/linearSet.hpp"
// #include "./sample/cc.h"
#include "./sample/sampleBasedOnColors.h"
// #include "./sample/shadow.h"
// #include "./sample/ccpath.h"
// #include "../tools/multiThreads.hpp"
#include "./sample/ccpathParallel.h"
#include "./sample/cccpathParallel.h"

#include <cassert>
#include <tuple>
#include <random>
#include <vector>
#include <algorithm>
#include <random>

using Pair = std::pair<v_size, v_size>;
using std::vector;
#define pP first
#define pR second

class samplePlusExactParallel {
private:
    Graph * g;
    v_size k;
    double ** cnt;
    c_size ** C = nullptr;
    hopstotchHash * hashTable;
    v_size len = 1000000;
    v_size edges;
    LinearSet ** S;

    v_size maxDeepth = 0;
    v_size threads = 1;

public:
    samplePlusExactParallel(Graph * g_, v_size threads_ = 1)
        {g = g_; threads = threads_;}
    ~samplePlusExactParallel() {
        delete [] hashTable;
        
        for(v_size i = 0; i < threads; i++) {
            delete [] cnt[i];
            delete S[i];
        }
        delete [] cnt;
        delete [] S;
        for(v_size i = 0; i <= g->degeneracy; i++) {
            delete [] C[i];
        }
        delete [] C;
    }

    void previousWork();
  
    bool cc(v_size u, v_size v) {
        // for(v_size i = g->pIdx2[u]; i < g->pIdx[u + 1]; i++) {
        //     if(g->pEdge[i] == v) return true;
        // }
        // return false;
        return hashTable[u].contain(v);
    }

    void run(v_size k_, v_size deb, double exactN, double a, e_size N=50000000);
    void shadowFinder(v_size l);
    void search(v_size hNum, v_size pNum, Pair section, LinearSet * S, v_size threadId);
    void searchForSpecificK(v_size h, v_size p, Pair section, LinearSet * S, v_size k, v_size d = 0);

};

void samplePlusExactParallel::run(v_size k_, v_size deb, double exCnt, 
double alpha, e_size N) {
    // printf("begin run, dengenerarcy %u, vCnt %u\n", g->degeneracy, g->vCnt);
    k = k_;
    printf("||%u| %.1f| %llu", k, alpha, N); fflush(stdout);

     //pre work
    g->colorG();
    previousWork();

    // vector<v_size> *nodesPerThread;
    vector<v_size> nodes(g->vCnt);
    // nodesPerThread = new vector<v_size>[threads];
    // for(v_size i = 0; i < threads; i++) {
    //     nodesPerThread[i].resize(g->vCnt/threads);
    // }

    // vector<v_size> *nid;
    // nid = new vector<v_size>[threads];
    // for(v_size i = 0; i < threads; i++) {
    //     nid[i].resize(128);
    // }
    // vector<v_size> nid(threads);

    v_size sz = 0;
    
    omp_set_num_threads(threads);
    double tS = omp_get_wtime();
    
    #pragma omp parallel 
    {
        int threadId = omp_get_thread_num();
        // int sz = 0;
        #pragma omp for schedule(dynamic, 64)
        for(v_size u = 0; u < g->vCnt; u++) {
            v_size st = g->pIdx2[u], ed = g->pIdx[u+1], l = 0;
            if(ed - st + 1 < k) continue;
            
            for(v_size i = st; i < ed; i++) {
                v_size v = g->pEdge[i];

                for(v_size j = i + 1; j < ed; j++) {
                    v_size w = g->pEdge[j];
                    if(cc(v, w)) l++;
                }
            }
            // eCnt[u] = l;
            if(l > alpha* (k-1) * (ed - st)) {
                // nodesPerThread[threadId][sz++] = u;
                v_size t = __sync_fetch_and_add(&sz, 1);
                nodes[t] = u;
            }
            else {
                Pair section = S[threadId]->sort(u, {0, g->vCnt});
                search(1, 0, section, S[threadId], threadId);
            }
        }

        // #pragma omp for schedule(dynamic, 64)
        // for(v_size u = 0; u < g->vCnt; u++) {
        //     v_size st = g->pIdx2[u], ed = g->pIdx[u+1];
        //     if(ed - st + 1 < k) continue;
        //     if(eCnt[u] > alpha* (k-1) * (ed - st)) continue;
            
        //     Pair section = S[threadId]->sort(u, {0, g->vCnt});
        //     search(1, 0, section, S[threadId], threadId);
        // }

        // nid[threadId] = sz;
    }

    
    double t1 = omp_get_wtime();
    // printf("timesExact %fs\n", (t1 - t)/CLOCKS_PER_SEC);
    printf("| %.2f", t1 - tS);
    
    // v_size totalNodes = 0;
    // for(v_size i = 0; i < threads; i++) {
    //     totalNodes += nid[i];
    // }
    // nodes.resize(totalNodes);
    // totalNodes = 0;
    // for(v_size i = 0; i < threads; i++) {
    //     nodes.insert(nodes.begin() + totalNodes, 
    //         nodesPerThread[i].begin(), nodesPerThread[i].end());
    //     totalNodes += nid[i];
    // }
    // sz = totalNodes;
    printf("| %u", sz);
    
    double exactCnt = 0;
    for(v_size tId = 0; tId < threads; tId++) {
        exactCnt += cnt[tId][k];
    }

    printf("| %.0f", exactCnt);

//     if(sz > 0) {//ccpathParallel
//         ccpathParallel * scpObj = new ccpathParallel();

//         scpObj->initForSingleNode(k-1, g, hashTable, threads);

//         scpObj->init(sz, nodes);

//         double t = omp_get_wtime();
//         double sampleCnt_ = scpObj->sample(nodes, N);
//         // double sampleCntMultiLayer = scpObj->sampleMultilayer(nodes, N);

//         double t1 = omp_get_wtime();
//         // printf("ccpath sample times %fs\n", (t1 - t)/CLOCKS_PER_SEC);
//         double totalCnt = exactCnt + sampleCnt_;
//         printf("| %.0f", totalCnt);
//         printf("| %.2f%%", (exactCnt / totalCnt)*100);
//         // double totalCntML = exactCnt + sampleCntMultiLayer;
//         // printf("| %.2f%%", (exactCnt / totalCntML)*100);
//         printf("| %.2f", t1 - t);
//         printf("| %.2f%%\n", (abs(totalCnt - exCnt) / exCnt)*100);
//         printf("tm:%u %.2f\n", k, t1 - tS);
//         // printf("| %.2f%%", (abs(totalCntML - exCnt) / exCnt)*100);
// // scpObj->print();
//         delete scpObj;
//     }

    if(sz > 0) {//ccpathParallel
        double t = omp_get_wtime();

        cccpathParallel * scpObj = new cccpathParallel();
        scpObj->initForSingleNode(k-1, g, hashTable, threads);
        scpObj->init(sz, nodes);
        double sampleCnt_ = scpObj->sample(nodes, N);

        double t1 = omp_get_wtime();
        // printf("ccpath sample times %fs\n", (t1 - t)/CLOCKS_PER_SEC);
        printf(" %.6f%%", ((exCnt - exactCnt) / scpObj->sumW)*100);
        double totalCnt = exactCnt + sampleCnt_;
        printf("| %.0f", totalCnt);
        printf("| %.2f%%", (exactCnt / totalCnt)*100);
        // double totalCntML = exactCnt + sampleCntMultiLayer;
        // printf("| %.2f%%", (exactCnt / totalCntML)*100);
        printf("| %.2f", t1 - t);
        printf("| %.2f%%\n", (abs(totalCnt - exCnt) / exCnt)*100);
        printf("tm:%u %.6f\n", k, t1 - tS);
        // printf("| %.2f%%", (abs(totalCntML - exCnt) / exCnt)*100);
// scpObj->print();
        delete scpObj;
    }

    fflush(stdout);
};

void samplePlusExactParallel::previousWork() {
    hashTable = new hopstotchHash[g->vCnt];
    for(v_size u = 0; u < g->vCnt; u++) {
        if(g->pIdx[u + 1] == g->pIdx[u]) continue;
        hashTable[u].build(g->pEdge + g->pIdx[u], 
            g->pIdx[u + 1] - g->pIdx[u]);
    }

    S = new LinearSet*[threads];
    for(v_size i = 0; i < threads; i++) {
        S[i] = new LinearSet(g, hashTable);
    }

    cnt = new double*[threads];
    // double * memoryPool = new double[threads*(k+1)];
    // cnt[0] = memoryPool;
    for(v_size i = 0; i < threads; i++) {
        // cnt[i] = cnt[i - 1] + k + 1;
        cnt[i] = new double[k + 1]();
    }

    C = new c_size*[g->degeneracy + 1];
    for(v_size i = 0; i <= g->degeneracy; i++) {
        C[i] = new c_size[k + 1]();
    }
    C[0][0] = 1;
    C[1][0] = 1;
    C[1][1] = 1;
    for(v_size i = 2; i <= g->degeneracy; i++) {
        C[i][0] = 1;
        if(i < k + 1) C[i][i] = 1;
        for(v_size j = 1; j < i && j < k + 1; j++) {
            C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
        }
    }

    // eCnt = new v_size[g->vCnt]();
}


void samplePlusExactParallel::search(v_size hNum, v_size pNum, Pair section, LinearSet * S, v_size threadId) {
    if(hNum > k) return;

    if(section.pP == section.pR) {
        for(v_size j = 0; j <= pNum && hNum + j <= k; j++) {
            cnt[threadId][hNum + j] += C[pNum][j];
        }
// calls[threadId]++;
        return;
    }
// den[section.pR - section.pP]++;
    if(section.pP + 1 == section.pR) {
        for(v_size j = 0; j <= pNum + 1 && hNum + j <= k; j++) {
            cnt[threadId][hNum + j] += C[pNum + 1][j];
        }
// calls[threadId]++;
        return;
    }
    if(section.pP + 2 == section.pR) {
        if(hashTable[(*S)[0]].contain((*S)[1])) {
            for(v_size j = 0; j <= pNum + 2 && hNum + j <= k; j++) {
                cnt[threadId][hNum + j] += C[pNum + 2][j];
            }
        }
        else {
            for(v_size j = 0; j <= pNum + 1 && hNum + j <= k; j++) {
                cnt[threadId][hNum + j] += C[pNum + 1][j];
            }
            for(v_size j = 0; j <= pNum && hNum + 1 + j <= k; j++) {
                cnt[threadId][hNum + 1 + j] += C[pNum][j];
            }
        }
// calls[threadId]++;
        return;
    }
    v_size pivot, pivotDeg;
    v_size * tmpP = new v_size[section.pR - section.pP];
    S->findPivotAndCopyToTmpMem(section, tmpP, pivot, pivotDeg);
    section.pR--;

    search(hNum, pNum + 1, {section.pP, section.pP + pivotDeg}, S, threadId);
    
    v_size ed = section.pR;
    for(v_size i = pivotDeg; i < ed; i++) {
        v_size v = tmpP[i];
        
        section.pR--;
        S->changeTo(v, section.pR);
        Pair sec = S->sort2(v, section);
        search(hNum + 1, pNum, sec, S, threadId);
    }

    delete [] tmpP;
};


#endif
