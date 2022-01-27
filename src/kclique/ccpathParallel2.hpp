#ifndef samplePlusExactParallelPARALLEL2_HPP
#define samplePlusExactParallelPARALLEL2_HPP

#include "../graph/graph.hpp"
#include "../tools/type.hpp"
#include "../tools/hopstotchHash.hpp"
#include "../tools/bitmap.hpp"
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
    
    v_size maxV = 0;
    Graph ** subG;
    ListLinearHeap ** lheap = nullptr;
    v_size ** rank = nullptr;
    v_size ** ids = nullptr;
    v_size ** keys = nullptr;

public:
    samplePlusExactParallel(Graph * g_, v_size threads_ = 1)
        {g = g_; threads = threads_;}
    ~samplePlusExactParallel() {
        delete [] hashTable;
        
        for(v_size i = 0; i < threads; i++) {
            delete [] cnt[i];
            delete S[i];
            lheap[i]->~ListLinearHeap();
            delete [] rank[i];
            delete [] ids[i];
            delete [] keys[i];
        }
        delete [] cnt;
        delete [] S;
        delete [] lheap;
        delete [] rank;
        delete [] ids;
        delete [] keys;

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
    void computeDegeneracy(v_size threadId, v_size & mV);
};

void samplePlusExactParallel::run(v_size k_, v_size deb, double exCnt, 
double alpha, e_size N) {
    // printf("begin run, dengenerarcy %u, vCnt %u\n", g->degeneracy, g->vCnt);
    k = k_;
    printf("||%u| %.1f| %lu| %u", k, alpha, N, threads); fflush(stdout);

     //pre work
    g->colorG();
    previousWork();
    double tS = omp_get_wtime();

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
    //cpu_set_t cpu_mask;                    //Allocate mask
      //  CPU_ZERO(     &cpu_mask);           //Set mask to zero
        //for(int i = 0; i < threads; i++) {
          //  CPU_SET(i, &cpu_mask);           //Set mask with thread #
       // }
        //int err = sched_setaffinity( (pid_t)0,  //Set the affinity
          //                          sizeof(cpu_mask),
            //                        &cpu_mask );
    
    #pragma omp parallel
    {
        int threadId = omp_get_thread_num();
        // int sz = 0;
        v_size mV = 0;
        #pragma omp for schedule(static, 16) nowait
        for(v_size u = 0; u < g->vCnt; u++) {
            v_size st = g->pIdx2[u], ed = g->pIdx[u+1], l = 0;
            if(ed - st + 1 < k) continue;
            
            // for(v_size i = st; i < ed; i++) {
            //     v_size v = g->pEdge[i];

            //     for(v_size j = i + 1; j < ed; j++) {
            //         v_size w = g->pEdge[j];
            //         if(cc(v, w)) l++;
            //     }
            // }

            v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];
            subG[threadId]->pIdx[0] = 0;
            subG[threadId]->pIdx2[0] = 0;
            for(v_size i = st; i < g->pIdx[u+1]; i++) {
                v_size v = g->pEdge[i];
                subG[threadId]->pIdx[i -st + 1] = subG[threadId]->pIdx[i - st];

                for(v_size j = st; j < g->pIdx[u+1]; j++) {
                    v_size w = g->pEdge[j];
                    if(i == j) {
                        subG[threadId]->pIdx2[i - st] = l;
                    }
                    else if(cc(v, w)) {
                        subG[threadId]->pEdge[l++] = j - st;
                        subG[threadId]->pIdx[i - st + 1]++;
                    }
                }
            }
            subG[threadId]->vCnt = outDegree;
            subG[threadId]->eCnt = l;

            if(l/2 > alpha* (k-1) * (ed - st)) {
                // nodesPerThread[threadId][sz++] = u;
                v_size t = __sync_fetch_and_add(&sz, 1);
                nodes[t] = u;
            }
            else {
                // Pair section = S[threadId]->sort(u, {0, g->vCnt});
                // search(1, 0, section, S[threadId], threadId);
                computeDegeneracy(threadId, mV);

                v_size l = 0;
                
                for(v_size i = 0; i < subG[threadId]->vCnt; i++) {
                    v_size deg = subG[threadId]->pIdx[i + 1] - subG[threadId]->pIdx[i];
                    if(deg + 2 < k) continue;
                    v_size rki = rank[threadId][i];
                    
                    l = 0;
                    for(v_size j = subG[threadId]->pIdx[i]; j < subG[threadId]->pIdx[i + 1]; j++) {
                        v_size p = subG[threadId]->pEdge[j];
                        if(rki > rank[threadId][p]) continue;
                        S[threadId]->changeTo(g->pEdge[g->pIdx2[u] + p], l++);
                    }
                    
                    if(l + 2 < k) continue;
    // printf("%u %u %u\n", maxK, g->degeneracy, l);
                    // cnt[1][k-1] = 0;
                    search(2, 0, {0, l}, S[threadId], threadId);
                }
            }
        }
        
        maxV = std::max(maxV, mV);
    }

    double t1 = omp_get_wtime();
    // printf("timesExact %fs\n", (t1 - t)/CLOCKS_PER_SEC);
    printf("| %.2f", t1 - tS);fflush(stdout);
    
    printf("| %u", sz);
    
    double exactCnt = 0;
    for(v_size tId = 0; tId < threads; tId++) {
        exactCnt += cnt[tId][k];
    }

    printf("| %.1f", exactCnt);
    
    printf("| d%u", maxV);
    printf("\n");

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
        printf("| %.4f%%\n", (abs(totalCnt - exCnt) / exCnt)*100);
        printf("tm:%u %.6f\n", k, t1 - tS);
        // printf("| %.2f%%", (abs(totalCntML - exCnt) / exCnt)*100);
// scpObj->print();
        delete scpObj;
    }

    fflush(stdout);
};

void samplePlusExactParallel::computeDegeneracy(v_size threadId, v_size &mV) {
    for(v_size i = 0; i < subG[threadId]->vCnt; i++) {
        ids[threadId][i] = i;
        keys[threadId][i] = subG[threadId]->pIdx[i + 1] - subG[threadId]->pIdx[i] + 1;
    }
    
    lheap[threadId]->init(subG[threadId]->vCnt, subG[threadId]->vCnt, 
        ids[threadId], keys[threadId]);

    for(v_size i = 0; i < subG[threadId]->vCnt; i++) {
        v_size v, degV;

        if(!lheap[threadId]->pop_min(v, degV)) printf("error\n");
// printf("%u %u\n", v, degV-1);
        rank[threadId][v] = i;
        mV = std::max(mV, degV);

        for(v_size j = subG[threadId]->pIdx[v]; j < subG[threadId]->pIdx[v + 1]; j++) {
            lheap[threadId]->decrement(subG[threadId]->pEdge[j]);
        }
    }
}

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

    subG = new Graph*[threads]();
    lheap = new ListLinearHeap*[threads];
    rank = new v_size*[threads];
    ids = new v_size*[g->degeneracy];
    keys = new v_size*[g->degeneracy + 1];

    for(int i = 0; i < threads; i++) {
        subG[i] = new Graph();
        subG[i]->pIdx = new v_size[g->degeneracy];
        subG[i]->pIdx2 = new v_size[g->degeneracy];
        subG[i]->pEdge = new v_size[g->degeneracy*g->degeneracy];

        lheap[i] = new ListLinearHeap(g->degeneracy, g->degeneracy*g->degeneracy);
        rank[i] = new v_size[g->degeneracy];
        ids[i] = new v_size[g->degeneracy];
        keys[i] = new v_size[g->degeneracy + 1];
    }   
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
