#ifndef SAMPLEPLUSEXACT_HPP
#define SAMPLEPLUSEXACT_HPP

#include "../graph/graph.hpp"
#include "../tools/type.hpp"
#include "../tools/hopstotchHash.hpp"
#include "../tools/linearSet.hpp"
#include "./sample/cc.h"
#include "./sample/sampleBasedOnColors.h"
//#include "./sample/shadow.h"
#include "./sample/ccpath.h"
#include "./sample/cccpath.h"

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

class samplePlusExact {
private:
    Graph * g, *subG;
    v_size k;
    double ** cnt;
    c_size ** C = nullptr;
    hopstotchHash * hashTable;
    v_size len = 1000000;
    v_size edges;
    LinearSet * S;
    v_size ** tmp;
    v_size maxDeepth = 0;

    ListLinearHeap *lheap=nullptr;
    v_size * rank = nullptr;
    v_size * ids = nullptr;
    v_size * keys = nullptr;

    //multiRun:
    v_size maxK = 0;
    v_size * eCnt = nullptr;
    v_size * vCnt = nullptr;
    double ** ansTmp = nullptr;
    double * ansTmpMemoryPool = nullptr;

public:
    samplePlusExact(Graph * g_) {g = g_;}
    ~samplePlusExact() {
        for(v_size i = 0; i <= g->degeneracy; i++) {
            delete [] C[i];
            // delete [] tmp[i];
        }
        delete S;
        delete [] tmp[0];
        delete [] tmp;
        delete [] C;
        delete [] cnt[0];

        delete [] cnt;
        delete [] hashTable;
        // delete subG;
        if(eCnt != nullptr) delete [] eCnt;
        if(vCnt != nullptr) delete [] vCnt;
        if(ansTmp != nullptr) delete [] ansTmp;
        if(ansTmpMemoryPool != nullptr) delete [] ansTmpMemoryPool;
        // if(lheap != nullptr) delete lheap;
        lheap->~ListLinearHeap();
 
        if(rank != nullptr) delete [] rank;	
        if(ids != nullptr) delete [] ids;	
        if(keys != nullptr) delete [] keys;	
    }
    void previousWork();
    void buildSubGraphFalse(v_size u);
    void buildSubGraphTrue(v_size u);
    bool cc(v_size u, v_size v) {
        return hashTable[u].contain(v);
    }

    void run(v_size k_, v_size deb, double exactN, double a, e_size N=50000000);
    void shadowFinder(v_size l);
    void search(Pair section, LinearSet * S, v_size d = 0, v_size h=1);
    void searchForSpecificK(v_size h, v_size p, Pair section, LinearSet * S, v_size k, v_size d = 0);

    void multiRunInit(v_size maxK, double maxA);
    void multiRun(v_size k_, double a, double exCnt, e_size N);
    void multiRun2(v_size k_, double a, double exCnt);
    // double sampleColor();

    void runCC(v_size k_, v_size deb, double exactN, double a, e_size N=50000000);

    void runCCPath(v_size k_, v_size deb, double exactN, double a, e_size N=50000000);

    void computeDegeneracy();
};

void samplePlusExact::multiRunInit(v_size maxK_, double maxA) {
    double t = clock();
    k = maxK = maxK_;
    // g->colorG();
    
    previousWork();
    g->colorG();

    ansTmp = new double*[g->vCnt];
    ansTmpMemoryPool = new double[g->vCnt * (k+1)];
    std::fill(ansTmpMemoryPool, 
        ansTmpMemoryPool + g->vCnt * (k+1), 0.0);
    double * p = ansTmpMemoryPool;
    for(v_size i = 0; i < g->vCnt; i++) {
        ansTmp[i] = p;
        p += k + 1;
    }
    // std::fill(ansTmp, ansTmp + g->vCnt, -1.0);
    eCnt = new v_size[g->vCnt];
    vCnt = new v_size[g->vCnt];

    std::fill(cnt[0], cnt[0] + k + 1, 0.0);
    cnt[0][1] = g->vCnt;

    for(v_size u = 0; u < g->vCnt; u++) {
        // v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];

        // if(outDegree + 1 < 6) continue;
        buildSubGraphFalse(u);
        
        eCnt[u] = subG->eCnt;
        vCnt[u] = subG->vCnt;

        // if(eCnt[u] / 2 > 13 * vCnt[u]) continue;

        Pair section = S->sort(u, {0, g->vCnt});
        cnt[1][k-1] = 0.0;

        search(section, S, 1);

        for(v_size i = 2; i <= k && i <= section.pR + 1; i++) {
            ansTmp[u][i] += cnt[1][i-1];
            if(cnt[1][i-1] == 0.0) break;
        }
    }

    printf("previous Work : %.2fs\n", (clock() - t) / CLOCKS_PER_SEC);
}

void samplePlusExact::multiRun2(v_size k_, double alpha, double exCnt) {
    k = k_;
    printf("%u ", k_);
    double t = clock();
    v_size st = 0, ed = g->vCnt;

    std::fill(cnt[0], cnt[0] + k+1, 0.0);
    cnt[0][1] = g->vCnt;

    vector<v_size> nodes;
    v_size sz = 0;
    nodes.resize(g->vCnt/4);

    v_size a[1000];
    for(v_size u = st; u < ed; u++) {
        v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];

        if(outDegree + 1 < k) continue;
        
        // if(vCnt[u] == g->vCnt) {
        //     buildSubGraphFalse(u);
        //     eCnt[u] = subG->eCnt;
        //     vCnt[u] = subG->vCnt;
        // }
        // else {
        //     subG->eCnt = eCnt[u];
        //     subG->vCnt = vCnt[u];
        // }
        for(v_size i = 0; i < g->cc; i++) {
            a[i] = 0;
        }
        v_size st = g->pIdx2[u], l = 0;
        for(v_size i = st; i < g->pIdx[u+1]; i++) {
            v_size v = g->pEdge[i];

            for(v_size j = i + 1; j < g->pIdx[u+1]; j++) {
                v_size w = g->pEdge[j];
                if(cc(w, v)) {
                    l++;
                    a[g->color[w]]++;
                    a[g->color[v]]++;
                }
            }
        }
        v_size ret = 0;
        for(v_size i = 0; i < g->cc; i++) {
            if(a[i] > k - 1) ret++;
        }

        // if(subG->eCnt / 2 > alpha*k* subG->vCnt) {
        if(ret > k-1) {
            if(sz == nodes.size()) {
                nodes.resize(sz + ed - u);
            }
            nodes[sz] = u;
            sz++;
        }
        else {

            Pair section = S->sort(u, {0, g->vCnt});
            search(section, S, 1);

            // ansTmp[u] = cnt[1][k-1];
            for(v_size i = 2; i <= k && i <= section.pR + 1; i++) {
                cnt[0][i] += cnt[1][i-1];
                if(cnt[1][i-1] == 0.0) break;
            }
        }
    }

    // printf("timesExact %fs\n", (clock()-t)/CLOCKS_PER_SEC);
    // printf(" %u ", sz);
    
    double exactCnt = cnt[0][k];
   
    // t = clock();
    
    // printf("color graph time %.2f\n", (clock()-t)/CLOCKS_PER_SEC);
    if(sz > 0) {//c+t
        // t1 = clock();
        shadowPlusColor * spcObj = new shadowPlusColor();
        spcObj->initForSingleNode(k-1, g, hashTable);
        spcObj->init(sz, nodes);

        double sampleCnt_ = spcObj->sampleAuto(nodes);

        double sampleT_ = clock();
        printf("%.3f\n", (sampleT_-t)/CLOCKS_PER_SEC);
        
    double totalCnt = exactCnt + sampleCnt_;
        printf("ex/total %.2f%% ", (exactCnt / exCnt)*100);
        // double exCnt = 59636674809701648;
    printf(" %.2f\n", ((totalCnt - exCnt) / exCnt)*100);
        // double exCnt = 59636674809701648;
        printf("%f\n", (totalCnt - exCnt) / exCnt);

        delete spcObj;
    }
    else printf("n n\n");
    

    fflush(stdout);
}

void samplePlusExact::multiRun(v_size k_, double alpha, double exCnt, e_size N) {
    k = k_;
    printf("k:%u alpha:%.1f N:%lu", k_, alpha, N);fflush(stdout);
    
    double t = clock(), t1;
    
    double exactCnt = 0.0;

    vector<v_size> nodes;
    v_size sz = 0;
    nodes.resize(g->vCnt/4);

    for(v_size u = 0; u < g->vCnt; u++) {
        v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];

        if(outDegree + 1 < k) continue;
        // buildSubGraphFalse(u);

        // if(subG->eCnt / 2 > alpha*k* subG->vCnt) {
        if(eCnt[u] / 2 > alpha * (k-1) * vCnt[u]) {
            if(sz == nodes.size()) {
                nodes.resize(g->vCnt - u);
            }
            nodes[sz] = u;
            sz++;
        }
        else {
            exactCnt += ansTmp[u][k];
        }
    }

    t1 = clock();
    // printf(" exT:%.2fs", (t1-t)/CLOCKS_PER_SEC);
    t = t1;
    printf(" sz:%u", sz);
   
    // t = clock();
    
    // printf("color graph time %.2f\n", (clock()-t)/CLOCKS_PER_SEC);
    if(sz > 0) {//c+t
        // ccpath * scpObj = new ccpath();
        // scpObj->initForSingleNode(k-1, g, hashTable);
        // scpObj->init(sz, nodes);

        // double sampleCnt_ = scpObj->sample(nodes, N);

        // double t1 = clock();
        // printf(" spT:%.2fs", (t1 - t)/CLOCKS_PER_SEC);
        
        // double totalCnt = exactCnt + sampleCnt_;
        // printf(" ex/total:%.2f%%", (exactCnt / totalCnt)*100);
        // printf(" ex:%.0f", exactCnt);
        // printf(" totalT:%.2fs", (t1 - tS)/CLOCKS_PER_SEC);
        // printf(" error:%.2f%%\n", ((totalCnt - exCnt) / exCnt)*100);

        // delete scpObj;

        shadowPlusColor * spcObj = new shadowPlusColor();
        spcObj->initForSingleNode(k-1, g, hashTable);
        spcObj->init(sz, nodes);

        t = clock();
        double sampleCnt_ = spcObj->sample(nodes, N);
        double t1 = clock();
        double totalCnt = exactCnt + sampleCnt_;

        if(totalCnt > 0)
        printf("| %.2f%%", (exactCnt / totalCnt)*100);
        printf("| %.2f", (t1 - t)/CLOCKS_PER_SEC);

        printf("| %.2f%%|\n", (abs(totalCnt - exCnt) / exCnt)*100);

        delete spcObj;
    }
    else printf("\n");
    

    fflush(stdout);
}

void samplePlusExact::run(v_size k_, v_size deb, double exCnt, double alpha, e_size N) {
    // printf("begin run, dengenerarcy %u, vCnt %u\n", g->degeneracy, g->vCnt);
    k = k_;
printf("|%u| %.1f| %lu", k, alpha, N);
    double t = clock(), t1, tS;
    
     //pre work
    g->colorG();
    previousWork();

    t1 = clock();
    double tS2 = omp_get_wtime();
    tS = clock();
    // printf("previousWork %.2fs\n", (t1 - t) / CLOCKS_PER_SEC);
    t = t1;

    v_size st = 0, ed = g->vCnt;

    std::fill(cnt[0], cnt[0] + k+1, 0.0);
    cnt[0][1] = g->vCnt;

    vector<v_size> nodes;
    v_size sz = 0;
    nodes.resize(g->vCnt/4);

    for(v_size u = st; u < ed; u++) {
        v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];

        if(outDegree + 1 < k) continue;
        // buildSubGraphFalse(u);
        buildSubGraphTrue(u);

        if(subG->eCnt / 2 > alpha* (k-1) * subG->vCnt) {
            if(sz == nodes.size()) {
                nodes.resize(sz + ed - u);
            }
            nodes[sz] = u;
            sz++;
        }
        else {
            // for(v_size i = g->pIdx2[u]; i < g->pIdx[u+1]; i++) {
            //     v_size v = g->pEdge[i];
            //     v_size l = 0;
            //     for(v_size j = i+1; j < g->pIdx[u+1]; j++) {
            //         v_size w = g->pEdge[j];
            //         if(cc(v, w)) {
            //             S->changeTo(w, l++);
            //         }
            //     }
            //     cnt[1][k-2] = 0;
            //     search({0, l}, S, 1, 2);
            //     cnt[0][k] += cnt[1][k-2];
            // }
// printf("%.0f\n", cnt[0][k]);
            computeDegeneracy();

            v_size l = 0;

            for(v_size i = 0; i < subG->vCnt; i++) {
                v_size deg = subG->pIdx[i + 1] - subG->pIdx[i];
                if(deg + 2 < k) continue;
                v_size rki = rank[i];
                
                l = 0;
                for(v_size j = subG->pIdx[i]; j < subG->pIdx[i + 1]; j++) {
                    v_size p = subG->pEdge[j];
                    if(rki > rank[p]) continue;
                    S->changeTo(g->pEdge[g->pIdx2[u] + p], l++);
                }
                
                if(l + 2 < k) continue;
// printf("%u %u %u\n", maxK, g->degeneracy, l);
                // cnt[1][k-1] = 0;
                // search({0, l}, S, 1, 2);
                // cnt[0][k] += cnt[1][k-2];
                searchForSpecificK(2, 0, {0, l}, S, k, 0);
            }
            

            // Pair section = S->sort(u, {0, g->vCnt});
            // cnt[1][k-1] = 0;
            // search(section, S, 1, 1);
            // cnt[0][k] += cnt[1][k-1];

            // Pair section = S->sort(u, {0, g->vCnt});
            // searchForSpecificK(1, 0, section, S, k, 0);
        }
    }

    t1 = clock();
    // printf("timesExact %fs\n", (t1 - t)/CLOCKS_PER_SEC);
    double tS3 = omp_get_wtime();
    printf("| %.2f", tS3 - tS2);
    printf("| %.2f", (t1 - t)/CLOCKS_PER_SEC);fflush(stdout);
    t = t1;

    printf("| %u", sz);
    printf("| d%u %u", maxDeepth, maxK);
    
    double exactCnt = cnt[0][k];

    printf("| %.0f", exactCnt);

    printf("\n");
   
    if(sz > 0) {//cccpath
        cccpath * scpObj = new cccpath();
        scpObj->initForSingleNode(k-1, g, hashTable);
        scpObj->init(sz, nodes);

        double sampleCnt_ = scpObj->sample(nodes, N, exCnt - exactCnt);
        // double sampleCntMultiLayer = scpObj->sampleMultilayer(nodes, N);
        
        printf("|s %f", (exCnt - exactCnt) / scpObj->sumW);
        double t1 = clock();
        // printf("ccpath sample times %fs\n", (t1 - t)/CLOCKS_PER_SEC);
        double totalCnt = exactCnt + sampleCnt_;
        printf("| %.0f", totalCnt);
        printf("| %.2f%%", (exactCnt / totalCnt)*100);
        // double totalCntML = exactCnt + sampleCntMultiLayer;
        // printf("| %.2f%%", (exactCnt / totalCntML)*100);
        printf("| %.2f", (t1 - t)/CLOCKS_PER_SEC);
        
        
        // printf("%d-clique: %.2f\n", k, totalCnt);
        // printf("ex/total %.2f%%\n", (exactCnt / totalCnt)*100);
        
        
        // printf("error %.2f%%\n", ((totalCnt - exCnt) / exCnt)*100);
        // printf("total run time %.2fs\n", (t1 - tS)/CLOCKS_PER_SEC);
        
        printf("| %.2f", (t1 - tS)/CLOCKS_PER_SEC);
        printf("| %.9f%%\n", (abs(totalCnt - exCnt) / exCnt)*100);
        // printf("| %.2f%%", (abs(totalCntML - exCnt) / exCnt)*100);
// scpObj->print();
        delete scpObj;

// for(v_size i = 0; i < sz; i++) {
//     v_size u = nodes[i];
//     v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];
// printf("%u:", outDegree);
//     v_size st = g->pIdx2[u], l = 0;
    
// double avg = 0, D = 0;
//     for(v_size i = st; i < g->pIdx[u+1]; i++) {
//         v_size v = g->pEdge[i];
//         v_size d = 0;
//         for(v_size j = st; j < g->pIdx[u+1]; j++) {
//             v_size w = g->pEdge[j];
//             if(cc(w, v)) {
//                 d++;
//             }
//         }
//         avg += d;
//     }
//     avg /= outDegree;
//     for(v_size i = st; i < g->pIdx[u+1]; i++) {
//         v_size v = g->pEdge[i];
//         v_size d = 0;
//         for(v_size j = st; j < g->pIdx[u+1]; j++) {
//             v_size w = g->pEdge[j];
//             if(cc(w, v)) {
//                 d++;
//             }
//         }
// D += (avg-d)*(avg-d);
//     }
// printf("%.2f", sqrt(D/outDegree));
// printf("\n");
// }
    }

    fflush(stdout);
};

void samplePlusExact::computeDegeneracy() {
    for(v_size i = 0; i < subG->vCnt; i++) {
        ids[i] = i;
        keys[i] = subG->pIdx[i + 1] - subG->pIdx[i] + 1;
    }
    lheap->init(subG->vCnt, subG->vCnt, ids, keys);

    //maxK = 0;
    for(v_size i = 0; i < subG->vCnt; i++) {
        v_size v, degV;

        if(!lheap->pop_min(v, degV)) printf("error\n");
// printf("%u %u\n", v, degV-1);
        rank[v] = i;
        maxK = std::max(maxK, degV);
        for(v_size j = subG->pIdx[v]; j < subG->pIdx[v + 1]; j++) {
            lheap->decrement(subG->pEdge[j]);
        }
    }
}

void samplePlusExact::runCC(v_size k_, v_size deb, double exCnt, double alpha, e_size N) {
    printf("begin run, dengenerarcy %u, vCnt %u\n", g->degeneracy, g->vCnt);
    k = k_;
printf("|%u| %.1f| %lu", k, alpha, N);
    previousWork();
    double t = clock(), t1, tS;
    tS = t;
     //pre work
    
    g->colorG();

    t1 = clock();
    t = t1;

    v_size st = 0, ed = g->vCnt;

    shadowPlusColor * spcObj = new shadowPlusColor();
    spcObj->initForSingleNode(k-1, g, hashTable);

    std::fill(cnt[0], cnt[0] + k+1, 0.0);
    cnt[0][1] = g->vCnt;

    vector<v_size> nodes;
    v_size sz = 0;
    nodes.resize(g->vCnt/4);

    for(v_size u = st; u < ed; u++) {
        v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];

        if(outDegree + 1 < k) continue;
        buildSubGraphFalse(u);
        
        if(subG->eCnt / 2 > alpha* (k-1) * subG->vCnt) {
            if(sz == nodes.size()) {
                nodes.resize(sz + ed - u);
            }
            nodes[sz] = u;
            sz++;
        }
        else {
            Pair section = S->sort(u, {0, g->vCnt});
            cnt[1][k-1] = 0;
            search(section, S, 1);

            // for(v_size i = 2; i <= k && i <= section.pR + 1; i++) {
            //     cnt[0][i] += cnt[1][i-1];
            //     if(cnt[1][i-1] == 0.0) break;
            // }
            cnt[0][k] += cnt[1][k-1];
        }
    }

    t1 = clock();
    // printf("timesExact %fs\n", (t1 - t)/CLOCKS_PER_SEC);
    printf("| %.2f", (t1 - t)/CLOCKS_PER_SEC);
    t = t1;

    printf("| %u", sz);
    
    double exactCnt = cnt[0][k];

    printf("| %.0f", exactCnt);
   
    if(sz > 0) {//ccpath
        spcObj->init(sz, nodes);

        double sampleCnt_ = spcObj->sample(nodes, N);
        // double sampleCnt_ = scpObj->sampleAuto(nodes);
        printf("| %f", (exCnt - exactCnt) / spcObj->sumW);

        double t1 = clock();
        // printf("ccpath sample times %fs\n", (t1 - t)/CLOCKS_PER_SEC);
        double totalCnt = exactCnt + sampleCnt_;
        printf("| %.2f%% %.2f%%", (exactCnt / totalCnt)*100, (exactCnt / exCnt)*100);
        printf("| %.2f", (t1 - t)/CLOCKS_PER_SEC);
        
        // printf("%d-clique: %.2f\n", k, totalCnt);
        // printf("ex/total %.2f%%\n", (exactCnt / totalCnt)*100);
        
        // printf("error %.2f%%\n", ((totalCnt - exCnt) / exCnt)*100);
        // printf("total run time %.2fs\n", (t1 - tS)/CLOCKS_PER_SEC);
        
        printf("| %.2f", (t1 - tS)/CLOCKS_PER_SEC);
        printf("| %.6f%%|\n", (abs(totalCnt - exCnt) / exCnt)*100);
// scpObj->print();
        delete spcObj;
    }

    fflush(stdout);
};

void samplePlusExact::runCCPath(v_size k_, v_size deb, double exCnt, double alpha, e_size N) {
    // printf("begin run, dengenerarcy %u, vCnt %u\n", g->degeneracy, g->vCnt);
    k = k_;
    printf("|%u| %.1f| %lu", k, alpha, N);
    previousWork();
    g->colorG();
    
    double t = clock(), t1, tS;
    tS = t;

    v_size st = 0, ed = g->vCnt;
    
    
    cccpath * scpObj = new cccpath();

    std::fill(cnt[0], cnt[0] + k+1, 0.0);
    cnt[0][1] = g->vCnt;

    vector<v_size> nodes;
    v_size sz = 0;
    nodes.resize(g->vCnt/4);

    for(v_size u = st; u < ed; u++) {
        v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];

        if(outDegree + 1 < k) continue;
        buildSubGraphFalse(u);

        if(subG->eCnt / 2 > alpha* (k-1) * subG->vCnt) {
            if(sz == nodes.size()) {
                nodes.resize(sz + ed - u);
            }
            nodes[sz] = u;
            sz++;
        }
        else {
            Pair section = S->sort(u, {0, g->vCnt});
            cnt[1][k-1] = 0;
            search(section, S, 1, 1);
            cnt[0][k] += cnt[1][k-1];
// printf("%u:%.0f\n", u, cnt[0][k]);
            // Pair section = S->sort(u, {0, g->vCnt});
            // searchForSpecificK(1, 0, section, S, k, 0);
        }
    }

    t1 = clock();
    // printf("timesExact %fs\n", (t1 - t)/CLOCKS_PER_SEC);
    printf("| %.2f", (t1 - t)/CLOCKS_PER_SEC);fflush(stdout);
    t = t1;

    printf("| %u", sz);
    printf("| d%u", maxDeepth);
    
    double exactCnt = cnt[0][k];

    printf("| %.1f", exactCnt);
   
    if(sz > 0) {//ccpath
        t = clock();
        scpObj->initForSingleNode(k-1, g, hashTable);
        scpObj->init(sz, nodes);

        double sampleCnt_ = scpObj->sample(nodes, N, exCnt - exactCnt);

        double t1 = clock();
        // printf("ccpath sample times %fs\n", (t1 - t)/CLOCKS_PER_SEC);
        printf(" %.8f%%", ((exCnt - exactCnt) / exCnt)*100);
        double totalCnt = exactCnt + sampleCnt_;
        
        printf("| %.0f", totalCnt);
        printf("| %.2f%%", (exactCnt / totalCnt)*100);
        // double totalCntML = exactCnt + sampleCntMultiLayer;
        // printf("| %.2f%%", (exactCnt / totalCntML)*100);
        printf("| %.2f", (t1 - t)/CLOCKS_PER_SEC);
        printf("| %.2f", (t1 - tS)/CLOCKS_PER_SEC);
        printf("| %.2f%%\n", (abs(totalCnt - exCnt) / exCnt)*100);
        // printf("| %.2f%%", (abs(totalCntML - exCnt) / exCnt)*100);
// scpObj->print();
        delete scpObj;
    }

    fflush(stdout);
};

void samplePlusExact::previousWork() {
    hashTable = new hopstotchHash[g->vCnt];
    for(v_size u = 0; u < g->vCnt; u++) {
        if(g->pIdx[u + 1] == g->pIdx[u]) continue;
        hashTable[u].build(g->pEdge + g->pIdx[u], 
            g->pIdx[u + 1] - g->pIdx[u]);
    }

    S = new LinearSet(g, hashTable);

    cnt = new double*[g->degeneracy + 1];
    double * memoryPool2 = new double[(g->degeneracy + 1)*(k+1)];
    cnt[0] = memoryPool2;
    for(v_size i = 1; i <= g->degeneracy; i++) {
        cnt[i] = cnt[i-1] + (k + 1);
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

    tmp = new v_size*[g->degeneracy];
    v_size * memoryPool = new v_size[g->degeneracy*g->degeneracy];
    v_size p = 0;
    for(v_size i = 0; i < g->degeneracy; i++) {
        // tmp[i] = new v_size[g->degeneracy]();
        tmp[i] = memoryPool + p;
        p += g->degeneracy;
    }

    subG = new Graph();
    subG->pIdx = new v_size[g->degeneracy];
    subG->pIdx2 = new v_size[g->degeneracy];
    subG->pEdge = new v_size[g->degeneracy*g->degeneracy];

    lheap = new ListLinearHeap(g->degeneracy, g->degeneracy*g->degeneracy);
    rank = new v_size[g->degeneracy];
    ids = new v_size[g->degeneracy];
    keys = new v_size[g->degeneracy + 1];
}

void samplePlusExact::buildSubGraphFalse(v_size u) {
    v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];

    v_size st = g->pIdx2[u], l = 0;
    for(v_size i = st; i < g->pIdx[u+1]; i++) {
        v_size v = g->pEdge[i];

        for(v_size j = st; j < g->pIdx[u+1]; j++) {
            v_size w = g->pEdge[j];
            if(cc(v, w)) l++;
        }
    }

    subG->vCnt = outDegree;
    subG->eCnt = l;
}

void samplePlusExact::buildSubGraphTrue(v_size u) {
    v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];

    v_size st = g->pIdx2[u], l = 0;
    subG->pIdx[0] = 0;
    subG->pIdx2[0] = 0;

    for(v_size i = st; i < g->pIdx[u+1]; i++) {
        v_size v = g->pEdge[i];
        subG->pIdx[i -st + 1] = subG->pIdx[i - st];

        for(v_size j = st; j < g->pIdx[u+1]; j++) {
            v_size w = g->pEdge[j];
            if(i == j) {
                subG->pIdx2[i - st] = l;
            }
            else if(cc(v, w)) {
                subG->pEdge[l++] = j - st;
                subG->pIdx[i - st + 1]++;
            }
        }
    }

    subG->vCnt = outDegree;
    subG->eCnt = l;
}


void samplePlusExact::search(Pair section, LinearSet * S, v_size d, v_size h) {
    std::fill(cnt[d], cnt[d] + std::min(k + 1, section.pR + 2), 0.0);
    cnt[d][1] = section.pR - section.pP;
    maxDeepth = std::max(maxDeepth, d);

    if(h > k) return;

    if(section.pP == section.pR) return;
    if(section.pP + 1 == section.pR) return;
    if(section.pP + 2 == section.pR) {
        if(hashTable[(*S)[section.pP]].contain((*S)[section.pP+1])) cnt[d][2] = 1.0;
        return;
    }
    if(section.pP + 3 == section.pR) {
        int f1 = hashTable[(*S)[section.pP]].contain((*S)[section.pP+1]);
        int f2 = hashTable[(*S)[section.pP]].contain((*S)[section.pP+2]);
        int f3 = hashTable[(*S)[section.pP+1]].contain((*S)[section.pP+2]);
        cnt[d][2] = double(f1 + f2 + f3);
        if(f1 && f2 && f3) cnt[d][3] = 1.0;
        return;
    }

    v_size pivot, pivotDeg;
    int numMax = S->findPivotsAndClique(section, tmp[d], pivot, pivotDeg);

    section.pR--;

    if((v_size)numMax == pivotDeg && section.pR == pivotDeg) {
        for(v_size i = 2; i <= pivotDeg + 1 && i <= k; i++) {
            cnt[d][i] = C[pivotDeg + 1][i];
        }
        return;
    }

    search({section.pP, section.pP + pivotDeg}, S, d+1, h);
    
    for(v_size i = 2; i <= k && i <= pivotDeg + 1; i++) {
        cnt[d][i] = cnt[d+1][i-1] + cnt[d+1][i];
    }

    v_size ed = section.pR;
    for(v_size i = section.pP + pivotDeg; i < ed; i++) {
        v_size v = tmp[d][i];
        
        section.pR--;
        S->changeTo(v, section.pR);

        Pair sec = S->sort2(v, section);
        if(sec.pR == pivotDeg && ed == section.pP + pivotDeg + 1) {
            for(v_size i = 2; i <= k && i <= pivotDeg + 1; i++) {
                cnt[d][i] += cnt[d+1][i-1];
            }
            continue;
        }

        search(sec, S, d+1, h + 1);
        for(v_size i = 2; i <= k && i <= sec.pR - sec.pP + 1; i++) {
            cnt[d][i] += cnt[d+1][i-1];
        }
    }
};

void samplePlusExact::searchForSpecificK(v_size hNum, 
    v_size pNum, Pair section, LinearSet * S, v_size k, v_size d) {
    if(hNum > k) return;

    maxDeepth = std::max(maxDeepth, d);

    if(section.pP == section.pR) {
        for(v_size j = 0; j <= pNum && hNum + j <= k; j++) {
            cnt[0][hNum + j] += C[pNum][j];
        }
        return;
    }
    if(section.pP + 1 == section.pR) {
        for(v_size j = 0; j <= pNum + 1 && hNum + j <= k; j++) {
            cnt[0][hNum + j] += C[pNum + 1][j];
        }
        return;
    }
    if(section.pP + 2 == section.pR) {
        if(hashTable[(*S)[section.pP]].contain((*S)[section.pP+1])) {
            for(v_size j = 0; j <= pNum + 2 && hNum + j <= k; j++) {
                cnt[0][hNum + j] += C[pNum + 2][j];
            }
        }
        else {
            for(v_size j = 0; j <= pNum + 1 && hNum + j <= k; j++) {
                cnt[0][hNum + j] += C[pNum + 1][j];
            }
            for(v_size j = 0; j <= pNum && hNum + 1 + j <= k; j++) {
                cnt[0][hNum + 1 + j] += C[pNum][j];
            }
        }
        return;
    }

    v_size pivot, pivotDeg;
    // v_size * tmpP = new v_size[section.pR - section.pP];
    S->findPivotAndCopyToTmpMem(section, tmp[d], pivot, pivotDeg);
    section.pR--;
// printf("pivot, left %u,right %u, num %d\n", pivotDeg, section.pR - pivotDeg, numMax);

    if(pivotDeg == 0) {
        for(v_size j = 0; j <= pNum && hNum + 1 + j <= k; j++) {
            cnt[0][hNum + 1 + j] += C[pNum][j]*section.pR;
        }
        for(v_size j = 0; j <= pNum+1 && hNum + j <= k; j++) {
            cnt[0][hNum + j] += C[pNum+1][j];
        }
        return;
    }

    // if(pivotDeg == (v_size)numMax && pivotDeg == section.pR) {
    //     for(v_size j = 0; j <= pNum+pivotDeg+1 && hNum + j <= k; j++) {
    //         cnt[hNum + j] += C[pNum+pivotDeg+1][j];
    //     }
    //     return;
    // }

    searchForSpecificK(hNum, pNum + 1, {section.pP, section.pP + pivotDeg}, S, k, d+1);

    // if(pivotDeg == 1) {
    //     for(v_size j = 0; j <= pNum && hNum + 1 + j <= k; j++) {
    //         cnt[hNum + 1 + j] += C[pNum][j]*(section.pR - pivotDeg - numMax);
    //     }
    //     if(numMax > 0)
    //     for(v_size j = 0; j <= pNum+1 && hNum + 1 + j <= k; j++) {
    //         cnt[hNum + 1 + j] += C[pNum+1][j]*numMax;
    //     }
    //     return;
    // }
// if(d == 0) printf("%u:", section.pR);
    v_size ed = section.pR;
    for(v_size i = pivotDeg; i < ed; i++) {
        v_size v = tmp[d][i];
        
        section.pR--;

        S->changeTo(v, section.pR);
// if(section.pR == pivotDeg) printf("1\n");
// printf("deep%u:search %u %u-%u,h %u\n", d, ed+1, hNum, pNum, section.pR);
        Pair sec = S->sort2(v, section);
        searchForSpecificK(hNum + 1, pNum, sec, S, k, d+1);
    }

}

#endif
// bin\pivoter -edge data\GowallaEdge.bindeg.bin -idx data\Gowallaidx196591.bindeg.bin -v 196591 -p shadow -debug 1 -k 10 -exCnt 106264724
// bin\pivoter.exe -edge data\Stanfordedge.bindeg.bin -idx data\Stanfordidx281903.bindeg.bin -v 281903 -p shadow -k 10 -exCnt 5833322749668 -a 1
// bin\pivoter.exe -edge data\ljedge.bindeg.bin -idx data\ljidx4036538.bindeg.bin -v 4036538 -p shadow -k 10 -a 1
// bin\pivoter.exe -edge data\Ama0601edge.bindeg.bin -idx data\Ama0601idx403394.bindeg.bin -v 403394 -p shadow -k 5 -exCnt 3606466
// bin\pivoter.exe -edge data\Googleedge.bindeg.bin -idx data\Googleidx916428.bindeg.bin -v 916428 -p shadow -k 10 -exCnt 13006798813
// bin\pivoter.exe -idx data\as-skitter_idx_1696415.bindeg.bin -edge data\as-skitter_edge_1696415.bindeg.bin -v 1696415 -p shadow -exCnt 14217188170569 -k 10
// bin\pivoter.exe -edge data\BerkStanedge.bindeg.bin -idx data\BerkStanidx685230.bindeg.bin -v 685230 -p shadow -k 10 -exCnt 59636674809701648
// bin\pivoter.exe -idx data\as-skitter_idx_1696415.bindeg.bin -edge data\as-skitter_edge_1696415.bindeg.bin -v 1696415 -p shadow -debug 111111111 -k 10 
// bin\pivoter.exe -edge data\BerkStanedge.bindeg.bin -idx data\BerkStanidx685230.bindeg.bin -v 685230 -p shadow -k 10 -debug 111111111
// void samplePlusExact::sample() {
//     for(v_size i = 0; i < subG->vCnt; i++) {
//         rng[i] = i;
//     }

//     v_size Xr = 0;
//     v_size sampleTimes = std::max(subG->vCnt * subG->vCnt, 10000u);
//     for(v_size i = 0; i < sampleTimes; i++) {
//         std::random_shuffle(rng, rng+subG->vCnt);
//         bool f = true;

//         for(v_size j = 1; j < k - 1; j++) {
//             for(v_size l = 0; l < j; l++) {
//                 if(!G[rng[j]][rng[l]]) {
//                     f = false; break;
//                 }
//             }
//             if(!f) break;
//         }

//         if(f) Xr++;
//     }

//     printf("shadow %f\n", Xr * C[subG->vCnt][k - 1] / sampleTimes);
// };

// auto sample2 = [&](vector<v_size> & nodes, v_size sz)->double {
//     v_size Xr = 0;
//     v_size * rng = new v_size[sz];
//     for(v_size i = 0; i < sz; i++) {
//         rng[i] = i;
//     }

//     v_size sampleTimes = sz;

//     for(v_size i = 0; i < sampleTimes; i++) {
//         std::random_shuffle(rng, rng + sz);
//         bool f = true;

//         for(v_size j = 1; j < k-1; j++) {
//             v_size u = nodes[rng[j]];
//             for(v_size m = 0; m < j; m++) {
//                 if(!cc(u, nodes[rng[m]])) {
//                     f = false; break;
//                 }
//             }
//             if(!f) break;
//         }

//         if(f) Xr++;
//     }

//     delete [] rng;

//     double c = 1.0;//sz!/k-1!(sz-k+1)!
//     for(v_size i = sz; i >= sz - k + 2; i--) {
//         c *= 1.0 * i / (i - sz + k - 1);
//     }

//     return 1.0*Xr/sampleTimes*c;
// };

//lj
// begin run, dengenerarcy 360, vCnt 4036538
// 2462
// sptimes 1724076
// times 799.320000s
// shadow size 2462
// 10-clique: 19000867068328607744.0

// D:\codes\SCANwindows>bin\pivoter.exe -edge data\BerkStanedge.bindeg.bin -idx data\BerkStanidx685230.bindeg.bin -v 685230 -p msk -k 10 -debug 111111111
// 5.098000 s
// 2-clique: 6649470.00
// 3-clique: 64690980.00
// 4-clique: 1065796916.00
// 5-clique: 21870178738.00
// 6-clique: 460155286971.00
// 7-clique: 9398610960254.00
// 8-clique: 183727787254533.00
// 9-clique: 3410179718341183.00
// 10-clique: 59636674809701648.00

// shadowsObj.init(sz, g->degeneracy, nodes, g);

        // for(v_size i = 0; i < sz; i++) {
        //     v_size u = nodes[i];

        //     v_size l = k - 1;
        //     e_size eCnt;
        //     v_size vCnt = g->pIdx[u] - g->pIdx2[u];
        //     double experiments;

        //     buildSubGraphTrue(u, shadowsObj.graphs[i], l, eCnt, vCnt, experiments);

        //     shadowsObj.set(i, l, eCnt, vCnt, experiments);
        // }

        // double sampleCnt = shadowsObj.sample();

// double samplePlusExact::sampleColor() {
//     //dp
//     double **dp = new double*[g->cc + 1];
//     for(v_size i = 0; i < g->cc + 1; i++) {
//         dp[i] = new double[k + 1]();
//     }
//     v_size * a = new v_size[g->cc]();
//     for(v_size i = 0; i < g->vCnt; i++) {
//         a[g->color[i]]++;
//     }

//     dp[0][0] = 1.0;
//     for(v_size i = 1; i <= g->cc; i++) {
//         dp[i][0] = 1.0;
//         for(v_size j = 1; j <= k && j <= i; j++) {
//             dp[i][j] = dp[i - 1][j - 1]*a[i - 1] + dp[i - 1][j];
//         }
//     }

//     std::default_random_engine e;
//     std::uniform_real_distribution<double> u(0, 1);
//     v_size sampleTimes = 1000, s = 0, c = 0;
//     v_size * clique = new v_size[k];

//     auto rdC = [&](v_size cl)->v_size {
//         std::uniform_int_distribution<unsigned> u(1, a[cl]);
//         std::default_random_engine e;
//         v_size idx = u(e);

//     };

//     while(s < sampleTimes) {
//         v_size i = g->cc, j = k;
//         v_size l = 0;
        
//         while(j > 0) {
//             if(i < j) break;

//             double total = dp[i][j];
//             double choose = dp[i - 1][j - 1] * a[i - 1];
//             // v_size noChoose = dp[i - 1][j];
//             double p = 1.0 * choose / total;
//             double rd = u(e);

//             if(rd < p) {
//                 i--; j--;
//                 clique[l++] = rdC(i);
//             } 
//             else {
//                 i--;
//             }
//         }

//         if(j == 0) {
//             s++;
//             if(isClique(clique, k)) c++;
//         }
//     }

//     delete [] clique;
//     for(v_size i = 0; i < g->cc + 1; i++) {
//         delete [] dp[i];
//     }
//     delete [] dp;

//     return 1.0 * c / sampleTimes * dp[g->cc][k];
// }


// v_size preT = 0, preTotals = 0;
//         v_size cnt = 0;
        
//         v_size buffer[batchSize];
//         v_size preDP[batchSize];
//         bool vis[batchSize];

//         std::default_random_engine generator;
// 		std::discrete_distribution<int> 
//             distribution(experiments, experiments + sz);

//         while(true) {
//             v_size t = 0, totals = 0;

//             for(v_size i = 0; i < batchSize; i++) {
//                 buffer[i] = distribution(generator);
//                 vis[i] = false;
//             }
//             std::sort(buffer, buffer + batchSize);

//             for(v_size i = 0, j = 0; i < batchSize; i++) {
//                 while(j < batchSize && preDP[j] < buffer[i]) j++;
                
//                 if(j == batchSize) break;
//                 else if(preDP[j] == buffer[i]) {
//                     vis[i] = true;
//                 }
//             }


//             for(v_size i = 0; i < batchSize; i++) {
//                 int ret = sampleByNode(nodes[buffer[i]]);
//                 t += ret;
//                 totals++;
//             }

//             if(abs(1.0*t/totals - 1.0*preT/preTotals) < 0.01) {
//                 cnt++;
//                 if(cnt == 2) break;
//             }
//             else {
//                 preT += t;
//                 preTotals += totals;
//                 cnt = 0;
//             }
//         }


// void samplePlusExact::buildSubGraphTrue(v_size u, Bitmap*subGraph, 
//         v_size &l, e_size & eCnt, 
//         v_size & vCnt, double & experiments) {
//     v_size outDegree = g->pIdx[u+1]-g->pIdx2[u];
//     v_size edges = 0;

//     v_size st = g->pIdx2[u];
//     for(v_size i = st; i < g->pIdx[u+1]; i++) {
//         v_size v = g->pEdge[i];
//         subGraph[i-st].setNowSize(outDegree);
//         subGraph[i-st].Clear();
//         for(v_size j = st; j < g->pIdx[u+1]; j++) {
//             v_size w = g->pEdge[j];

//             if(cc(w, v)) {
//                 subGraph[i-st].SetBit(j - st);
//                 edges++;
//             }
//         }
//     }

//     vCnt = outDegree;
//     eCnt = edges;
//     experiments = C[outDegree][l];
// }

