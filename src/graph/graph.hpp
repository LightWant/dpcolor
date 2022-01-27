#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <algorithm>
#include <string>
#include <vector>
#include <thread>

//#include "../tools/constants.hpp"
// #include "../tools/filesystem.hpp"
// #include "../tools/mpi.hpp"
// #include "../tools/time.hpp"
#include "../tools/type.hpp"
#include "../tools/fastIO.hpp"
// #include "../tools/atomic.hpp"
#include "../tools/listLinearHeap.hpp"

class Graph
{
public:
    v_size *pEdge = nullptr;
    v_size *pIdx = nullptr, *pIdx2 = nullptr;
    e_size *pSumDeg = nullptr;
    // v_size *partition_offset;
    v_size vCnt;
    v_size eCnt;
    v_size maxDeg;
    v_size degeneracy = 0;

    v_size *degNum = nullptr;
    v_size * color = nullptr, cc = 0;
    // multiThreads * threadController;
    v_size *partitionOffset = nullptr;

    v_size *pDEdge = nullptr;
    v_size *pDIdx = nullptr;
    v_size *mp = nullptr, *mp2 = nullptr;
    e_size outDegSum;

public:
    Graph() {
        // threadController = new multiThreads();
        // partitionOffset = new v_size[threadController->totalThreads + 1];
    }

    void load(std::string pEdgePath, std::string pIdxPath, v_size vCnt_) {
        fastIO<v_size> * readEdge = new fastIO<v_size>(pEdgePath, "rb");
        fastIO<v_size> * readIdx = new fastIO<v_size>(pIdxPath, "rb");

        vCnt = vCnt_;
        pIdx = new v_size[vCnt + 1];
        // pSumDeg = new e_size[vCnt + 1];
        pIdx2 = new v_size[vCnt + 1];
        
        //load pIdx
        for(v_size i = 0; i <= vCnt; i++) {
            pIdx2[i] = pIdx[i] = readIdx->getFromBin();
// if(i==7 || i == 8)
// printf("%u \n", pIdx[i]);
        }
        delete readIdx;
        eCnt = pIdx[vCnt];
        // printf("edges %u\n", eCnt);

        pEdge = new v_size[eCnt];

        for(e_size i = 0; i < eCnt; i++) {
            pEdge[i] = readEdge->getFromBin();
            // if(pEdge[i] == (unsigned int)EOF) {
            //     printf("EOF %llu\n", i);
            //     exit(0);
            // }
        }
        delete readEdge;

        degeneracy = 0;
        outDegSum = 0;
        for(v_size i = 0; i < vCnt; i++) {
            while(pIdx2[i] < pIdx[i + 1] && pEdge[pIdx2[i]] < i) pIdx2[i]++;
            if(pIdx[i + 1] < pIdx2[i]) {
                printf("error pIdx2 %u\n", i);
                continue;
            }
            degeneracy = std::max(degeneracy, pIdx[i + 1] - pIdx2[i]);
            outDegSum += pIdx[i + 1] - pIdx2[i];
            maxDeg = std::max(maxDeg, pIdx[i + 1] - pIdx[i]);
        }

        // e_size sumdeg = 0;
        // e_size maxDegTmp = 0;
        //compute pSumDeg
        // #pragma omp parallel for reduction(+:sumdeg) reduction(max:maxDegTmp)
        // for(v_size i = 0; i < vCnt; i++) {
        //     e_size cnt = 0;
        //     for(v_size j = pIdx[i]; j < pIdx[i + 1]; j++) {
        //         v_size u = pEdge[j];
        //         cnt += pIdx[u+1] - pIdx[u];
        //     }
        //     // if(pIdx[i + 1] - pIdx[i] > maxDeg) write_max<v_size>((v_size*)(&maxDeg), pIdx[i + 1] - pIdx[i]) ;
        //     if(pIdx[i + 1] - pIdx[i] > maxDegTmp) maxDegTmp = pIdx[i + 1] - pIdx[i];
        //     // pSumDeg[i] = cnt;
        //     sumdeg += cnt;
        // }
        // maxDeg = maxDegTmp;

        // degNum = new v_size[maxDeg + 1]();

        // // #pragma omp parallel for reduction(+:degNum) 
        // for(v_size u = 0; u < vCnt; u++) {
        //     // ++degNum[pIdx[u + 1] - pIdx[u]];
        //     // write_add<v_size>(&degNum[pIdx[u + 1] - pIdx[u]], 1);
        //     degNum[pIdx[u + 1] - pIdx[u]]++;
        // }
        // for(v_size u = maxDeg - 1; u >= 2; u--) {
        //     degNum[u] += degNum[u + 1];
        // }

        // printf("threads %d\n", threadController->totalThreads);
        // partitionOffset = new v_size[threadController->totalThreads + 1];
        // partitionOffset[0] = 0;
        // v_size u = 0;
        // e_size partedEdges = 0;
        // for(int i = 1; i < threadController->totalThreads; i++) {
        //     e_size expected = partedEdges + 
        //         (sumdeg - partedEdges) / (threadController->totalThreads - i + 1);
        //     if(expected > sumdeg) expected = sumdeg;

        //     while(partedEdges < expected) {
        //         partedEdges += pSumDeg[u++];
        //     }
        //     partitionOffset[i] = u;
        //     // printf("%d %lu %u\n", i, partedEdges, u);
        // }
        
        // partitionOffset[threadController->totalThreads] = vCnt;
        // printf("%d %lu %u\n", threadController->totalThreads, sumdeg, vCnt);
    }

    void changeToDegeneracy(std::string pEdgePath, std::string pIdxPath) {
// printf("%u %lu\n", vCnt, eCnt);
        ListLinearHeap lheap(vCnt, eCnt);
        v_size * ids = new v_size[vCnt];
        v_size * keys = new v_size[vCnt + 1];
        for(v_size i = 0; i < vCnt; i++) {
            ids[i] = i;
            keys[i] = pIdx[i + 1] - pIdx[i] + 1;
        }
        lheap.init(vCnt, vCnt, ids, keys);
        mp = ids; //v to v
        mp2 = new v_size[vCnt];

        pDIdx = keys;
        pDEdge = new v_size[eCnt];

        for(v_size i = 0; i < vCnt; i++) {
            v_size v, degV;

            if(!lheap.pop_min(v, degV)) printf("error\n");
// printf("%u %u\n", v, degV-1);
            mp[i] = v; mp2[v] = i;
            for(v_size j = pIdx[v]; j < pIdx[v + 1]; j++) {
                lheap.decrement(pEdge[j]);
            }
        }
        pDIdx[0] = 0;
        for(v_size i = 1; i <= vCnt; i++) {
            v_size v = mp[i - 1];
            pDIdx[i] = pDIdx[i - 1] + pIdx[v + 1] - pIdx[v];
        }

        for(v_size i = 0; i < vCnt; i++) {
            v_size k = pIdx[mp[i]];
            for(v_size j = pDIdx[i]; j < pDIdx[i + 1]; j++) {
                pDEdge[j] = mp2[pEdge[k++]];
            }
            std::sort(pDEdge + pDIdx[i], pDEdge + pDIdx[i + 1]);
        }

        delete [] mp;
        delete [] mp2;

        FILE * fEdge = fopen(pEdgePath.c_str(), "wb");
        FILE * fIdx = fopen(pIdxPath.c_str(), "wb");
        fwrite(pDEdge, 4, eCnt, fEdge);
        fwrite(pDIdx, 4, vCnt + 1, fIdx);
        fclose(fEdge);
        fclose(fIdx);
    }

    ~Graph() {
        // printf("delete graph\n");
        if(pEdge != nullptr) delete [] pEdge;
        if(pIdx != nullptr) delete [] pIdx;
        if(pIdx2 != nullptr) delete [] pIdx2;
        // delete [] partitionOffset;
        if(pSumDeg != nullptr) delete [] pSumDeg;
        if(degNum != nullptr) delete [] degNum;
        // delete threadController;

        if(pDEdge != nullptr) {
            delete [] pDEdge;
            delete [] pDIdx;
        }
        
        if(color != nullptr) {
            delete [] color;
        }
    }

    void colorG() {
        color = new v_size[vCnt]();//0-cc
        v_size *f = new v_size[maxDeg]();
        
        cc = 0;
        for(v_size u = vCnt - 1; u >= 0; u--) {
            for(v_size j = pIdx2[u]; j < pIdx[u+1]; j++) {
                v_size v = pEdge[j];
                f[color[v]] = u;   
            }
            
            v_size c = 0;
            while(f[c] == u) c++;
            cc = std::max(cc, c);
            color[u] = c;

            if(u == 0) break;
        }
        cc++;

        delete [] f;

        // for(v_size i = 0; i < vCnt; i++) {
        //     for(v_size j = pIdx2[i]; j < pIdx[i+1]; j++) {
        //         assert(color[i] != color[pEdge[j]]);  
        //     }
        // }
    }

    void colorG2() {
        if(color == nullptr) color = new v_size[vCnt]();
        else
        for(v_size u = 0; u < vCnt; u++) {
            color[u] = 0;
        }

        cc = 0;
        for(v_size u = vCnt - 1; u >= 0; u--) {

            if(pIdx[u+1] == pIdx2[u]) {
                if(u == 0) break;
                continue;
            }

            v_size maxC = 0;
            for(v_size j = pIdx2[u]; j < pIdx[u+1]; j++) {
                v_size v = pEdge[j];
                maxC = std::max(maxC, color[v]);
            }

            cc = std::max(cc, maxC + 1);
            color[u] = maxC + 1;

            if(u == 0) break;
        }

        cc++;
    }
};



#endif
