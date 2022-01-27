
#ifndef LINEARSET
#define LINEARSET
#include <tuple>
#include <utility>
#include "type.hpp"
#include "hopstotchHash.hpp"
#include "../graph/graph.hpp"

using Pair = std::pair<v_size, v_size>;

#define pP first
#define pR second


class LinearSet {
private:
    v_size * vSet;
    v_size * fIndex;
    Graph * g;
    hopstotchHash * hashTable;
    bool * vis;
    
public:
    LinearSet(Graph * g_, hopstotchHash * hashTable_) {
        g = g_;
        hashTable = hashTable_;
        vSet = new v_size[g->vCnt];
        fIndex = new v_size[g->vCnt];
        for(v_size i = 0; i < g->vCnt; i++) {
            vSet[i] = fIndex[i] = i;
        }
        vis = new bool[g->degeneracy];
    }

    ~LinearSet() { delete [] fIndex; delete [] vSet; delete [] vis;}
    v_size * begin() {
        return vSet;
    }
    bool isIn(v_size v, v_size l, v_size r) {
        return l <= fIndex[v] && fIndex[v] < r;
    }

    v_size operator [] (v_size i) {
        if(i >= g->vCnt) {
            printf("error index\n"); return -1;
        }
        return vSet[i];
    }

    void changeTo(v_size u, v_size p) {
        v_size pU = fIndex[u];
        std::swap(fIndex[u], fIndex[vSet[p]]);
        std::swap(vSet[pU], vSet[p]);
    }

    void changeToByPos(v_size pU, v_size p) {
        std::swap(fIndex[vSet[pU]], fIndex[vSet[p]]);
        std::swap(vSet[pU], vSet[p]);
    }

    Pair sort(v_size u, const Pair & sec) {
        v_size l, r;
        l = r = sec.pP;

        for(v_size i = g->pIdx2[u]; i < g->pIdx[u + 1]; i++) {
            v_size v = g->pEdge[i];
            // if(sec.pP <= fIndex[v] && fIndex[v] < sec.pR) {
            changeTo(v, r++);
            // }
        }

        return {l, r};
    }

    Pair sort2(v_size u, const Pair & sec) {
        v_size i = sec.pP, j = sec.pP;
        while(i < sec.pR && j < sec.pR) {
            if(hashTable[u].contain(vSet[i])) {
                changeTo(vSet[i], j++);
            }
            i++;
        }

        return {sec.pP, j};
    }

    Pair sort3(v_size u, const Pair & sec, v_size p, v_size & newP_) {
        //clique
        v_size l = p;
        for(v_size i = p - 1; i >= sec.pP; i--) {
            if(hashTable[u].contain(vSet[i])) {
                changeToByPos(i, --l);
            }
        }

        v_size i = p, j = p;
        while(i < sec.pR) {
            if(hashTable[u].contain(vSet[i])) {
                changeTo(vSet[i], j++);
            }
            i++;
        }

        v_size newP = p;
        for(v_size i = p; i < j; i++) {
            bool f = true;
            for(v_size j = l; j < newP; j++) {
                if(!hashTable[vSet[i]].contain(vSet[j])) {
                    f = false; break;
                }
            }
            if(f) {
                changeToByPos(i, newP++);
            }
        }
        newP_ = newP;
        return {l, j};
    }

    v_size findClique(const Pair & sec) {
        v_size pivot = vSet[sec.pP];
        v_size pivotDeg = 0;

        for(auto i = sec.pP; i < sec.pR; i++) {
            v_size v = vSet[i];
            if(g->pIdx[v+1] - g->pIdx[v] > pivotDeg) {
                v_size tmp = 0;

                for(auto j = sec.pP; j < sec.pR; j++) {
                    if(hashTable[v].contain(vSet[j]) ) tmp++;
                }

                if(tmp > pivotDeg) {
                    pivot = v; pivotDeg = tmp;
                }
            }
        }

        changeTo(pivot, sec.pP);

        v_size l = sec.pP, r = sec.pR;
        do {
            for(v_size i = l + 1; i < r; ) {
                if(!hashTable[vSet[l]].contain(vSet[i])) {
                    changeToByPos(i, --r);
                }
                else i++;
            }
            
            l++;
        }while(l < r);
        return r;
    }

    Pair mtMaxClique(const Pair & sec, v_size p) {
        //p >= 2
        v_size pivot = vSet[sec.pP];
        // changeToByPos(p - 1, sec.pR - 1);

        v_size j = p;
        for(v_size i = p; i < sec.pR; i++) {
            if(hashTable[pivot].contain(vSet[i])) {
                changeToByPos(i, j++);
            }
        }
        
        v_size newP = p;
        if(p == sec.pP + 1) newP += 1;
        for(v_size i = newP; i < j; i++) {
            bool f = true;
            for(v_size j = sec.pP + 1; j < newP; j++) {
                if(!hashTable[vSet[i]].contain(vSet[j])) {
                    f = false; break;
                }
            }
            if(f) changeToByPos(i, newP++);
        }

        return {newP, j};
    }

    void findPivot(const Pair & sec, v_size & pivot_, v_size & pivotDeg_) {
        v_size pivot = vSet[sec.pP];
        v_size pivotDeg = 0;

        for(auto i = sec.pP; i < sec.pR; i++) {
            v_size v = vSet[i];
            if(g->pIdx[v+1] - g->pIdx[v] > pivotDeg) {
                v_size tmp = 0;

                for(auto j = sec.pP; j < sec.pR; j++) {
                    if(hashTable[v].contain(vSet[j]) ) tmp++;
                }

                if(tmp > pivotDeg) {
                    pivot = v; pivotDeg = tmp;
                }
            }
        }

        changeTo(pivot, sec.pR - 1);

        v_size i = sec.pP, j = sec.pP;
        while(j < sec.pP + pivotDeg) {
            if(hashTable[pivot].contain(vSet[i])) {
                changeTo(vSet[i], j++);
            }
            i++;
        }
// printf("%u %u %u %u\n", j, i, sec.pR, vSet[sec.pR - 1]);
        pivot_ = pivot;
        pivotDeg_ = pivotDeg;
// if(sec.pP > g->vCnt) {
//     printf("%u\n", sec.pP);fflush(stdout);
// }
        // v_size pivot = vSet[sec.pP];
        // v_size pivotDeg = sec.pR - sec.pP;
        // v_size mid = g->degeneracy;
        // memset(vis, false, sizeof(bool)*(sec.pR - sec.pP));

        // for(auto i = sec.pP; i < sec.pR; i++) {
        //     v_size v = vSet[i];
        //     v_size tmp = sec.pP, tmp2 = 0;

        //     for(auto j = sec.pP; j < sec.pR; j++) {
        //         v_size u = vSet[j];
        //         if(hashTable[u].contain(v)) {
        //             // changeToByPos(j, tmp++);
        //             tmp++;
        //             for(auto k = j + 1; k < sec.pR; k++)  {
        //                 v_size w = vSet[k];
        //                 if(!vis[k - sec.pP] && hashTable[u].contain(w)) {
        //                     vis[k - sec.pP] = true;
        //                     tmp2++;
        //                 }
        //             }
        //         }
        //     }

        //     if(tmp2 < mid) {
        //         mid = tmp2;
        //         pivot = v;
        //         pivotDeg = tmp;
        //     }
        // }

        // changeTo(pivot, sec.pR - 1);

        // v_size i = sec.pP, j = sec.pP;
        // while(j < sec.pP + pivotDeg) {
        //     if(hashTable[pivot].contain(vSet[i])) {
        //         changeTo(vSet[i], j++);
        //     }
        //     i++;
        // }

        // pivot_ = pivot;
        // pivotDeg_ = pivotDeg;
    }
    
    int findPivotAndCopyToTmpMem(const Pair & sec, v_size * tmpP, v_size & pivot_, v_size & pivotDeg_) {
        v_size pivot = vSet[sec.pP];
        v_size pivotDeg = 0;
        int num = 0;

        for(auto i = sec.pP; i < sec.pR; i++) {
            v_size v = vSet[i];
            if(g->pIdx[v+1] - g->pIdx[v] > pivotDeg) {
                v_size tmp = 0;

                for(auto j = sec.pP; j < sec.pR; j++) {
                    if(hashTable[v].contain(vSet[j]) ) tmp++;
                }

                if(tmp > pivotDeg) {
                    pivot = v; pivotDeg = tmp; num = 0;
                }
                else if(tmp == pivotDeg){
                    num++;
                }
            }
        }

        changeTo(pivot, sec.pR - 1);

        v_size i = sec.pP, j = sec.pP;
        while(j < sec.pP + pivotDeg) {
            if(hashTable[pivot].contain(vSet[i])) {
                changeTo(vSet[i], j++);
            }
            i++;
        }

        pivot_ = pivot;
        pivotDeg_ = pivotDeg;

        memcpy(tmpP, vSet + sec.pP, sizeof(v_size) * (sec.pR - sec.pP));

        return num;
    }

    v_size findDensity(const Pair & section) {
        v_size n = 0;
        for(v_size i = section.pP; i < section.pR; i++) {
            for(v_size j = i + 1; j < section.pR; j++) {
                if(hashTable[vSet[i]].contain(vSet[j])) {
                    n++;
                }
            } 
        }
        return 2*n/section.pR;
    }

    int findPivotsAndClique(const Pair & sec, v_size * tmpP, v_size & pivot_, v_size & pivotDeg_) {
        v_size pivot = vSet[sec.pP];
        v_size pivotDeg = 0;
        int num = 0;

        // int a[20];
        // bool b[20];
        for(auto i = sec.pP; i < sec.pR; i++) {
            v_size v = vSet[i];
            if(g->pIdx[v+1] - g->pIdx[v] > pivotDeg) {
                v_size tmp = 0;

                for(auto j = sec.pP; j < sec.pR; j++) {
                    if(hashTable[v].contain(vSet[j]) ) tmp++;
                }

                if(tmp > pivotDeg) {
                    pivot = v; pivotDeg = tmp; num = 0;
                }
                else if(tmp == pivotDeg){
                    // a[num++] = v;
                    num++;
                }
            }
        }

        // for(v_size i = 0; i < num; i++) {
        //     if(hashTable[pivot].contain(a[i])) {
        //         b[i] = true;
        //     }
        //     else b[i] = false;
        // }
        // for(v_size i = 0; i < num; i++) {
        //     if(!b[i]) continue;
        //     for(v_size j = i + 1; j < num; j++) {
        //         if(!b[i]) continue;
        //         if(!hashTable[a[i]].contain(a[j])) b[j] = false;
        //     }
        // }
        // int cnt = 0;
        // for(v_size i = 0; i < num; i++) {
        //     if(b[i]) cnt++;
        // }
        // printf("c %d ", cnt);

        changeTo(pivot, sec.pR - 1);

        v_size i = sec.pP, j = sec.pP;
        while(j < sec.pP + pivotDeg) {
            if(hashTable[pivot].contain(vSet[i])) {
                changeTo(vSet[i], j++);
            }
            i++;
        }

        pivot_ = pivot;
        pivotDeg_ = pivotDeg;

        memcpy(tmpP, vSet + sec.pP, sizeof(v_size) * (sec.pR - sec.pP));

        return num;
    }

    void copy(v_size * tmpP, const Pair & sec) {
        for(v_size i = 0; i < sec.pR - sec.pP; i++) {
            v_size v = tmpP[i];
            changeTo(v, i);
        }
    }
};

#endif