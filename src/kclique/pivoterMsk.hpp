#ifndef PIVOTERMSK_HPP
#define PIVOTERMSK_HPP

#include "../graph/graph.hpp"
#include "../tools/type.hpp"
#include "../tools/hopstotchHash.hpp"
#include "../tools/linearSet.hpp"
#include <cassert>
#include <assert.h>
#include <tuple>
#include <vector>
#include <chrono> 
using std::tuple;
using Pair = std::pair<v_size, v_size>;
using std::vector;

#define pP first
#define pR second

class PivoterMsk {
private:
    Graph * g;
    v_size k;
    double ** cnt = nullptr;
    double ** C = nullptr;
    v_size ** tmp;
    hopstotchHash * hashTable;
    v_size nodeNum = 0;

public:
    PivoterMsk(Graph * g_) { g = g_;}
    ~PivoterMsk() {
        for(v_size i = 0; i <= g->degeneracy; i++) {
            delete [] C[i];
            delete [] cnt[i];
            delete [] tmp[i];
        }
        delete [] tmp;
        delete [] cnt;
        delete [] C;
        delete [] hashTable;
    }

    void runV(v_size k_, v_size deb=111111111) {
    // printf("begin run, dengenerarcy %u, vCnt %u\n", g->degeneracy, g->vCnt);
        hashTable = new hopstotchHash[g->vCnt];
        for(v_size u = 0; u < g->vCnt; u++) {
            // printf("%u\n", u);fflush(stdout);
            if(g->pIdx[u + 1] == g->pIdx[u]) continue;
            hashTable[u].build(g->pEdge + g->pIdx[u], g->pIdx[u + 1] - g->pIdx[u]);
        }
    // printf("build hash tables\n");fflush(stdout);

        k = k_;
        cnt = new double*[g->degeneracy + 1];
        for(v_size i = 0; i <= g->degeneracy; i++)
            cnt[i] = new double[k + 1]();
        tmp = new v_size*[g->degeneracy];
        for(v_size i = 0; i <= g->degeneracy; i++)
            tmp[i] = new v_size[g->degeneracy]();
        C = new double*[g->degeneracy + 1];
        for(v_size i = 0; i <= g->degeneracy; i++) {
            C[i] = new double[k + 1]();
        }
        C[0][0] = 1.0;
        C[1][0] = 1.0;
        C[1][1] = 1.0;
        for(v_size i = 2; i <= g->degeneracy; i++) {
            C[i][0] = 1.0;
            if(i < k + 1) C[i][i] = 1.0;
            for(v_size j = 1; j < i && j < k + 1; j++) {
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
            }
        }
        // printf("Combine\n");
// printf("%u\n", g->pIdx[deb+1]-g->pIdx2[deb]);
        LinearSet * S = new LinearSet(g, hashTable);
    // printf("build linear set\n");

        double timeStape = clock();
        std::fill(cnt[0], cnt[0] + k+1, 0.0);
        cnt[0][1] = g->vCnt;

        if(deb != 111111111u)
        for(v_size u = deb; u < deb+1; u++) {
    // double timeStape = clock();
            Pair section = S->sort(u, {0, g->vCnt});
    v_size deg = g->pIdx[u+1]-g->pIdx2[u];
    printf("deg %u\n", deg);
            search(section, S, 1);
            // searchWithoutPivot(section, S, 1);
    // v_size d = 0;
    // gao(section, S, d);
    // printf("%u %u\n", u, d);
            cnt[0][1] = section.pR;

            for(v_size i = 2; i <= k && i <= section.pR + 1; i++) {
                cnt[0][i] = cnt[1][i-1];
            }
    
    printf("%u:%.10f\n", u, (clock() - timeStape) / CLOCKS_PER_SEC);
    v_size n = 0;
    for(v_size i = 0; i < section.pR; i++) {
        for(v_size j = i + 1; j < section.pR; j++) {
            if(hashTable[(*S)[i]].contain((*S)[j])) {
                n++;
                // printf("%u %u\n", (*S)[i], (*S)[j]);
            }
        } 
    }
    printf("edges %u\n", n);
    printf(" edge density %fs\n", 2.0*n/section.pR );
    fflush(stdout);
        }
        else {
    // printf("all\n");
            for(v_size u = 0; u < g->vCnt; u++) {
    // double timeStape = clock();
    // auto start = std::chrono::high_resolution_clock::now();
                Pair section = S->sort(u, {0, g->vCnt});
                search(section, S, 1);
                // searchWithoutPivot(section, S, 1);
                for(v_size i = 2; i <= k && i <= section.pR + 1; i++) {
                    cnt[0][i] += cnt[1][i-1];
                }
    // auto finish = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed = finish - start;
    // printf("%u:%.10f s\n", u, elapsed.count() );
    // printf("%u:%.10f s\n", u, (clock() - timeStape) / CLOCKS_PER_SEC );
            }
            // fflush(stdout);
        }

        printf("%f s\n", (clock() - timeStape) / CLOCKS_PER_SEC );
        for(v_size i = 2; i <= k; i++)
            printf("%d-clique: %.2f\n", i, cnt[0][i]);
        fflush(stdout);
    }

    void runE(v_size k_, v_size deb=111111111) {
    printf("begin run, dengenerarcy %u, vCnt %u\n", g->degeneracy, g->vCnt);
        hashTable = new hopstotchHash[g->vCnt];
        for(v_size u = 0; u < g->vCnt; u++) {
            hashTable[u].build(g->pEdge + g->pIdx[u], g->pIdx[u + 1] - g->pIdx[u]);
        }
    printf("build hash tables\n");fflush(stdout);

        k = k_;
        cnt = new double*[g->degeneracy + 1];
        for(v_size i = 0; i <= g->degeneracy; i++)
            cnt[i] = new double[k + 1]();
        tmp = new v_size*[g->degeneracy];
        for(v_size i = 0; i <= g->degeneracy; i++)
            tmp[i] = new v_size[g->degeneracy]();
        C = new double*[g->degeneracy + 1];
        for(v_size i = 0; i <= g->degeneracy; i++) {
            C[i] = new double[k + 1]();
        }
        C[0][0] = 1.0;
        C[1][0] = 1.0;
        C[1][1] = 1.0;
        for(v_size i = 2; i <= g->degeneracy; i++) {
            C[i][0] = 1.0;
            if(i < k + 1) C[i][i] = 1.0;
            for(v_size j = 1; j < i && j < k + 1; j++) {
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
            }
        }
        printf("Combine\n");

        LinearSet * S = new LinearSet(g, hashTable);
    printf("build linear set\n");

        e_size outDegs = 0;
        for(v_size i = 0; i < g->vCnt; i++) {
            outDegs += g->pIdx[i + 1] - g->pIdx2[i];
        }    
        Pair * outNei = new Pair[outDegs];
        for(v_size i = 0, j = 0; i < g->vCnt; i++) {
            for(v_size k = g->pIdx2[i]; k < g->pIdx[i + 1]; k++) {
                outNei[j].first = i;
                outNei[j++].second = g->pEdge[k];
            }
        }

        double timeStape = clock();
        std::fill(cnt[0], cnt[0] + k+1, 0.0);
        cnt[0][1] = g->vCnt;

        //#pragma omp parallel for 
        for(v_size p = 0; p < outDegs; p++) {
            v_size u = outNei[p].first;
            v_size v = outNei[p].second;
            v_size j = 0;
            for(v_size i = g->pIdx2[u]; i < g->pIdx[u + 1]; i++) {
                v_size w = g->pEdge[i];
                if(w > v && hashTable[v].contain(w)) {
                    S->changeTo(w, j++);
                } 
            }

            search({0, j}, S, 1);
            for(v_size i = 3; i <= k; i++) {
                cnt[0][i] += cnt[1][i-2];
            }
        }

        printf("%f s\n", (clock() - timeStape) / CLOCKS_PER_SEC );
        for(v_size i = 3; i <= k; i++)
            printf("%d-clique: %.2f\n", i, cnt[0][i]);
    }

    void gao(Pair section, LinearSet * S, v_size &d) {
        if(section.pP == section.pR) return;
        if(d>=3) {
            printf("%u\n", section.pR - section.pP);
            return;
        }
        auto sec = section;
        
        while(true) {
            v_size pivot = (*S)[sec.pP];
            v_size pivotDeg = 0;
            int num = 0;

            for(auto i = sec.pP; i < sec.pR; i++) {
                v_size v = (*S)[i];
                if(g->pIdx[v+1] - g->pIdx[v] > pivotDeg) {
                    v_size tmp = 0;

                    for(auto j = sec.pP; j < sec.pR; j++) {
                        if(hashTable[v].contain((*S)[j]) ) tmp++;
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

            S->changeTo(pivot, sec.pP);

            v_size i = sec.pP+1, j = sec.pP+1;
            while(j < sec.pP + pivotDeg+1) {
                if(hashTable[pivot].contain((*S)[i])) {
                    S->changeTo((*S)[i], j++);
                }
                i++;
            }

            sec.pP++;
            sec.pR = j;

            if(sec.pP == sec.pR) break;
        }
        if(section.pP + 2 < sec.pR) d++;
        for(v_size i = section.pP; i < sec.pR; i++) {
            printf("%u ", (*S)[i]);
        }
        printf("\n");

        gao({sec.pR, section.pR}, S, d);
    }

    void search(Pair section, LinearSet * S, v_size d = 0, bool h=false) {
// {
// v_size n = 0, m =0 ;
// for(v_size i = section.pP; i < section.pR; i++) {
//     for(v_size j = i + 1; j < section.pR; j++) {
//         if(hashTable[(*S)[i]].contain((*S)[j])) {
//             n++;
//             // printf("%u %u\n", (*S)[i], (*S)[j]);
//         }
//         m++;
//     } 
// }
// double edgeDen = 1.0*n/m;
// double threshold = (1.0*k-2-d)/(k-1-d);
// // if(edgeDen < threshold)
// // printf("%u density %f, t %f\n", d, edgeDen, threshold);
// if(edgeDen >= threshold) {
//     printf("%u ", d);
//     return;
// }
// }

        std::fill(cnt[d], cnt[d] + std::min(k + 1, section.pR + 2), 0.0);
        cnt[d][1] = section.pR - section.pP;
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
// printf("%u %u\n", d, pivot);
// printf("pivot%u, left %u,right %u, num %d\n", pivotDeg, section.pR - pivotDeg, numMax);
// fflush(stdout);
        section.pR--;

        if((v_size)numMax == pivotDeg && section.pR == pivotDeg) {
            for(v_size i = 2; i <= pivotDeg + 1 && i <= k; i++) {
                cnt[d][i] = C[pivotDeg + 1][i];
            }
            return;
        }

        // if(pivotDeg == 0) return;
        // if(pivotDeg == (v_size)numMax && pivotDeg == section.pR) {
        //     for(v_size i = 2; i <= k && i <= pivotDeg; i++) {
        //         // cnt[d][i] = cnt[d+1][i-1] + cnt[d+1][i];
        //         cnt[d][i] = C[pivotDeg][i-1] + C[pivotDeg][i];
        //     }
        //     cnt[d][pivotDeg + 1] = cnt[d+1][pivotDeg];
        //     return;
        // }

        search({section.pP, section.pP + pivotDeg}, S, d+1, false);
        
        for(v_size i = 2; i <= k && i <= pivotDeg + 1; i++) {
            cnt[d][i] = cnt[d+1][i-1] + cnt[d+1][i];
        }


        v_size ed = section.pR;
        for(v_size i = section.pP + pivotDeg; i < ed; i++) {
            v_size v = tmp[d][i];
            
            section.pR--;
            S->changeTo(v, section.pR);
// if(section.pR == pivotDeg) printf("1\n");
// printf("deep%u:search %u %u-%u,h %u\n", d, ed+1, hNum, pNum, section.pR);
            Pair sec = S->sort2(v, section);
            if(sec.pR == pivotDeg && ed == section.pP + pivotDeg + 1) {
                for(v_size i = 2; i <= k && i <= pivotDeg + 1; i++) {
                    cnt[d][i] += cnt[d+1][i-1];
                }
                continue;
            }

            search(sec, S, d+1, true);
            for(v_size i = 2; i <= k && i <= sec.pR - sec.pP + 1; i++) {
                cnt[d][i] += cnt[d+1][i-1];
            }

// printf("deep%u:hsearch %u\n", d+1, sec.pR);
// for(v_size i = 0; i < sec.pR; i++) {
//     printf("%u ", (*S)[i]);
// }printf("\n");
// for(v_size i = 1; i <= sec.pR+1; i++) {
//     printf("%.0f ", cnt[d+1][i]);
// }printf("\n");
        }

        // delete [] tmpP;
    }

    void searchWithoutPivot(Pair section, LinearSet * S, v_size d = 0, bool h=false) {

        std::fill(cnt[d], cnt[d] + std::min(k + 1, section.pR + 2), 0.0);
        cnt[d][1] = section.pR - section.pP;
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

        memcpy(tmp[d], S->begin()+section.pP, sizeof(v_size)*(section.pR-section.pP));
        v_size ed = section.pR;
        for(v_size i = section.pP; i < ed; i++) {
            v_size v = tmp[d][i];
            
            section.pR--;
            S->changeTo(v, section.pR);
// if(section.pR == pivotDeg) printf("1\n");
// printf("deep%u:search %u %u-%u,h %u\n", d, ed+1, hNum, pNum, section.pR);
            Pair sec = S->sort2(v, section);
            // if(sec.pR == pivotDeg && ed == section.pP + pivotDeg + 1) {
            //     for(v_size i = 2; i <= k && i <= pivotDeg + 1; i++) {
            //         cnt[d][i] += cnt[d+1][i-1];
            //     }
            //     continue;
            // }

            search(sec, S, d+1, true);
            for(v_size i = 2; i <= k && i <= sec.pR - sec.pP + 1; i++) {
                cnt[d][i] += cnt[d+1][i-1];
            }

// printf("deep%u:hsearch %u\n", d+1, sec.pR);
// for(v_size i = 0; i < sec.pR; i++) {
//     printf("%u ", (*S)[i]);
// }printf("\n");
// for(v_size i = 1; i <= sec.pR+1; i++) {
//     printf("%.0f ", cnt[d+1][i]);
// }printf("\n");
        }

        // delete [] tmpP;
    }


};

#undef pP
#undef pR

#endif

//             v_size pR = section.pR;
//             vector<v_size>* v = new vector<v_size>[pR];
//             bool * f = new bool[pR]();
//             v_size c = 0, i = 0, p = 0;
//             while(true) {
//                 while(i < pR && f[i]) i++;
//                 if(i == pR) break;
// // printf("%u ", i);
//                 v[p].push_back((*S)[i]);
//                 f[i] = true;

//                 for(v_size j = i + 1; j < pR; j++) {
//                     if(f[j]) continue;
//                     v_size w = (*S)[j];
//                     bool ok = true;
//                     for(v_size l = 0; l < v[p].size(); l++) {
//                         v_size u = v[p][l];
//                         if(!hashTable[w].contain(u)) {
//                             ok = false;
//                             break;
//                         }
//                     }
//                     if(ok) {
//                         f[j] = true;
//                         v[p].push_back(w);
//                     }
//                 }

//                 p++;
//             }

//             printf("%u\n", p);
//             delete [] v;
//             delete [] f;