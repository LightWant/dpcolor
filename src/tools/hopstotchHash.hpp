#ifndef HOPSTOTCH
#define HOPSTOTCH

#include "../tools/type.hpp"
#ifndef _WIN32
#include <x86intrin.h>
#else
#include <intrin.h>
#endif

#include <immintrin.h>

#include <utility>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

#define EMPTY 0xffffffff

constexpr v_size L2_CACHE_LINE_SIZE = 64;
constexpr v_size H = L2_CACHE_LINE_SIZE / sizeof(v_size);

#ifdef __AVX512__
#define setFunction _mm512_set1_epi32
#define cmpFunction _mm512_cmpeq_epi32_mask
#define loadFunction _mm512_loadu_si512
#endif

#ifdef __AVX512__
const __m256i eightEMPTY = _mm256_set1_epi32(EMPTY);
const __m512i sixteenEMPTY = _mm512_set1_epi32(EMPTY);
#endif

class hopstotchHash {
public:
    v_size * v;
    v_size n; //邻接点数量
    v_size roundN, preZero;
    v_size hashKey = 1e9 + 7;
    
public:
#ifndef __SSE__
#define __SSE__
#endif

#ifdef 	__SSE__
#ifndef __AVX512__
    bool containSIMD(v_size u) {
        v_size hashV = hash(u);
        v_size p = hashV;

        __m256i eightU = _mm256_set1_epi32(u);
        
        auto address = reinterpret_cast<const __m256i *>(v + p);
        // while(true) {
        __m256i eightNextV = _mm256_loadu_si256(address);
        // __mmask8 cmpRes = _mm256_cmpeq_epi32_mask(eightU, eightNextV);
        // if(cmpRes) return true;
        __m256i cmpRes = _mm256_cmpeq_epi32(eightU, eightNextV);
        auto msk = _mm256_movemask_epi8(cmpRes);
#ifdef _WIN32
        if(__popcnt(msk) > 0) return true;
#else
        if(_popcnt32(msk) > 0) return true;
#endif
        eightNextV = _mm256_loadu_si256(address + 1);
        // cmpRes = _mm256_cmpeq_epi32_mask(eightU, eightNextV);
        // if(cmpRes) return true;
        cmpRes = _mm256_cmpeq_epi32(eightU, eightNextV);
        msk = _mm256_movemask_epi8(cmpRes);
#ifdef _WIN32
        if(__popcnt(msk) > 0) return true;
#else
        if(_popcnt32(msk) > 0) return true;
#endif 
    
        return false;
    }
#else
    bool containSIMD(v_size u) {
        v_size hashV = hash(u);
        v_size p = hashV;

        __m256i eightU = _mm256_set1_epi32(u);
        
        // while(true) {
        __m256i eightNextV = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(v + p));
        __mmask8 cmpRes = _mm256_cmpeq_epi32_mask(eightU, eightNextV);
        if(cmpRes) return true;    
        
        eightNextV = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(v + p + 8));
        cmpRes = _mm256_cmpeq_epi32_mask(eightU, eightNextV);
        if(cmpRes) return true;    

        return false;
    }
#endif
#endif

    ~hopstotchHash() {
        delete [] v;
    }

    void build(v_size * adjList, v_size n_) {
        n = n_;
        // preZero = __builtin_clz(2 * n - 1);
        // if(n == 0) return;

        if(n > 0) {
            preZero = 0;
            while((1u<<preZero) <= 2*n-1) preZero++;
            preZero = 32 - preZero;

            roundN = 1u << (32 - preZero);
            if(roundN < H) roundN = H;
            v = new v_size[roundN + H];
        }
        else {
            preZero = 32;
            roundN = 0;
            v = new v_size[H];
        }
        
// printf("zero:%u roundN%u, n%u\n", preZero, roundN, n);fflush(stdout);
        memset(v, EMPTY, sizeof(v_size) * (roundN + H));
        for(v_size i = 0; i < n; i++) {
            v_size u = adjList[i];
            if(!insert(u)) {
                bool f = reBuild();
                if(!f) {
                    printf("error build hash table\n");
                    fflush(stdout);
                    exit(-1);
                }
                else insert(u);
            }
        }
        // memcpy(v + roundN, v, sizeof(v_size)*(H-1));
// printf("zero:%u roundN%u\n", preZero, roundN);
        for(v_size i = 0; i < H - 1; i++) {
            v[i+roundN] = v[i];
        }

        for(v_size i = 0; i < n; i++) {
            v_size u = adjList[i];
            if(!contain(u)) {
                printf("error build hash table 2\n");
                fflush(stdout);
                exit(-1);
            }
        }
        // for(v_size i = 0; i < H - 1; i++) {
        //     if(v[i] != v[i + roundN]) {
        //         printf("error st and ed of hashTable\n");fflush(stdout);
        //     }
        // }
    }

    bool insert(v_size u) {
        v_size hashV = hash(u);
        v_size p = hashV;
        v_size i = 0;
        for(; i < roundN; i++) {
            if(v[p] == EMPTY) break;
            p = (p + 1) % roundN;
            // if(v[p] == u) return true;
        }

        while((p - hashV + roundN) % roundN >= H) {
            v_size t = (p - H + 1 + roundN) % roundN;
            bool f = false;

            for(v_size i = t; i != p; i = (i + 1) % roundN) {
                if((p - hash(v[i]) + roundN) % roundN < H) {
                    v[p] = v[i];
                    v[i] = EMPTY;
                    p = i;
                    f = true;
                    break;
                }
            }

            if(!f) return false;
        }

        v[p] = u;
        return true;
    }

    bool containNormal(v_size u) {
        v_size hashV = hash(u);

        for(v_size i = hashV; i < hashV + H; i++) {
            // v_size j = i % roundN;
            if(v[i] == EMPTY) return false;
            if(v[i] == u) return true;
        }

        return false;
    }

    bool contain(v_size u) {

#ifdef __AVX512__
        return containSIMD512(u);
#else
        return containSIMD(u);
        // return containNormal(u);
#endif
    }
#ifdef __AVX512__
    bool containSIMD512(v_size u) {
        v_size hashV = hash(u);
        v_size p = hashV;

        __m512i sixteenU = _mm512_set1_epi32(u);
        
        // while(true) {
        __m512i sixteenNextV = _mm512_loadu_si512(reinterpret_cast<const __m512i *>(v + p));
        __mmask16 cmpRes = _mm512_cmpeq_epi32_mask(sixteenU, sixteenNextV);
        if(cmpRes) return true;    
        // if(_mm512_cmpeq_epi32_mask(sixteenNextV, sixteenEMPTY)) return false;

        //     p += parallelism;
        // }
        
        return false;
    }
#endif
    bool reBuild() {
        v_size * tmpV = v;
        v = new v_size[roundN];
        memset(v, 0xff, sizeof(v_size) * roundN);

        v_size rebuildTimes = 0;
        while(rebuildTimes < 100) {
            bool f = true;
            for(v_size i = rebuildTimes; i < roundN + rebuildTimes; i++) {
                v_size j = i % roundN;
                if(tmpV[j] != EMPTY) {
                    if(!insert(tmpV[j])) {
                        f = false;
                        break;
                    }
                }
            }

            if(f) break;
            else {
                rebuildTimes++;
                memset(v, 0xff, sizeof(v_size) * roundN);
            }
        }
        
        delete [] tmpV;

        if(rebuildTimes == 100) return false;

        return true;
    }
//h(k) = (A∗k mod 2^w) >> (w − r), w = 32, A常数，r=log(roundN)
    v_size hash(v_size t) {
        v_size tmp = hashKey * t;
        v_size ans = tmp >> preZero;
        // assert(ans < roundN);
        return ans;
    }
};

#endif