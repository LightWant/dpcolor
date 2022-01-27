#include "../graph/graph.hpp"
// #include "./maxCliques/cliques.hpp"
// #include "./maxCliques/BK.hpp"
#include "../tools/getArgs.hpp"
#include <string>
using std::string;

int main(int argc, char *argv[])
{
    argsController * aC = new argsController(argc, argv);

    Graph * g = new Graph();
    g->load(aC->get("-edge"), aC->get("-idx"), atoi(aC->get("-v").c_str()));
    // printf("load\n");

    bool *f = new bool[g->vCnt]();
    v_size * reId = new v_size[g->vCnt]();

    for(v_size i = 0; i < g->vCnt; i++) {
        for(v_size j = g->pIdx2[i]; j < g->pIdx[i+1]; j++) {
            f[g->pEdge[j]] = true;
        }
        if(g->pIdx[i+1] > g->pIdx[i]) f[i] = true;
    }

    v_size id = 0;
    for(v_size i = 0; i < g->vCnt; i++) {
        if(!f[i]) continue;
        reId[i] = id++;
    }
    
    printf("%u %u\n", id, g->eCnt / 2);
    for(v_size i = 0; i < g->vCnt; i++) {
        if(!f[i]) continue;
        v_size u = reId[i];
        for(v_size j = g->pIdx2[i]; j < g->pIdx[i+1]; j++) {
            printf("%u %u\n", u, reId[g->pEdge[j]]);
        }
    }

    delete g;

    return 0;
}