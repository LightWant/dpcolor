#include "../tools/getArgs.hpp"

#include "../graph/graph.hpp"

#include <cassert>
#include <string>
#include <iostream>
using std::string;
#include <cassert>

int main(int argc, char * argv[])
{
    argsController * aC = new argsController(argc, argv);
    Graph * g = new Graph();

    string filePath = aC->get("-f");
    v_size n;

    FILE * f = fopen((filePath+"s.txt").c_str(), "r");
    fscanf(f, "%u", &n);
    fclose(f);
    
    g->load(filePath+"edge.bin", filePath+"idx.bin", n);

    f = fopen((filePath+"data.txt").c_str(), "w");
    fprintf(f, "%u %u\n", g->vCnt, g->eCnt/2);
    for(v_size i = 0; i < g->vCnt; i++) {
        for(v_size j = g->pIdx2[i]; j < g->pIdx[i+1]; j++) {
            fprintf(f, "%u %u\n", i, g->pEdge[j]);
        }
    }
    fclose(f);

    delete g;

    return 0;
}
