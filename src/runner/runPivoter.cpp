#include "../tools/getArgs.hpp"

#include "../graph/graph.hpp"

#include "../kClique/ccpath.hpp"
#include "../kClique/pivoterMsk.hpp"

#include <cassert>
#include <string>
#include <iostream>
using std::string;
#include <cassert>
/**
 * nvcc src\runPivoter.cpp -rdc=true --gpu-architecture=sm_70 -I "D:\cuda\NIVIDIA  GPU computing Toolklt_v10.1\CUDA\v10.1\include" -I "D:\cuda\NVIDIA Corporation_v10.1\common\inc" -o bin\runPivoter src\kClique\pivoterGPU.cu  -Xcompiler "/wd 4819"
 * **/
int main(int argc, char * argv[])
{
    argsController * aC = new argsController(argc, argv);
    Graph * g = new Graph();
    g->load(aC->get("-edge"), aC->get("-idx"), atoi(aC->get("-v").c_str()));

    PivoterMsk * pt = new PivoterMsk(g);
    v_size deb = atoi(aC->get("-debug").c_str());
    // if(aC->get("-debug").c_str() != "") deb = ;
    pt->runV(atoi(aC->get("-k").c_str()), deb);

    delete pt;
    delete g;

    return 0;
}