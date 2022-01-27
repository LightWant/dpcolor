#include "../graph/graph.hpp"
// #include "./maxCliques/cliques.hpp"
// #include "./maxCliques/BK.hpp"
#include "../tools/getArgs.hpp"

int main(int argc, char *argv[])
{
    argsController * aC = new argsController(argc, argv);

    Graph * g = new Graph();
    g->load(aC->get("-edge"), aC->get("-idx"), atoi(aC->get("-v").c_str()));
    printf("load\n");

    g->changeToDegeneracy(aC->get("-edge")+"deg.bin" , aC->get("-idx")+"deg.bin");

    delete g;

    return 0;
}
