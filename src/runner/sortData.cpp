#include <vector>
#include <algorithm>
#include "../tools/fastIO.hpp"
using namespace std;

typedef pair<unsigned, unsigned> P;

int main(int argc, char ** argv)
{
    // FILE *fp = fopen(argv[1], "r");
    FILE *fp2 = fopen(argv[2], "w");

    fastIO<char> red(argv[1], "r");

    unsigned n, m, u, v;
    red.getUInt(n);

    red.getUInt(m);
printf("%u %u", n, m);

    vector<P> edges;

    edges.resize(m);
printf("there\n");fflush(stdout);

    n = 0;
    for(unsigned i = 0; i < m; i++) {
        red.getUInt(u);
        red.getUInt(v);
        if(u > v) swap(u, v);

        edges[i].first = u;
        edges[i].second = v;

        n = max(n, v);
    }
    n++;

    sort(edges.begin(), edges.end());
printf("there\n");fflush(stdout);
    unsigned l = 0;
    for(unsigned i = 1; i < m; i++) {
        if(edges[i] == edges[l]) continue;
        l++;
        swap(edges[i], edges[l]);
    }
printf("%u %u %u\n", n, m, l);
    fprintf(fp2, "%u %u\n", n, l + 1);
    for(unsigned i = 0; i <= l; i++) {
        fprintf(fp2, "%u %u\n", edges[i].first, edges[i].second);
    }

    fclose(fp2);

    return 0;
}