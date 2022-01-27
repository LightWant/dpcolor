#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <atomic>
#include <fcntl.h>
// #include <unistd.h>
#include <bits/stdc++.h>
using namespace std;
#define PAGESIZE (1<<20)
typedef unsigned ui;
typedef unsigned long long ull;

void gao(char * f1, char * f2, char * f3) {
	FILE *fpEdge = fopen(f2, "wb");
	FILE *fpIdx = fopen(f3, "wb");
	FILE *fp = fopen(f1, "r");

	ui n, u, v;
	ui m;
	ui *pEdge = nullptr, *pIdx = nullptr, *pDeg = nullptr;
	ui i = 0, preU = 0, preV = 0;

	fscanf(fp, "%u %u", &n, &m);
	m *= 2;

	pDeg = new ui[n]();
	pEdge = new ui[m];
	pIdx = new ui[n + 1]();

	while(fscanf(fp, "%u %u", &u, &v) != EOF) {
		pDeg[u]++;
		pDeg[v]++;
	}

	pIdx[0] = 0;
	for(ui u = 1; u <= n; u++) {
		pIdx[u] = pIdx[u - 1] + pDeg[u - 1];
	}
	
	fseek(fp, 0, SEEK_SET);
	fscanf(fp, "%u%u", &u, &v);
	printf("%u %u\n", n, m);

	while(fscanf(fp, "%u %u", &u, &v) != EOF) {
		pEdge[pIdx[u]++] = v;
		pEdge[pIdx[v]++] = u;
	}

	pIdx[0] = 0;
	for(ui u = 1; u <= n; u++) {
		pIdx[u] = pIdx[u - 1] + pDeg[u - 1];
	}

	fclose(fp);

	fwrite(pIdx, 4, n + 1, fpIdx);
	fclose(fpIdx);

	fwrite(pEdge, 4, m, fpEdge);
	fclose(fpEdge);

	delete [] pIdx;
	delete [] pEdge;
	delete [] pDeg;	
}

int main(int argc, char ** argv)
{
	// gaoOrkut();
	// gaoTwitter();
	// testTwitter();
	// gaomy(argv);
	// gaotmp();
	// gaoAm();
	// gaoOrkut();
	gao(argv[1], argv[2], argv[3]);

    return 0;
}
