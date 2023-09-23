# The work has been extended to a journal paper "Efficient k-Clique Counting on Large Graphs: The Power of Color-Based Sampling Approaches". https://github.com/LightWant/tripath_tkde.

The source code of the work "Lightning Fast and Space Efficient k-clique Counting, www'22, https://doi.org/10.1145/3485447.3512167".

# kclique counting

## Build

makedir bin

make

make changeToD

make makeCSR

### Data format

use com-DBLP as the example:

At first, we have text file 

​	./data/dblp/dblp.txt:

n,m

0 1

0 2

2 3

...



run: bin/makeCSR ./data/dblp/dblp.txt ./data/dblp/tmpedge.bin ./data/dblp/tmpidx.bin, and we get:

​	./data/dblp/tmpedge.bin ./data/dblp/tmpidx.bin

run bin/changeToD -edge ./data/dblp/tmpedge.bin  -idx ./data/dblp/tmpidx.bin -v 425957

-v is the node num.

and we get ./data/dblp/tmpedge.bindeg.bin ./data/dblp/tmpidx.bindeg.bin

rename ./data/dblp/tmpedge.bindeg.bin->./data/dblp/edge.bin

rename ./data/dblp/tmpidx.bindeg.bin->./data/dblp/idx.bin

add file s.txt under the ./data/dblp/:

425957
30
2-clique: 1049866.00
3-clique: 2224385.00
4-clique: 16713192.00

...  

30-clique: 2952522344434315298898706432.00

The next is the s.txt for com-lj:

4036538
8
3-clique: 177820130.00
4-clique: 5216918441.00
5-clique: 246378629120.00
6-clique: 10990740312954.00
7-clique: 449022426169164.00
8-clique: 16890998195437619.00

At last,  ./data/dblp/edge.bin, ./data/dblp/idx.bin and ./data/dblp/s.txt are what we need.



If you don't have the exact count of the cliques, just write the node number:

./data/dblp/s.txt:

425957



### run

dpcolorpath: bin/run -f data folder -k clique size -N number of samples -cccpath

bin/run -f data/dblp/ -k 8 -N 100000 -cccpath

output:|k| alpha| count of samples|time of exact part| size of exact part| the depth of search| count in exact part|not expected:the lack count | appDensity hittedcount real_sample_times  the percentage of the cliques in the sparse part| approximate count | the percentage of the cliques in the dense part among the approximate count | sampling running time| total running time| error

example:|8| 1.0| 100000| 0.03| 1739| d13| 5008853.0|not expected 53 | 1.000000 100000 100000 99.99935555%| 777243510729| 0.00%| 0.07| 0.10| 0.00%



dpcolor:bin/run -f data folder -k clique size -N number of samples -cc

bin/run -f data/dblp/ -k 8 -N 100000 -cc



parallel dpcolorpath: bin/run -f data folder -k clique size -N number of samples -deb 3 -p threads

bin/run -f data/dblp/ -k 8 -N 100000 -deb 3 -p 10



parallel dpcolor:bin/run -f data folder -k clique size -N number of samples -cc -deb 3 -p threads

bin/run -f data/dblp/ -k 8 -N 100000 -deb 3 -cc -p 10

### remark that the Randomness lacks in this version.
