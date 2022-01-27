CXXFLAGS= -O3 -Wall -std=c++11 -g -fopenmp -march=native

all:bin/run

bin/run:src/runner/run.cpp $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $<

bin/makeTxt:src/runner/makeTxt.cpp
	g++ $(CXXFLAGS) -o $@ $<

bin/makeCSR:src/runner/makeCSR.cpp
	g++ $(CXXFLAGS) -o $@ $<
bin/sortData:src/runner/sortData.cpp
	g++ $(CXXFLAGS) -o $@ $<

bin/changeToD:src/runner/changeToD.cpp $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $<
bin/tmp:src/runner/tmp.cpp $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $<


clean: 
	rm bin/*
