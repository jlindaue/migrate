// #include <stdio.h>
#define WITHOUT_NUMPY

#include <string>
#include <iostream>
#include "matplotlibcpp.h"
#include <fstream>

namespace plt = matplotlibcpp;
//using namespace std;


int main(int argc, char **argv) {

    std::string homepath = std::getenv("HOME");
    std::ifstream infile(homepath + "/jlindaue/gradient_timer.txt");

    std::vector<int> Ns;
    std::vector<long> ts;
    int N;
    long t;
    float dist;
    while (infile >> N >> dist >> t){
	Ns.push_back(N);
	ts.push_back(t);
    }
    plt::plot(Ns, ts, "rx");
    plt::ylabel("Time[s]");
    plt::xlabel("N");
    plt::xlim(0,200);
    plt::ylim(0,200000);
    //plt::axis("equal");
    plt::legend();
    plt::show();
    //plt::save(homepath + "/jlindaue/bezier.pdf");
    return 0;
}

