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
    std::ifstream infile(homepath + "/jlindaue/migrate/data/gradient_timer.txt");

    std::vector<int> Ns;
    std::vector<long> ts;
    int N;
    long t;
    float dist;
    while (infile >> N >> dist >> t){
	Ns.push_back(N);
	ts.push_back(float(t)/1000);
    }
    plt::plot(Ns, ts, "rx");
    plt::ylabel("Time[ms]");
    plt::xlabel("N");
    //plt::xlim(0,300);
    //plt::ylim(0,500);
    //plt::axis("equal");
    plt::legend();
    plt::show();
    //plt::save(homepath + "/jlindaue/bezier.pdf");
    return 0;
}

