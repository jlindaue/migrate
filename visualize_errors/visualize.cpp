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
    std::string filename = argv[1];
    int plot_col = std::stod(argv[2]);;
    std::ifstream infile(homepath + "/jlindaue/organizovana_data/" + filename);

    std::vector<float> d_es;
    std::vector<float> yaw_es;
    std::vector<float> speeds;
    std::vector<float> steers;
    std::vector<float> accs;
    std::vector<float> speedrs;
    std::vector<float> steerrs;
    std::vector<long> ts;
    std::vector<int> is;
    float d_e;
    float yaw_e;
    float speed;
    float steer;
    float acc;
    float speedr;
    float steerr;
    long t;
    int i;
    while (infile >> d_e >> yaw_e >> speed >> steer >> acc >> speedr >> steerr >> t >> i){
	d_es.push_back(d_e);
	yaw_es.push_back(yaw_e);
	speeds.push_back(speed);
	steers.push_back(steer);
	accs.push_back(acc);
	speedrs.push_back(speedr);
	steerrs.push_back(steerr);
	ts.push_back(float(t)/1000);
	is.push_back(i);
    }
    if (plot_col==1) plt::plot(is,d_es);
    if (plot_col==2) plt::plot(is,yaw_es);
    if (plot_col==3) plt::plot(is,speeds);
    if (plot_col==4) plt::plot(is,steers);
    if (plot_col==5) plt::plot(is,accs);
    if (plot_col==6) plt::plot(is,speedrs);
    if (plot_col==7) plt::plot(is,steerrs);
    plt::ylabel("Time[ms]");
    //plt::xlabel("N");
    //plt::xlim(0,300);
    //plt::ylim(0,500);
    //plt::axis("equal");
    plt::legend();
    plt::show();
    //plt::save(homepath + "/jlindaue/bezier.pdf");
    return 0;
}

