// #include <stdio.h>
#define WITHOUT_NUMPY

#include "gradient_planner_lib3.h"
#include <string>
#include <iostream>
#include <fstream>
#include "matplotlibcpp.h"
#include "timer.h"


namespace plt = matplotlibcpp;
using namespace std;
using namespace Eigen;

Vector4d x0_backup, x_backup, flex;
float var = 0.15;


void compute_trajectory_caller(params &p, std::vector<Vector4d> &my_wpts) {
    MatrixXd part1;
    std::cout << "DONE1\n";
    {
	Timer timer;
    	part1 = compute_trajectory(p, my_wpts);
    }
    std::cout << "DONE1\n";
    std::vector<double> x(part1.cols());
    std::vector<double> y(part1.cols());
    for (size_t i = 0; i < x.size(); i++) {
        x[i] = part1(0, i);
        y[i] = part1(1, i);
    }
    std::cout << "DONE1\n";
    
    for (size_t i = 0; i < 4; i++) {
        flex(i) = part1(i, x.size()-1);
    }
    std::cout << "DONE1\n";


    plt::named_plot("Gradient", x, y);
    //cout << part1;
}


int main(int argc, char **argv) {
    params p;
    Vector4d start_pt;
    Vector4d end_pt;
    MatrixXd ret;

    p.T = 0.05;
    p.v = 1;
    p.Wye = 0;
    p.Wyaw = 0;
    p.Wc = 1;
    p.Wdc = 0.1;
    p.WB = 0.32;
    p.wr = 0.5;

    std::vector <float> xs;
    std::vector <float> ys;
    std::vector <float> us;
    std::vector <float> vs;

    std::ifstream mystream("wpts.txt");
    std::vector <Vector4d> my_wpts;
    Vector4d vec;
    std::string line;
    while(std::getline(mystream,line)){
	std::stringstream ss(line);
        for (int i=0;i<4;i++){
	    std::string str;
	    std::getline(ss,str,',');
	    vec(i) = std::stod(str);
	}
    	xs.push_back(vec(0));    
    	ys.push_back(vec(1));    
    	us.push_back(0);    
    	vs.push_back(0);    
    	my_wpts.push_back(vec);
	std::cout << vec << "\n";
    }
    std::cout << "DONE1\n";
    //plt::named_plot("WPTS", xs, ys, "rx");

    us[0]=std::cos(my_wpts[0][2]);
    vs[0]=std::sin(my_wpts[0][2]);

    int last=my_wpts.size()-1;
    us[last]=std::cos(my_wpts[last][2]);
    vs[last]=std::sin(my_wpts[last][2]);
    plt::quiver(xs,ys,us,vs);


    compute_trajectory_caller(p, my_wpts);
    plt::axis("equal");
    plt::legend();
    plt::grid(true);
    plt::xlabel("x coordinate [m]");
    plt::ylabel("y coordinate [m]");
    plt::show();
    //plt::save("bezier.pdf");
    return 0;
}

