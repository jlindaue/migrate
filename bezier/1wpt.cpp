// #include <stdio.h>
#define WITHOUT_NUMPY
#include "gradient_planner_lib3.h"
#include <string>
#include <iostream>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;
using namespace Eigen;

Vector4d x0_backup, x_backup, flex;

float k1=0.4;
float k4=0;
float c1=0.4;
float c4=0;
float t1=0.2;
float t4=0;

float side(Vector4d& a, Vector4d& b){
	float mul=sin(a(2))*(b(0)-a(0))+cos(a(2))*(b(1)-a(1));
	float curv=a(3);
	std::cout << mul << "\n";
	std::cout << curv << "\n";
	std::cout << mul*curv << "\n";
	if (mul*curv>=0) return fabs(curv);
	return -fabs(curv);
}

void bezier(Vector4d& start_point, Vector4d& end_point){
    Vector4d P2, P3;
    float d=sqrt(pow(start_point(0)-end_point(0),2)+pow(start_point(1)-end_point(1),2));
    float s1=side(start_point,end_point);
    P2(0)=start_point(0)+cos(start_point(2)+start_point(3)*t1)*d*(k1-c1*s1);
    P2(1)=start_point(1)+sin(start_point(2)+start_point(3)*t1)*d*(k1-c1*s1);
    float s4=-side(end_point,start_point);
    P3(0)=end_point(0)-cos(end_point(2)-end_point(3)*t4)*d*(k4-c4*s4);
    P3(1)=end_point(1)-sin(end_point(2)-end_point(3)*t4)*d*(k4-c4*s4);

    int steps=d/0.15;

    std::vector<double> x2(steps+1);
    std::vector<double> y2(steps+1);

    for(int i=0;i<=steps;i++) {
        float t=float(i)/float(steps);
        x2[i] = (1 - t) * (1 - t) * (1 - t) * start_point(0) + 3 * t * (1 - t) * (1 - t) * P2(0) + 3 * t * t * (1 - t) * P3(0) +
                       t * t * t * end_point(0);
        y2[i] = (1 - t) * (1 - t) * (1 - t) * start_point(1) + 3 * t * (1 - t) * (1 - t) * P2(1) + 3 * t * t * (1 - t) * P3(1) +
                       t * t * t * end_point(1);

    }

    plt::named_plot("Bezier", x2,y2);
}


void compute_trajectory_caller(params& p, Vector4d& start_point, Vector4d& end_point){
    MatrixXd part1;
    x0_backup = start_point;
    x_backup = end_point;
    std::vector<Vector4d> my_wpts;
    my_wpts.push_back(x0_backup);
    my_wpts.push_back(flex);
    my_wpts.push_back(x_backup);
    //cout << x2[steps]; 
    //cout << x_backup; 
    //cout << p.v; 
    part1 = compute_trajectory(p, my_wpts);

    std::vector<double> x(part1.cols());
    std::vector<double> y(part1.cols());
    for(size_t i = 0; i < x.size(); i++){
	x[i]=part1(0, i);
	y[i]=part1(1, i);
    }

    plt::named_plot("Gradient", x,y);
    //cout << part1;
}

int main (int argc, char** argv){
        params p;
        Vector4d start_pt;
        Vector4d end_pt;
        MatrixXd ret;

        p.T=0.1;
        p.v=2.5;
        c1=stod(argv[1]);
        c4=stod(argv[2]);
        k1=stod(argv[3]);
        k4=stod(argv[4]);
	p.Wye=0;
    	p.Wyaw=0;
    	p.Wc=1;
    	p.Wdc=0.1;
        start_pt(0)=stod(argv[5]);
        start_pt(1)=stod(argv[6]);
        start_pt(2)=stod(argv[7]);
        start_pt(3)=stod(argv[8]);
        end_pt(0)=stod(argv[9]);
        end_pt(1)=stod(argv[10]);
        end_pt(2)=stod(argv[11]);
        end_pt(3)=stod(argv[12]);
        flex(0)=stod(argv[13]);
        flex(1)=stod(argv[14]);
        flex(2)=end_pt(2);
        flex(3)=end_pt(3);
	
        bezier(start_pt, flex);
        compute_trajectory_caller(p, start_pt, end_pt);
	plt::axis("equal");
	plt::legend();
        plt::show();
	//plt::save("bezier.pdf");
        return 0;
}

