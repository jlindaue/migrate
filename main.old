// #include <stdio.h>
#include "gradient_planner_lib.h"
#include <string>
#include <iostream>


using namespace std;
using namespace Eigen;


int main (int argc, char** argv){
        params p;
        Vector4d start_pt;
        Vector4d end_pt;
        MatrixXd ret;

        p.T=stod(argv[1]);
        p.v=stod(argv[2]);
        p.WB=stod(argv[3]);
        p.wr=stod(argv[4]);
        start_pt(0)=stod(argv[5]);
        start_pt(1)=stod(argv[6]);
        start_pt(2)=stod(argv[7]);
        start_pt(3)=stod(argv[8]);
        end_pt(0)=stod(argv[9]);
        end_pt(1)=stod(argv[10]);
        end_pt(2)=stod(argv[11]);
        end_pt(3)=stod(argv[12]);

        ret=compute_trajectory(p, start_pt, end_pt);
        cout << ret;
        return 0;
}

