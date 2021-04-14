// #include <stdio.h>
#define WITHOUT_NUMPY

#include "gradient_planner_lib3.h"
#include <string>
#include <iostream>
#include "matplotlibcpp.h"
#include "timer.h"


namespace plt = matplotlibcpp;
using namespace std;
using namespace Eigen;

Vector4d x0_backup, x_backup, flex;
float var = 0.15;

struct WEIGHTS {
    float k1, k4, c1, c4, t1, t4, dif, dif2;
};

class Bezier {
    float side(Vector4d &a, Vector4d &b) {
        float mul = sin(a(2)) * (b(0) - a(0)) + cos(a(2)) * (b(1) - a(1));
        float curv = a(3);
        std::cout << mul << "\n";
        std::cout << curv << "\n";
        std::cout << mul * curv << "\n";
        if (mul * curv >= 0) return fabs(curv);
        return -fabs(curv);
    }

    float angle(Vector4d &a, Vector4d &b) {
        return atan2(b(1) - a(1), b(0) - a(0));
    }

    WEIGHTS w;
    Vector4d P2, P3;
    float alpha, beta;
    std::vector<double> x2;
    std::vector<double> y2;

    void comp_curve(Vector4d &start_point, Vector4d &end_point) {
        float d = sqrt(pow(start_point(0) - end_point(0), 2) + pow(start_point(1) - end_point(1), 2));
        int steps = d / 0.15;

        for (int i = 0; i <= steps; i++) {
            float t = float(i) / float(steps);
            float x = (1 - t) * (1 - t) * (1 - t) * start_point(0) + 3 * t * (1 - t) * (1 - t) * P2(0) +
                      3 * t * t * (1 - t) * P3(0) +
                      t * t * t * end_point(0);
            float y = (1 - t) * (1 - t) * (1 - t) * start_point(1) + 3 * t * (1 - t) * (1 - t) * P2(1) +
                      3 * t * t * (1 - t) * P3(1) +
                      t * t * t * end_point(1);
            x2.push_back(x);
            y2.push_back(y);

        }
    }

    void comp_ctrl_pts(Vector4d &start_point, Vector4d &end_point) {
        float d = sqrt(pow(start_point(0) - end_point(0), 2) + pow(start_point(1) - end_point(1), 2));
        float s1 = side(start_point, end_point);

        //in pi range
        float thetad = fabs(in_pi_range(alpha - beta));

        P2(0) = start_point(0) + cos(alpha) * d * max(0.1f, w.k1 - w.c1 * s1 + thetad * w.dif);
        P2(1) = start_point(1) + sin(alpha) * d * max(0.1f, w.k1 - w.c1 * s1 + thetad * w.dif);
        float s4 = -side(end_point, start_point);
        P3(0) = end_point(0) - cos(beta) * d * max(0.1f, w.k4 - w.c4 * s4 * s4 + thetad * w.dif2);
        P3(1) = end_point(1) - sin(beta) * d * max(0.1f, w.k4 - w.c4 * s4 * s4 + thetad * w.dif2);
    }

public:
    void bezier(Vector4d &start_point, Vector4d &end_point) {
        w = {0.5, 0.4, 0., 0, 0.0, 0.0, 0, -0};
        alpha = start_point(2) + start_point(3) * w.t1;
        beta = end_point(2) - end_point(3) * w.t4;
        comp_ctrl_pts(start_point, end_point);
        comp_curve(start_point, end_point);
        plt::named_plot("Bezier", x2, y2);
    }

    void bezier(Vector4d &start_point, Vector4d &end_point, Vector4d &flex) {
        w = {0.55, 0.4, 0.2, 0, 0.25, 0.25, 0., -0.};
        alpha = start_point(2) + start_point(3) * w.t1;
        beta = angle(start_point, end_point) - start_point(3) * w.t1 + end_point(3) * w.t4;
        beta += var * in_pi_range(alpha - beta);
        comp_ctrl_pts(start_point, flex);
        P3 = P2;
        comp_curve(start_point, flex);
        w = {0.4, 0.4, 0.01, 0, 0.2, 0.2, -0., 0.};
        alpha = beta;
        beta = end_point(2) - end_point(3) * w.t4;
        comp_ctrl_pts(flex, end_point);
        P2 = P3;
        comp_curve(flex, end_point);
        plt::named_plot("Bezier", x2, y2);
    }
    /*
    void bezier(Vector4d &start_point, Vector4d &end_point, Vector4d &flex, bool my_bool) {
        w = {0.55, 0.4, 0.2, 0, 0.25, 0.25, 0., -0.};
        alpha = start_point(2) + start_point(3) * w.t1;
        P2(0) = start_point(0) + cos(alpha) * d * max(0.1f, w.k1 - w.c1 * s1);
        P2(1) = start_point(1) + sin(alpha) * d * max(0.1f, w.k1 - w.c1 * s1);
        float s4 = -side(end_point, start_point);

        Vector4d PN;
        w = {0.4, 0.4, 0.01, 0, 0.2, 0.2, -0., 0.};
        PN(0) = end_point(0) - cos(beta) * d * max(0.1f, w.k4 - w.c4 * s4 * s4);
        PN(1) = end_point(1) - sin(beta) * d * max(0.1f, w.k4 - w.c4 * s4 * s4);

        float result1=in_pi_range(2*(in_pi_range(-angle(flex,start_point) + angle(flex,P2))));
        result1=in_pi_range(angle(flex,start_point)+result1);
        float result2=in_pi_range(2*(in_pi_range(-angle(flex,end_point) + angle(flex,PN))));
        result2=in_pi_range(angle(flex,end_point)+result2);

        float result=result1+result2


        beta = angle(start_point, end_point) - start_point(3) * w.t1 + end_point(3) * w.t4;
        beta += var * in_pi_range(alpha - beta);
        comp_ctrl_pts(start_point, flex);
        P3 = P2;
        comp_curve(start_point, flex);

        alpha = beta;
        comp_ctrl_pts(flex, end_point);
        P2 = P3;
        comp_curve(flex, end_point);
        plt::named_plot("Bezier", x2, y2);
    }*/
};

void compute_trajectory_caller(params &p, Vector4d &start_point, Vector4d &end_point) {
    MatrixXd part1;
    x0_backup = start_point;
    x_backup = end_point;
    std::vector <Vector4d> my_wpts;
    my_wpts.push_back(x0_backup);
    //my_wpts.push_back(flex);
    my_wpts.push_back(x_backup);
    cout << "\nWPT1:\n" << x0_backup; 
    cout << "\nWPT2:\n" << x_backup << "\n"; 
    //cout << x_backup; 
    //cout << p.v; 
    {
	Timer timer;
    	part1 = compute_trajectory(p, my_wpts);
    }
    std::vector<double> x(part1.cols());
    std::vector<double> y(part1.cols());
    for (size_t i = 0; i < x.size(); i++) {
        x[i] = part1(0, i);
        y[i] = part1(1, i);
    }
    
    for (size_t i = 0; i < 4; i++) {
        flex(i) = part1(i, x.size()-1);
    }


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
    start_pt(0) = stod(argv[1]);
    start_pt(1) = stod(argv[2]);
    start_pt(2) = stod(argv[3]);
    start_pt(3) = stod(argv[4]);
    end_pt(0) = stod(argv[5]);
    end_pt(1) = stod(argv[6]);
    end_pt(2) = stod(argv[7]);
    end_pt(3) = stod(argv[8]);
    flex(0) = stod(argv[9]);
    flex(1) = stod(argv[10]);
    flex(2) = stod(argv[11]);
    flex(3) = stod(argv[12]);
    //p.v = stod(argv[13])/p.T;
    std::cout << "Set to: " << p.v*p.T << "\n";

    flex(2)=in_pi_range((start_pt(2)+end_pt(2))/2);

    Bezier b;
    b.bezier(start_pt,flex);

    Bezier b2;
    b2.bezier(flex, end_pt);


    compute_trajectory_caller(p, start_pt, flex);
    compute_trajectory_caller(p, flex, end_pt);
    plt::axis("equal");
    plt::legend();
    plt::show();
    //plt::save("bezier.pdf");
    return 0;
}

