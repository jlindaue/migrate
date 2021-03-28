#include "gradient_planner_lib3.h"

MatrixXd compute_trajectory(params& p, std::vector<Vector4d>& wpts)
{
    std::vector<Vector4d> wpts_t;
    transform_to_car_frame(wpts, wpts_t);
    Vector3d x0 = wpts_t[0].tail(3);
    Vector3d xf = wpts_t.back().tail(3);
    
    // get inner waypoints
    std::vector<double> Yoff;
    for (size_t i = 1; i < wpts_t.size()-1; i++){
        Yoff.push_back(wpts_t[i](1));
    }
    int Yoff_cnt = Yoff.size();
    
    // get number of samples
    std::vector<double> d_init;
    std::vector<int> Ns_init;
    std::vector<int> Ns;
    int N = 0;
    Ns.push_back(0);
    for (size_t i = 0; i < wpts_t.size()-1; i++){
        double d = (wpts_t.at(i).head(2) - wpts_t.at(i+1).head(2)).norm();
        int n = round(d / p.v / p.T);
        N += n;
        d_init.push_back(d);
        Ns_init.push_back(n);
        Ns.push_back(N);
    }
    
    
    // 3*N states + 1*N inputs + 3*(N+1) constraints
    VectorXd X = VectorXd::Zero(4*N + 3*(N+1) + Yoff_cnt); // states, inputs, costates
    MatrixXd Xx(4,N+1);
    VectorXd u;
    
    // compute optimal trajectory
    for(int j=0;j<2;j++){
        // recompute N
        if(j == 1){
            std::vector<int> Ns_new;
            Ns_new.push_back(0);
            double N_new = 0;
            for (size_t i = 0; i < Ns.size()-1; i++){
                double d = (Xx.col(Ns.at(i)).head(2) - Xx.col(Ns.at(i+1)).head(2)).norm();
                int n = round(Ns_init.at(i) * d_init.at(i) / d);
                N_new += n;
                Ns_new.push_back(N_new);
            }
            N = N_new;
            Ns = Ns_new;
            X = VectorXd::Zero(4*N + 3*(N+1) + Yoff_cnt);
            Xx.resize(4,N+1);
        }
        
        SparseMatrix<double> A( 4*N + 3*(N+1) + Yoff_cnt, 4*N + 3*(N+1) + Yoff_cnt);
        std::vector<T> A_t;
        A_t.resize(N+2*(10*N+6)); // TODO size
        VectorXd b = VectorXd::Zero(4*N + 3*(N+1) + Yoff_cnt);
        
        // create Bk, gf
        for(int i=0;i<N;i++){
            // shaped input cost function
            // just diagonal matrix, directly filled in the matrix A
            double Wye = p.Wye * (0.0001 + fabs(X(3*i)));
            double Wyaw = p.Wyaw * (0.0001 + fabs(X(3*i+1)));
            double Wc = p.Wc * (0.0001 + fabs(X(3*i+2)));
            double Wdc = p.Wdc * (0.0001 + fabs(X(3*N+i)));
            A_t.push_back( T( 3*i, 3*i, Wye));
            A_t.push_back( T( 3*i+1, 3*i+1, Wyaw) );
            A_t.push_back( T( 3*i+2, 3*i+2, Wc) );
            A_t.push_back( T( 3*N+i, 3*N+i, Wdc) );
            // gradient of cost functions are just cost weights in the b vector
            b(3*i) = X(3*i) * p.Wye;
            b(3*i+1) = X(3*i+1) * p.Wyaw;
            b(3*i+2) = X(3*i+2) * p.Wc;
            b(3*N+i) = X(3*N+i) * p.Wdc;
        }
        
        G(X, x0, xf, Yoff, Ns, p, b);    
        grad_G(X, Yoff, Ns, p, A_t);
        A.setFromTriplets(A_t.begin(), A_t.end());
        A.makeCompressed();
        
        // solve
        // call some magic solver from eigen lib and hope that it will works
        SparseLU<SparseMatrix<double> > solver;
        solver.compute(A);
        if(solver.info()!=Success) {
            cout<<"decomposition failed"<<endl;
        }
        VectorXd sol = solver.solve(b);
        if(solver.info()!=Success) {
            cout<<"solving failed"<<endl;
        }
        
        X = X - sol;

        // recompute states
        u = X.segment(3*N,N);
        // by starting in initial state without transformation,
        // we reconstruct trajectory in the original coord. frame
        Xx.col(0) = wpts.at(0);
        
        for(int i=0;i<N;i++){
            Vector4d tmp = Xx.col(i);
            Xx.col(i+1) = F(tmp, u(i), p);
        }
    }
    
    VectorXd t = VectorXd::LinSpaced(N+1,0,N*p.T);
    
    // fill in the output matrix
    MatrixXd ret(6,N+1);
    ret << Xx,
           u.transpose(), 0,
           t.transpose();
    
    //transform_to_input_frame(ret, init_pt);
    
    return ret;
}

// state equation x[k+1] = f(x[k],u[k])
Vector4d F(Vector4d& X, double u, params& p)
{
    Vector4d ret;
    double x = X(0);
    double y = X(1);
    double yaw = X(2);
    double curv = X(3);
    
    ret(0) = p.v * cos(yaw + curv*p.WB*p.wr) * p.T + x;
    ret(1) = p.v * sin(yaw + curv*p.WB*p.wr) * p.T + y;
    ret(2) = p.v * curv * p.T + yaw;
    ret(3) = u * p.T + curv;
    
    return ret;
}

// constraints g(X) = 0
void G(VectorXd& X, Vector3d& x0, Vector3d& xf, std::vector<double> Yoff, std::vector<int> Ns, params& p, VectorXd& b)
{
    int N = Ns.back();
    const int start = 4*N;
    // reduced state space eq.
    b[4*N] = p.v * sin( x0(1)+x0(2)*p.WB*p.wr ) * p.T + x0(0) - X(0);
    b[4*N+1] = p.v * x0(2) * p.T + x0(1) - X(1);
    b[start+2] = X(3*N) * p.T + x0(2) - X(2);
    
    for(int i=1;i<N-1;i++){
        b[4*N+3*i] = p.v * sin( X(3*(i-1)+1)+X(3*(i-1)+2)*p.WB*p.wr ) * p.T + X(3*(i-1)) - X(3*i);
        b[4*N+3*i+1] = p.v * X(3*(i-1)+2) * p.T + X(3*(i-1)+1) - X(3*i+1);
        b[4*N+3*i+2] = X(3*N+i) * p.T + X(3*(i-1)+2) - X(3*i+2);
    }
    
    b[4*N+3*N] = X(3*N-3)-xf(0);
    b[4*N+3*N+1] = X(3*N-2)-xf(1);
    b[4*N+3*N+2] = X(3*N-1)-xf(2);
    
    // y offset constraints
    for (size_t i = 0; i < Yoff.size(); i++){
        b[4*N+3*N+2 + i+1] = Yoff[i] - X[3*(Ns[i+1]-1)+1];
    }
}

// gradient of constraints grad(g(x))
void grad_G(VectorXd& X, std::vector<double> Yoff, std::vector<int> Ns, params& p, std::vector<T>& A_t)
{
    const int N = Ns.back();
    const int start = 4*N;
    A_t.push_back( T(start+0,0, -1));
    A_t.push_back( T(start+1,1, -1));
    A_t.push_back( T(start+2,2, -1));
    A_t.push_back( T(start+2,3*N, p.T));
    
    // gg'
    A_t.push_back( T(0,start+0, -1));
    A_t.push_back( T(1,start+1, -1));
    A_t.push_back( T(2,start+2, -1));
    A_t.push_back( T(3*N,start+2, p.T));
    
    for(int i=0;i<N-1;i++){
        double arg = X(3*i+1)+p.WB*p.wr*X(3*i+2);
        A_t.push_back( T(start+ 3*(i+1), 3*i, 1));
        A_t.push_back( T(start+ 3*(i+1), 3*i+1, p.v*p.T*cos(arg) ));
        A_t.push_back( T(start+ 3*(i+1), 3*i+2, p.v*p.T*p.WB*p.wr*cos(arg) ));
        A_t.push_back( T(start+ 3*(i+1), 3*i+3, -1 ));
        
        A_t.push_back( T(start+ 3*(i+1)+1, 3*i+1, 1 ));
        A_t.push_back( T(start+ 3*(i+1)+1, 3*i+2, p.v*p.T ));
        A_t.push_back( T(start+ 3*(i+1)+1, 3*i+4, -1 ));
        
        A_t.push_back( T(start+ 3*(i+1)+2, 3*i+2, 1 ));
        A_t.push_back( T(start+ 3*(i+1)+2, 3*i+5, -1 ));
        
        A_t.push_back( T(start+ 3*(i+1)+2, 3*N+i+1, p.T ));
        
        // gg'
        A_t.push_back( T(3*i, start+ 3*(i+1), 1));
        A_t.push_back( T(3*i+1, start+ 3*(i+1), p.v*p.T*cos(arg) ));
        A_t.push_back( T(3*i+2, start+ 3*(i+1), p.v*p.T*p.WB*p.wr*cos(arg) ));
        A_t.push_back( T(3*i+3, start+ 3*(i+1), -1 ));
        
        A_t.push_back( T(3*i+1, start+ 3*(i+1)+1, 1 ));
        A_t.push_back( T(3*i+2, start+ 3*(i+1)+1, p.v*p.T ));
        A_t.push_back( T(3*i+4, start+ 3*(i+1)+1, -1 ));
        
        A_t.push_back( T(3*i+2, start+ 3*(i+1)+2, 1 ));
        A_t.push_back( T(3*i+5, start+ 3*(i+1)+2, -1 ));
        
        A_t.push_back( T(3*N+i+1, start+ 3*(i+1)+2, p.T ));
        
    }
    // final state
    A_t.push_back( T(start+ 3*N, 3*(N-1), 1));
    A_t.push_back( T(start+ 3*N+1, 3*(N-1)+1, 1));
    A_t.push_back( T(start+ 3*N+2, 3*(N-1)+2, 1));
    
    // gg'
    A_t.push_back( T(3*(N-1), start+ 3*N, 1));
    A_t.push_back( T(3*(N-1)+1, start+ 3*N+1, 1));
    A_t.push_back( T(3*(N-1)+2, start+ 3*N+2, 1));
    
    // y offset constraints
    for (size_t i = 0; i < Yoff.size(); i++){
        A_t.push_back( T(start+ 3*(N+1)+i, 3*(Ns[i+1]-1), -1));
        // gg'
        A_t.push_back( T( 3*(Ns[i+1]-1), start+ 3*(N+1)+i, -1));
    }
}

// phi in range <-PI,PI>
double in_pi_range(double ang)
{
    double ret = ang;
    if(ang > M_PI) ret -= 2*M_PI;
    if(ang < -M_PI) ret += 2*M_PI;
    return ret;
}

void transform_to_car_frame(std::vector<Vector4d>& wpts, std::vector<Vector4d>& wpts_t)
{
    std::vector<Vector4d> ret;
    auto x0 = -wpts.at(0);
    double ang = -in_pi_range(atan2(wpts.back()(1) - wpts[0](1), wpts.back()(0) - wpts[0](0)));
    
    Matrix2d rot;
    rot << cos(ang), -sin(ang), sin(ang), cos(ang);
    
    for(size_t i = 0; i < wpts.size(); i++) {
        Vector4d v;
        v.head(2) = rot * (wpts[i].head(2) + x0.head(2));
        v(2) = in_pi_range(wpts[i](2) + ang);
        v(3) = wpts[i](3);
        wpts_t.push_back(v);
    }
}
