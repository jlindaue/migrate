#include "gradient_planner_lib.h"
#include <string>

// return clothoid trajectory, parameter, starting point, ending point
MatrixXd compute_trajectory(params& p, Vector4d& start_pt, Vector4d& end_pt)
{
    Vector4d init_pt = start_pt;
    transform_to_car_frame(start_pt, end_pt);
    Vector3d x0 = start_pt.tail(3);
    Vector3d xf = end_pt.tail(3);
    
    SparseMatrix<double> A( 4*p.N + 3*(p.N+1), 4*p.N + 3*(p.N+1) );
        std::vector<T> A_t;
        //         Bk   gg+gg'
        A_t.resize(p.N+2*(10*p.N+6));
    VectorXd b = VectorXd::Zero(4*p.N + 3*(p.N+1));
    
    // 3*N states + 1*N inputs + 3*(N+1) constraints
    VectorXd X = VectorXd::Zero(4*p.N + 3*(p.N+1)); // states, inputs, costates
    
    // compute optimal trajectory
    for(int i=0;i<2;i++){
        
        // create Bk, gf
        for(int i=0;i<p.N;i++){
            // shaped input cost function
            // just diagonal matrix, directli filled in the matrix A
            A_t[i] = T( 3*p.N+i, 3*p.N+i, 0.01 + fabs(X(3*p.N+i)) );
            // gradient of cost functions are just cost weights in the b vector
            b(3*p.N+i) = X(3*p.N+i) * (0.01 + fabs(X(3*p.N+i)));
        }
        
        G(X, x0, xf, p, b);    
        grad_G(X, p, A_t);
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
    }
    
    
    //cout<<X<<endl;
    //exit(0);
    
    // recompute states
    VectorXd u = X.segment(3*p.N,p.N);
    MatrixXd Xx(4,p.N+1);
    // by starting in initial state without transformation,
    // we reconstruct trajectory in the original coord. frame
    Xx.col(0) = init_pt;//start_pt;
    
    for(int i=0;i<p.N;i++){
        Vector4d tmp = Xx.col(i);
        Xx.col(i+1) = F(tmp, u(i), p);
    }
    
    VectorXd t = VectorXd::LinSpaced(p.N+1,0,p.N*p.T);
    
    // fill in the output matrix
    MatrixXd ret(6,p.N+1);
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
void G(VectorXd& X, Vector3d& x0, Vector3d& xf, params& p, VectorXd& b)
{
    const int N = p.N;
    const int start = 4*N;
    // reduced state space eq.
    b[start+0] = p.v * sin( x0(1)+x0(2)*p.WB*p.wr ) * p.T + x0(0) - X(0);
    b[start+1] = p.v * x0(2) * p.T + x0(1) - X(1);
    b[start+2] = X(3*p.N) * p.T + x0(2) - X(2);
    
    for(int i=1;i<p.N-1;i++){
        b[start+3*i] = p.v * sin( X(3*(i-1)+1)+X(3*(i-1)+2)*p.WB*p.wr ) * p.T + X(3*(i-1)) - X(3*i);
        b[start+3*i+1] = p.v * X(3*(i-1)+2) * p.T + X(3*(i-1)+1) - X(3*i+1);
        b[start+3*i+2] = X(3*N+i) * p.T + X(3*(i-1)+2) - X(3*i+2);
    }
    
    int size = 3*(N+1);
    b[start+3*N] = X(3*p.N-3)-xf(0);
    b[start+3*N+1] = X(3*p.N-2)-xf(1);
    b[start+3*N+2] = X(3*p.N-1)-xf(2);
}

// gradient of constraints grad(g(x))
void grad_G(VectorXd& X, params& p, std::vector<T>& A_t)
{
    const int N = p.N;
    const int start = 4*N;
    A_t[N+0] = T(start+0,0, -1);
    A_t[N+1] = T(start+1,1, -1);
    A_t[N+2] = T(start+2,2, -1);
    A_t[N+3] = T(start+2,3*N, p.T);
    
    // gg'
    A_t[N+4] = T(0,start+0, -1);
    A_t[N+5] = T(1,start+1, -1);
    A_t[N+6] = T(2,start+2, -1);
    A_t[N+7] = T(2,start+3*N, p.T);
    
    for(int i=0;i<p.N-1;i++){
        double arg = X(3*i+1)+p.WB*p.wr*X(3*i+2);
        A_t[N+7+20*i+1] = T(start+ 3*(i+1), 3*i, 1);
        A_t[N+7+20*i+2] = T(start+ 3*(i+1), 3*i+1, p.v*p.T*cos(arg) );
        A_t[N+7+20*i+3] = T(start+ 3*(i+1), 3*i+2, p.v*p.T*p.WB*p.wr*cos(arg) );
        A_t[N+7+20*i+4] = T(start+ 3*(i+1), 3*i+3, -1 );
        
        A_t[N+7+20*i+5] = T(start+ 3*(i+1)+1, 3*i+1, 1 );
        A_t[N+7+20*i+6] = T(start+ 3*(i+1)+1, 3*i+2, p.v*p.T );
        A_t[N+7+20*i+7] = T(start+ 3*(i+1)+1, 3*i+4, -1 );
        
        A_t[N+7+20*i+8] = T(start+ 3*(i+1)+2, 3*i+2, 1 );
        A_t[N+7+20*i+9] = T(start+ 3*(i+1)+2, 3*i+5, -1 );
        
        A_t[N+7+20*i+10] = T(start+ 3*(i+1)+2, 3*p.N+i+1, p.T );
        
        // gg'
        A_t[N+7+20*i+11] = T(3*i, start+ 3*(i+1), 1);
        A_t[N+7+20*i+12] = T(3*i+1, start+ 3*(i+1), p.v*p.T*cos(arg) );
        A_t[N+7+20*i+13] = T(3*i+2, start+ 3*(i+1), p.v*p.T*p.WB*p.wr*cos(arg) );
        A_t[N+7+20*i+14] = T(3*i+3, start+ 3*(i+1), -1 );
        
        A_t[N+7+20*i+15] = T(3*i+1, start+ 3*(i+1)+1, 1 );
        A_t[N+7+20*i+16] = T(3*i+2, start+ 3*(i+1)+1, p.v*p.T );
        A_t[N+7+20*i+17] = T(3*i+4, start+ 3*(i+1)+1, -1 );
        
        A_t[N+7+20*i+18] = T(3*i+2, start+ 3*(i+1)+2, 1 );
        A_t[N+7+20*i+19] = T(3*i+5, start+ 3*(i+1)+2, -1 );
        
        A_t[N+7+20*i+20] = T(3*N+i+1, start+ 3*(i+1)+2, p.T );
        
        //gg.block(3*(i+1),3*i,3,6) <<
            //1, p.v*p.T*cos(arg), p.v*p.T*p.WB*p.wr*cos(arg), -1,  0,  0,
            //0, 1,                p.v*p.T,                  0, -1,  0,
            //0, 0,                1,                        0,  0, -1;
        //gg(3*(i+1)+2,3*p.N+i+1) = p.T;
    }
    // final state
    A_t[21*N-12] = T(start+ 3*N, 3*(N-1), 1);
    A_t[21*N-11] = T(start+ 3*N+1, 3*(N-1)+1, 1);
    A_t[21*N-10] = T(start+ 3*N+2, 3*(N-1)+2, 1);
    
    // gg'
    A_t[21*N-9] = T(3*(N-1), start+ 3*N, 1);
    A_t[21*N-8] = T(3*(N-1)+1, start+ 3*N+1, 1);
    A_t[21*N-7] = T(3*(N-1)+2, start+ 3*N+2, 1);
    
    //gg.block(3*p.N,3*(p.N-1),3,3) << 1, 0, 0,
                                     //0, 1, 0,
                                     //0, 0, 1;
}

// phi in range <-PI,PI>
void in_pi_range(double *phi)
{
    if(*phi > M_PI) *phi -= 2*M_PI;
    if(*phi < -M_PI) *phi += 2*M_PI;
}

void transform_to_car_frame(Vector4d& x0, Vector4d& xf)
{
    double phi = x0(2);
    Vector2d DX = xf.head(2) - x0.head(2);
    
    Matrix2d rot;
    rot << cos(phi), sin(phi), -sin(phi), cos(phi);
    Vector2d Xrot = rot * DX;
    double yaw_rot = xf(2) - phi;
    in_pi_range(&yaw_rot);
    
    x0 << 0, 0, 0, x0(3);
    xf << Xrot(0), Xrot(1), yaw_rot, xf(3);
}

//void transform_to_input_frame(MatrixXd& traj, Vector4d& xf)
//{
    //double phi = x0(2);
    //Vector2d DX = xf.head(2) - x0.head(2);
    
    //Matrix2d rot;
    //rot << cos(phi), sin(phi), -sin(phi), cos(phi);
    //Vector2d Xrot = rot * DX;
    //double yaw_rot = xf(2) - phi;
    //in_pi_range(&yaw_rot);
    
    //x0 << 0, 0, 0, x0(3);
    //xf << Xrot(0), Xrot(1), yaw_rot, xf(3);
//}

int main (int argc, char** argv){
	params p;
	Vector4d start_pt;
	Vector4d end_pt;

    	p.N=stoi(argv[1]);
	p.T=stod(argv[2]);
	p.v=stod(argv[3]);
    	p.WB=stod(argv[4]);
	p.wr=stod(argv[5]);
	start_pt(0)=stod(argv[6]);
    	start_pt(1)=stod(argv[7]);
	start_pt(2)=stod(argv[8]);
	start_pt(3)=stod(argv[9]);
    	end_pt(0)=stod(argv[10]);
	end_pt(1)=stod(argv[11]);
	end_pt(2)=stod(argv[12]);
	end_pt(3)=stod(argv[13]);

	MatrixXd ret=compute_trajectory(p, start_pt, end_pt);
	cout << ret;
	// MatrixXd ret=compute_trajectory(p, start_pt, end_pt);
	return 0; 
}

