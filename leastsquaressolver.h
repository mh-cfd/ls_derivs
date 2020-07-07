#ifndef LEASTSQUARESSOLVER_H
#define LEASTSQUARESSOLVER_H
#include <vector>


//number of variables
#define VAR_NUM 10

//maximum number of points
#define MAX_EQNS 10000


typedef struct
{double x ,y,z,f,ax,ay,az;
int is_boundary; //if<0 than it's not the boundary oterwise it's boundary plane number
double f_bound;//value of the gradient at the boundary
double rhs;// value of the rhs for a poisson equation at the boundary

double INV[VAR_NUM][VAR_NUM];
} node3d;


typedef struct
{double d[VAR_NUM];} deriv3D;

typedef struct
{double nx,ny,nz;
 double d;//normal distance from  coordinate origin  i.e.  nx*x+ny*y+nz*z=d
} boundary_plane;





class leastSquaresSolver
{
private:
    double M_[MAX_EQNS][50],
           M_0[MAX_EQNS][50],
           MWM[MAX_EQNS][50],
           x_m[MAX_EQNS],
           b_m[MAX_EQNS],
           w[MAX_EQNS],
           mwb[MAX_EQNS],
           LU[MAX_EQNS][50],
           Inv[MAX_EQNS][50];
    int ps[MAX_EQNS];

public:
    enum deriv_order
    {
        F,
        FX,
        FY,
        FZ,
        FXX,
        FXY,
        FXZ,
        FYY,
        FYZ,
        FZZ
    };

    leastSquaresSolver();

    std::vector<node3d> m_p; //points cloud
    std::vector<deriv3D> m_d,m_d0; //derivs cloud
    std::vector<boundary_plane> walls; //all boundaries are here for now

    double pww;
    double rad;
    void init_with_points();
    void init_for_benchmark();
    void LU_decompose(void);
    void m_solve(void);
    void m_invert(void);

    void LU_decompose2(int num);
    void m_solve2(int num);
    void m_invert2(int num);

    void LU_decompose_krig(int num);
    void m_solve_krig(int num);
    void m_invert_krig(int num);

    void nullify_m(double m[MAX_EQNS][50]);
    void nullify_v(double v[MAX_EQNS]);
    void get_derivs(node3d &p, deriv3D &res, double delta);
    void get_derivs_fast(node3d &p, deriv3D &res, double delta);

            void interpKrig(node3d &p);
                 void interpKrigDx(node3d &p);
            double dist(node3d &p1,node3d &p2);
            double distDx(node3d &p1,node3d &p2);
                 double distLapl(node3d &p1, node3d &p2);

     void get_poisson_internal(node3d &p, deriv3D &res, double delta);
     void get_poisson_boundary(node3d &p, deriv3D &res, double delta);

     void get_poisson_combo(node3d &p, deriv3D &res, double delta);

     void get_poisson_matr(node3d &p, deriv3D &res, double delta);// only gets the inverse matirx
     void get_poisson_nomatr(node3d &p, deriv3D &res, double delta); //solves with given inverse matrix;
     void draw_points(double sc);
     void get_derivs_bench(node3d &p, deriv3D &res);
     void getKrigInv();
     void interpKrig_nomatr(node3d &p);

     void solvePoisson_krig();
     double distAx(node3d &p1, node3d &p2);
     double distAy(node3d &p1, node3d &p2);
     double distAz(node3d &p1, node3d &p2);
     void solveKrig_grad(int itn, double delta);
};

#endif // LEASTSQUARESSOLVER_H
