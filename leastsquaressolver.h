#ifndef LEASTSQUARESSOLVER_H
#define LEASTSQUARESSOLVER_H
#include <vector>


//number of variables
#define VAR_NUM 10

//maximum number of points
#define MAX_EQNS 100

typedef struct
{double x ,y,z,f;} node3d;


typedef struct
{double d[VAR_NUM];} deriv3D;



class leastSquaresSolver
{
private:
    double M_[MAX_EQNS][MAX_EQNS],
           M_0[MAX_EQNS][MAX_EQNS],
           MWM[MAX_EQNS][MAX_EQNS],
           x_m[MAX_EQNS],
           b_m[MAX_EQNS],
           w[MAX_EQNS],
           mwb[MAX_EQNS],
           LU[MAX_EQNS][MAX_EQNS],
           Inv[MAX_EQNS][MAX_EQNS];
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

    void init_with_points();
    void init_for_benchmark();
    void LU_decompose(void);
    void m_solve(void);
    void m_invert(void);

    void nullify_m(double m[MAX_EQNS][MAX_EQNS]);
    void nullify_v(double v[VAR_NUM]);
    void get_derivs(node3d &p, deriv3D &res, double delta);

    void draw_points(double sc);
    void get_derivs_bench(node3d &p, deriv3D &res);
};

#endif // LEASTSQUARESSOLVER_H
