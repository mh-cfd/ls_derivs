#include "leastsquaressolver.h"
#include <math.h>
#include <stdio.h>
//#include <my_include/gl.h>
//#include <my_include/glu.h>
//#include <my_include/glut.h>
#include <GL/glut.h>


double LU_krig[3000][3000],M_krig[3000][3000],Inv_krig[3000][3000],M0_krig[3000][3000],Wx_krig[3000][3000],Wy_krig[3000][3000],Wz_krig[3000][3000];
int ps_krig[3000];

leastSquaresSolver::leastSquaresSolver()
{
    pww=0.5;
    rad=0.125;
}



void leastSquaresSolver::init_for_benchmark()
{
    node3d n;
    deriv3D d,d0;
    n.x=0.0;
    n.y=0.0;
    n.z=0.0;
    get_derivs_bench(n,d0);
    n.f=d0.d[F];

    m_p.clear();
    m_d0.clear();
    m_d.clear();



    int i=0;
    m_p.push_back(n);
    m_d0.push_back(d0);
    m_d.push_back(d);

    double sc=0.5;
    for ( i=1;i<39;i++)
    {
        n.x=sc*(rand()*1.0/RAND_MAX-0.5);
        n.y=sc*(rand()*1.0/RAND_MAX-0.5);
        n.z=sc*(rand()*1.0/RAND_MAX-0.5);
        get_derivs_bench(n,d0);
        n.f=d0.d[F];
        m_p.push_back(n);
        m_d0.push_back(d0);
        m_d.push_back(d);
    }


    for (int i=0;i<39;i++)
    {
        get_derivs(m_p[i],m_d[i],0.5);
        //printf("i=%d x=%f y=%f z=%f f=%f fx=%f fy=%f fz=%f fnew=%f \n",
        //      i,m_p[i].x,m_p[i].y,m_p[i].z,m_p[i].f,m_d[i].d[FX],m_d[i].d[FY],m_d[i].d[FZ],m_d[i].d[F]);
    }
}

void leastSquaresSolver::LU_decompose(void)
{
    int i,j,k,pivotindex;
    double scales[50];
    double normrow,pivot,size,biggest,mult;

    for (i=0;i<VAR_NUM;i++) //заполнение начальными данными
    {
        ps[i]=i;//маппинг изначального порядка на переставленный.
        normrow=0;//максимум в итой строке

        for (j=0;j<VAR_NUM;j++)
        {
            LU[i][j]=M_[i][j];
            if (normrow<fabs(LU[i][j]))
                normrow=fabs(LU[i][j]);
        }
        if (normrow!=0)
            scales[i]=1.0/normrow;//для общих множителей
        else
        {
            scales[i]=0.0;
            //     err_code(DIV_ZERO);
        }
    }
    //метод гаусса с частичным упорядочиванием

    for (k=0;k<VAR_NUM-1;k++)
    {
        biggest=0;
        for (i=k; i<VAR_NUM;i++)
        {
            size=fabs(LU[ps[i]][k])*scales[ps[i]];
            if (biggest<size)
            {
                biggest=size;
                pivotindex=i;
            }
        }

        if (biggest==0)
        {
            //	err_code(1);
            pivotindex=0;
        }

        if (pivotindex!=k)
        {
            j=ps[k];
            ps[k]=ps[pivotindex];
            ps[pivotindex]=j;
        }

        pivot=LU[ps[k]][k];

        for (i=k+1;i<VAR_NUM;i++)
        {
            mult=LU[ps[i]][k]/pivot;
            LU[ps[i]][k]=mult;

            if (mult!=0.0)
            {
                for (j=k+1; j<VAR_NUM;j++)
                    LU[ps[i]][j]-=mult*LU[ps[k]][j];
            }
        }
    }
    //      if (LU[ps[VAR_NUM-1]][VAR_NUM-1]==0.0) err_code(1);
}


void leastSquaresSolver::LU_decompose2(int num)
{
    int i,j,k,pivotindex;
    double scales[50];
    double normrow,pivot,size,biggest,mult;

    for (i=0;i<num;i++) //заполнение начальными данными
    {
        ps[i]=i;//маппинг изначального порядка на переставленный.
        normrow=0;//максимум в итой строке

        for (j=0;j<num;j++)
        {
            LU[i][j]=M_[i][j];
            if (normrow<fabs(LU[i][j]))
                normrow=fabs(LU[i][j]);
        }
        if (normrow!=0)
            scales[i]=1.0/normrow;//для общих множителей
        else
        {
            scales[i]=0.0;
            //     err_code(DIV_ZERO);
        }
    }
    //метод гаусса с частичным упорядочиванием

    for (k=0;k<num-1;k++)
    {
        biggest=0;
        for (i=k; i<num;i++)
        {
            size=fabs(LU[ps[i]][k])*scales[ps[i]];
            if (biggest<size)
            {
                biggest=size;
                pivotindex=i;
            }
        }

        if (biggest==0)
        {
            //	err_code(1);
            pivotindex=0;
        }

        if (pivotindex!=k)
        {
            j=ps[k];
            ps[k]=ps[pivotindex];
            ps[pivotindex]=j;
        }

        pivot=LU[ps[k]][k];

        for (i=k+1;i<num;i++)
        {
            mult=LU[ps[i]][k]/pivot;
            LU[ps[i]][k]=mult;

            if (mult!=0.0)
            {
                for (j=k+1; j<num;j++)
                    LU[ps[i]][j]-=mult*LU[ps[k]][j];
            }
        }
    }
    //      if (LU[ps[VAR_NUM-1]][VAR_NUM-1]==0.0) err_code(1);
}

void leastSquaresSolver::m_solve2(int num)
{
    int i,j;
    double dot;

    for (i=0;i<num;i++)
    {
        dot=0;
        for (j=0;j<i;j++)
            dot+=LU[ps[i]][j]*x_m[j];

        x_m[i]=b_m[ps[i]]-dot;
    }

    for (i=num-1; i>=0;i--)
    {
        dot=0.0;

        for (j=i+1;j<num;j++)
            dot+=LU[ps[i]][j]*x_m[j];

        x_m[i]=(x_m[i]-dot)/LU[ps[i]][i];
    }
}


void leastSquaresSolver::m_solve(void)
{
    int i,j;
    double dot;

    for (i=0;i<VAR_NUM;i++)
    {
        dot=0;
        for (j=0;j<i;j++)
            dot+=LU[ps[i]][j]*x_m[j];

        x_m[i]=b_m[ps[i]]-dot;
    }

    for (i=VAR_NUM-1; i>=0;i--)
    {
        dot=0.0;

        for (j=i+1;j<VAR_NUM;j++)
            dot+=LU[ps[i]][j]*x_m[j];

        x_m[i]=(x_m[i]-dot)/LU[ps[i]][i];
    }
}

void leastSquaresSolver::m_invert(void)
{

    int i,j;
    //err_code(1);
    LU_decompose();

    for (j=0;j<VAR_NUM;j++)
    {

        for (i=0;i<VAR_NUM;i++)
        {
            if (i==j)
                b_m[i]=1;
            else
                b_m[i]=0;
        }

        m_solve();

        for (i=0;i<VAR_NUM;i++)
            Inv[i][j]=x_m[i];
    }
}

void leastSquaresSolver::m_invert2(int num)
{

    int i,j;
    //err_code(1);
    LU_decompose2(num);

    for (j=0;j<num;j++)
    {

        for (i=0;i<num;i++)
        {
            if (i==j)
                b_m[i]=1;
            else
                b_m[i]=0;
        }

        m_solve2(num);

        for (i=0;i<num;i++)
            Inv[i][j]=x_m[i];
    }
}


////////////////Kriging below

void leastSquaresSolver::LU_decompose_krig(int num)
{
    int i,j,k,pivotindex;
    static double scales[3000];
    double normrow,pivot,size,biggest,mult;

    for (i=0;i<num;i++) //заполнение начальными данными
    {
        ps_krig[i]=i;//маппинг изначального порядка на переставленный.
        normrow=0;//максимум в итой строке

        for (j=0;j<num;j++)
        {
            LU_krig[i][j]=M_krig[i][j];
            if (normrow<fabs(LU_krig[i][j]))
                normrow=fabs(LU_krig[i][j]);
        }
        if (normrow!=0)
            scales[i]=1.0/normrow;//для общих множителей
        else
        {
            scales[i]=0.0;
            //     err_code(DIV_ZERO);
        }
    }
    //метод гаусса с частичным упорядочиванием

    for (k=0;k<num-1;k++)
    {
        biggest=0;
        for (i=k; i<num;i++)
        {
            size=fabs(LU_krig[ps_krig[i]][k])*scales[ps_krig[i]];
            if (biggest<size)
            {
                biggest=size;
                pivotindex=i;
            }
        }

        if (biggest==0)
        {
            //	err_code(1);
            pivotindex=0;
        }

        if (pivotindex!=k)
        {
            j=ps_krig[k];
            ps_krig[k]=ps_krig[pivotindex];
            ps_krig[pivotindex]=j;
        }

        pivot=LU_krig[ps_krig[k]][k];

        for (i=k+1;i<num;i++)
        {
            mult=LU_krig[ps_krig[i]][k]/pivot;
            LU_krig[ps_krig[i]][k]=mult;

            if (mult!=0.0)
            {
                for (j=k+1; j<num;j++)
                    LU_krig[ps_krig[i]][j]-=mult*LU_krig[ps_krig[k]][j];
            }
        }
    }
    //      if (LU[ps[VAR_NUM-1]][VAR_NUM-1]==0.0) err_code(1);
}

void leastSquaresSolver::m_solve_krig(int num)
{
    int i,j;
    double dot;

    for (i=0;i<num;i++)
    {
        dot=0;
        for (j=0;j<i;j++)
            dot+=LU_krig[ps_krig[i]][j]*x_m[j];

        x_m[i]=b_m[ps_krig[i]]-dot;
    }

    for (i=num-1; i>=0;i--)
    {
        dot=0.0;

        for (j=i+1;j<num;j++)
            dot+=LU_krig[ps_krig[i]][j]*x_m[j];

        x_m[i]=(x_m[i]-dot)/LU_krig[ps_krig[i]][i];
    }
}



void leastSquaresSolver::m_invert_krig(int num)
{

    int i,j;
    //err_code(1);
    LU_decompose_krig(num);

    for (j=0;j<num;j++)
    {

        for (i=0;i<num;i++)
        {
            if (i==j)
                b_m[i]=1;
            else
                b_m[i]=0;
        }

        m_solve_krig(num);

        for (i=0;i<num;i++)
            Inv_krig[i][j]=x_m[i];
    }
}

///////////////Kriging above

void leastSquaresSolver::nullify_m(double m[MAX_EQNS][50])
{
    for (int i=0;i<40;i++)
        for(int j=0;j<VAR_NUM;j++)
            m[i][j]=0;
}

void leastSquaresSolver::nullify_v(double v[MAX_EQNS])
{
    for (int i=0;i<100;i++)
        v[i]=0;
}
void leastSquaresSolver::get_derivs(node3d &p, deriv3D &res,double delta)
{
    int var_num=VAR_NUM;
    int eq_num=m_p.size();
    // printf("eq_num=%d \n",eq_num);
    //first index is a row number
    //second index is a column number  M[eq_num][var_num]
    // get Mt*W*M
    nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dx,dy,dz;
        dx=m_p[i].x-p.x;
        dy=m_p[i].y-p.y;
        dz=m_p[i].z-p.z;
        double r2=dx*dx+dy*dy+dz*dz;
        M_0[i][0]=1.0;
        M_0[i][1]=dx;             M_0[i][2]=dy;       M_0[i][3]=dz;
        M_0[i][4]=0.5*dx*dx;     M_0[i][5]=dx*dy;   M_0[i][6]=dx*dz;
        M_0[i][7]=0.5*dy*dy;     M_0[i][8]=dy*dz;   M_0[i][9]=0.5*dz*dz;

        w[i]=exp(-r2/(delta*delta));//1.0;

        b_m[i]=0.0;
    }

    nullify_m(M_);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }
    m_invert();//now mvm inverted is in Inv; and its[var_num][var_num];

    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=m_p[i].f;
    }

    nullify_v(mwb);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    nullify_v(res.d); //now to the solution

    nullify_m(MWM);

    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int k=0;k<var_num;k++)
            {
                MWM[i][j]+=M_[i][k]*Inv[k][j];

            }
            //    printf("%f ",MWM[i][j]);

        }
        //  printf("\n ");
    }

    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            res.d[i]+=Inv[i][j]*mwb[j];
        }
    }
}


void leastSquaresSolver::get_derivs_fast(node3d &p, deriv3D &res,double delta)
{
    int var_num=VAR_NUM;
    int eq_num=m_p.size();
    // printf("eq_num=%d \n",eq_num);
    //first index is a row number
    //second index is a column number  M[eq_num][var_num]
    // get Mt*W*M
    nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dx,dy,dz;
        dx=m_p[i].x-p.x;
        dy=m_p[i].y-p.y;
        dz=m_p[i].z-p.z;
        double r2=dx*dx+dy*dy+dz*dz;
        M_0[i][0]=1.0;
        M_0[i][1]=dx;             M_0[i][2]=dy;       M_0[i][3]=dz;
        M_0[i][4]=0.5*dx*dx;     M_0[i][5]=dx*dy;   M_0[i][6]=dx*dz;
        M_0[i][7]=0.5*dy*dy;     M_0[i][8]=dy*dz;   M_0[i][9]=0.5*dz*dz;

        w[i]=exp(-r2/(delta*delta));//1.0;

        b_m[i]=0.0;
    }

    nullify_m(M_);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }


    for (int i=0;i<eq_num;i++)
    {
        b_m[i]=m_p[i].f;
    }

    nullify_v(mwb);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose();
    m_solve();

    for (int i=0;i<var_num;i++)
    {

        res.d[i]=x_m[i];

    }
}






double leastSquaresSolver::dist(node3d &p1,node3d &p2)
{

    return pow((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z),pww);//sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

double leastSquaresSolver::distDx(node3d &p1,node3d &p2)
{
    // return (p1.x - p2.x)/sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
    //  2.0*(p1.x - p2.x);
    // pow((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z),0.25);


    //first deriv
    // return pww*2.0*(p1.x - p2.x)/pow((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z),1.0-pww);


    double dx=(p1.x - p2.x);
    double r2=dx*dx + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z)+rad*rad;

    //zero deriv
    //return pow(r2,pww-1);
    //first driv
    //return pww*2.0*dx*pow(r2,pww-1.0);
    //second deriv
    return pww*2.0*pow(r2,pww-1) + 2.0*pww*dx*2.0*(pww-1.0)*dx*pow(r2,pww-2.0);
}


double leastSquaresSolver::distAx(node3d &p1,node3d &p2)
{
    // return (p1.x - p2.x)/sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
    //  2.0*(p1.x - p2.x);
    // pow((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z),0.25);


    //first deriv
    // return pww*2.0*(p1.x - p2.x)/pow((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z),1.0-pww);


    double dx=(p1.x - p2.x);
    double dy=(p1.y - p2.y);
    double dz=(p1.z - p2.z);
    double r2=dx*dx + dy*dy + dz*dz + rad*rad;

    //zero deriv
    //return pow(r2,pww-1);
    //first driv
    return pww*2.0*dx*pow(r2,pww-1.0);
    //second deriv

}

double leastSquaresSolver::distAy(node3d &p1,node3d &p2)
{
    // return (p1.x - p2.x)/sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
    //  2.0*(p1.x - p2.x);
    // pow((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z),0.25);


    //first deriv
    // return pww*2.0*(p1.x - p2.x)/pow((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z),1.0-pww);


    double dx=(p1.x - p2.x);
    double dy=(p1.y - p2.y);
    double dz=(p1.z - p2.z);
    double r2=dx*dx + dy*dy + dz*dz + rad*rad;

    //zero deriv
    //return pow(r2,pww-1);
    //first driv
    return pww*2.0*dy*pow(r2,pww-1.0);
    //second deriv

}

double leastSquaresSolver::distAz(node3d &p1,node3d &p2)
{
    // return (p1.x - p2.x)/sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
    //  2.0*(p1.x - p2.x);
    // pow((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z),0.25);


    //first deriv
    // return pww*2.0*(p1.x - p2.x)/pow((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z),1.0-pww);


    double dx=(p1.x - p2.x);
    double dy=(p1.y - p2.y);
    double dz=(p1.z - p2.z);
    double r2=dx*dx + dy*dy + dz*dz + rad*rad;

    //zero deriv
    //return pow(r2,pww-1);
    //first driv
    return pww*2.0*dz*pow(r2,pww-1.0);
    //second deriv

}

double leastSquaresSolver::distLapl(node3d &p1,node3d &p2)
{
    // return (p1.x - p2.x)/sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
    //  2.0*(p1.x - p2.x);
    // pow((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z),0.25);


    //first deriv
    // return pww*2.0*(p1.x - p2.x)/pow((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z),1.0-pww);


    double dx=(p1.x - p2.x);
    double dy=(p1.y - p2.y);
    double dz=(p1.z - p2.z);
    double r2=dx*dx + dy*dy + dz*dz + rad*rad;

    //zero deriv
    //return pow(r2,pww-1);
    //first driv
    //return pww*2.0*dx*pow(r2,pww-1.0);
    //second deriv
    double d2x=pww*2.0*pow(r2,pww-1) + 2.0*pww*dx*2.0*(pww-1.0)*dx*pow(r2,pww-2.0);
    double d2y=pww*2.0*pow(r2,pww-1) + 2.0*pww*dy*2.0*(pww-1.0)*dy*pow(r2,pww-2.0);
    double d2z=pww*2.0*pow(r2,pww-1) + 2.0*pww*dz*2.0*(pww-1.0)*dz*pow(r2,pww-2.0);

    return d2x+d2y+d2z;
}

void leastSquaresSolver::interpKrig(node3d &p)
{
    int size=m_p.size()+1;

    for (int i=0;i<size-1;i++)
    {
        M_[i][i] = 0;
        for (int j=0;j<i;j++)
        {
            M_[i][j] = M_[j][i] = dist(m_p.at(i), m_p.at(j));
        }
    }
    for (int i=0;i<size;i++)
    {
        M_[i][size-1] = M_[size-1][i] = 1.0;
    }
    M_[size-1][size-1] = 0;

    //ax:
    for (int i=0;i<size-1;i++)
    {
        b_m[i]=dist(p, m_p.at(i));
        //printf("b_m=%f \n", b_m[i]);
    }
    b_m[size-1] = 1.0;

    LU_decompose2(size);
    m_solve2(size);


    //ay:
    p.f_bound = 0;
    for (int i=0;i<size-1;i++)
    {
        p.f_bound+=m_p.at(i).f * x_m[i];
    }
    //printf("pf=%f \n", p.f);
}


void leastSquaresSolver::interpKrigDx(node3d &p)
{
    int size=m_p.size()+1;

    for (int i=0;i<size-1;i++)
    {
        M_[i][i] = 0;
        for (int j=0;j<i;j++)
        {
            M_[i][j] = M_[j][i] = dist(m_p.at(i), m_p.at(j));
        }
    }
    for (int i=0;i<size;i++)
    {
        M_[i][size-1] = M_[size-1][i] = 1.0;
    }
    M_[size-1][size-1] = 0;

    //ax:
    for (int i=0;i<size-1;i++)
    {
        b_m[i]=distDx(p, m_p.at(i));
        //printf("b_m=%f \n", b_m[i]);
    }
    b_m[size-1] = 1.0;

    LU_decompose2(size);
    m_solve2(size);


    //ay:
    p.f_bound = 0;
    for (int i=0;i<size-1;i++)
    {
        p.f_bound+=m_p.at(i).f * x_m[i];
    }
    //printf("pf=%f \n", p.f);
}


void leastSquaresSolver::getKrigInv()
{
    int size=m_p.size()+1;

    for (int i=0;i<size-1;i++)
    {
        M_krig[i][i] = 0;
        for (int j=0;j<i;j++)
        {
            M_krig[i][j] = M_krig[j][i] = dist(m_p.at(i), m_p.at(j));
        }
    }
    for (int i=0;i<size;i++)
    {
        M_krig[i][size-1] = M_krig[size-1][i] = 1.0;
    }
    M_krig[size-1][size-1] = 0;


    LU_decompose_krig(size);
    printf("got krig LU \n");
    for (int j=0;j<size-1;j++)
    {
        for (int i=0;i<size-1;i++)
        {
            b_m[i]=distAx(m_p.at(j), m_p.at(i));//distAx(p, m_p.at(i));
            //printf("b_m=%f \n", b_m[i]);
        }
        b_m[size-1] = 1.0;
        m_solve_krig(size);
        for (int i=0;i<size-1;i++)
        {
            Wx_krig[j][i]=x_m[i];
        }

        for (int i=0;i<size-1;i++)
        {
            b_m[i]=distAy(m_p.at(j), m_p.at(i));//distAx(p, m_p.at(i));
            //printf("b_m=%f \n", b_m[i]);
        }
        b_m[size-1] = 1.0;
        m_solve_krig(size);
        for (int i=0;i<size-1;i++)
        {
            Wy_krig[j][i]=x_m[i];
        }

        for (int i=0;i<size-1;i++)
        {
            b_m[i]=distAz(m_p.at(j), m_p.at(i));//distAx(p, m_p.at(i));
            //printf("b_m=%f \n", b_m[i]);
        }
        b_m[size-1] = 1.0;
        m_solve_krig(size);
        for (int i=0;i<size-1;i++)
        {
            Wz_krig[j][i]=x_m[i];
        }
    }
    printf("got krig W \n");


}


void leastSquaresSolver::interpKrig_nomatr(node3d &p)
{

    /* int size=m_p.size()+1;

    //ax:
    for (int i=0;i<size-1;i++)
    {
        b_m[i]=dist(p, m_p.at(i));
        //printf("b_m=%f \n", b_m[i]);
    }
    b_m[size-1] = 1.0;


    m_solve_krig(size);


    //ay:
    p.f_bound = 0;
    for (int i=0;i<size-1;i++)
    {
        p.f_bound+=m_p.at(i).f * x_m[i];
    }
*/


    int size=m_p.size()+1;

    /*  for (int i=0;i<size-1;i++)
    {
        M_krig[i][i] = 0;
        for (int j=0;j<i;j++)
        {
            M_krig[i][j] = M_krig[j][i] = dist(m_p.at(i), m_p.at(j));
        }
    }
    for (int i=0;i<size;i++)
    {
        M_krig[i][size-1] = M_krig[size-1][i] = 1.0;
    }
    M_krig[size-1][size-1] = 0;*/

    //getKrigInv();

    //LU_decompose_krig(size);
    //ax:
    for (int i=0;i<size-1;i++)
    {
        b_m[i]=dist(p, m_p.at(i));
        //printf("b_m=%f \n", b_m[i]);
    }
    b_m[size-1] = 1.0;

    //getKrigInv();
    m_solve_krig(size);


    //ay:
    p.f_bound = 0;
    for (int i=0;i<size-1;i++)
    {
        p.f_bound+=m_p.at(i).f * x_m[i];
    }
}


double f_krig[3000], df_krig[3000],df1_krig[3000];
void leastSquaresSolver::solveKrig_grad(int itn,double d)
{
    int size=m_p.size()+1;
    for (int n=0;n<size-1; n++)
    {
        f_krig[n]=m_p.at(n).f;
        df_krig[n]=0.0;
    }
    double delta=d;
    for (int nn=0;nn<itn;nn++)
    {
        //delta=delta*1.01;
        for (int j=0;j<size-1;j++)
        {
            double ax,ay,az;
            ax=m_p.at(j).ax;
            ay=m_p.at(j).ay;
            az=m_p.at(j).az;

            double fwx,fwy,fwz;
            fwx=0.0; fwy=0.0;fwz=0.0;
            for (int k=0;k<size-1;k++)
            {
                fwx+=f_krig[k]*Wx_krig[j][k];
                fwy+=f_krig[k]*Wy_krig[j][k];
                fwz+=f_krig[k]*Wz_krig[j][k];
            }
            fwx-=ax;
            fwy-=ay;
            fwz-=az;

            printf("fwx = %e,fwy=%e,fwz=%e \n",fwx,fwy,fwz);

           for (int i=0;i<size-1;i++)
             {
               //if (i!=j)
                //df_krig[i]+=fwx*Wx_krig[j][i] + fwy*Wy_krig[j][i]  + fwz*Wz_krig[j][i] ;
               df_krig[i]+=fwx*Wx_krig[j][i] + fwy*Wy_krig[j][i]  + fwz*Wz_krig[j][i] ;

               //double df=fwx*Wx_krig[j][i] + fwy*Wy_krig[j][i]  + fwz*Wz_krig[j][i] ;
               //f_krig[i]-=df*delta;
            }
        }
        for (int i=0;i<size-1;i++)
        {
           f_krig[i]-=df_krig[i]*delta;
           df_krig[i]=0.0;
        }
    }

    for (int i=0;i<size-1;i++)
    {
        m_p.at(i).f=f_krig[i];
    }
}

double x_lapl[3000], x_dx[3000],x_dy[3000],x_dz[3000];

void leastSquaresSolver::solvePoisson_krig()
{


    int size=m_p.size()+1;


    double lapl,f_lapl;
    double ax,f_ax;
    double ay,f_ay;
    double az,f_az;

    for (int n=0;n<size-1; n++)
    {
        for (int i=0;i<size-1;i++)
        {
            b_m[i]=distLapl(m_p.at(n), m_p.at(i));
        }
        b_m[size-1] = 1.0;
        m_solve_krig(size);
        lapl = 0.0;

        for (int i=0;i<size-1;i++)
        {
            if (i!=n)
                lapl+=m_p.at(i).f * x_m[i];
        }
        //lapl+m_p.at(n).f*x_m[n]=m_p.at(n).rhs;
        f_lapl=(m_p.at(n).rhs-lapl)/x_m[n];


        for (int i=0;i<size-1;i++)
        {
            b_m[i]=distAx(m_p.at(n), m_p.at(i));
        }
        b_m[size-1] = 1.0;
        m_solve_krig(size);
        ax = 0.0;

        for (int i=0;i<size-1;i++)
        {
            if (i!=n)
                ax+=m_p.at(i).f * x_m[i];
        }
        //lapl+m_p.at(n).f*x_m[n]=m_p.at(n).rhs;
        f_ax=(m_p.at(n).ax-ax)/x_m[n];



        for (int i=0;i<size-1;i++)
        {
            b_m[i]=distAy(m_p.at(n), m_p.at(i));
        }
        b_m[size-1] = 1.0;
        m_solve_krig(size);
        ay = 0.0;

        for (int i=0;i<size-1;i++)
        {
            if (i!=n)
                ay+=m_p.at(i).f * x_m[i];
        }
        //lapl+m_p.at(n).f*x_m[n]=m_p.at(n).rhs;
        f_ay=(m_p.at(n).ay-ay)/x_m[n];




        for (int i=0;i<size-1;i++)
        {
            b_m[i]=distAz(m_p.at(n), m_p.at(i));
        }
        b_m[size-1] = 1.0;
        m_solve_krig(size);
        az = 0.0;

        for (int i=0;i<size-1;i++)
        {
            if (i!=n)
                az+=m_p.at(i).f * x_m[i];
        }
        //lapl+m_p.at(n).f*x_m[n]=m_p.at(n).rhs;
        f_az=(m_p.at(n).az-az)/x_m[n];


        m_p.at(n).f=m_p.at(n).f*0.5+ 0.5*(f_lapl);//*0.5+(0.5/3.0)*(f_ax+f_ay+f_az));
    }
}




/////////////////////

void leastSquaresSolver::get_poisson_internal(node3d &p, deriv3D &res, double delta)
{
    int var_num=VAR_NUM;
    int eq_num=m_p.size();

    nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dx,dy,dz;
        dx=m_p[i].x-p.x;
        dy=m_p[i].y-p.y;
        dz=m_p[i].z-p.z;
        double r2=dx*dx+dy*dy+dz*dz;
        M_0[i][0]=1.0;
        M_0[i][1]=dx;             M_0[i][2]=dy;       M_0[i][3]=dz;
        M_0[i][4]=0.5*dx*dx;     M_0[i][5]=dx*dy;   M_0[i][6]=dx*dz;
        M_0[i][7]=0.5*dy*dy;     M_0[i][8]=dy*dz;   M_0[i][9]=0.5*dz*dz;

        w[i]=exp(-r2/(delta*delta));//or just 1.0;

        b_m[i]=0.0;
    }
    int i=eq_num; //constraint for laplacian
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0;             M_0[i][2]=0;       M_0[i][3]=0;
    M_0[i][4]=1;             M_0[i][5]=0;       M_0[i][6]=0;
    M_0[i][7]=1;             M_0[i][8]=0;       M_0[i][9]=1;

    w[i]=1.0;

    b_m[i]=0.0;

    nullify_m(M_);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }


    for (int i=0;i<eq_num-1;i++)
    {
        b_m[i]=m_p[i].f;//from prev interation
    }
    b_m[eq_num-1]=p.rhs;//this is the rhs value


    nullify_v(mwb);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose();
    m_solve();

    for (int i=0;i<var_num;i++)
    {

        res.d[i]=x_m[i];

    }
    p.f=res.d[F];
}

void leastSquaresSolver::get_poisson_boundary(node3d &p, deriv3D &res, double delta)
{
    //we do not chack if p is a boundary point this should be done before calling this function

    int var_num=VAR_NUM;
    int eq_num=m_p.size();

    nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dx,dy,dz;
        dx=m_p[i].x-p.x;
        dy=m_p[i].y-p.y;
        dz=m_p[i].z-p.z;
        double r2=dx*dx+dy*dy+dz*dz;
        M_0[i][0]=1.0;
        M_0[i][1]=dx;             M_0[i][2]=dy;       M_0[i][3]=dz;
        M_0[i][4]=0.5*dx*dx;     M_0[i][5]=dx*dy;   M_0[i][6]=dx*dz;
        M_0[i][7]=0.5*dy*dy;     M_0[i][8]=dy*dz;   M_0[i][9]=0.5*dz*dz;

        w[i]=exp(-r2/(delta*delta));//1.0;  todo: make w anisotrpopic in the direction normal to the wall

        b_m[i]=0.0;
    }
    int i=eq_num; //constraint for laplacian
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0;             M_0[i][2]=0;       M_0[i][3]=0;
    M_0[i][4]=1;             M_0[i][5]=0;       M_0[i][6]=0;
    M_0[i][7]=1;             M_0[i][8]=0;       M_0[i][9]=1;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;
    /////////////////////

    i=eq_num; //constraint for gradient at the wall
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=walls[p.is_boundary].nx;             M_0[i][2]=walls[p.is_boundary].ny;       M_0[i][3]=walls[p.is_boundary].nz;
    M_0[i][4]=0.0;     M_0[i][5]=0.0;   M_0[i][6]=0.0;
    M_0[i][7]=0.0;     M_0[i][8]=0.0;   M_0[i][9]=0.0;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;

    nullify_m(M_);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }


    for (int i=0;i<eq_num-2;i++)
    {
        b_m[i]=m_p[i].f;//from prev interation
    }
    b_m[eq_num-2]=p.rhs;//from prev interation
    b_m[eq_num-1]=p.f_bound;//from prev interation



    nullify_v(mwb);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose();
    m_solve();

    for (int i=0;i<var_num;i++)
    {

        res.d[i]=x_m[i];

    }
    p.f=res.d[F];
}


void leastSquaresSolver::get_poisson_combo(node3d &p, deriv3D &res, double delta)
{
    //we do not chack if p is a boundary point this should be done before calling this function

    int var_num=VAR_NUM;
    int eq_num=m_p.size();

    //nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=m_p.at(i).x-p.x;
        dyl=m_p.at(i).y-p.y;
        dzl=m_p.at(i).z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=exp(-r2/(delta*delta));//1.0;  todo: make w anisotrpopic in the direction normal to the wall

        b_m[i]=0.0;
    }
    int i=eq_num; //constraint for laplacian
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0;             M_0[i][2]=0;       M_0[i][3]=0;
    M_0[i][4]=1;             M_0[i][5]=0;       M_0[i][6]=0;
    M_0[i][7]=1;             M_0[i][8]=0;       M_0[i][9]=1;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;
    /////////////////////

    i=eq_num; //constraint for gradient at the dx
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=1.0;     M_0[i][2]=0.0;       M_0[i][3]=0.0;
    M_0[i][4]=0.0;     M_0[i][5]=0.0;   M_0[i][6]=0.0;
    M_0[i][7]=0.0;     M_0[i][8]=0.0;   M_0[i][9]=0.0;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;


    i=eq_num; //constraint for gradient at the dy
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0.0;     M_0[i][2]=1.0;   M_0[i][3]=0.0;
    M_0[i][4]=0.0;     M_0[i][5]=0.0;   M_0[i][6]=0.0;
    M_0[i][7]=0.0;     M_0[i][8]=0.0;   M_0[i][9]=0.0;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;


    i=eq_num; //constraint for gradient at the dz
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0.0;     M_0[i][2]=0.0;   M_0[i][3]=1.0;
    M_0[i][4]=0.0;     M_0[i][5]=0.0;   M_0[i][6]=0.0;
    M_0[i][7]=0.0;     M_0[i][8]=0.0;   M_0[i][9]=0.0;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;

    nullify_m(M_);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }


    for (int i=0;i<eq_num-4;i++)
    {
        b_m[i]=m_p.at(i).f;//from prev interation
    }
    b_m[eq_num-4]=p.rhs;//from prev interation
    b_m[eq_num-3]=p.ax;//from prev interation
    b_m[eq_num-2]=p.ay;//from prev interation
    b_m[eq_num-1]=p.az;//from prev interation


    nullify_v(mwb);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }

    for (int i=0;i<var_num;i++)
    {
        b_m[i]=mwb[i];  //its mvb

    }
    LU_decompose();
    m_solve();

    /*for (int i=0;i<var_num;i++)
    {

            res.d[i]=x_m[i];

    }*/
    p.f=x_m[0];//res.d[F];
}




void leastSquaresSolver::get_poisson_matr(node3d &p, deriv3D &res, double delta)
{

    int var_num=VAR_NUM;
    int eq_num=m_p.size();

    //nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=m_p.at(i).x-p.x;
        dyl=m_p.at(i).y-p.y;
        dzl=m_p.at(i).z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=exp(-r2/(delta*delta));//1.0;  todo: make w anisotrpopic in the direction normal to the wall

        b_m[i]=0.0;
    }
    int i=eq_num; //constraint for laplacian
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0;             M_0[i][2]=0;       M_0[i][3]=0;
    M_0[i][4]=1;             M_0[i][5]=0;       M_0[i][6]=0;
    M_0[i][7]=1;             M_0[i][8]=0;       M_0[i][9]=1;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;
    /////////////////////

    i=eq_num; //constraint for gradient at the dx
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=1.0;     M_0[i][2]=0.0;       M_0[i][3]=0.0;
    M_0[i][4]=0.0;     M_0[i][5]=0.0;   M_0[i][6]=0.0;
    M_0[i][7]=0.0;     M_0[i][8]=0.0;   M_0[i][9]=0.0;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;


    i=eq_num; //constraint for gradient at the dy
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0.0;     M_0[i][2]=1.0;   M_0[i][3]=0.0;
    M_0[i][4]=0.0;     M_0[i][5]=0.0;   M_0[i][6]=0.0;
    M_0[i][7]=0.0;     M_0[i][8]=0.0;   M_0[i][9]=0.0;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;


    i=eq_num; //constraint for gradient at the dz
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0.0;     M_0[i][2]=0.0;   M_0[i][3]=1.0;
    M_0[i][4]=0.0;     M_0[i][5]=0.0;   M_0[i][6]=0.0;
    M_0[i][7]=0.0;     M_0[i][8]=0.0;   M_0[i][9]=0.0;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;

    nullify_m(M_);
    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {
            for (int n=0;n<eq_num;n++)
            {
                M_[i][j]+=M_0[n][i]*M_0[n][j]*w[n];  //its mvm
            }
        }
    }

    m_invert();

    for (int i=0;i<var_num;i++)
    {
        for (int j=0;j<var_num;j++)
        {

            p.INV[i][j]=Inv[i][j];


        }
    }


}

void leastSquaresSolver::get_poisson_nomatr(node3d &p, deriv3D &res, double delta)
{
    //we do not chack if p is a boundary point this should be done before calling this function

    int var_num=VAR_NUM;
    int eq_num=m_p.size();

    //nullify_m(M_0);
    for (int i=0;i<eq_num;i++)
    {
        double dxl,dyl,dzl;
        dxl=m_p.at(i).x-p.x;
        dyl=m_p.at(i).y-p.y;
        dzl=m_p.at(i).z-p.z;
        double r2=dxl*dxl+dyl*dyl+dzl*dzl;
        M_0[i][0]=1.0;
        M_0[i][1]=dxl;             M_0[i][2]=dyl;       M_0[i][3]=dzl;
        M_0[i][4]=0.5*dxl*dxl;     M_0[i][5]=dxl*dyl;   M_0[i][6]=dxl*dzl;
        M_0[i][7]=0.5*dyl*dyl;     M_0[i][8]=dyl*dzl;   M_0[i][9]=0.5*dzl*dzl;

        w[i]=exp(-r2/(delta*delta));//1.0;  todo: make w anisotrpopic in the direction normal to the wall

        b_m[i]=0.0;
    }
    int i=eq_num; //constraint for laplacian
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0;             M_0[i][2]=0;       M_0[i][3]=0;
    M_0[i][4]=1;             M_0[i][5]=0;       M_0[i][6]=0;
    M_0[i][7]=1;             M_0[i][8]=0;       M_0[i][9]=1;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;
    /////////////////////

    i=eq_num; //constraint for gradient at the dx
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=1.0;     M_0[i][2]=0.0;       M_0[i][3]=0.0;
    M_0[i][4]=0.0;     M_0[i][5]=0.0;   M_0[i][6]=0.0;
    M_0[i][7]=0.0;     M_0[i][8]=0.0;   M_0[i][9]=0.0;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;


    i=eq_num; //constraint for gradient at the dy
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0.0;     M_0[i][2]=1.0;   M_0[i][3]=0.0;
    M_0[i][4]=0.0;     M_0[i][5]=0.0;   M_0[i][6]=0.0;
    M_0[i][7]=0.0;     M_0[i][8]=0.0;   M_0[i][9]=0.0;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;


    i=eq_num; //constraint for gradient at the dz
    eq_num++;
    M_0[i][0]=0.0;
    M_0[i][1]=0.0;     M_0[i][2]=0.0;   M_0[i][3]=1.0;
    M_0[i][4]=0.0;     M_0[i][5]=0.0;   M_0[i][6]=0.0;
    M_0[i][7]=0.0;     M_0[i][8]=0.0;   M_0[i][9]=0.0;

    w[i]=1.0;//1.0;

    b_m[i]=0.0;


    for (int i=0;i<eq_num-4;i++)
    {
        b_m[i]=m_p.at(i).f;//from prev interation
    }
    b_m[eq_num-4]=p.rhs;//from prev interation
    b_m[eq_num-3]=p.ax;//from prev interation
    b_m[eq_num-2]=p.ay;//from prev interation
    b_m[eq_num-1]=p.az;//from prev interation


    nullify_v(mwb);
    for (int i=0;i<var_num;i++)
    {
        for (int n=0;n<eq_num;n++)
        {
            mwb[i]+=M_0[n][i]*w[n]*b_m[n];  //its mvb
        }
    }


    for (int i=0;i<var_num;i++)
    {
        res.d[i]=0.0;
        for (int j=0;j<var_num;j++)
        {
            res.d[i]+=p.INV[i][j]*mwb[j];
        }
    }

    p.f=res.d[F];
}

void leastSquaresSolver::draw_points(double sc)
{
    glPointSize(4);
    /*glBegin(GL_POINTS);
    for (int i=0;i<m_p.size();i++)
    {
        glColor3f(m_p[i].f,m_p[i].f,-m_p[i].f);
        glVertex3f(m_p[i].x,m_p[i].y,m_p[i].z);
    }
    glEnd();*/

    glPointSize(4);
    glBegin(GL_POINTS);
    for (int i=0;i<m_p.size();i++)
    {
        //glColor3f(10.0*m_d0[i].d[F],10.0*m_d0[i].d[F],-10.0*m_d0[i].d[F]);
        glColor3f(sc*10.0*m_p[i].f,sc*10.0*m_p[i].f,-sc*10.0*m_p[i].f);
        //glColor3f(1,1,0);
        glVertex3f(m_p[i].x,m_p[i].y,m_p[i].z);
    }
    glEnd();

    /* glBegin(GL_LINES);
    for (int i=0;i<m_p.size();i++)
    {
        glColor3f(1,1,1);
        glVertex3f(m_p[i].x,m_p[i].y,m_p[i].z);
        glVertex3f(m_p[i].x + sc*m_d[i].d[FX] ,
                   m_p[i].y + sc*m_d[i].d[FY] ,
                   m_p[i].z + sc*m_d[i].d[FX]);
    }
    glEnd();

    glBegin(GL_LINES);
    for (int i=0;i<m_p.size();i++)
    {
        glColor3f(1,0,0);
        glVertex3f(m_p[i].x,m_p[i].y,m_p[i].z);
        glVertex3f(m_p[i].x + sc*m_d0[i].d[FX] ,
                   m_p[i].y + sc*m_d0[i].d[FY] ,
                   m_p[i].z + sc*m_d0[i].d[FX]);
    }
    glEnd();*/
}

double bench_f(double x,double y,double z)//just a gauss function for benchmark
{
    double r2=x*x+y*y+z*z;
    return exp(-r2/0.1);
}
double bench_dfdx(double x,double y,double z)
{
    double r2=x*x+y*y+z*z;
    return -2*x*exp(-r2/0.1)/0.1;
}

double bench_dfdy(double x,double y,double z)
{
    double r2=x*x+y*y+z*z;
    return -2*y*exp(-r2);
}

double bench_dfdz(double x,double y,double z)
{
    double r2=x*x+y*y+z*z;
    return -2*z*exp(-r2);
}

double bench_d2fdx(double x,double y,double z)
{
    double r2=x*x+y*y+z*z;
    double dd=0.1;
    return -(2/dd)*exp(-r2/dd)*(1.0-2*x*x/dd);
}

double bench_d2fdy(double x,double y,double z)
{
    double r2=x*x+y*y+z*z;
    return (4*y*y - 2)*exp(-r2);
}

double bench_d2fdz(double x,double y,double z)
{
    double r2=x*x+y*y+z*z;
    return (4*z*z - 2)*exp(-r2);
}

void leastSquaresSolver::get_derivs_bench(node3d &p, deriv3D &res)
{
    double x,y,z;

    x=p.x;
    y=p.y;
    z=p.z;

    res.d[F]=bench_f(x,y,z);
    res.d[FX]=bench_dfdx(x,y,z);
    res.d[FY]=bench_dfdy(x,y,z);
    res.d[FZ]=bench_dfdz(x,y,z);

    res.d[FXX]=bench_d2fdx(x,y,z);
    res.d[FYY]=bench_d2fdy(x,y,z);
    res.d[FZZ]=bench_d2fdz(x,y,z);
}

