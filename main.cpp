
#include <stdio.h>
#include <stdlib.h>

//#include <my_include/gl.h>
//#include <my_include/glu.h>
//#include <my_include/glut.h>
#include <GL/glut.h>
#include  <math.h>
#include <time.h>

#include "leastsquaressolver.h"

double sc=0.5;

void display(void);
void init();

leastSquaresSolver ls;
node3d node; //for testing
deriv3D deriv,deriv_true; //for testing

double ck=1.0;//0.250;
double nu0=0.001;
int j_curr=0;

void display(void)
{
    int i,j;
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();

    glLineWidth(1.0);
    glPointSize(2.0);

    glBegin(GL_LINES);
        glColor3f(1,0,0);
        glVertex3f(node.x,node.y,node.z);
        glColor3f(0,0,0);
        glVertex3f(node.x + sc*deriv.d[leastSquaresSolver::FX] ,
                   node.y + sc*deriv.d[leastSquaresSolver::FY] ,
                   node.z + sc*deriv.d[leastSquaresSolver::FZ]);

        glColor3f(1,1,1);
        glVertex3f(node.x,node.y,node.z);
        glColor3f(0,0,0);
        glVertex3f(node.x + sc*deriv_true.d[leastSquaresSolver::FX] ,
                   node.y + sc*deriv_true.d[leastSquaresSolver::FY] ,
                   node.z + sc*deriv_true.d[leastSquaresSolver::FZ]);

    glEnd();

    ls.draw_points(sc);
    glutSwapBuffers();
}

double amp=1.0;
void solve_poisson()
{
    for(int i=0;i<ls.m_p.size();i++)
    {
        if (ls.m_p[i].is_boundary>=0)
        {
            ls.get_poisson_boundary(ls.m_p[i],ls.m_d[i],0.125);

        }else
        {
            ls.get_poisson_internal(ls.m_p[i],ls.m_d[i],0.125);
        }
        if (i%10==0)
        {printf("%f \n",i*100.0/ls.m_p.size());}
    }
    double f0=ls.m_p[0].f;
    for(int i=0;i<ls.m_p.size();i++)
    {

        ls.m_p[i].f-=f0;
    }

}

void kb(unsigned char key, int x, int y)
{


    if (key=='q')
    {
        sc*=1.01;
    }
    if (key=='e')
    {
        sc/=1.01;
    }
    if (key=='w')
    {
        node.y+=0.001;
        ls.get_derivs_fast(node,deriv,0.25);
        ls.get_derivs_bench(node,deriv_true);        
        printf("f_true=%f %f %f\n",deriv_true.d[leastSquaresSolver::F],deriv.d[leastSquaresSolver::F],deriv_true.d[leastSquaresSolver::F]-deriv.d[leastSquaresSolver::F]);

    }
    if (key=='s')
    {
        node.y-=0.001;
        ls.get_derivs_fast(node,deriv,0.25);
        ls.get_derivs_bench(node,deriv_true);
        printf("f_true=%f %f %f\n",deriv_true.d[leastSquaresSolver::F],deriv.d[leastSquaresSolver::F],deriv_true.d[leastSquaresSolver::F]-deriv.d[leastSquaresSolver::F]);

    }
    if (key=='a')
    {
        node.x-=0.001;
        ls.get_derivs_fast(node,deriv,0.25);
        ls.get_derivs_bench(node,deriv_true);
        printf("f_true=%f %f %f\n",deriv_true.d[leastSquaresSolver::F],deriv.d[leastSquaresSolver::F],deriv_true.d[leastSquaresSolver::F]-deriv.d[leastSquaresSolver::F]);

    }
    if (key=='d')
    {
        node.x+=0.001;
        ls.get_derivs_fast(node,deriv,0.25);
        ls.get_derivs_bench(node,deriv_true);
        printf("f_true=%f %f %f\n",deriv_true.d[leastSquaresSolver::F],deriv.d[leastSquaresSolver::F],deriv_true.d[leastSquaresSolver::F]-deriv.d[leastSquaresSolver::F]);

    }
    if (key==' ')
    {
       // printf(" sweep \n");
      /*  for (int i=0;i<50;i++)
            sweep(1,0.001);*/
        // go =! go;
       // ls.init();
        solve_poisson();
    }

    glutPostRedisplay();
}




void init()
{ 
  //нужно просто заполнить ккординаты точек (m_p[i].x .y. z) и значение функции (то есть ускорение в нашем случве) m_p[i].f
    ls.m_p.clear();
     ls.m_d.clear();
      ls.m_d0.clear();
    node3d n;
    deriv3D d0,d;
    //for derivatives bench:
/*    for (int i=0;i<39;i++)
    {
        n.x=0.5*(rand()*1.0/RAND_MAX-0.5);
        n.y=0.5*(rand()*1.0/RAND_MAX-0.5);
        n.z=0.5*(rand()*1.0/RAND_MAX-0.5);
        ls.get_derivs_bench(n,d0);//тут просто тестовый гаусс вычисляется, со всеми производными. Этого делать в обзем случае не надо, порсто завполнить m_p[i].f
        n.f=d0.d[leastSquaresSolver::F];
        ls.m_p.push_back(n);
    }
    //все. теперрь для получения производных нужно вызвать get_derivs(точка, массив проивзодных)
    //вычисляются сразу первые и вторые производные, см объявление класса
*/

    //for poisson bench

    for (int i=0;i<1000;i++)
    {
        n.x=0.0125+(rand()*0.975/RAND_MAX);
        n.y=0.0125+(rand()*0.975/RAND_MAX);
        n.z=0.1*(0.0125+(rand()*0.975/RAND_MAX));
        n.f=0.0+0.1*(0.0125+(rand()*0.975/RAND_MAX));
        n.is_boundary=-1;
        n.f_bound=0.0;
        n.rhs=-cos(M_PI*n.x);


        d0.d[leastSquaresSolver::F]=-cos(M_PI*n.x)/(M_PI*M_PI);
        ls.m_p.push_back(n);
        ls.m_d.push_back(d);
        ls.m_d0.push_back(d0);
    }
    boundary_plane bp;

    bp.nx=1.0; bp.ny=0.0; bp.nz=0.0; bp.d=0.0; //0 x0
    ls.walls.push_back(bp);

    bp.nx=1.0; bp.ny=0.0; bp.nz=0.0; bp.d=1.0; //1 x1
    ls.walls.push_back(bp);

    bp.nx=0.0; bp.ny=1.0; bp.nz=0.0; bp.d=0.0; //2 y0
    ls.walls.push_back(bp);

    bp.nx=0.0; bp.ny=1.0; bp.nz=0.0; bp.d=1.0; //3 y2
    ls.walls.push_back(bp);

    bp.nx=0.0; bp.ny=0.0; bp.nz=1.0; bp.d=0.0; //4 z0
    ls.walls.push_back(bp);

    bp.nx=0.0; bp.ny=0.0; bp.nz=1.0; bp.d=0.1; //5 z0
    ls.walls.push_back(bp);

    for (int i=0;i<70;i++)
    {
        n.x=0;
        n.y=(rand()*1.0/RAND_MAX);
        n.z=0.1*((rand()*1.0/RAND_MAX));
             n.f=0.0+0.2*(0.0125+(rand()*0.975/RAND_MAX));
        n.is_boundary=0;
        n.f_bound=0.0;
        n.rhs=-cos(M_PI*n.x);
        d0.d[leastSquaresSolver::F]=-cos(M_PI*n.x)/(M_PI*M_PI);
        ls.m_p.push_back(n);
        ls.m_d.push_back(d);
        ls.m_d0.push_back(d0);
       //
        n.x=1.0;
        n.y=(rand()*1.0/RAND_MAX);
        n.z=0.1*((rand()*1.0/RAND_MAX));
             n.f=0.0+0.2*(0.0125+(rand()*0.975/RAND_MAX));
        n.is_boundary=1;
        n.f_bound=0.0;
        n.rhs=-cos(M_PI*n.x);
        d0.d[leastSquaresSolver::F]=-cos(M_PI*n.x)/(M_PI*M_PI);
        ls.m_p.push_back(n);
        ls.m_d.push_back(d);
        ls.m_d0.push_back(d0);

        //
        n.x=(rand()*1.0/RAND_MAX);
        n.y=0.0;
        n.z=0.1*((rand()*1.0/RAND_MAX));
             n.f=0.0+0.2*(0.0125+(rand()*0.975/RAND_MAX));
        n.is_boundary=2;
        n.f_bound=0.0;
        n.rhs=-cos(M_PI*n.x);
         d0.d[leastSquaresSolver::F]=-cos(M_PI*n.x)/(M_PI*M_PI);
        ls.m_p.push_back(n);
        ls.m_d.push_back(d);
        ls.m_d0.push_back(d0);
       //
        n.x=(rand()*1.0/RAND_MAX);
        n.y=1.0;
        n.z=0.1*((rand()*1.0/RAND_MAX));
             n.f=0.0+0.2*(0.0125+(rand()*0.975/RAND_MAX));
        n.is_boundary=3;
        n.f_bound=0.0;
        n.rhs=-cos(M_PI*n.x);
        d0.d[leastSquaresSolver::F]=-cos(M_PI*n.x)/(M_PI*M_PI);
        ls.m_p.push_back(n);
        ls.m_d.push_back(d);
        ls.m_d0.push_back(d0);

    }

    for (int i=0;i<200;i++)
    {
        n.x=(rand()*1.0/RAND_MAX);
        n.y=(rand()*1.0/RAND_MAX);
        n.z=0.0;
             n.f=0.0+0.2*(0.0125+(rand()*0.975/RAND_MAX));
        n.is_boundary=4;
        n.f_bound=0.0;
        n.rhs=-cos(M_PI*n.x);
        d0.d[leastSquaresSolver::F]=-cos(M_PI*n.x)/(M_PI*M_PI);
        ls.m_p.push_back(n);
        ls.m_d.push_back(d);
        ls.m_d0.push_back(d0);
       //
        n.x=(rand()*1.0/RAND_MAX);
        n.y=(rand()*1.0/RAND_MAX);
        n.z=0.1;
             n.f=0.0+0.2*(0.0125+(rand()*0.975/RAND_MAX));
        n.is_boundary=5;
        n.f_bound=0.0;
        n.rhs=-cos(M_PI*n.x);
        d0.d[leastSquaresSolver::F]=-cos(M_PI*n.x)/(M_PI*M_PI);
        ls.m_p.push_back(n);
        ls.m_d.push_back(d);
        ls.m_d0.push_back(d0);

    }




    glClearColor (0.0, 0.1, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho(-1.1,  1.1, -1.1, 1.1, -10.0, 10.0);
    glMatrixMode (GL_MODELVIEW);
}
int main(int argc, char** argv)
{
    srand(time(NULL));
    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    //  glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(800,800);
    glutInitWindowPosition(0,0);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutKeyboardFunc(kb);
    init();
    glutMainLoop();
}
