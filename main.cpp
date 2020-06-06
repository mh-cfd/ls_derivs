
#include <stdio.h>
#include <stdlib.h>

#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
//#include <GL/glut.h>
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
        ls.get_derivs(node,deriv,0.25);
        ls.get_derivs_bench(node,deriv_true);

    }
    if (key=='s')
    {
        node.y-=0.001;
        ls.get_derivs(node,deriv,0.25);
        ls.get_derivs_bench(node,deriv_true);
    }
    if (key=='a')
    {
        node.x-=0.001;
        ls.get_derivs(node,deriv,0.25);
        ls.get_derivs_bench(node,deriv_true);
    }
    if (key=='d')
    {
        node.x+=0.001;
        ls.get_derivs(node,deriv,0.25);
        ls.get_derivs_bench(node,deriv_true);
    }
    if (key==' ')
    {
       // printf(" sweep \n");
      /*  for (int i=0;i<50;i++)
            sweep(1,0.001);*/
        // go =! go;
       // ls.init();
    }

    glutPostRedisplay();
}



void init()
{ 
  //нужно просто заполнить ккординаты точек (m_p[i].x .y. z) и значение функции (то есть ускорение в нашем случве) m_p[i].f
    ls.m_p.clear();
    node3d n;
    deriv3D d0;
    for (int i=0;i<39;i++)
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
