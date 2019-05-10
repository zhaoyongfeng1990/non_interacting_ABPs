#ifndef PARTICLES_H
#define PARTICLES_H

#include <cstdint>
#include <cmath>
#include "mt64.h"
#include <iostream>
#include <fstream>
#include <mpi.h>
#include "config.h"
#include "./vectorclass/vectorclass.h"
#include "./vectorclass/vectormath_exp.h"
#include "./vectorclass/vectormath_trig.h"
#include "./vectorclass/vectormath_hyp.h"
#include "./vectorclass/ranvec1.h"

#define PERIODIC_IN_X

using namespace std;
using namespace vcl;

const double E = 2.7182818284590452354;
const double PI = 3.14159265358979323846264338327950288419716939937510;

typedef struct
{
    Vec4d x;
    Vec4d y;
    Vec4d theta;
} particle;

static Ranvec1 ran(3);
static MPI_Datatype MPI_PARTICLE;

double Force(Vec4d &__restrict__ x, Vec4d &__restrict__ y, Vec4d &__restrict__ fx, Vec4d &__restrict__ fy)
{
    const Vec4d vec0(0);
    const Vec4d vecLx(Lx);
    const Vec4d vecLy(Ly);
    fx=0;
    fy=0;
    Vec4d xsup=max(x-Lx,0);
    Vec4d ysup=max(y-Ly,0);
    Vec4d xinf=min(x,0);
    Vec4d yinf=min(y,0);
    fy += -4.0*k*pow(ysup,3)*(Lambda1+(Lambda2-Lambda1)*0.5*(tanh((x-x1)/w1)-tanh((x-x2)/w2)));
    fx += -k*pow(ysup,4)*(Lambda2-Lambda1)*0.5*(1.0/(w1*square(cosh((x-x1)/w1)))-1.0/(w2*square(cosh((x-x2)/w2))));
    fy += -4.0*k*pow(yinf,3);

#ifndef PERIODIC_IN_X
    fx += -4.0*k*pow(xinf,3);
    fx += -4.0*k*pow(xsup,3);
#endif

}

void randn(Vec4d &__restrict__ result)
{
    double u1,v1,u2,v2,w1,w2;
    w1=2.0;
    while(w1>=1.0)
    {
        Vec2d temp=ran.random2d();
        u1=2.0*temp[0]-1.0;
        v1=2.0*temp[1]-1.0;
        w1=u1*u1+v1*v1;
    }
    w2=2.0;
    while(w2>=1.0)
    {
        Vec2d temp=ran.random2d();
        u2=2.0*temp[0]-1.0;
        v2=2.0*temp[1]-1.0;
        w2=u2*u2+v2*v2;
    }

    Vec2d z={w1,w2};
    z=sqrt(-2.0*log(z)/z);
    result.insert(0, u1*z[0]);
    result.insert(1, v1*z[0]);
    result.insert(2, u2*z[1]);
    result.insert(3, v2*z[1]);
}

void step(particle *__restrict__ p)
{
    // static ofstream debug("debug.txt");
    Vec4d fx, fy, fxstar, fystar;
    Force(p->x, p->y, fx, fy);
    Vec4d noise;
    randn(noise);
    noise *= dev;

    //debug << noise << endl;
    Vec4d vx=V*cos(p->theta);
    Vec4d vy=V*sin(p->theta);

    Vec4d xstar = p->x + (vx + fx) * dt;
    Vec4d ystar = p->y + (vy + fy) * dt;
    Force(xstar, ystar, fxstar, fystar);
    p->x += (vx + (fx + fxstar) * 0.5) * dt;
    p->y += (vy + (fy + fystar) * 0.5) * dt;
    p->theta +=noise;
    
#ifdef PERIODIC_IN_X
    p->x=p->x-Lx*floor(p->x/Lx);
#endif
}

void collectingResult(particle *__restrict__ pList)
{
    int pid = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    int NP=0;
    MPI_Comm_size(MPI_COMM_WORLD, &NP);
    uint64_t localN=NumSample/NP;

    MPI_Status status;
    if (pid != 0)
    {
        MPI_Ssend(pList, localN, MPI_PARTICLE, 0, pid, MPI_COMM_WORLD);
    }
    else
    {
        for (int ic = 1; ic < NP; ++ic)
        {
            MPI_Recv(pList+ic*localN, localN, MPI_PARTICLE, ic, ic, MPI_COMM_WORLD, &status);
        }
    }
}

#endif