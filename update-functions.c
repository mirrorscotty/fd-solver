#include <stdlib.h>
#include "fd-solver.h"
#include "can.h"

double UpdateTest(struct Node1D *node, int Var)
{
    int t = node->TimeIndex;
    /* t is t */
//    node->Value[t] = t;
    return (double) t;
}

/* Heat conduction in a slab */
double UpdateSubdomain(struct Node1D *node, int Var)
{
    double T, Ta, Tb, Tc;
    int t = node->TimeIndex;
    double Cond, Dens, Cap, alpha, M;

    Ta = node->Prev->Value[Var][t];
    Tb = node->Value[Var][t];
    Tc = node->Next->Value[Var][t];

    /* Calculate the thermal properties based on the temperature of the node at
    the pervious time step. */
    Cond = k(Tb); //0.5985; // W/(m K)
    Dens = rho(Tb); //1000; // kg/m^3
    Cap = Cp(Tb); //1865; // J/(kg K)

    alpha = Cond/(Dens*Cap);
    M = (node->dx)*(node->dx)/(alpha*node->dt);

    T = 1/M * (Ta + (M-2)*Tb + Tc);
//    node->Value[t+1] = T;
    
    return T;
}

double UpdateSubdomainRxn1(struct Node1D *node, int Var)
{
    double k; /* Rate Constant */
    double T;
    int t = node->TimeIndex;
    double dt = node->dt;
    double c = node->Value[Var][t];

    T = NodeGetValue(node, 'T', t);

    k = reaction_rate1(T, c);

    return c * (1-k*dt);
}

double UpdateSubdomainRxn2(struct Node1D *node, int Var)
{
    double k; /* Rate Constant */
    double T;
    int t = node->TimeIndex;
    double dt = node->dt;
    double c = node->Value[Var][t];

    T = NodeGetValue(node, 'T', t);

    k = reaction_rate2(T, c);

    return c * (1-k*dt);
}

double UpdateConvectiveBoundary(struct Node1D *node, int Var)
{
    double T, Text, Ta, Tb;
    int t = node->TimeIndex;
    double Cond, Dens, Cap, alpha, M, N, h;

    Tb = node->Value[Var][t];
    Text = T_ext(t*node->dt);

    /* Figure out which side the boundary is on */
    if(node->Next == NULL)
        Ta = node->Prev->Value[Var][t];
    else
        Ta = node->Next->Value[Var][t];

    h = 3000; //W/(m^2 K)
    Cond = k(Tb); //0.5985; // W/(m K)
    Dens = rho(Tb); //1000; // kg/m^3
    Cap = Cp(Tb); //1865; // J/(kg K)

    alpha = Cond/(Dens*Cap);
    M = (node->dx)*(node->dx)/(alpha*node->dt);
    N = h*node->dx/Cond;

    T = 1/M * (2*N*Text + (M-(2*N+2))*Tb + 2*Ta);

//    node->Value[t+1] = T;
    return T;
}

double UpdateInsulatedBoundary(struct Node1D *node, int Var)
{
    double T, Ta, Tb;
    int t = node->TimeIndex;
    double Cond, Dens, Cap, alpha, M;

    Tb = node->Value[Var][t];
    if(node->Next == NULL)
        Ta = node->Prev->Value[Var][t];
    else
        Ta = node->Next->Value[Var][t];

    Cond = k(Tb); //0.5985; // W/(m K)
    Dens = rho(Tb); //1000; // kg/m^3
    Cap = Cp(Tb); //1865; // J/(kg K)

    alpha = Cond/(Dens*Cap);
    M = (node->dx)*(node->dx)/(alpha*node->dt);

    T = 1/M * ((M-2)*Tb + 2*Ta);
//    node->Value[t+1] = T;

    return T;
}
