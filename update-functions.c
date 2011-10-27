#include <stdlib.h>
#include "fd-solver.h"

double UpdateTest(struct Node1D *node)
{
    int t = node->TimeIndex;
    /* t is t */
    node->Value[t] = t;
    return node->Value[t];
}

/* Heat conduction in a slab */
double UpdateSubdomain(struct Node1D *node)
{
    double T, Ta, Tb, Tc;
    int t = node->TimeIndex;
    double k, rho, Cp, alpha, M;

    k = 0.5985; // W/(m K)
    rho = 1000; // kg/m^3
    Cp = 1865; // J/(kg K)

    alpha = k/(rho*Cp);
    M = (node->dx)*(node->dx)/(alpha*node->dt);

    Ta = node->Prev->Value[t];
    Tb = node->Value[t];
    Tc = node->Next->Value[t];

    T = 1/M * (Ta + (M-2)*Tb + Tc);
    node->Value[t+1] = T;
    
    return T;
}

double UpdateConvectiveBoundary(struct Node1D *node)
{
    double T, Text, Ta, Tb;
    int t = node->TimeIndex;
    double k, rho, Cp, alpha, M, N, h;

    h = 3000; //W/(m^2 K)
    k = 0.5985; // W/(m K)
    rho = 1000; // kg/m^3
    Cp = 1865; // J/(kg K)

    alpha = k/(rho*Cp);
    M = (node->dx)*(node->dx)/(alpha*node->dt);
    N = h*node->dx/k;

    Tb = node->Value[t];
    Text = 273.15;

    /* Figure out which side the boundary is on */
    if(node->Next == NULL)
        Ta = node->Prev->Value[t];
    else
        Ta = node->Next->Value[t];

    T = 1/M * (2*N*Text + (M-(2*N+2))*Tb + 2*Ta);

    node->Value[t+1] = T;
    return T;
}

double UpdateInsulatedBoundary(struct Node1D *node)
{
    double T, Ta, Tb;
    int t = node->TimeIndex;
    double k, rho, Cp, alpha, M;

    k = 0.5985; // W/(m K)
    rho = 1000; // kg/m^3
    Cp = 1865; // J/(kg K)

    alpha = k/(rho*Cp);
    M = (node->dx)*(node->dx)/(alpha*node->dt);

    Tb = node->Value[t];
    if(node->Next == NULL)
        Ta = node->Prev->Value[t];
    else
        Ta = node->Next->Value[t];

    T = 1/M * ((M-2)*Tb + 2*Ta);
    node->Value[t+1] = T;

    return T;
}
