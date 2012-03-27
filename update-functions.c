#include <stdlib.h>
#include "fd-solver.h"
#include "can.h"

double h = 0;

// Yay for meeting deadlines!
void seth(double x)
{
    h = x;
    printh();
    return;
}

void printh()
{
    printf("h = %g\n", h);
    return;
}

double UpdateTest(struct Node1D *node, int Var)
{
    int t = Var; // Get rid of the annoying warning about not using "Var"
    t = node->TimeIndex;
    /* t is t */
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
    
    return T;
}

/* Heat conduction in a cylinder */
double UpdateSubdomainCyl(struct Node1D *node, int Var)
{
    double T, Ta, Tb, Tc;
    int t = node->TimeIndex;
    double n = node->NodeNum; // Node number where n=0 at the center
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

    T = 1/M * ( (2*n+1)/(2*n) * Tc + (M-2) * Tb + (2*n-1)/(2*n) * Ta );
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
    double Cond, Dens, Cap, alpha, M, N;
    //double h;

    Tb = node->Value[Var][t];
    Text = T_ext(t*node->dt);

    /* Figure out which side the boundary is on */
    if(node->Next == NULL)
        Ta = node->Prev->Value[Var][t];
    else
        Ta = node->Next->Value[Var][t];

    //h = 3000; //W/(m^2 K)
    //h = find_val("HConv", node->varlst);
    printf("%g\n", h);

    Cond = k(Tb); //0.5985; // W/(m K)
    Dens = rho(Tb); //1000; // kg/m^3
    Cap = Cp(Tb); //1865; // J/(kg K)

    alpha = Cond/(Dens*Cap);
    M = (node->dx)*(node->dx)/(alpha*node->dt);
    N = h*node->dx/Cond;

    T = 1/M * (2*N*Text + (M-(2*N+2))*Tb + 2*Ta);

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

    return T;
}

double UpdateRadialSymmetryBoundary(struct Node1D *node, int Var)
{
    int t = node->TimeIndex;
    double Ta, Tb, T, Cond, Dens, Cap, alpha, M;

    Tb = node->Value[Var][t];

    if(node->Next == NULL)
        Ta = node->Prev->Value[Var][t];
    else
        Ta = node->Next->Value[Var][t];

    /* Calculate the thermal properties based on the temperature of the node at
    the pervious time step. */
    Cond = k(Tb); //0.5985; // W/(m K)
    Dens = rho(Tb); //1000; // kg/m^3
    Cap = Cp(Tb); //1865; // J/(kg K)

    alpha = Cond/(Dens*Cap);
    M = (node->dx)*(node->dx)/(alpha*node->dt);

    T = 4/M * Ta + (M-4)/M * Tb;

    return T;
}

double UpdateConvectiveBoundaryCyl(struct Node1D *node, int Var)
{
    double T, Text, Ta, Tb;
    int t = node->TimeIndex;
    int n = node->NodeNum;
    double Cond, Dens, Cap, alpha, M, N;
    //double h;

    Tb = node->Value[Var][t];
    Text = T_ext(t*node->dt);

    /* Figure out which side the boundary is on */
    if(node->Next == NULL)
        Ta = node->Prev->Value[Var][t];
    else
        Ta = node->Next->Value[Var][t];

    //h = 3000; //W/(m^2 K)
    //h = find_val("HConv", node->varlst);
    Cond = k(Tb); //0.5985; // W/(m K)
    Dens = rho(Tb); //1000; // kg/m^3
    Cap = Cp(Tb); //1865; // J/(kg K)

    alpha = Cond/(Dens*Cap);
    M = (node->dx)*(node->dx)/(alpha*node->dt);
    N = h*node->dx/Cond;

    T = n*N/((2*n-1)/2 + n*N) * Text + (2*n-1)/2/((2*n-1)/2 + n*N) * Ta;

    return T;
}
