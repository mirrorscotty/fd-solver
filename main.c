#include <stdio.h>
#include "fd-solver.h"
#include "can.h"

/* Loaded from the data file. */
extern double NNodes, Deltax, Deltat, NTimeSteps;
extern double To;

int main(int argc, char *argv[])
{
    struct Domain1D *domain;
    struct Node1D *node;
    //int i;

    /* Load in the material property parameters from the data file */
    initialize_variables();
    init("can_data.dat");

    domain = CreateDomain1D("x", (int) NNodes, Deltax, Deltat, 3, (int) NTimeSteps);
    DomainApplyInitialCondition(domain, 0, To);
    DomainApplyInitialCondition(domain, 1, 1.0);
    DomainApplyInitialCondition(domain, 2, 1.0);
    node = domain->Nodes[0];
    while(node) {
        NodeSetUpdateFunction(node, 0, 'T', &UpdateSubdomain);
        NodeSetUpdateFunction(node, 1, 'c', &UpdateSubdomainRxn1);
        NodeSetUpdateFunction(node, 2, 'd', &UpdateSubdomainRxn2);
        node = node->Next;
    }
    NodeSetUpdateFunction(domain->Nodes[0], 0, 'T', &UpdateInsulatedBoundary);
    NodeSetUpdateFunction(domain->Nodes[5], 0, 'T', &UpdateConvectiveBoundary);

    while(UpdateDomain(domain) != 1);

    node = domain->Nodes[3];
    //for(i = 0; i < node->TimeIndex; i++) {
    //    fprintf(stdout, "%3.0f, ", node->Value[0][i]);
    //}
    //fprintf(stdout, "\n\n");
    //while(node) {
        //fprintf(stdout, "%d -> %f\n", node->nid, node->Value[0][3000]);
        
        //fprintf(stdout, "%d -> %f\n", node->nid, NodeGetValue(node, 'd', 100));
        //node = node->Next;
    //}
    //
    CSVOutFixedTime(domain, 4000, "test.csv");

    DestroyDomain1D(domain);

    return 0;
}

