#include <stdio.h>
#include "fd-solver.h"

int main(int argc, char *argv[])
{
    struct Domain1D *domain;
    struct Node1D *node;
    //int i;

    domain = CreateDomain1D("x", 6, .001, .01, 2, 5000);
    DomainApplyInitialCondition(domain, 0, 373.15);
    DomainApplyInitialCondition(domain, 1, 1.0);
    node = domain->Nodes[0];
    while(node) {
        NodeSetUpdateFunction(node, 0, 'T', &UpdateSubdomain);
        NodeSetUpdateFunction(node, 1, 'c', &UpdateSubdomainRxn);
        node = node->Next;
    }
    NodeSetUpdateFunction(domain->Nodes[0], 0, 'T', &UpdateInsulatedBoundary);
    NodeSetUpdateFunction(domain->Nodes[5], 0, 'T', &UpdateConvectiveBoundary);

    while(UpdateDomain(domain) != 1);

    node = domain->Nodes[0];
    //for(i = 0; i < node->TimeIndex; i++) {
    //    fprintf(stdout, "%3.0f, ", node->Value[0][i]);
    //}
    //fprintf(stdout, "\n\n");
    while(node) {
        //fprintf(stdout, "%d -> %f\n", node->nid, node->Value[0][3000]);
        
        fprintf(stdout, "%d -> %f\n", node->nid, NodeGetValue(node, 'c', 100));
        node = node->Next;
    }

    DestroyDomain1D(domain);

    return 0;
}

