#include <stdio.h>
#include "fd-solver.h"

int main(int argc, char *argv[])
{
    struct Domain1D *domain;
    struct Node1D *node;
    int i;

    domain = CreateDomain1D("x", 6, .001, .01, 5000);
    DomainApplyInitialCondition(domain, 373.15);
    node = domain->Nodes[0];
    while(node) {
        node->Update = &UpdateSubdomain;
        node = node->Next;
    }
    domain->Nodes[0]->Update = &UpdateInsulatedBoundary;
    domain->Nodes[5]->Update = &UpdateConvectiveBoundary;

    while(UpdateDomain(domain) != 1);

    node = domain->Nodes[3];
    for(i = 0; i < node->TimeIndex; i++) {
        fprintf(stdout, "%3.0f, ", node->Value[i]);
    }
    fprintf(stdout, "\n\n");
    while(node) {
        fprintf(stdout, "%d -> %f\n", node->nid, node->Value[0]);
        node = node->Next;
    }

    DestroyDomain1D(domain);

    return 0;
}

