#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fd-solver.h"

struct Node1D* CreateNode1D(double dx, int NumValues)
{
    struct Node1D *node;
    static int id = 0;
    node = (struct Node1D*) calloc(1, sizeof(struct Node1D));

    node->nid = id;
    id++;

    node->Update = NULL;
    node->dx = dx;
    node->Type = SUBDOMAIN;
    node->Next = NULL;
    node->Prev = NULL;
    node->NextIsDep = 0;
    node->PrevIsDep = 0;
    node->Value = NULL;
    node->Value = (double*) calloc(NumValues, sizeof(double));
    node->TimeIndex = 0;

    return node;
}

void DestroyNode1D(struct Node1D *node)
{
    if(node) {
        if(node->Value)
            free(node->Value);
        free(node);
    }
}

struct Domain1D* CreateDomain1D(char *name,
        int width,
        double dx,
        double dt,
        int times)
{
    struct Domain1D *domain;
    struct Node1D *PrevNode;
    int i;

    PrevNode = NULL;

    domain = (struct Domain1D*) calloc(1, sizeof(struct Domain1D));
    
    domain->Name = (char*) calloc(80, sizeof(char));
    strncpy(domain->Name, name, 80);

    domain->NumNodes = width;

    domain->dt = dt;

    domain->Nodes = (struct Node1D**) calloc(times, sizeof(struct Node1D*));
    for(i=0; i<width; i++) {
        domain->Nodes[i] = CreateNode1D(dx, times);
        domain->Nodes[i]->Prev = PrevNode;
        if(PrevNode)
            PrevNode->Next = domain->Nodes[i];
        PrevNode = domain->Nodes[i];
    }
    domain->NumTimes = times;

    return domain;
}

void DestroyDomain1D(struct Domain1D *domain)
{
    int i;
    if(domain) {
        if(domain->Name)
            free(domain->Name);
        if(domain->Nodes) {
            for(i=0; i<domain->NumNodes; i++)
                DestroyNode1D(domain->Nodes[i]);
            free(domain->Nodes);
        }
        free(domain);
    }
}

int main(int argc, char *argv[])
{
    struct Domain1D *domain;
    struct Node1D *node;

    fprintf(stdout, "BLAH!\n");
    
    domain = CreateDomain1D("x", 5, .5, .1, 500);
    
    node = domain->Nodes[0];
    while(node) {
        fprintf(stdout, "%d -> %f\n", node->nid, node->Value[0]);
        node = node->Next;
    }

    DestroyDomain1D(domain);

    return 0;
}

