#include <stdlib.h>
#include <string.h>
#include "fd-solver.h"

struct Node1D* CreateNode1D(double dx, double dt, int NumValues)
{
    struct Node1D *node;
    static int id = 0;
    node = (struct Node1D*) calloc(1, sizeof(struct Node1D));

    node->nid = id;
    id++;

    node->Update = NULL;
    node->dx = dx;
    node->dt = dt;
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

int UpdateNode(struct Node1D *node)
{
    /* Return 1 if the node depends on values of adjacent nodes at the current
     * time step. */
    if(node->NextIsDep) {
        return 1;
    } else if(node->PrevIsDep) {
        return 1;
    }

    node->Update(node);
    node->TimeIndex++;
    
    return 0;
}

void NodeApplyInitialCondition(struct Node1D *node, double value)
{
    node->Value[node->TimeIndex] = value;
}
