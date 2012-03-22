#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "fd-solver.h"

struct Node1D* CreateNode1D(double dx, double dt, int NumVars, int NumValues)
{
    int i;
    struct Node1D *node;
    static int id = 0;
    node = (struct Node1D*) calloc(1, sizeof(struct Node1D));

    node->nid = id;
    id++;

    node->Update = (double (**)(struct Node1D*, int)) calloc(NumVars, sizeof(double (*)(struct Node1D*, int)));
    node->dx = dx;
    node->dt = dt;
    node->Type = SUBDOMAIN;
    node->Next = NULL;
    node->Prev = NULL;
    node->NextIsDep = 0;
    node->PrevIsDep = 0;
    node->Value = NULL;
    node->Value = (double**) calloc(NumVars, sizeof(double*));
    for(i=0; i<NumVars; i++) {
        node->Value[i] = (double*) calloc(NumValues, sizeof(double));
    }
    node->Vars = NULL;
    node->Vars = (char*) calloc(NumVars, sizeof(char));
    node->NumVars = NumVars;
    node->NumValues = NumValues;
    node->NodeNum = 0;

    node->TimeIndex = 0;

    return node;
}

void DestroyNode1D(struct Node1D *node)
{
    int i;
    if(node) {
        if(node->Value) {
            for(i=0; i<node->NumVars; i++)
                free(node->Value[i]);
        }
        free(node->Value);
        if(node->Vars) {
            free(node->Vars);
        }
        if(node->Update) {
            free(node->Update);
        }
        free(node);
    }
}

int UpdateNode(struct Node1D *node)
{
    int i;
    /* Return 1 if the node depends on values of adjacent nodes at the current
     * time step. */
    if(node->NextIsDep) {
        return 1;
    } else if(node->PrevIsDep) {
        return 1;
    }

    if(node->TimeIndex+1 == node->NumValues)
        return -1;

    for(i=0; i<node->NumVars; i++) {
        node->Value[i][node->TimeIndex+1] = node->Update[i](node, i);
    }
    node->TimeIndex++;
    
    return 0;
}

double NodeGetValue(struct Node1D *node, char Var, int TimeIndex)
{
    int i;
    if(TimeIndex >= node->NumValues)
        return 0;
    for(i=0; i<node->NumVars; i++) {
        if(Var == node->Vars[i])
            break;
    }
    return node->Value[i][TimeIndex];
}

void NodeSetUpdateFunction(struct Node1D *node, int VarNumber, char VarName, double (*f)(struct Node1D*, int))
{
    node->Update[VarNumber] = f;
    node->Vars[VarNumber] = VarName;
}

void NodeApplyInitialCondition(struct Node1D *node, int VarNumber, double value)
{
    node->Value[VarNumber][node->TimeIndex] = value;
}

