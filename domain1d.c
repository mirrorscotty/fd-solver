#include <stdlib.h>
#include <string.h>
#include "fd-solver.h"

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
        domain->Nodes[i] = CreateNode1D(dx, dt, times);
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

int UpdateDomain(struct Domain1D *domain)
{
    int i;
    struct NodeStack *UpdateLater;
    struct Node1D *node;

    UpdateLater = NULL;

    if(domain->TimeIndex == domain->NumTimes) {
        return 1;
    }
    
    for(i = 0; i < domain->NumNodes; i++) {
        if(UpdateNode(domain->Nodes[i]) == 1)
            UpdateLater = Push(UpdateLater, domain->Nodes[i]);
    }

    while(UpdateLater) {
        node = UpdateLater->Node;
        node->Update(node);
        node->TimeIndex++;
        UpdateLater = Pop(UpdateLater);
    }

    domain->TimeIndex++;

    return 0;
}

void DomainApplyInitialCondition(struct Domain1D *domain, double value)
{
    int i;
    for(i = 0; i < domain->NumNodes; i++) {
        NodeApplyInitialCondition(domain->Nodes[i], value);
    }
}

struct NodeStack* Push(struct NodeStack *stack, struct Node1D *node)
{
    struct NodeStack *tmp;
    tmp = (struct NodeStack*) calloc(1, sizeof(struct NodeStack));
    stack->Node = node;
    stack->Next = stack;
    return tmp;
}

struct NodeStack* Pop(struct NodeStack *stack)
{
    struct NodeStack *tmp;
    tmp = stack->Next;
    free(stack);
    return tmp;
}
