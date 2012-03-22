#ifndef FD_SOLVER_H
#define FD_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

/* Relavent data structures to make the simulation work. */

enum NodeType { BOUNDARY, SUBDOMAIN };

struct Node1D {
    /* Node ID (mostly for debugging purposes */
    int nid;
    
    /* Update function */
    double (**Update)(struct Node1D*, int);

    /* Node index in the domain (center = 0) */
    int NodeNum;

    /* Node width */
    double dx;

    /* Time step size */
    double dt;

    /* Type of node */
    enum NodeType Type;

    /* Node that's next in the positive X direction */
    struct Node1D *Next;

    /* Previous node */
    struct Node1D *Prev;

    /* True if the Next node requires a value for the current node a the current
     * time step
     */
    int NextIsDep;

    /* Same as above, only for the previous node */
    int PrevIsDep;

    /* List of values by time index. First value is the initial value. First
     * index determines which variable the values are for. */
    double **Value;

    /* Names of each of the variables */
    char *Vars;

    int NumVars;
    int NumValues;

    /* Current time index. */
    int TimeIndex;
};

struct Domain1D {
    /* Domain name (in case anyone wants to name their domains for fun. */
    char *Name;

    /* Number of nodes in the domain. */
    int NumNodes;

    /* Time step to be used when calculating stuff. */
    double dt;

    /* Simulation length */
    double NumTimes;

    /* Current time index. */
    int TimeIndex;

    /* An array of all the nodes in the domain. */
    struct Node1D **Nodes;
};

struct NodeStack {
    struct Node1D *Node;
    struct NodeStack *Next;
};

struct Node1D* CreateNode1D(double, double, int, int);
void DestroyNode1D(struct Node1D*);
struct Domain1D* CreateDomain1D(char*, int, double, double, int, int);
void DestroyDomain1D(struct Domain1D*);

int UpdateNode(struct Node1D*);
int UpdateDomain(struct Domain1D*);
void NodeApplyInitialCondition(struct Node1D*, int, double);
void DomainApplyInitialCondition(struct Domain1D*, int, double);
void NodeSetUpdateFunction(struct Node1D*, int, char, double (*)(struct Node1D*, int));
double NodeGetValue(struct Node1D*, char, int);
struct NodeStack* Push(struct NodeStack*, struct Node1D*);
struct NodeStack* Pop(struct NodeStack*);

double UpdateTest(struct Node1D*, int);
double UpdateSubdomain(struct Node1D*, int);
double UpdateSubdomainCyl(struct Node1D*, int);
double UpdateInsulatedBoundary(struct Node1D*, int);
double UpdateConvectiveBoundary(struct Node1D*, int);
double UpdateRadialSymmetryBoundary(struct Node1D*, int);
double UpdateConvectiveBoundaryCyl(struct Node1D*, int);
double UpdateSubdomainRxn1(struct Node1D*, int);
double UpdateSubdomainRxn2(struct Node1D*, int);

void CSVOutFixedNode(struct Node1D*, char*);
void CSVOutFixedTime(struct Domain1D*, int, char*);

#ifdef __cplusplus
}
#endif

#endif

