/* Relavent data structures to make the simulation work. */

enum NodeType { BOUNDARY, SUBDOMAIN };

struct Node1D {
    /* Node ID (mostly for debugging purposes */
    int nid;
    
    /* Update function */
    double (*Update)(struct Node1D*);

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

    /* List of values by time index. First value is the initial value. */
    double *Value;

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

struct Node1D* CreateNode1D(double, double, int);
void DestroyNode1D(struct Node1D*);
struct Domain1D* CreateDomain1D(char*, int, double, double, int);
void DestroyDomain1D(struct Domain1D*);

int UpdateNode(struct Node1D*);
int UpdateDomain(struct Domain1D*);
void NodeApplyInitialCondition(struct Node1D*, double);
void DomainApplyInitialCondition(struct Domain1D*, double);
struct NodeStack* Push(struct NodeStack*, struct Node1D*);
struct NodeStack* Pop(struct NodeStack*);

double UpdateTest(struct Node1D*);
double UpdateSubdomain(struct Node1D*);
double UpdateInsulatedBoundary(struct Node1D*);
double UpdateConvectiveBoundary(struct Node1D*);

