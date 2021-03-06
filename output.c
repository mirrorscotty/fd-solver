#include <stdio.h>

#include "fd-solver.h"
#include "can.h"

void CSVOutFixedNode(struct Node1D *node, char *filename)
{
    int i, j;
    FILE *fp;
    fp = fopen(filename, "w+");

    if(!fp) {
        fprintf(stderr, "Unable to open file for writing.\n");
        return;
    }

    /* Print out a column for time. */
    fprintf(fp, "t,");

    /* Print out the variable names */
    for(i=0; i<node->NumVars; i++) {
        fprintf(fp, "%c", node->Vars[i]);
        //if(i==node->NumVars-1)
        //    fprintf(fp, "\n");
        //else
            fprintf(fp, ",");
    }
    fprintf(fp, "%s,%s,%s\n", "rho", "Cp", "k");

    /* Print out the values */
    for(i=0; i<node->NumValues; i++) {
        fprintf(fp, "%g,", node->dt*i);
        for(j=0; j<node->NumVars; j++) {
            fprintf(fp, "%g", node->Value[j][i]);
            //if(j==node->NumVars-1)
            //    fprintf(fp, "\n");
            //else
                fprintf(fp, ",");
        }
        // Hack to print out some extra stuff that wasn't explicitely solved for
        fprintf(fp, "%g,%g,%g\n", rho(node->Value[0][i]),
                                    Cp(node->Value[0][i]),
                                    k(node->Value[0][i]));
    }

    /* Put an extra new line at the end of the file for fun. */
    fprintf(fp, "\n");

    fclose(fp);

    return;
}

void CSVOutFixedTime(struct Domain1D *domain, int tindex, char *filename)
{
    int i, j;
    FILE *fp;
    struct Node1D *n;

    fp = fopen(filename, "w+");

    if(!fp) {
        fprintf(stderr, "Unable to open file for writing.\n");
        return;
    }

    n = *(domain->Nodes);

    fprintf(fp, "%s,", domain->Name);
    for(i=0; i<n->NumVars; i++) {
        fprintf(fp, "%c", n->Vars[i]);
        //if(i==n->NumVars-1)
        //    fprintf(fp, "\n");
        //else
            fprintf(fp, ",");
    }
    fprintf(fp, "%s,%s,%s\n", "rho", "Cp", "k");

    for(j=0; j<domain->NumNodes; j++) {
        n = domain->Nodes[j];
        fprintf(fp, "%g,", j*n->dx);
        for(i=0; i<n->NumVars; i++) {
            fprintf(fp, "%g", n->Value[i][tindex]);
            //if(i==n->NumVars-1)
            //    fprintf(fp, "\n");
            //else
                fprintf(fp, ",");
        }
        // Same as above.
        fprintf(fp, "%g,%g,%g\n", rho(n->Value[0][tindex]),
                                    Cp(n->Value[0][tindex]),
                                    k(n->Value[0][tindex]));
    }

    fprintf(fp, "\n");
    fclose(fp);

    return;
}
