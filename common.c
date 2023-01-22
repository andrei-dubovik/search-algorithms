#include <stdio.h>
#include "common.h"

const int M = 10;
const int N = 2*M;
const int T = 320;
const int tree_size = (T + 1)*(M*T + 1);


/** Compute the objective function at a given point in the tree */
double objf(const int i, const int j, const int uj) {
    const double dt = 4.0/T;
    const double du = 1.0/M;
    double u = uj*du - 1.0;
    double x = -i*dt + j*dt*du;
    return (u*u + x)*dt;
}


/** Compute a flat index for a tree node */
int idx(int i, int j) {
    return i*(M*i - M + 1) + j;
}


/** Traverse the optimal path */
void dump(const int *path, const char *name) {
    const double dt = 4.0/T;
    const double du = 1.0/M;

    FILE *file = fopen(name, "w");
    fprintf(file, "t u x v\n");

    int j = 0;
    double x = .0;
    double v = .0;
    for (int i = 0; i <= T; i++) {
        double t = i*dt;
        double x = (j - i*M)*dt*du;
        int uj = path[i];
        double u = uj*du - 1.0;
        v += objf(i, j, uj);
        fprintf(file, "%f %f %f %f\n", t, u, x, v);
        j += uj;
    }

    fclose(file);
    printf("Dumped solution to \"%s\"\n", name);
}


#ifdef TRACK_VISITS
void print_visits(int *stats) {
    int count = 0;
    long total = 0;
    for(int i = 0; i < tree_size; i++) {
        count += stats[i] > 0;
        total += stats[i];
    }
    printf("Visited states: %'i (%.1f\%)\n", count, 100.0*count/tree_size);
    printf("Total visits: %'li (%.1f per visited state)\n", total, (double)total/count);
}
#endif
