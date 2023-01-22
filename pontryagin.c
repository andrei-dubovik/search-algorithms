#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "common.h"
#include "monitor.h"


/** Convergence threshold on lambda */
const double eps = 1e-8;

#ifdef TRACK_VISITS
int *visit_stats;
#endif


/** Swap two numbers */
void swap(double *x, double *y) {
    double z = *x;
    *x = *y;
    *y = z;
}


/** Compute the instantaneous objective function */
double F(const double t, const double y, const double u) {
    return -(u*u + y);
}


/** Compute dy (motion equation) */
double f(const double t, const double y, const double u) {
    return u;
}


/** Compute the Hamiltonian */
double H(const double t, const double y, const double u, const double lambda) {
    return F(t, y, u) + lambda*f(t, y, u);
}


/** Compute the derivative of H w.r.t. y */
double H_y(const double t, const double y, const double u, const double lambda) {
    return -1.0;
}


/** Maximize the Hamiltonian in the control variable */
double max_H(const double t, const double y, const double lambda) {
    double maxh = -INFINITY;
    double u_star = 0.0;
    for(int k = 0; k < N + 1; k++) {
        double u = 2.0*k/N - 1.0;
        double h = H(t, y, u, lambda);
        if(h > maxh) {
            maxh = h;
            u_star = u;
        }
    }
    return u_star;
}


/** Compute a control path that satisfies all the necessary conditions
    exluding the transversality condition */
double draw_path(double lambda, double *v, int **path) {
    double y = 0.0;
    const double dt = 4.0/T;
    const double du = 1.0/M;

    if(v != NULL) *v = 0.0;
    if(path != NULL) *path = calloc(T + 1, sizeof(int));

    for(int i = 0; i < T + 1; i++) {
        #ifdef TRACK_VISITS
        int j = lround((y + i*dt)/(dt*du));
        if(j < 0) j = 0;
        if(j > 2*M*i) j = 2*M*i;
        visit_stats[idx(i, j)]++;
        #endif
        double t = i*dt;
        double u = max_H(t, y, lambda);
        double dy = f(t, y, u);
        double dl = -H_y(t, y, u, lambda);
        if(v != NULL) *v += -F(t, y, u)*dt;
        if(path != NULL) {
            int uj = lround(u/du) + M;
            if(uj < 0) uj = 0;
            if(uj > 2*M) uj = 2*M;
            (*path)[i] = uj;
        }
        if(i < T) {
            y += dy*dt;
            lambda += dl*dt;
        }
    }
    return lambda;
}


/** Find the path satisfying lambda(T) = 0 using the bisection method */
double solve_lambda_T(double l0a, double l0b) {
    double lTa = draw_path(l0a, NULL, NULL);
    double lTb = draw_path(l0b, NULL, NULL);
    if(lTb < lTa) {
        swap(&l0a, &l0b);
        swap(&lTa, &lTb);
    }
    if(lTa > 0 || lTb < 0) {
        printf("Error: lambda_T_0 and lambda_T_1 have the same sign.\n");
        exit(1);
    }

    int count = 0;
    while(lTb - lTa > eps) {
        count++;
        double lm = (l0a + l0b)/2.0;
        // Alternatively, a secant can be used but in this case it finds a
        // solution too fast to make for an interesting visualization
        // double lm = (l0a + l0b - (l0b - l0a)/(lTb - lTa)*(lTa + lTb))/2.0;
        double lTm = draw_path(lm, NULL, NULL);
        if(lTm < 0) {
            l0a = lm;
            lTa = lTm;
        } else {
            l0b = lm;
            lTb = lTm;
        }
    }

    printf("Convergence achieved in %i steps.\n", count);
    return (l0a + l0b)/2.0;
}


/** Solve the problem using the Pontryagin's maximum principle */
void pontryagin(double *value, int **path) {
    #ifdef TRACK_VISITS
    visit_stats = calloc(tree_size, sizeof(int));
    #endif

    // Find the path satisfying the transversality conditions
    double l0 = solve_lambda_T(-10.0, +10.0);

    // Report various statistics
    #ifdef TRACK_VISITS
    print_visits(visit_stats);
    #endif

    // Save the solution
    draw_path(l0, value, path);

    #ifdef TRACK_VISITS
    free(visit_stats);
    #endif
}


int main() {
    setlocale(LC_NUMERIC, "");

    // Solve the problem
    double V, timer;
    int *path;
    TIMEIT(pontryagin(&V, &path), timer);

    // Save the solution
    dump(path, "pontryagin.dat");
    free(path);

    // Report the optimal objective function
    printf("V = %f\n", V);

    // Report resource usage
    printf("Running time: %.2f s\n", timer);
    printf("Peak dynamic memory usage: %'i bytes\n", peak_memory);

    return 0;
}
