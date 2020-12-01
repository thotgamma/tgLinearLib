#ifndef SOLVER_H
#define SOLVER_H

#include "linear.h"


void incomp_chol_decomp(const mat_csr_t* A, mat_csr_t* L, vector_t* d);
void solve_tri(vector_t* ans, const mat_csr_t* L, const mat_csr_t* U, const vector_t* b);
double cg_solver(vector_t* ans, const mat_csr_t* A, const vector_t* b, const double e, const int max_itr);
double iccg_solver(vector_t* ans, const mat_csr_t* A, const vector_t* b, const double e, const int max_itr);
double irls_solver(vector_t* ans, const mat_csr_t* A, const vector_t* b, const double p, const int max_itr);


#endif

