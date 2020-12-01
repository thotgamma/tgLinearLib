#ifndef UTILITY_H
#define UTILITY_H

#include <stdbool.h>


bool check_sorted_acend(int* target, int n);
int intersect_i(const int* A, const int nA, const int* B, const int nB, int* X);
void print_array_i(const int* array, int n);
void print_array_d(const double* array, int n);

#endif
