
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <stdbool.h>

#include "utility.h"
#include "minitrace/minitrace.h"


bool check_sorted_acend(int* target, int n) {
	for (int i = 0; i < n-1; i++) {
		if (!(target[i] <= target[i+1])) return false;
	}
	return true;
}

int intersect_i(const int* A, const int nA, const int* B, const int nB, int* X) {
	//MTR_BEGIN("utility", "intersect_i");
	/*
	int n = (nA > nB) ? nA : nB;
	*X = (int*)malloc(sizeof(int) * n);
	*/

	int itr_x = 0;

	for (int i = 0; i < nA; i++){
		for (int j = 0; j < nB; j++) {

			if (B[j] > A[i]) continue; //AとBはそれぞれソートされているという前提

			if (A[i] == B[j]) {
				X[itr_x++] = A[i];
			}

		}
	}
	//MTR_END("utility", "intersect_i");
	return itr_x;

}

void print_array_i(const int* array, int n) {
	for (int i = 0; i < n; i++) {
		printf("%d, ", array[i]);
	}
	printf("\n");
}

void print_array_d(const double* array, int n) {
	for (int i = 0; i < n; i++) {
		printf("%f, ", array[i]);
	}
	printf("\n");
}
