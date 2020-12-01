
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include "utility.h"
#include "linear.h"
#include "minitrace/minitrace.h"


mat_dense_t* mat_dense(int cols, int rows, double* data) {
	mat_dense_t* dense = (mat_dense_t*)malloc(sizeof(mat_dense_t));
	dense->cols = cols;
	dense->rows = rows;
	dense->data = data;

	return dense;
}

mat_dense_t* new_mat_dense(int cols, int rows){
	mat_dense_t* dense = (mat_dense_t*)malloc(sizeof(mat_dense_t));
	dense->cols = cols;
	dense->rows = rows;
	dense->data = (double*)malloc(sizeof(double) * cols * rows);
	return dense;
}

void destroy_mat_dense(mat_dense_t* mat) {
	free(mat->data);
	free(mat);
}

mat_dense_t* create_dense_from_coo(const mat_coo_t* input) {
	mat_csr_t* csr = create_csr_from_coo(input);
	mat_dense_t* dense = create_dense_from_csr(csr);
	destroy_mat_csr(csr);
	return dense;
}

mat_dense_t* create_dense_from_csr(const mat_csr_t* input) {
	int cols     = input->cols;
	int rows     = input->rows;
	double* data = input->data;
	int* colind  = input->colind;
	int* rowind  = input->rowind;

	mat_dense_t* dense = new_mat_dense(cols, rows);

	int whole_itr = 0;
	int data_itr = 0;
	int row_itr = 1;
	int current_row = rowind[0];

	for (int j = 0; j < rows; j++) {
		for (int i = 0; i < cols; i++) {
			if (i == colind[data_itr] && current_row == j) {
				dense->data[whole_itr] = data[data_itr];
				data_itr++;
			} else {
				dense->data[whole_itr] = 0.0;
			}

			if (data_itr == rowind[row_itr]){
				current_row++;
				row_itr++;
			}

			whole_itr++;
		}
	}

	return dense;

}

mat_coo_t* mat_coo(int cols, int rows, int nnz, double* data, int* colind, int* rowind) {
	mat_coo_t* coo = (mat_coo_t*)malloc(sizeof(mat_coo_t));
	coo->cols = cols;
	coo->rows = rows;
	coo->nnz = nnz;
	coo->data = data;
	coo->colind = colind;
	coo->rowind = rowind;

	detect_coo_axis(coo);

	return coo;
}

mat_coo_t* new_mat_coo(int cols, int rows, int nnz) {
	mat_coo_t* coo = (mat_coo_t*)malloc(sizeof(mat_coo_t));
	coo->cols = cols;
	coo->rows = rows;
	coo->nnz = nnz;
	coo->data = (double*)malloc(sizeof(double)*nnz);
	coo->colind = (int*)malloc(sizeof(int)*nnz);
	coo->rowind = (int*)malloc(sizeof(int)*nnz);

	return coo;
}

void detect_coo_axis(mat_coo_t* mat) {
	if (check_sorted_acend(mat->colind, mat->nnz)) {
		mat->sorted_axis = COO_COL_SORTED;
		return;
	}
	if (check_sorted_acend(mat->rowind, mat->nnz)) {
		mat->sorted_axis = COO_ROW_SORTED;
		return;
	}

	mat->sorted_axis = COO_RANDOM;
	return;
}

void destroy_mat_coo(mat_coo_t* mat) {
	free(mat->data);
	free(mat->colind);
	free(mat->rowind);
	free(mat);
}

mat_csr_t* mat_csr(int cols, int rows, double* data, int* colind, int* rowind) {
	mat_csr_t* csr = (mat_csr_t*)malloc(sizeof(mat_csr_t));
	csr->cols      = cols;
	csr->rows      = rows;
	csr->data      = data;
	csr->colind    = colind;
	csr->rowind    = rowind;

	return csr;
}

mat_csr_t* new_mat_csr(int cols, int rows, int nnz) {
	mat_csr_t* csr = (mat_csr_t*)malloc(sizeof(mat_csr_t));
	csr->cols = cols;
	csr->rows = rows;
	csr->data = (double*)malloc(sizeof(double) * nnz);
	csr->colind = (int*)malloc(sizeof(int) * nnz);
	csr->rowind = (int*)malloc(sizeof(int) * (cols + 1));

	return csr;
}

void destroy_mat_csr(mat_csr_t* mat) {
	free(mat->data);
	free(mat->colind);
	free(mat->rowind);
	free(mat);
}

mat_csc_t* mat_csc(int cols, int rows, double* data, int* colind, int* rowind) {
	mat_csc_t* csc = (mat_csc_t*)malloc(sizeof(mat_csc_t));
	csc->cols      = cols;
	csc->rows      = rows;
	csc->data      = data;
	csc->colind    = colind;
	csc->rowind    = rowind;

	return csc;
}

mat_csc_t* new_mat_csc(int cols, int rows, int nnz) {
	mat_csc_t* csc = (mat_csc_t*)malloc(sizeof(mat_csc_t));
	csc->cols = cols;
	csc->rows = rows;
	csc->data = (double*)malloc(sizeof(double) * nnz);
	csc->colind = (int*)malloc(sizeof(int) * (rows + 1));
	csc->rowind = (int*)malloc(sizeof(int) * nnz);

	return csc;
}

void destroy_mat_csc(mat_csc_t* mat) {
	free(mat->data);
	free(mat->colind);
	free(mat->rowind);
	free(mat);
}

vector_t* vector(int size, double* data) {
	vector_t* vec = (vector_t*)malloc(sizeof(vector_t));
	vec->size  = size;
	vec->data = data;

	return vec;
}

vector_t* new_vector(int size) {
	vector_t* vec = (vector_t*)malloc(sizeof(vector_t));
	double* data = (double*)calloc(size, sizeof(double));

	vec->size= size;
	vec->data = data;

	return vec;
}

void destroy_vector(vector_t* vec) {
	free(vec->data);
	free(vec);
}

// -------------------------------------------------------------

void mul_mat_vec(vector_t* ans, const mat_csr_t* A, const vector_t* x) {
	MTR_BEGIN_FUNC();
	for (int j = 0; j < A->cols; j++) {
		ans->data[j] = 0;
		for (int i = A->rowind[j]; i < A->rowind[j + 1]; i++) {
			ans->data[j] += A->data[i] * x->data[A->colind[i]];
		}
	}
	MTR_END_FUNC();
}

void mul_double_vec(vector_t* ans, const double a, const vector_t* x) {
	MTR_BEGIN_FUNC();
	for (int i = 0; i < x->size; i++) {
		ans->data[i] = a * x->data[i];
	}
	MTR_END_FUNC();
}


double dot_vec_vec(const vector_t* x, const vector_t* y) {
	MTR_BEGIN_FUNC();
	double ans = 0;

	for (int i = 0; i < x->size; i++) {
		ans += x->data[i] * y->data[i];
	}
	MTR_END_FUNC();
	return ans;
}


void add_vec_vec(vector_t* ans, const vector_t* x, const vector_t* y) {
	MTR_BEGIN_FUNC();
	assert(x->size == y->size);

	for (int i = 0; i < x->size; i++) {
		ans->data[i] = x->data[i] + y->data[i];
	}
	MTR_END_FUNC();
}


void sub_vec_vec(vector_t* ans, const vector_t* x, const vector_t* y) {
	MTR_BEGIN_FUNC();
	assert(x->size == y->size);

	for (int i = 0; i < x->size; i++) {
		ans->data[i] = x->data[i] - y->data[i];
	}
	MTR_END_FUNC();
}

void div_vec_double(vector_t* ans, const vector_t* x, const double a) {
	MTR_BEGIN_FUNC();
	for (int i = 0; i < x->size; i++) {
		ans->data[i] = x->data[i]/a;
	}
	MTR_END_FUNC();
}

void pow_vec_double(vector_t* ans, const vector_t* x, const double a) {
	MTR_BEGIN_FUNC();
	for (int i = 0; i < x->size; i++) {
		ans->data[i] = pow(x->data[i], a);
	}
	MTR_END_FUNC();
}

bool eq_vector(const vector_t* x, const vector_t* y) {
	if (x->size != y->size) return false;

	for (int i = 0; i < x->size; i++) {
		if (!(fabs(x->data[i] - y->data[i]) < 10e-10)) return false;
	}

	return true;
}


vector_t* clone_vec(const vector_t* x) {
	MTR_BEGIN_FUNC();
	vector_t* ans = new_vector(x->size);

	for (int i = 0; i < x->size; i++) {
		ans->data[i] = x->data[i];
	}
	MTR_END_FUNC();
	return ans;
}


void copy_vec(const vector_t* x, vector_t* y) {
	MTR_BEGIN_FUNC();
	assert(x->size == y->size);

	for (int i = 0; i < x->size; i++) {
		x->data[i] = y->data[i];
	}
	MTR_END_FUNC();
}


double norm_vec(const vector_t* x) {
	MTR_BEGIN_FUNC();
	double ans = 0;

	for (int i = 0; i < x->size; i++) {
		ans += x->data[i] * x->data[i];
	}
	ans = sqrt(ans);
	MTR_END_FUNC();
	return ans;
}

bool eq_mat_csr(const mat_csr_t* A_csr, const mat_csr_t* B_csr) {

	mat_dense_t* A = create_dense_from_csr(A_csr);
	mat_dense_t* B = create_dense_from_csr(B_csr);

	bool result = true;

	if (A->cols != B->cols) result = false;
	if (A->rows != B->rows) result = false;

	int elems = A->cols * A->rows;

	for (int i = 0; i < elems; i++) {
		if (!(fabs(A->data[i] - B->data[i]) < 10e-10)) result = false;
	}

	destroy_mat_dense(A);
	destroy_mat_dense(B);

	return result;
}

int nnz_csr(const mat_csr_t* A) {
	return A->rowind[A->cols];
}

int nnz_csc(const mat_csc_t* A) {
	return A->colind[A->rows];
}

void abs_vec(vector_t* ans, const vector_t* x) {
	MTR_BEGIN_FUNC();
	ans->size = x->size;

	for (int i = 0; i < x->size; i++) {
		ans->data[i] = fabs(x->data[i]);
	}
	MTR_END_FUNC();
}

double sum_vec(const vector_t* x) {
	MTR_BEGIN_FUNC();
	double ans = 0;

	for (int i = 0; i < x->size; i++) {
		ans += x->data[i];
	}
	MTR_END_FUNC();
	return ans;
}


void print_vector(const vector_t* vec) {
	MTR_BEGIN_FUNC();
	printf("--- vector(%d) ---\n", vec->size);

	for (int i = 0; i < vec->size; i++) {
		printf("%f, ", vec->data[i]);
	}
	printf("\n");
	printf("-----------------\n");
	MTR_END_FUNC();
}


//TODO doubleに対して==0を使わない
mat_csr_t* create_csr_from_dense(const mat_dense_t* input) {
	MTR_BEGIN_FUNC();
	int cols     = input->cols;
	int rows     = input->rows;
	double* data = input->data;


	//STEP1. count nnz

	int ind = 0;
	int nnz = 0;

	for (int j = 0; j < cols; j++) {
		for (int i = 0; i < rows; i++) {

			if (data[ind] != 0) {
				nnz++;
			}

			ind++;
		}
	}

	//STEP2. create csr

	mat_csr_t* csr = (mat_csr_t*)malloc(sizeof(mat_csr_t));
	csr->cols      = cols;
	csr->rows      = rows;
	csr->data      = (double*)malloc(sizeof(double) * nnz);
	csr->colind    = (int*)malloc(sizeof(int) * nnz);
	csr->rowind    = (int*)malloc(sizeof(int) * (rows + 1));


	//STEP3. copy data

	ind = 0;
	int data_itr = 0;
	int row_itr= 0;
	int last_row = 0;
	bool first = true;

	for (int j = 0; j < rows; j++) {
		for (int i = 0; i < cols; i++) {

			if (data[ind] != 0) {
				csr->data[data_itr] = data[ind];
				csr->colind[data_itr] = i;

				if (first == true) {
					first = false;
					last_row = j;
					csr->rowind[0] = j;
					row_itr++;
				} else if (last_row != j) {
					last_row = j;
					csr->rowind[row_itr] = data_itr;
					row_itr++;
				}

				data_itr++;

			}

			ind++;
		}
	}

	csr->rowind[row_itr] = data_itr;
	MTR_END_FUNC();
	return csr;
}


void print_mat_csr(const mat_csr_t* input) {
	MTR_BEGIN_FUNC();
	int cols     = input->cols;
	int rows     = input->rows;
	double* data = input->data;
	int* colind  = input->colind;
	int* rowind  = input->rowind;

	printf("-- matrix(%dx%d) --\n", cols, rows);

	int data_itr = 0;
	int row_itr = 1;
	int current_row = rowind[0];

	for (int j = 0; j < rows; j++) {
		for (int i = 0; i < cols; i++) {
			if (i == j) printf("\x1b[33m");
			if (i == colind[data_itr] && current_row == j) {
				printf("%f\t", data[data_itr]);
				data_itr++;
			} else {
				printf("%f\t", 0.0);
			}

			if (data_itr == rowind[row_itr]){
				current_row++;
				row_itr++;
			}
			if (i == j) printf("\x1b[39m");
		}
		printf("\n");
	}

	printf("-----------------\n");
	MTR_END_FUNC();
}

void swap_int(int *x, int *y) {
	/*
	*x ^= *y;
	*y = *x ^ *y;
	*x ^= *y;
	*/
	int tmp = *x;
	*x = *y;
	*y = tmp;
}

void swap_double(double *x, double *y) {
	double tmp = *x;
	*x = *y;
	*y = tmp;
}


int min(int x, int y) {
	return x > y ? y : x;
}

bool compare(int ref, int subref, int axis, int sub) {
	if (ref == axis) return sub > subref;
	return axis > ref;
}

void insert_sort(int* axis, int left, int right, double* data) {

	for (int i = left +1; i <= right; i++) {

		int ref = axis[i];
		double tmp = data[i];

		if (ref < axis[i-1]) {
			int j = i;
			do {
				axis[j] = axis[j-1];
				data[j] = data[j-1];
				--j;
			} while ( j > left && ref < axis[j-1]);
			axis[j] = ref;
			data[j] = tmp;
		}
	}
}


void insert_sort_with_sub(int* axis, int left, int right, int* sub, double* data) {

	for (int i = left +1; i <= right; i++) {

		int ref = axis[i];
		int subref = sub[i];
		double tmp = data[i];

		if (compare(ref, subref, axis[i-1], sub[i-1])) {
			int j = i;
			do {
				axis[j] = axis[j-1];
				sub[j] = sub[j-1];
				data[j] = data[j-1];
				--j;
			} while ( j > left && compare(ref, subref, axis[j-1], sub[j-1]));
			axis[j] = ref;
			sub[j] = subref;
			data[j] = tmp;
		}
	}
}

void quick_sort(int* axis, int left, int right, double* data) {

	if (left >= right) return;

	if ((right - left) < _QUICK_SORT_LIMIT_) {
		insert_sort(axis, left, right, data);
		return;
	}

	int i, j, pivot;
	i = left;
	j = right + 1;
	pivot = left;

	while (i < j) {

		i++;
		while (axis[i] < axis[pivot]) i++;

		j--;
		while (axis[pivot] < axis[j]) j--;

		if (i < j) {
			swap_int(&axis[i], &axis[j]);
			swap_double(&data[i], &data[j]);
		}

	}

	swap_int(&axis[pivot], &axis[j]);
	swap_double(&data[pivot], &data[j]);

	quick_sort(axis, left, j-1, data);
	quick_sort(axis, j+1, right, data);


}

void mkey_quick_sort(int* axis, int left, int right, int* sub, double* data) {

	if (left >= right) return;

	if ((right - left) < _QUICK_SORT_LIMIT_) {
		insert_sort_with_sub(axis, left, right, sub, data);
		return;
	}

	int pivot = axis[left];

	int i = left;
	int m = left;

	int j = right;
	int n = right;


	while (true) {

		while (i <= j) {
			if (pivot < axis[i]) break;
			if (axis[i] == pivot) {
				swap_int(&axis[i], &axis[m]);
				swap_int(&sub[i], &sub[m]);
				swap_double(&data[i], &data[m]);
				m++;
			}
			i++;
		}

		while (i <= j) {
			if (axis[j] < pivot) break;
			if (axis[j] == pivot) {
				swap_int(&axis[j], &axis[n]);
				swap_int(&sub[j], &sub[n]);
				swap_double(&data[j], &data[n]);
				n--;
			}
			j--;
		}

		if (i > j) break;

		swap_int(&axis[i], &axis[j]);
		swap_int(&sub[i], &sub[j]);
		swap_double(&data[i], &data[j]);
		i++;
		j--;

	}

	int trans_left = min(m - left, i - m);
	for (int k = 0; k < trans_left; k++) {
		swap_int(&axis[left+k], &axis[j-k]);
		swap_int(&sub[left+k], &sub[j-k]);
		swap_double(&data[left+k], &data[j-k]);
	}

	int trans_right = min(right - n, n - j);
	for (int k = 0; k < trans_right; k++) {
		swap_int(&axis[i+k], &axis[right-k]);
		swap_int(&sub[i+k], &sub[right-k]);
		swap_double(&data[i+k], &data[right-k]);
	}

	int mid_left  = left  + (i-m);
	int mid_right = right - (n-j);

	mkey_quick_sort(axis, left, mid_left-1, sub, data); //left
	quick_sort(sub, mid_left, mid_right, data); //mid
	mkey_quick_sort(axis, mid_right+1, right, sub, data); //right

}


void sort_for_spmat(int* axis, int left, int right, int* sub, double* data) {
	MTR_BEGIN_FUNC();
	mkey_quick_sort(axis, left, right, sub, data);
	MTR_END_FUNC();
}


int compress_array(int* target, int* sub, double* data, int length) {
	MTR_BEGIN_FUNC();
	int lastvalue = target[0];
	int idx = 1;

	for (int i = 1; i < length; i++) {
		if (lastvalue != target[i]) {
			target[idx++] = i;
			lastvalue = target[i];
		}
	}

	target[idx] = length;

	MTR_END_FUNC();

	return idx;
}

mat_csr_t* create_csr_from_coo(const mat_coo_t* input) {
	MTR_BEGIN_FUNC();

	mat_csr_t* csr = NULL;

	double* data = (double*)malloc(sizeof(double)*input->nnz);
	memcpy(data, input->data, sizeof(double)*input->nnz);
	int* colind = (int*)malloc(sizeof(int)*input->nnz);
	memcpy(colind, input->colind, sizeof(int)*input->nnz);
	int* rowind_tmp = (int*)malloc(sizeof(int)*input->nnz);
	memcpy(rowind_tmp, input->rowind, sizeof(int)*input->nnz);

	if (input->sorted_axis != COO_ROW_SORTED) sort_for_spmat(rowind_tmp, 0, input->nnz-1, colind, data);
	compress_array(rowind_tmp, colind, data, input->nnz);

	csr =  mat_csr(input->cols, input->rows, data, colind, rowind_tmp);

	MTR_END_FUNC();
	return csr;
}

void create_csr_from_coo_with_dest(mat_csr_t* dest, const mat_coo_t* input) {
	MTR_BEGIN_FUNC();

	memcpy(dest->data, input->data, sizeof(double)*input->nnz);
	memcpy(dest->colind, input->colind, sizeof(int)*input->nnz);
	int* rowind_tmp = (int*)malloc(sizeof(int)*input->nnz);
	memcpy(rowind_tmp, input->rowind, sizeof(int)*input->nnz);

	if (input->sorted_axis != COO_ROW_SORTED) sort_for_spmat(rowind_tmp, 0, input->nnz-1, dest->colind, dest->data);
	compress_array(rowind_tmp, dest->colind, dest->data, input->nnz);

	dest->cols = input->cols;
	dest->rows = input->rows;

	memcpy(dest->rowind, rowind_tmp, sizeof(int) * (input->cols + 1));
	free(rowind_tmp);

	MTR_END_FUNC();
}


mat_csc_t* create_csc_from_coo(const mat_coo_t* input) {
	MTR_BEGIN_FUNC();

	mat_csc_t* csc = NULL;

	double* data = (double*)malloc(sizeof(double)*input->nnz);
	memcpy(data, input->data, sizeof(double)*input->nnz);
	int* colind_tmp = (int*)malloc(sizeof(int)*input->nnz);
	memcpy(colind_tmp, input->colind, sizeof(int)*input->nnz);
	int* rowind = (int*)malloc(sizeof(int)*input->nnz);
	memcpy(rowind, input->rowind, sizeof(int)*input->nnz);

	if (input->sorted_axis != COO_COL_SORTED) sort_for_spmat(colind_tmp, 0, input->nnz-1, rowind, data);
	compress_array(colind_tmp, rowind, data, input->nnz);

	csc =  mat_csc(input->cols, input->rows, data, colind_tmp, rowind);

	MTR_END_FUNC();
	return csc;

}


void trans_coo(mat_coo_t* ans, const mat_coo_t* input) {
	MTR_BEGIN_FUNC();
	ans->cols = input->cols;
	ans->rows = input->rows;
	ans->nnz = input->nnz;
	memcpy(ans->data, input->data, sizeof(double)*input->nnz);
	memcpy(ans->colind, input->rowind, sizeof(int)*input->nnz);
	memcpy(ans->rowind, input->colind, sizeof(int)*input->nnz);
	if (input->sorted_axis == COO_COL_SORTED) ans->sorted_axis = COO_ROW_SORTED;
	if (input->sorted_axis == COO_ROW_SORTED) ans->sorted_axis = COO_COL_SORTED;
	MTR_END_FUNC();
}

mat_coo_t* create_coo_from_csr(const mat_csr_t* input) {
	MTR_BEGIN_FUNC();
	int n = nnz_csr(input);

	mat_coo_t* coo = new_mat_coo(input->cols, input->rows, n);

	for (int j = 0; j < input->cols; j++) {
		for (int i = input->rowind[j]; i < input->rowind[j + 1]; i++) {
			coo->data[i] = input->data[i];
			coo->colind[i] = input->colind[i];
			coo->rowind[i] = j;
		}
	}
	
	coo->sorted_axis = COO_ROW_SORTED;

	MTR_END_FUNC();
	return coo;
}



void apply_diag_Ad(mat_csr_t* ans, const mat_csr_t* A, const vector_t* d) {
	MTR_BEGIN_FUNC();

	ans->cols = A->cols;
	ans->rows = A->rows;

	memcpy(ans->colind, A->colind, sizeof(int) * nnz_csr(A));
	memcpy(ans->rowind, A->rowind, sizeof(int) * (A->cols + 1));

	for (int j = 0; j < A->cols; j++) {
		for (int i = A->rowind[j]; i < A->rowind[j + 1]; i++) {
			ans->data[i] = A->data[i] * d->data[A->colind[i]];
		}
	}


	MTR_END_FUNC();
}

void apply_diag_dA(mat_csr_t* ans, const mat_csr_t* A, const vector_t* d) {
	MTR_BEGIN_FUNC();

	ans->cols = A->cols;
	ans->rows = A->rows;

	memcpy(ans->colind, A->colind, sizeof(int) * nnz_csr(A));
	memcpy(ans->rowind, A->rowind, sizeof(int) * (A->cols + 1));

	for (int j = 0; j < A->cols; j++) {
		for (int i = A->rowind[j]; i < A->rowind[j + 1]; i++) {
			ans->data[i] = A->data[i] * d->data[j];
		}
	}

	MTR_END_FUNC();
}

void trans_csr(mat_csr_t* ans, const mat_csr_t* input) {
	MTR_BEGIN_FUNC();
	mat_coo_t* A_coo = create_coo_from_csr(input);
	mat_coo_t* At_coo = new_mat_coo(A_coo->cols, A_coo->rows, A_coo->nnz);
	trans_coo(At_coo, A_coo);
	create_csr_from_coo_with_dest(ans, At_coo);

	destroy_mat_coo(A_coo);
	destroy_mat_coo(At_coo);
	MTR_END_FUNC();
}

vector_t* create_vec_from_mmfile(const char* filename) {
	MTR_BEGIN_FUNC();
	const int BufferSize = 64;

	FILE* fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("target mmfile not found: %s\n", filename);
		exit(-1);
	}


	char line[BufferSize];

	int state = 0; //0: size, else: content
	vector_t* vec;
	int itr = 0;

	while (fgets(line, BufferSize, fp) != NULL) {

		if (line[0] == '%') continue;

		if (state == 0) {
			int cols;
			int rows;
			sscanf(line, "%d %d", &cols, &rows);

			int size = 0;

			if (cols == 1) {
				size = rows;
			} else if (rows == 1) {
				size = cols;
			} else {
				assert(false);
			}

			vec = new_vector(size);

			state = 1;
		} else {
			double val;
			sscanf(line, "%lf", &val);

			vec->data[itr++] = val;

		}
	}


	fclose(fp);
	MTR_END_FUNC();

	return vec;
}

mat_coo_t* create_coo_from_mmfile(const char* filename) {
	MTR_BEGIN_FUNC();
	const int BufferSize = 1024;

	FILE* fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("target mmfile not found: %s\n", filename);
		exit(-1);
	}

	char line[BufferSize];

	int state = 0; //0: size, else: content
	mat_coo_t* coo;
	int itr = 0;

	enum {
		GENERAL,
		SYMMETRIC,
		SKEW,
		HERMITIAN
	} attr;

	if (fgets(line, BufferSize, fp) != NULL) {
		char object[16];
		char format[16];
		char field[16];
		char attribute[16];
		sscanf(line, "%%%%MatrixMarket %15s %15s %15s %15s", object, format, field, attribute);

		if (strcmp(object, "matrix") != 0) {
			printf("mmfile: not a matrix file(%s)\n", object);
			assert(false);
		}

		if (strcmp(format, "coordinate") != 0) {
			printf("mmfile: unsupported format(%s)\n", format);
			assert(false);
		}

		if (strcmp(field, "real") != 0) {
			printf("mmfile: unsupported field(%s)\n", field);
			assert(false);
		}

		if (strcmp(attribute, "general") == 0) {
			attr = GENERAL;
		} else if(strcmp(attribute, "symmetric") == 0) {
			attr = SYMMETRIC;
		} else {
			printf("mmfile: unsupported attribute (%s)\n", attribute);
			assert(false);
		}


	} else {
		printf("empty file\n");
		assert(false);
	}

	while (fgets(line, BufferSize, fp) != NULL) {

		if (line[0] == '%') continue;

		if (state == 0) {
			int cols;
			int rows;
			int nnz;
			sscanf(line, "%d %d %d", &cols, &rows, &nnz);

			if (attr == SYMMETRIC) {
				nnz -= cols;
				nnz *= 2;
				nnz += cols;
			}

			coo = new_mat_coo(cols, rows, nnz);

			state = 1;
		} else {
			int col;
			int row;
			double val;
			sscanf(line, "%d %d %lf", &col, &row, &val);

			coo->colind[itr] = col-1;
			coo->rowind[itr] = row-1;
			coo->data[itr] = val;
			itr++;

			if (attr == SYMMETRIC) {
				if (col != row) {
					coo->colind[itr] = row-1;
					coo->rowind[itr] = col-1;
					coo->data[itr] = val;
					itr++;
				}
			}

		}
	}


	fclose(fp);
	MTR_END_FUNC();

	return coo;
}

// TODO: 未実装
mat_csc_t* create_csc_from_csr(const mat_csr_t* input) {
	MTR_BEGIN_FUNC();

	int cols = input->cols;
	int rows = input->rows;
	int nnz = nnz_csr(input);
	mat_csc_t* csc = new_mat_csc(cols, rows, nnz);

	MTR_END_FUNC();
	return csc;
}

// TODO: ゲロ遅い
mat_csr_t* create_csr_from_csc(const mat_csc_t* input) {
	MTR_BEGIN_FUNC();

	int cols = input->cols;
	int rows = input->rows;
	int nnz = nnz_csc(input);
	mat_csr_t* csr = new_mat_csr(cols, rows, nnz);

	MTR_END_FUNC();
	return csr;
}

