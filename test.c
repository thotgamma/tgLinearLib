
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include "test.h"
#include "utility.h"
#include "linear.h"
#include "solver.h"
#include "minitrace/minitrace.h"

void test() {
	MTR_BEGIN_FUNC();

	// ========================== :: UTLITY_H :: ==========================
	// bool check_sorted_acend(int* target, int n);
	int list_acend[] = {1, 2, 3, 3, 4, 5, 6};
	int list_random[] = {1, 3, 3, 3, 4, 0, 2};

	assert(check_sorted_acend(list_acend, 7) == true);
	assert(check_sorted_acend(list_random, 7) == false);

	// int intersect_i(const int* A, const int nA, const int* B, const int nB, int* X);

	// ========================== :: LINEAR_H:: ==========================

	// :: function :: generator
	// ## create dense ##
	// mat_dense_t* mat_dense(int cols, int rows, double* data);
	// void destroy_mat_dense(mat_dense_t* mat);
	// mat_dense_t* create_dense_from_coo(const mat_coo_t* input);
	// mat_dense_t* create_dense_from_csr(const mat_csr_t* input);
	// mat_dense_t* create_dense_from_csc(const mat_csc_t* input);

	// ## create coo ##
	// mat_coo_t* mat_coo(int cols, int rows, int nnz, double* data, int* colind, int* rowind);
	// mat_coo_t* new_mat_coo(int cols, int rows, int nnz);
	// void detect_coo_axis(mat_coo_t* mat);
	// void destroy_mat_coo(mat_coo_t* mat);
	// mat_coo_t* create_coo_from_csr(const mat_csr_t* input);
	// mat_coo_t* create_coo_from_mmfile(const char* filename);

	// ## create csr ##
	// mat_csr_t* mat_csr(int cols, int rows, double* data, int* colind, int* rowind);
	// mat_csr_t* new_mat_csr(int cols, int rows, int nnz);
	// void destroy_mat_csr(mat_csr_t* mat);
	// mat_csr_t* create_csr_from_dense(const mat_dense_t* input);
	// mat_csr_t* create_csr_from_coo(const mat_coo_t* input);
	// mat_csr_t* create_csr_from_csc(const mat_csc_t* input);
	// void create_csr_from_coo_with_dest(mat_csr_t* dest, const mat_coo_t* input);

	// ## create csc ##
	// mat_csc_t* mat_csc(int cols, int rows, double* data, int* colind, int* rowind);
	// mat_csc_t* new_mat_csc(int cols, int rows, int nnz);
	// void destroy_mat_csc(mat_csc_t* mat);
	// mat_csc_t* create_csc_from_coo(const mat_coo_t* input);
	// mat_csc_t* create_csc_from_csr(const mat_csr_t* input);

	// ## create vec ##
	// vector_t* vector(int size, double* data);
	// vector_t* new_vector(int size);
	// void destroy_vector(vector_t* vec);
	// vector_t* create_vec_from_mmfile(const char* filename);

	// :: function :: operate
	// ## vector ##
	// vector_t* clone_vec(const vector_t* x);
	{
	double x_data[] = {1, 2, 3, 4};
	vector_t* x = vector(4, x_data);
	vector_t* y = clone_vec(x);
	assert(eq_vector(x, y));
	}

	// void copy_vec(const vector_t* x, vector_t* y);
	{
	double x_data[] = {1, 2, 3, 4};
	vector_t* x = vector(4, x_data);
	vector_t* y = new_vector(4);
	copy_vec(y, x);
	assert(eq_vector(x, y));
	destroy_vector(y);
	}

	// double norm_vec(const vector_t* x);
	{
	double x_data[] = {1, 2, 3, 4};
	vector_t* x = vector(4, x_data);
	double ans = norm_vec(x);
	assert(fabs(ans - 5.4722) < 10e-3);
	}

	// void abs_vec(vector_t* ans, const vector_t* x);
	{
	double x_data[] = {1, -2, 3, -4};
	vector_t* x = vector(4, x_data);
	double ans_data[] = {1, 2, 3, 4};
	vector_t* ans = vector(4, ans_data);
	vector_t* y = new_vector(4);
	abs_vec(y, x);
	assert(eq_vector(y, ans));
	destroy_vector(y);
	}

	// double sum_vec(const vector_t* x);
	{
	double x_data[] = {1, 2, 3, 4};
	vector_t* x = vector(4, x_data);
	double ans = sum_vec(x);
	assert(fabs(ans - 10.0) < 10e-3);
	}

	// void mul_double_vec(vector_t* ans, const double a, const vector_t* x);
	{
	double x_data[] = {1, 2, 3, 4};
	vector_t* x = vector(4, x_data);
	double ans_data[] = {4, 8, 12, 16};
	vector_t* ans = vector(4, ans_data);
	vector_t* y = new_vector(4);
	mul_double_vec(y, 4.0, x);
	assert(eq_vector(y, ans));
	destroy_vector(y);
	}

	// void div_vec_double(vector_t* ans, const vector_t* x, const double a);
	{
	double x_data[] = {4, 8, 12, 16};
	vector_t* x = vector(4, x_data);
	double ans_data[] = {1, 2, 3, 4};
	vector_t* ans = vector(4, ans_data);
	vector_t* y = new_vector(4);
	div_vec_double(y, x, 4);
	assert(eq_vector(y, ans));
	destroy_vector(y);
	}

	// void pow_vec_double(vector_t* ans, const vector_t* x, const double a);
	{
	double x_data[] = {1, 2, 3, 4};
	vector_t* x = vector(4, x_data);
	double ans_data[] = {1, 4, 9, 16};
	vector_t* ans = vector(4, ans_data);
	vector_t* y = new_vector(4);
	pow_vec_double(y, x, 2.0);
	assert(eq_vector(y, ans));
	destroy_vector(y);
	}

	// void add_vec_vec(vector_t* ans, const vector_t* x, const vector_t* y);
	{
	double x_data[] = {1, 2, 3, 4};
	vector_t* x = vector(4, x_data);
	double y_data[] = {4, 8, 12, 16};
	vector_t* y = vector(4, y_data);

	double ans_data[] = {5, 10, 15, 20};
	vector_t* ans = vector(4, ans_data);

	vector_t* z = new_vector(4);
	add_vec_vec(z, x, y);
	assert(eq_vector(z, ans));
	destroy_vector(z);
	}

	// void sub_vec_vec(vector_t* ans, const vector_t* x, const vector_t* y);
	{
	double x_data[] = {4, 8, 12, 16};
	vector_t* x = vector(4, x_data);
	double y_data[] = {4, 5, 20, 15};
	vector_t* y = vector(4, y_data);

	double ans_data[] = {0, 3, -8, 1};
	vector_t* ans = vector(4, ans_data);

	vector_t* z = new_vector(4);
	sub_vec_vec(z, x, y);
	assert(eq_vector(z, ans));
	destroy_vector(z);
	}

	// double dot_vec_vec(const vector_t* x, const vector_t* y);
	{
	double x_data[] = {1, 2, 3, 4};
	vector_t* x = vector(4, x_data);
	double y_data[] = {4, 8, 12, 16};
	vector_t* y = vector(4, y_data);
	double a = dot_vec_vec(x, y);
	assert(a == 4 + 16 + 36 + 64);
	}

	// ## matrix ##
	// void mul_mat_vec(vector_t* ans, const mat_csr_t* A, const vector_t* x);
	{
	mat_dense_t A_dense;
	A_dense.cols = 4;
	A_dense.rows = 4;
	double A_data[] = {
					 1, 0, 0, 0,
					 0, 1, 0, 0,
					 0, 0, 1, 0,
					 0, 0, 0, 1};
	A_dense.data = A_data;

	mat_csr_t* A = create_csr_from_dense(&A_dense);

	double b_data[] = {1, 2, 3, 4};
	vector_t* b = vector(4, b_data);

	vector_t* x = new_vector(4);
	mul_mat_vec(x, A, b);

	assert(eq_vector(x, b));

	destroy_vector(x);
	}

	// int nnz_csr(const mat_csr_t* A);
	// int nnz_csc(const mat_csc_t* A);

	// void trans_coo(mat_coo_t* ans, const mat_coo_t* input);
	// void trans_csr(mat_csr_t* ans, const mat_csr_t* input);
	{
	mat_dense_t A_dense;
	A_dense.cols = 4;
	A_dense.rows = 4;
	double A_data[] = {
					 1, 0, 0, 7,
					 2, 1, 0, 0,
					 3, 2, 1, 0,
					 4, 3, 2, 1};
	A_dense.data = A_data;
	mat_csr_t* A = create_csr_from_dense(&A_dense);

	mat_dense_t B_dense;
	B_dense.cols = 4;
	B_dense.rows = 4;
	double B_data[] = {
					 1, 2, 3, 4,
					 0, 1, 2, 3,
					 0, 0, 1, 2,
					 7, 0, 0, 1};
	B_dense.data = B_data;
	mat_csr_t* B = create_csr_from_dense(&B_dense);

	mat_csr_t* At = new_mat_csr(4, 4, 11);
	trans_csr(At, A);

	assert(eq_mat_csr(At, B));

	destroy_mat_csr(At);
	}

	// trans_csc
	// void apply_diag_Ad(mat_csr_t* ans, const mat_csr_t* A, const vector_t* d);
	// void apply_diag_dA(mat_csr_t* ans, const mat_csr_t* A, const vector_t* d);


	// -- linear::detect_axis
	double coo0_data[] = {1, 3, 5, 8, 2, 4, 6, 9, 7};
	int coo0_col[]     = {1, 0, 1, 0, 3, 2, 2, 2, 3};
	int coo0_row[]     = {0, 1, 2, 3, 0, 1, 2, 3, 2};
	mat_coo_t* coo0 = mat_coo(4, 4, 9, coo0_data, coo0_col, coo0_row);

	assert(coo0->sorted_axis == COO_RANDOM);

	mat_csr_t* csr0 = create_csr_from_coo(coo0);
	mat_coo_t* coo1 = create_coo_from_csr(csr0);

	assert(coo1->sorted_axis == COO_ROW_SORTED);
	detect_coo_axis(coo1);
	assert(coo1->sorted_axis == COO_ROW_SORTED);

	// ========================== :: SOLVER_H :: ==========================
	// void incomp_chol_decomp(const mat_csr_t* A, mat_csr_t* L, vector_t* d);
	// void solve_tri(vector_t* ans, const mat_csr_t* L, const mat_csr_t* U, const vector_t* b);
	// double cg_solver(vector_t* ans, const mat_csr_t* A, const vector_t* b, const double e, const int max_itr); //TODO
	{
	mat_dense_t A_dense;
	A_dense.cols = 4;
	A_dense.rows = 4;
	double A_data[] = {4, 3, 2, 1,
					   3, 4, 3, 2,
					   2, 3, 4, 3,
					   1, 2, 3, 4};
	A_dense.data = A_data;

	mat_csr_t* A = create_csr_from_dense(&A_dense);


	double b_data[] = {1, 1, 1, 1};
	vector_t* b = vector(4, b_data);

	vector_t* x = new_vector(4);
	cg_solver(x, A, b, 10e-7, 64);

	double ans_data[] = {0.2, 0.0, 0.0, 0.2};
	vector_t* ans = vector(4, ans_data);

	assert(eq_vector(x, ans));

	}
	// double iccg_solver(vector_t* ans, const mat_csr_t* A, const vector_t* b, const double e, const int max_itr);
	// double irls_solver(vector_t* ans, const mat_csr_t* A, const vector_t* b, const double p, const int max_itr);

	printf("test passed.\n");
	MTR_END_FUNC();
}

