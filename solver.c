

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <math.h>
#include "utility.h"
#include "solver.h"
#include "minitrace/minitrace.h"



int arr_index(int* x, int target) {
	int ind = 0;
	while (1) {
		if (x[ind] == target) break;
		ind++;
	}
	return ind;
}



void incomp_chol_decomp(const mat_csr_t* A, mat_csr_t* L, vector_t* d) {
	MTR_BEGIN_FUNC();

	L->data[0] = A->data[0];
	L->colind[0] = A->colind[0];
	d->data[0] = 1.0/L->data[0];

	L->rowind[0] = 0;
	L->rowind[1] = 1;

	int* x = (int*)malloc(sizeof(int) * A->cols);

	int L_data_itr = 1;

	for (int j = 1; j < A->cols; j++) {

		double ii = 0;
		int datas_in_a_row = 0;

		for (int i = A->rowind[j]; i < A->rowind[j + 1]; i++) {

			int col = j;
			int row = A->colind[i];

			if (row > col) continue; // 最適化の余地あり

			double lld = A->data[i];

			int ncol = datas_in_a_row;
			int nrow = L->rowind[row+1] - L->rowind[row];
			if (col == row) nrow = datas_in_a_row;

			int mincol = (ncol > row+1) ? row+1 : ncol;
			int minrow = (nrow > row+1) ? row+1 : nrow;

			int n = intersect_i(&(L->colind[L->rowind[col]]), mincol, &(L->colind[L->rowind[row]]), minrow, x);

			for (int k = 0; k < n; k++) {
				int ind = arr_index(&(L->colind[L->rowind[col]]), x[k]);
				lld -= L->data[L->rowind[col] + ind]
					 * L->data[L->rowind[row] + ind]
					 * d->data[L->colind[L->rowind[col] + ind]];
			}

			assert(!(lld == 0));
			assert(!isnan(lld));
			assert(!isinf(lld));

			L->data[L_data_itr] = lld;
			L->colind[L_data_itr] = row;

			L_data_itr++;
			datas_in_a_row++;

			if (A->colind[i] == j) ii = lld;

		}

		L->rowind[j+1] = L_data_itr;
		d->data[j] = 1.0/ii;

	}

	free(x);
	MTR_END_FUNC();
}


void solve_tri(vector_t* ans, const mat_csr_t* L, const mat_csr_t* U, const vector_t* b) {
	MTR_BEGIN_FUNC();

	assert(L->rows == b->size);
	assert(L->rows == b->size);

	int n = L->rows;

	// TODO: xは0番目の要素を初期化したansでも良いはず
	double* x = (double*)calloc(n, sizeof(double));
	double* y = (double*)calloc(n, sizeof(double));

	// 前進代入
	for (int j = 0; j < n; j++) {

		double buff = b->data[j];

		for (int i = L->rowind[j]; i < L->rowind[j + 1]; i++) {
			buff -= L->data[i] * y[ L->colind[i] ];
		}

		y[j] = buff / L->data[ L->rowind[j + 1] - 1 ];
	}


	// 後進代入
	for (int j = n-1; j >= 0; j--) {

		double buff = y[j];

		for (int i = U->rowind[j + 1]-1; i >= U->rowind[j]; i--) {
			buff -= U->data[i] * x[ U->colind[i] ];
		}

		x[j] = buff / U->data[ U->rowind[j] ];
	}

	memcpy(ans->data, x, sizeof(double) * n);
	free(x);
	free(y);
	MTR_END_FUNC();
}


double cg_solver(vector_t* ans, const mat_csr_t* A, const vector_t* b, const double e, const int max_itr) {
	MTR_BEGIN_FUNC();

	int n = A->rows;

	vector_t* y = new_vector(n);
	vector_t* r = clone_vec(b);
	vector_t* r_= clone_vec(r);
	vector_t* p = clone_vec(r);

	vector_t* alpha_p = new_vector(n);
	vector_t* alpha_y = new_vector(n);
	vector_t* beta_p = new_vector(n);

	double alpha;
	double beta;

	double error = 0;

	MTR_BEGIN("cg_solver", "mainLoop");
	for (int i = 0; i < max_itr; i++) {
		mul_mat_vec(y, A, p);
		alpha = dot_vec_vec(r, r)/dot_vec_vec(p, y);
		mul_double_vec(alpha_p, alpha, p);
		add_vec_vec(ans, ans, alpha_p);
		mul_double_vec(alpha_y, alpha, y);
		sub_vec_vec(r_, r, alpha_y);

		error = norm_vec(r_);
		if (error < e) {
			printf("cg solved in %d\n", i);
			break;
		}

		beta = dot_vec_vec(r_, r_)/dot_vec_vec(r, r);
		mul_double_vec(beta_p, beta, p);
		add_vec_vec(p, r_, beta_p);
		copy_vec(r, r_);
	}
	MTR_END("cg_solver", "mainLoop");

	MTR_BEGIN("cg_solver", "destroy");
	destroy_vector(y);
	destroy_vector(r);
	destroy_vector(r_);
	destroy_vector(p);
	destroy_vector(alpha_p);
	destroy_vector(alpha_y);
	destroy_vector(beta_p);
	MTR_END("cg_solver", "destroy");
	MTR_END_FUNC();

	return error;
}


double iccg_solver(vector_t* ans, const mat_csr_t* A, const vector_t* b, const double e, const int max_itr) {
	MTR_BEGIN_FUNC();

	assert(A->cols == b->size);

	int n = A->rows;

	mat_csr_t* L = new_mat_csr(A->cols, A->rows, nnz_csr(A));
	vector_t* d = new_vector(n);

	incomp_chol_decomp(A, L, d);

	vector_t*  y  = new_vector(n);
	vector_t*  r = clone_vec(b);
	vector_t*  r_ = new_vector(n);

	mat_csr_t* Ld = new_mat_csr(L->cols, L->rows, nnz_csr(L)); apply_diag_Ad(Ld, L, d);
	mat_csr_t* Lt = new_mat_csr(L->cols, L->rows, nnz_csr(L)); trans_csr(Lt, L);

	destroy_mat_csr(L);
	free(d);

	vector_t*  p  = new_vector(n); solve_tri(p, Ld, Lt, r);
	vector_t*  rp = clone_vec(p);
	vector_t*  rp_= new_vector(n);

	vector_t* alpha_p = new_vector(n);
	vector_t* alpha_y = new_vector(n);
	vector_t* beta_p = new_vector(n);

	double alpha = 0;
	double beta = 0;

	double error = 0;

	MTR_BEGIN("iccg_solver", "mainLoop");
	for (int i = 0; i < max_itr; i++) {
		mul_mat_vec(y, A, p);
		alpha = dot_vec_vec(r, rp)/dot_vec_vec(p, y);
		assert(!isnan(alpha));
		assert(!isinf(alpha));
		mul_double_vec(alpha_p, alpha, p);
		add_vec_vec(ans, ans, alpha_p);
		mul_double_vec(alpha_y, alpha, y);
		sub_vec_vec(r_, r, alpha_y);

		error = norm_vec(r_);
		assert(!isnan(error));
		assert(!isinf(error));
		if (error < e) {
			printf("iccg solved in %d\n", i);
			break;
		}

		solve_tri(rp_, Ld, Lt, r_);
		beta = dot_vec_vec(r_, rp_)/dot_vec_vec(r, rp);
		mul_double_vec(beta_p, beta, p);
		add_vec_vec(p, rp_, beta_p);

		copy_vec(r, r_);
		copy_vec(rp, rp_);

	}
	MTR_END("iccg_solver", "mainLoop");

	MTR_BEGIN("iccg_solver", "destroy");
	destroy_vector(y);
	destroy_vector(r);
	destroy_vector(r_);
	destroy_mat_csr(Ld);
	destroy_mat_csr(Lt);
	destroy_vector(p);
	destroy_vector(rp);
	destroy_vector(rp_);
	destroy_vector(alpha_p);
	destroy_vector(alpha_y);
	destroy_vector(beta_p);
	MTR_END("iccg_solver", "destroy");
	MTR_END_FUNC();

	return error;
}


/*
void irls_solver(vector_t* ans, const mat_csr_t* A, const vector_t* b, const double p, const int max_itr) {

	for (int i = 0; i < max_itr; i++) {

		vector_t* A_ans = new_vector(A->rows); mul_mat_vec(A_ans, A, ans);
		vector_t* e = new_vector(A->rows); sub_vec_vec(e, A_ans, b);

		vector_t* e_abs = new_vector(A->row); vec_abs(e_abs, e);
		vector_t* w = new_vector(A->row); pow_vec_double(w, e_abs, (p-2)/2);

		double w_sum = vec_sum(w);
		vector_t* w_n = new_vector(A->row); div_vec_double(w_n, w, w_sum);

		mat_csr_t* wA = new_mat_csr(A->cols, A->rows, nnz(A)); apply_diag_dA(wA, A, w_n);

		mat_csr_t* wAt;
		mat_csr_t* wAt_wA;
		mat_csr_t* wAt_w;
		mat_csr_t* wAt_w_b;

		e = A*ans - b;
		w = abs(e).^((p-2)/2);
		W = diag(w/sum(w));
		WA = W*A;
		ans = myLUSolver(WA'*WA, (WA'*W)*b)';

	}
}
*/
