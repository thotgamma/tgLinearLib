
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "utility.h"
#include "linear.h"
#include "solver.h"
#include "test.h"
#include "minitrace/minitrace.h"


int main() {


	// -- prepare mini trace -----------------
	mtr_init("trace.json");

	mtr_register_sigint_handler();

	MTR_META_PROCESS_NAME("minitrace_test");
	MTR_META_THREAD_NAME("main thread");

	MTR_BEGIN("main", "main");
	// ---------------------------------------

	// -- run test --
		//test();
	// --------------



	mat_coo_t* A_coo = create_coo_from_mmfile("thermomech_dM.mtx"); //cols: 204316, rows: 204316, nnz: 1423116

	printf("matrix loaded. cols: %d, rows: %d, nnz: %d\n", A_coo->cols, A_coo->rows, A_coo->nnz);
	
	mat_csr_t* A = create_csr_from_coo(A_coo);
	//print_mat_csr(A);

	int n = A->cols;

	vector_t* b = new_vector(n);
	for (int i = 0; i < n; i++) b->data[i] = 1;

	vector_t* ans = new_vector(n);
	mul_mat_vec(ans, A, b);

	vector_t* cg_ans = new_vector(n);
	double cg_error = cg_solver(cg_ans, A, b, 1e-10, 512);
	printf("cg_error: %f\n", cg_error);
	//print_vector(cg_ans);

	vector_t* ic_ans = new_vector(n);
	double ic_error = iccg_solver(ic_ans, A, b, 1e-10, 512);
	printf("ic_error: %f\n", ic_error);
	//print_vector(ic_ans);

	// ---------------------------------------

	MTR_END("main", "main");

	mtr_shutdown();


	return 0;
}
