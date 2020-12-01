#ifndef LINEAR_H 
#define LINEAR_H

static const int _QUICK_SORT_LIMIT_ = 16;

// :: MODEL ::  mat_dense_t
typedef struct {
	int cols;
	int rows;
	double* data;
} mat_dense_t;

// :: MODEL :: mat_coo_t
typedef enum {
	COO_RANDOM,
	COO_COL_SORTED,
	COO_ROW_SORTED
}coo_sorted_axis_t;

typedef struct {
	int cols;
	int rows;
	int nnz;
	coo_sorted_axis_t sorted_axis;
	double* data;
	int* colind;
	int* rowind;
} mat_coo_t;

// :: MODEL :: mat_csr
typedef struct {
	int cols;
	int rows;
	double* data;
	int* colind;
	int* rowind; //(compressed)
} mat_csr_t;

// :: MODEL :: mat_csc
typedef struct {
	int cols;
	int rows;
	double* data;
	int* colind;
	int* rowind; //(compressed)
} mat_csc_t;

// :: MODEL :: vector
typedef struct {
	int size;
	double* data;
} vector_t;


// :: function :: generator
//create dense
mat_dense_t* mat_dense(int cols, int rows, double* data);
mat_dense_t* new_mat_dense(int cols, int rows);
void destroy_mat_dense(mat_dense_t* mat);
mat_dense_t* create_dense_from_coo(const mat_coo_t* input);
mat_dense_t* create_dense_from_csr(const mat_csr_t* input);
mat_dense_t* create_dense_from_csc(const mat_csc_t* input);

// create coo
mat_coo_t* mat_coo(int cols, int rows, int nnz, double* data, int* colind, int* rowind);
mat_coo_t* new_mat_coo(int cols, int rows, int nnz);
void detect_coo_axis(mat_coo_t* mat);
void destroy_mat_coo(mat_coo_t* mat);
mat_coo_t* create_coo_from_csr(const mat_csr_t* input);
mat_coo_t* create_coo_from_mmfile(const char* filename);

// create csr
mat_csr_t* mat_csr(int cols, int rows, double* data, int* colind, int* rowind);
mat_csr_t* new_mat_csr(int cols, int rows, int nnz);
void destroy_mat_csr(mat_csr_t* mat);
mat_csr_t* create_csr_from_dense(const mat_dense_t* input);
mat_csr_t* create_csr_from_coo(const mat_coo_t* input);
mat_csr_t* create_csr_from_csc(const mat_csc_t* input);
void create_csr_from_coo_with_dest(mat_csr_t* dest, const mat_coo_t* input);

// create csc
mat_csc_t* mat_csc(int cols, int rows, double* data, int* colind, int* rowind);
mat_csc_t* new_mat_csc(int cols, int rows, int nnz);
void destroy_mat_csc(mat_csc_t* mat);
mat_csc_t* create_csc_from_coo(const mat_coo_t* input);
mat_csc_t* create_csc_from_csr(const mat_csr_t* input);

// create vec
vector_t* vector(int size, double* data);
vector_t* new_vector(int size);
void destroy_vector(vector_t* vec);
vector_t* create_vec_from_mmfile(const char* filename);

// :: function :: pretty print
void print_mat_csr(const mat_csr_t* input);
void print_vector(const vector_t* vec);

// :: function :: operate
// vector
bool eq_vector(const vector_t* x, const vector_t* y);
vector_t* clone_vec(const vector_t* x);
void copy_vec(const vector_t* x, vector_t* y);
double norm_vec(const vector_t* x);
void abs_vec(vector_t* ans, const vector_t* x);
double sum_vec(const vector_t* x);
void mul_double_vec(vector_t* ans, const double a, const vector_t* x);
void div_vec_double(vector_t* ans, const vector_t* x, const double a);
void pow_vec_double(vector_t* ans, const vector_t* x, const double a);
void add_vec_vec(vector_t* ans, const vector_t* x, const vector_t* y);
void sub_vec_vec(vector_t* ans, const vector_t* x, const vector_t* y);
double dot_vec_vec(const vector_t* x, const vector_t* y);

// matrix
void mul_mat_vec(vector_t* ans, const mat_csr_t* A, const vector_t* x);

bool eq_mat_csr(const mat_csr_t* A, const mat_csr_t* B);
int nnz_csr(const mat_csr_t* A);
int nnz_csc(const mat_csc_t* A);

void trans_coo(mat_coo_t* ans, const mat_coo_t* input);
void trans_csr(mat_csr_t* ans, const mat_csr_t* input); // TODO
// TODO trans_csc
void apply_diag_Ad(mat_csr_t* ans, const mat_csr_t* A, const vector_t* d);
void apply_diag_dA(mat_csr_t* ans, const mat_csr_t* A, const vector_t* d);


#endif

