#include <fstream>
#include <stdio.h>
#include <iostream>

const int num_mat = 1; // total number of matrices = total number of threads const
const int N = 20; // square symmetric matrix dimension
const int nTPB = 256; // threads per block

// test symmetric matrices

double a1[N*N];

/* ---------------------------------------------------------------*/
//****************************************************************************80
//(n, a, d)
__host__ __device__
void r8mat_diag_get_vector(int n, double a[], double v[])
//****************************************************************************80
{
  int i;
  for (i = 0; i < n; i++) {
    v[i] = a[i + i * n];
  }

  return;
}
__host__ __device__
void r8mat_identity(int n, double a[]) {
  int i;
  int j;
  int k;
  k = 0;
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      if (i == j) {
        a[k] = 1.0;
      } else {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }
  return;
}
//****************************************************************************80

//(n, a + (idx * n * n), it_max, v + (idx * n * n), d + (idx * n), it_num, rot_num)
__host__ __device__
void jacobi_eigenvalue(int n, double a[], int it_max, double v[], double d[], int &it_num, int &rot_num) {
  double * bw;
  double c;
  double g;
  double gapq;
  double h;
  int i;
  int j;
  int k;
  int l;
  int m;
  int p;
  int q;
  double s;
  double t;
  double tau;
  double term;
  double termp;
  double termq;
  double theta;
  double thresh;
  double w;
  double * zw;
  r8mat_identity(n, v);
  r8mat_diag_get_vector(n, a, d);
  bw = new double[n];
  zw = new double[n];
  for (i = 0; i < n; i++) {

    bw[i] = d[i];
    zw[i] = 0.0;
  }
  it_num = 0;
  rot_num = 0;
  while (it_num < it_max) {
    it_num = it_num + 1;
    //
    // The convergence threshold is based on the size of the elements in
    // the strict upper triangle of the matrix.
    //
    thresh = 0.0;
    for (j = 0; j < n; j++) {
      for (i = 0; i < j; i++) {
        thresh = thresh + a[i + j * n] * a[i + j * n];
      }
    }
    thresh = sqrt(thresh) / (double)(4 * n);
    if (thresh == 0.0) {
      break;
    }
    for (p = 0; p < n; p++) {
      for (q = p + 1; q < n; q++) {
        gapq = 10.0 * fabs(a[p + q * n]);
        termp = gapq + fabs(d[p]);
        termq = gapq + fabs(d[q]);
        //
        // Annihilate tiny offdiagonal elements.
        //
        if (4 < it_num &&
          termp == fabs(d[p]) && termq == fabs(d[q])) {
          a[p + q * n] = 0.0;

        }
        //
        // Otherwise, apply a rotation.
        //
        else if (thresh <= fabs(a[p + q * n])) {
          h = d[q] - d[p];
          term = fabs(h) + gapq;
          if (term == fabs(h)) {
            t = a[p + q * n] / h;
          } else {
            theta = 0.5 * h / a[p + q * n];
            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0) {
              t = -t;
            }
          }
          c = 1.0 / sqrt(1.0 + t * t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * a[p + q * n];
          //
          // Accumulate corrections to diagonal elements.
          //
          zw[p] = zw[p] - h;
          zw[q] = zw[q] + h;
          d[p] = d[p] - h;
          d[q] = d[q] + h;
          a[p + q * n] = 0.0;
          //
          // Rotate, using information from the upper triangle of A only.
          //
          for (j = 0; j < p; j++) {
            g = a[j + p * n];
            h = a[j + q * n];
            a[j + p * n] = g - s * (h + g * tau);
            a[j + q * n] = h + s * (g - h * tau);
          }

          for (j = p + 1; j < q; j++) {
            g = a[p + j * n];
            h = a[j + q * n];
            a[p + j * n] = g - s * (h + g * tau);
            a[j + q * n] = h + s * (g - h * tau);
          }
          for (j = q + 1; j < n; j++) {
            g = a[p + j * n];
            h = a[q + j * n];
            a[p + j * n] = g - s * (h + g * tau);
            a[q + j * n] = h + s * (g - h * tau);
          }
          //
          // Accumulate information in the eigenvector matrix.
          //
          for (j = 0; j < n; j++) {
            g = v[j + p * n];
            h = v[j + q * n];
            v[j + p * n] = g - s * (h + g * tau);
            v[j + q * n] = h + s * (g - h * tau);
          }
          rot_num = rot_num + 1;
        }
      }
    }
    for (i = 0; i < n; i++) {
      bw[i] = bw[i] + zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }
  }
  //
  // Restore upper triangle of input matrix.
  //
  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      a[i + j * n] = a[j + i * n];
    }
  }

  //
  // Ascending sort the eigenvalues and eigenvectors.
  //
  for (k = 0; k < n - 1; k++) {
    m = k;
    for (l = k + 1; l < n; l++) {
      if (d[l] < d[m]) {
        m = l;
      }
    }
    if (m != k) {
      t = d[m];
      d[m] = d[k];
      d[k] = t;
      for (i = 0; i < n; i++) {
        w = v[i + m * n];
        v[i + m * n] = v[i + k * n];
        v[i + k * n] = w;
      }
    }
  }
  delete[] bw;
  delete[] zw;
  return;
}

// end of FSU code
/* ---------------------------------------------------------------- */

//(num_mat, N,d_a, max_iter, d_v, d_d)
__global__ void je(int num_matr, int n, double *a, int it_max, double *v, double *d) {
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  int it_num;
  int rot_num;
  if (idx < num_matr) {
    jacobi_eigenvalue(n, a + (idx * n * n), it_max, v + (idx * n * n), d + (idx * n), it_num, rot_num);
  }
}

//(0, N, a1, h_a)
void initialize_matrix(int mat_id, int n, double * mat, double * v) {
  for (int i = 0; i < n * n; i++) * (v + (mat_id * n * n) + i) = mat[i];
}
void print_vec(int vec_id, int n, double * d) {
  std::cout << "matrix " << vec_id << " Diagonal Values: " << std::endl;
  std::cout.precision(12);
  for (int i = 0; i < n; i++) std::cout << i << ": " << * (d + (n * vec_id) + i) << std::endl;
  std::cout << std::endl;
}
int main() {
  // make sure device heap has enough space for in-kernel new allocations
  const int heapsize = num_mat * N * sizeof(double) * 2;
  const int chunks = heapsize / (8192 * 1024) + 1;
  cudaError_t cudaStatus = cudaDeviceSetLimit(cudaLimitMallocHeapSize, (8192 * 1024) * chunks);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "set device heap limit failed!");
  }


  //reading matrix from file matrix.dat
  std::ifstream in;
  in.open("data/matrix.dat");
	if (!in)
	{
		std::cerr << "Unable to open file containing matrix. It should have the name matrix.dat. It should be present inside the data folder" << std::endl;
		exit(1);
	}
	else
	{
		for (int i = 0; i < N*N; i++)
		{
			in >> a1[i];
		}
 	  in.close();
  }

  const int max_iter = 1000;
  double *h_a, *d_a, *h_v, *d_v, *h_d, *d_d;

  //initialize_matrix
  h_a = (double * ) malloc(num_mat * N * N * sizeof(double));
  h_v = (double *) malloc(num_mat * N * N * sizeof(double));
  h_d = (double * ) malloc(num_mat * N * sizeof(double));

  //Allocating mem in GPU
  cudaMalloc( & d_a, num_mat * N * N * sizeof(double));
  cudaMalloc( & d_v, num_mat * N * N * sizeof(double));
  cudaMalloc( & d_d, num_mat * N * sizeof(double));

  //sets all the bytes to a specific value (0)
  memset(h_a, 0, num_mat * N * N * sizeof(double));
  memset(h_v, 0, num_mat * N * N * sizeof(double));
  memset(h_d, 0, num_mat * N * sizeof(double));

  for(int i = 0; i< num_mat ; i++)
  initialize_matrix(i, N, a1, h_a);

  //initialize_matrix(1, N, a2, h_a);
  //initialize_matrix(2, N, a3, h_a);

  //copys variable to GPU
  cudaMemcpy(d_a, h_a, num_mat * N * N * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_v,h_v, num_mat * N * N * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_d, h_d, num_mat *N * sizeof(double), cudaMemcpyHostToDevice);


  clock_t cpu_startTime, cpu_endTime;
  double cpu_ElapseTime=0;
  cpu_startTime = clock();

  //Launches Kernel
  //<<<N/BLOCK_SIZE,BLOCK_SIZE>>>
  je << < (num_mat + nTPB - 1) / nTPB, nTPB >>> (num_mat, N,d_a, max_iter, d_v, d_d);

  cpu_endTime = clock();
  cpu_ElapseTime = ((cpu_endTime - cpu_startTime)/(double)CLOCKS_PER_SEC);
  std::cout << "CPU Elapsed time: " <<cpu_ElapseTime<<"s\n";


  //Copy back to CPU
  cudaMemcpy(h_d, d_d, num_mat * N * sizeof(double), cudaMemcpyDeviceToHost);


  //for(int i = 0; i< num_mat ; i++)
  //  print_vec(i, N, h_d);

  print_vec(0, N, h_d);
  //print_vec(2, N, h_d);

  return 0;
}
