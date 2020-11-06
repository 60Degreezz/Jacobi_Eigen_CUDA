#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <time.h>
#include "utils.h"
#include "testMatrices.h"

using namespace std;

int main(int argc, char *argv[])
{
	int N, nRot;
	bool eivec;		   // a flag for whether we compute the eigen vectors or not
	double **u = NULL; // the matrix
	double **v = NULL;
	double **tmp = NULL;
	double sm, tresh, g, h, t, theta, tau, c, s;
	clock_t start, end;
	double cpu_time_used = 0;
	std::ofstream file;
	std::string mystring;
	std::string matrix_file;

	if (argc != 4)
	{
		cout << "wrong format" << endl;
		cout << "the input must be: " << endl;
		cout << "1) the size N of a square matrix; " << endl;
		cout << "2) the flag for eigen vectors computation: 1 - both eigen values and eigen vectors are computed; 0 - only eigen values are computed;" << endl;
		cout << "3) a parameter indicating whether it is a test or real computation. The former uses one of the test matrices, the latter a user matrix from a file." << endl;
		cout << "Use 'test' for test computations and 'real' for real ones." << endl;
		exit(1);
	}
	else
	{
		N = atoi(argv[1]);
		eivec = atoi(argv[2]);
		mystring = argv[3];
		matrix_file = "./data/matrix.dat";
	}

	MemAlloc(u, N, N);
	MemAlloc(v, N, N);
	MemAlloc(tmp, N, N);
	vector<double> d(N);
	vector<double> b(N);
	vector<double> z(N);
	vector<int> index(N);

	// fill the matrix
	if (mystring == "real")
	{
		ifstream in;
		in.open(matrix_file);
		if (!in)
		{
			cerr << "Unable to open file containing matrix. It should have the name matrix.dat and placed in the data folder. " << endl;
			exit(1);
		}
		else
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
					in >> u[i][j];
			}
			in.close();
		}
	}
	if (mystring == "test")
	{
		// test matrices. Comment this out when your tests are finished
		testMatrix_1(u, N); // use N=30 here
							//testMatrix_2(u, N); // use N=44 here
	}
	getchar();

	if (eivec)
	{
		for (int p = 0; p < N; p++)
		{
			for (int q = 0; q < N; q++)
			{
				if (p == q)
					v[p][q] = 1.0;
				else
					v[p][q] = 0.0;
			}
		}
	}
	for (int p = 0; p < N; p++)
	{
		d[p] = u[p][p];
		b[p] = d[p];
		z[p] = 0.0;
	}
	nRot = 0;
	// iterations
	start = clock();
	for (int i = 0; i < 1000; i++)
	{
		sm = 0.0;
		for (int p = 0; p < N - 1; p++)
			for (int q = p + 1; q < N; q++)
				sm += abs(u[p][q]);

		if (sm == 0)
			break;
		if (i < 4)
			tresh = 0.2 * sm / (N * N);
		else
			tresh = 0.0;

		for (int p = 0; p < N - 1; p++)
		{
			for (int q = p + 1; q < N; q++)
			{
				g = 100.0 * abs(u[p][q]);
				if (i > 4 && (abs(d[p]) + g == abs(d[p])) && (abs(d[q]) + g == abs(d[q])))
				{
					u[p][q] = 0.0;
				}
				else if (abs(u[p][q]) > tresh) // start rotations
				{
					h = d[q] - d[p];
					if ((abs(h) + g) == abs(h))
						t = u[p][q] / h;
					else
					{
						theta = 0.5 * h / u[p][q];
						t = 1.0 / (abs(theta) + sqrt(1. + theta * theta));
						if (theta < 0)
							t = -t;
					}

					c = 1. / sqrt(1. + t * t);
					s = t * c;
					tau = s / (1. + c);
					h = t * u[p][q];
					z[p] -= h;
					z[q] += h;
					d[p] -= h;
					d[q] += h;
					u[p][q] = 0;
					for (int j = 0; j < p; j++)
					{
						g = u[j][p];
						h = u[j][q];
						u[j][p] = g - s * (h + g * tau);
						u[j][q] = h + s * (g - h * tau);
					}
					for (int j = p + 1; j < q; j++)
					{
						g = u[p][j];
						h = u[j][q];
						u[p][j] = g - s * (h + g * tau);
						u[j][q] = h + s * (g - h * tau);
					}
					for (int j = q + 1; j < N; j++)
					{
						g = u[p][j];
						h = u[q][j];
						u[p][j] = g - s * (h + g * tau);
						u[q][j] = h + s * (g - h * tau);
					}

					if (eivec)
					{
						for (int j = 0; j < N; j++)
						{
							g = v[j][p];
							h = v[j][q];
							v[j][p] = g - s * (h + g * tau);
							v[j][q] = h + s * (g - h * tau);
						}
					}
					nRot++;
				}
				// end of rotations
			}
		}
		for (int p = 0; p < N; p++)
		{
			b[p] = b[p] + z[p];
			d[p] = b[p];
			z[p] = 0;
		}
	}
	// end of iterations

	end = clock();
	// sorting of eigen values and eigen vectors
	// indexes for eigen values
	sort_indexes(d, index);

	// use indexes to re-order the eigen vectors
	int k = 0;
	for (auto j : index)
	{
		for (int i = 0; i < N; i++)
			tmp[i][k] = v[i][j];
		k++;
	}

	// sorted eigen vectors
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			v[i][j] = tmp[i][j];
	// write results to files
	file.open("./data/eigenValues.dat");
	file << std::fixed << setprecision(10);
	for (auto j : index)
		file << d[j] << endl;
	file.close();

	if (eivec)
	{
		file.open("./data/eigenVectors.dat");
		file << std::fixed << setprecision(10);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
				file << v[i][j] << " ";
			file << endl;
		}
		file.close();
	}
	// -->
	cpu_time_used = ((end - start) / (double)CLOCKS_PER_SEC);

	std::cout << "CPU Elapsed time: " << cpu_time_used << "s\n";
	MemRelease(u, N);
	MemRelease(v, N);
	MemRelease(tmp, N);

	return 0;
}
