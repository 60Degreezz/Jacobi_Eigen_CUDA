# Jacobi_Eigen Serial and Parallel Comparison

## [Jacobi Eigen Parallel](/Jacobi_Eigen_Parallel.cu)

When the program starts, It takes in matrices or vectors as inputs from the users.  
Then the program will assign each individual vector to each of the threads.  
We then wrap the trivial kernel (je) around them to launch each thread operating on separate data sets (i.e. matrices), and as an output each thread produces its own set of eigenvalues (and eigenvectors - for this particular code base).

The code can handle Double-precision floating-point values in a matrix. For viweing the max number of digits edit the

```{ .cpp }
std::cout.precision(<<int value>>);
```

To execute the code:

### On Terminal

```{ .sh }
$ nvcc -o Jacobi_Eigen_Parallel Jacobi_Eigen_Parallel.cu
$ cuda-memcheck ./Jacobi_Eigen_Parallel
```

### On Notebook

```{ .sh }
!nvcc -o Jacobi_Eigen_Parallel Jacobi_Eigen_Parallel.cu
!cuda-memcheck ./Jacobi_Eigen_Parallel
```

user has to create a matrix for which (s)he wants to compute eigen values and eigen vectors manually. The matrix has to be written into the file matrix.dat which has to be put in the data folder. The matrix in the file matrix.dat has to be written column-wise in one long column whithout empty lines.  
For example, if your matrix is

```
        | 1 2 3 |
    A = | 2 1 2 |
        | 3 2 1 |
```

then this matrix has to be written into the file matrix.dat as follows

1  
2  
3  
2  
1  
2  
3  
2  
1

where first three numbers are the first column of A, second three numbers are the second column of A, etc.

## [Jacobi Eigen Serial](/JacobiEigen_Serial/Jacobi_Eigen_Serial.cpp)

1. the executable file of the code is supposed to be used in command line environment;

2. the command line has the form: ./Jacobi_Eigen_Serial N eivec tr_param
   where

N is the size of a symmetric square matrix;
eivec - a parameter equal to 0 or 1. 0 implies only eigen values computations, 1 includes both eigen vectors and eigen values;
tr-param - may have two values: 1) real; 2) test;
"test" uses one of the test matrices provided with the code; "real" implies a user's matrix would be used.

In the latter case a user has to create a matrix for which (s)he wants to compute eigen values and eigen vectors manually. The matrix has to be written into the file matrix.dat which has to be put in the data folder. The matrix in the file matrix.dat has to be written column-wise in one long column whithout empty lines.  
For example, if your matrix is

```
        | 1 2 3 |
    A = | 2 1 2 |
        | 3 2 1 |
```

then this matrix has to be written into the file matrix.dat as follows

1  
2  
3  
2  
1  
2  
3  
2  
1

where first three numbers are the first column of A, second three numbers are the second column of A, etc.

3. the output of the executable is written in one or two files depending on the parameter eivec. The names of the files are  
   eigenValues.dat  
   eigenVectors.dat

4. Windows users have to make minor changes in the code before trying to compile it. It has to do with paths representation (using back slashes, c:\ etc.)

To execute the code:

```{ .sh }
$ g++ JacobiEigen_Serial/Jacobi_Eigen_Serial.cpp -o JacobiEigen_Serial/Jacobi_Eigen_Serial
$ ./JacobiEigen_Serial/Jacobi_Eigen_Serial 20 0 real
```

> :warning: **Can execute for an infinite time.**

## [Jacobi Eigen Notebook](/Eigen_Value_parallel_cuda.ipynb)

Code is given to run both Jacobe Eigen Serial and parallel.  
> :space_invader: Place the Python Notebook in the same location as shown in the repository.  

The code was run in [Google Colab](https://colab.research.google.com/).  
The [matrix.dat](/data/matrix.dat) currently has a symmetric matrix of size 20 \* 20.  
This [link](https://onlinemathtools.com/generate-random-matrix) can be used to produce random matrices: <https://onlinemathtools.com/generate-random-matrix>
