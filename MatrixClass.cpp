#include <iostream>
#include <stdlib.h>

class Matrix{
private:
  double **M; int N_rows,N_cols;
public:
  /*Builder*/
  Matrix(int m,int n);
  /*Operator Overload*/
  Matrix operator+(const Matrix &A);
  Matrix operator-(const Matrix &A);
  Matrix operator*(const Matrix &A);
  /*Class Functions*/
  void Destroyer(void);
  void FillingRandom(int seed);
  void Transpose(void);
  void Print(void);
};
/*Build Matrix with parameters m (Number of rows), n (Number of columns)*/
Matrix::Matrix(int m,int n){
  N_rows=m; N_cols=n;
  M= new double*[N_rows];
  for(int i=0;i<N_rows;i++){
    M[i]=new double[N_cols];
  }
}
/*Define Matrix Sum*/
Matrix Matrix::operator+(const Matrix &A){
  Matrix R(N_rows,N_cols);
  for(int i=0;i<N_rows;i++){
    for(int j=0;j<N_cols;j++){
      R.M[i][j]=M[i][j]+A.M[i][j];
    }
  }
  return R;
}
/*Define Matrix Subtraction*/
Matrix Matrix::operator-(const Matrix &A){
  Matrix R(N_rows,N_cols);
  for(int i=0;i<N_rows;i++){
    for(int j=0;j<N_cols;j++){
      R.M[i][j]=M[i][j]-A.M[i][j];
    }
  }
  return R;
}
/*Define Matrix Multiplication*/
Matrix Matrix::operator*(const Matrix &A){
  Matrix R(N_rows,A.N_cols);double Sum;
  for(int i=0;i<N_rows;i++){
    for(int j=0;j<A.N_cols;j++){
      for(int k=0;k<N_cols;k++){
	Sum += M[i][k]*A.M[k][j];
      }
     R.M[i][j]=Sum;Sum=0;
    }
  }
  return R;
}
/*-------------------------------------------Matrix Functions-----------------------------------------------------------------*/
/*Destroy Matrix*/
void Matrix::Destroyer(void){
  for(int i=0;i<N_rows;i++){
    delete[] M[i];
  }
  delete[] M;
}
/*Fill Matrix with random numbers between 0 and 100*/
void Matrix::FillingRandom(int seed){
  srand(seed+1);
  for(int i=0;i<N_rows;i++){
    for(int j=0;j<N_cols;j++){
      M[i][j]=drand48()*100;
    }
  }
}
/*Calculate the transpose of the matrix*/
void Matrix::Transpose(void){
  Matrix R(N_cols,N_rows); 
  for(int i=0;i<N_cols;i++){
    for(int j=0;j<N_rows;j++){
      R.M[i][j]=M[j][i];
    }
  }
}
/*Print matrix*/
void Matrix::Print(void){
  for(int i=0;i<N_rows;i++){
    for(int j=0;j<N_cols;j++){
      std::cout<<M[i][j]<<"\t\t";
    }
    std::cout<<"\n";
  }
}

int main(void){
  int M=2,N=2;
  Matrix A(M,N),B(M,N),C(M,N),D(N,M),E(M,M);
  A.FillingRandom(1);
  B.FillingRandom(2);
  C=A+B; E=A*B;
  A.Print(); B.Print(); std::cout<<std::endl; E.Print();
  A.Destroyer();
  B.Destroyer();
  C.Destroyer();
  D.Destroyer();
  E.Destroyer();
  return 0;
}
