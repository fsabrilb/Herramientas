#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;

class Matrix{
private:
  double **M; int N_rows,N_cols;
public:
  Matrix(int m,int n);
  Matrix operator+(const Matrix &A);
  Matrix operator-(const Matrix &A);
  //Matrix operator*(const Matrix &A);
  void Destroyer(void);
  void FillingRandom(int k);
  void Transpose(void);
  void Print(void);
};

Matrix::Matrix(int m,int n){
  N_rows=m; N_cols=n;
  M= new double*[N_rows];
  for(int i=0;i<N_rows;i++){M[i]=new double[N_cols];}
}

Matrix Matrix::operator+(const Matrix &A){
  Matrix R(N_rows,N_cols);
  for(int i=0;i<N_rows;i++){for(int j=0;j<N_cols;j++){R.M[i][j]=M[i][j]+A.M[i][j];}}
  return R;
}

Matrix Matrix::operator-(const Matrix &A){
  Matrix R(N_rows,N_cols);
  for(int i=0;i<N_rows;i++){for(int j=0;j<N_cols;j++){R.M[i][j]=M[i][j]-A.M[i][j];}}
  return R;
}

void Matrix::Destroyer(void){for(int i=0;i<N_rows;i++){delete[] M[i];} delete[] M;}
void Matrix::FillingRandom(int k){
  srand(k+1); for(int j=0;j<N_cols;j++){for(int i=0;i<N_rows;i++){M[i][j]=drand48()*10;}}
}
void Matrix::Transpose(void){
  Matrix R(N_cols,N_rows); 
  for(int i=0;i<N_cols;i++){for(int j=0;j<N_rows;j++){R.M[i][j]=M[j][i];}}
  R.Print();
}
void Matrix::Print(void){
  for(int i=0;i<N_rows;i++){for(int j=0;j<N_cols;j++){
      cout<<M[i][j]<<"\t";}
    cout<<endl;
  }
}

int main(void){
  int M=2,N=2;
  Matrix A(M,N),B(M,N),C(M,N),D(M,N);
  A.FillingRandom(1);
  B.FillingRandom(2);
  C=A+B; D=A-B;
  A.Print(); B.Print(); C.Print(); D.Print();
  A.Transpose();
  A.Destroyer();
  B.Destroyer();
  C.Destroyer();
  return 0;
}
