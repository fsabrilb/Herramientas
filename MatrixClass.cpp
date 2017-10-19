#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <chrono>
/*------------------------------------------------CLASS MATRIX----------------------------------------------------------------*/
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
/*------------------------------------------------INSIDE CLASS MATRIX---------------------------------------------------------*/
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
/*------------------------------------------------MATRIX'S FUNCTIONS-----------------------------------------------------------*/
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
/*------------------------------------------------GLOBAL FUNCTIONS------------------------------------------------------------*/
double TimeTranspose(Matrix &A,int seed){
  std::chrono::time_point<std::chrono::system_clock> start, end;
  A.FillingRandom(seed);
  start = std::chrono::system_clock::now();  
  A.Transpose();
  end   = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  return elapsed_seconds.count();
}
/*------------------------------------------------MAIN PROGRAM----------------------------------------------------------------*/
int main(void){
  int i,N=20; /*Repetitions*/
  int j,M=10;  /*Matrix Size*/
  
  double sum=0,sum2=0,time=0;
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  for(j=0;j<M;j++){
    Matrix A(std::pow(2,j),std::pow(2,j));
    for(i=0;i<N;i++){ 
      time=TimeTranspose(A,i);
      sum+=time;
      sum2+=time*time;
    }
    A.Destroyer();
    double mean = sum/N;
    double sigma = std::sqrt(N*std::abs(sum2/N-mean*mean)/(N-1));
    std::cout<<std::pow(2,j)<<"\t"<<mean<<"\t"<<sigma<<std::endl;
  }
  return 0;
}
