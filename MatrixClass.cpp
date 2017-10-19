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
double TimeTransposeDirect(Matrix &A,int seed){
  std::chrono::time_point<std::chrono::system_clock> start, end;
  A.FillingRandom(seed);
  start = std::chrono::system_clock::now();  
  A.Transpose();
  end   = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  return elapsed_seconds.count();
}
double TimeTransposeEigen(Eigen::MatrixXd &A,int seed, int size){
  std::chrono::time_point<std::chrono::system_clock> start, end;
  A=Eigen::MatrixXd::Random(size, size);
  start = std::chrono::system_clock::now();  
  Eigen::MatrixXd AT = A.transpose().eval();
  end   = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  return elapsed_seconds.count();
}
/*------------------------------------------------MAIN PROGRAM----------------------------------------------------------------*/
int main(void){
  int i,N=20; /*Repetitions*/
  int j,M=11; /*Matrix Size*/
  
  double sum1=0,sum2=0,sum12=0,sum22=0,time1=0,time2=0;
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  for(j=0;j<M;j++){
    Matrix A1(std::pow(2,j),std::pow(2,j));
    Eigen::MatrixXd A2;
    for(i=0;i<N;i++){ 
      time1=TimeTransposeDirect(A1,i);
      sum1+=time1;
      sum12+=time1*time1;
      time2=TimeTransposeEigen(A2,i,std::pow(2,j));
      sum2+=time2;
      sum22+=time2*time2;
    }
    A1.Destroyer();
    double mean1 = sum1/N;
    double sigma1 = std::sqrt(N*std::abs(sum12/N-mean1*mean1)/(N-1));
    
    double mean2 = sum2/N;
    double sigma2 = std::sqrt(N*std::abs(sum22/N-mean2*mean2)/(N-1));
    
    std::cout<<std::pow(2,j)<<"\t"<<mean1<<"\t"<<sigma1<<"\t"<<mean2<<"\t"<<sigma2<<std::endl;
  }
  return 0;
}
