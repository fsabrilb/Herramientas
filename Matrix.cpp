#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <chrono>
/*------------------------------------------------CLASS MATRIX----------------------------------------------------------------*/
class Matrix{
private:
  double *M; int N_rows,N_cols,N_rowsBlock,N_colsBlock;
public:
  /*Builder*/
  Matrix(int m,int n);
  /*Operator Overload*/
  Matrix operator*(const Matrix &A);
  /*Class Functions*/
  void Destroyer(void);
  void FillingRandom(int seed);
  void Transpose(void);
  void Print(void);
  Matrix TransposeBlocking(int p,int q);
  Matrix MultiplicationBlocking(Matrix &A,int p,int q);
};
/*------------------------------------------------INSIDE CLASS MATRIX---------------------------------------------------------*/
/*Build Matrix with parameters m (Number of rows), n (Number of columns)*/
Matrix::Matrix(int m,int n){
  N_rows=m; N_cols=n;
  M= new double [N_rows*N_cols];
}
/*Define Matrix Multiplication*/
Matrix Matrix::operator*(const Matrix &A){
  Matrix R(N_rows,A.N_cols);double Sum=0;
  for(int i=0;i<N_rows;i++){
    for(int j=0;j<A.N_cols;j++){
      for(int k=0;k<N_cols;k++){
	Sum += M[i*A.N_cols+k]*A.M[k*A.N_cols+j];
      }
     R.M[i*A.N_cols+j]=Sum;Sum=0;
    }
  }
  return R;
}
/*------------------------------------------------MATRIX'S FUNCTIONS-----------------------------------------------------------*/
/*Destroy Matrix*/
void Matrix::Destroyer(void){
  delete[] M;
}
/*Fill Matrix with random numbers between 0 and 100*/
void Matrix::FillingRandom(int seed){
  srand(seed+1);
  for(int i=0;i<N_rows;i++){
    for(int j=0;j<N_cols;j++){
      M[i*N_cols+j]=drand48()*100;
    }
  }
}
/*Calculate the transpose of the matrix*/
void Matrix::Transpose(void){
  Matrix R(N_cols,N_rows); 
  for(int i=0;i<N_cols;i++){
    for(int j=0;j<N_rows;j++){
      R.M[i*N_cols+j]=M[j*N_rows+i];
    }
  }
}
/*Print matrix*/
void Matrix::Print(void){
  for(int i=0;i<N_rows;i++){
    for(int j=0;j<N_cols;j++){
      std::cout<<M[i*N_cols+j]<<"\t\t";
    }
    std::cout<<"\n";
  }
}
/*TransposeBlocking*/
Matrix Matrix::TransposeBlocking(int p,int q){
  int a=N_rows,b=N_cols; N_rowsBlock=p; N_colsBlock=q;
  Matrix R(N_cols,N_rows);
  for(int i=0;i<a/p;i++){
    for(int j=0;j<b/q;j++){
      for(int k=0;k<p;k++){
	for(int l=0;l<q;l++){
	  R.M[l+k*a+j*p*a+i*q]=M[l*b+k+j*q+i*p*b]; //[l+j*q][k+i*p]->[l*b+k+j*q+i*p*b];->l*(a or b)+k+j*q+i*p*(b or a) 
	}
      }
    }
  }
  return R;
}
/*MultiplicationBlocking*/
Matrix Matrix::MultiplicationBlocking(Matrix &A,int p,int q){
  int a=N_rows,b=N_cols,c=A.N_rows,d=A.N_cols; N_rowsBlock=p; N_colsBlock=q;
  Matrix R(a,d);double Sum=0;
  for(int i=0;i<a/p;i++){
    for(int j=0;j<d/q;j++){
      for(int k=0;k<d/q;k++){
	for(int l=0;l<p;l++){
	  for(int m=0;m<q;m++){
	    for(int n=0;n<q;n++){
	      Sum+= M[l*a+n+i*p*a+k*q]*A.M[n*a+m+j*q+k*p*a];
	    }
	    R.M[l*a+m+j*q+i*p*a]+=Sum; Sum=0;
	  }
	}
      }
    }
  }
  return R;
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
double TimeMultiplicationDirect(Matrix &A, Matrix &B, int seed){
  std::chrono::time_point<std::chrono::system_clock> start, end;
  A.FillingRandom(seed); B.FillingRandom(seed+1);
  start = std::chrono::system_clock::now();  
  A*B;
  end   = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  return elapsed_seconds.count();
}
double TimeTransposeEigen(Eigen::MatrixXd &A,int seed, int size){
  std::chrono::time_point<std::chrono::system_clock> start, end;
  srand(seed+1); A=Eigen::MatrixXd::Random(size, size);
  start = std::chrono::system_clock::now();  
  Eigen::MatrixXd AT = A.transpose().eval();
  end   = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  return elapsed_seconds.count();
}
double TimeMultiplicationEigen(Eigen::MatrixXd &A,Eigen::MatrixXd B, int seed, int size){
  std::chrono::time_point<std::chrono::system_clock> start, end;
  srand(seed+1); A=Eigen::MatrixXd::Random(size, size);
  srand(seed+2); B=Eigen::MatrixXd::Random(size, size);
  start = std::chrono::system_clock::now();  
  A*B;
  end   = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  return elapsed_seconds.count();
}

double TimeTransposeBlocking(Matrix &A,int seed,int s){
  std::chrono::time_point<std::chrono::system_clock> start, end;
  A.FillingRandom(seed);
  start = std::chrono::system_clock::now();  
  A.TransposeBlocking(s,s);
  end   = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  return elapsed_seconds.count();
}
double TimeMultiplicationBlocking(Matrix &A, Matrix &B, int seed,int s){
  std::chrono::time_point<std::chrono::system_clock> start, end;
  A.FillingRandom(seed); B.FillingRandom(seed+1);
  start = std::chrono::system_clock::now();  
  A.MultiplicationBlocking(B,s,s);
  end   = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  return elapsed_seconds.count();
}

/*------------------------------------------------MAIN PROGRAM----------------------------------------------------------------*/
int main(void){
  int i,N=20; /*Repetitions*/
  int j,M=9; /*Matrix Size*/
  int k,L=7; /*Blocking Size*/
  
  double sum1=0,sum2=0,sum3=0,sum4=0,sum5=0,sum6=0;
  double sum12=0,sum22=0,sum32=0,sum42=0,sum52=0,sum62=0;
  double time1=0,time2=0,time3=0,time4=0,time5=0,time6=0;
  std::cout.precision(6);
  std::cout.setf(std::ios::scientific);
  for(j=0;j<M;j++){
    Matrix A1(std::pow(2,j),std::pow(2,j));
    Matrix B1(std::pow(2,j),std::pow(2,j));
    Eigen::MatrixXd A2;
    Eigen::MatrixXd B2;
    for(i=0;i<N;i++){ 
      time1=TimeTransposeDirect(A1,i);
      sum1+=time1;
      sum12+=time1*time1;
      time2=TimeTransposeEigen(A2,i,std::pow(2,j));
      sum2+=time2;
      sum22+=time2*time2;
      time3=TimeMultiplicationDirect(A1,B1,i);
      sum3+=time3;
      sum32+=time3*time3;
      time4=TimeMultiplicationEigen(A2,B2,i,std::pow(2,j));
      sum4+=time4;
      sum42+=time4*time4;
    }
    A1.Destroyer(); B1.Destroyer();
    double mean1 = sum1/N;
    double sigma1 = std::sqrt(N*std::abs(sum12/N-mean1*mean1)/(N-1));
    
    double mean2 = sum2/N;
    double sigma2 = std::sqrt(N*std::abs(sum22/N-mean2*mean2)/(N-1));

    double mean3 = sum3/N;
    double sigma3 = std::sqrt(N*std::abs(sum32/N-mean3*mean3)/(N-1));

    double mean4 = sum4/N;
    double sigma4 = std::sqrt(N*std::abs(sum42/N-mean4*mean4)/(N-1));
    
    std::cout<<std::pow(2,j)<<"\t"<<mean1<<"\t"<<sigma1<<"\t"<<mean2<<"\t"<<sigma2<<"\t"<<mean3<<"\t"<<sigma3
	     <<"\t"<<mean4<<"\t"<<sigma4<<std::endl;
  }
  std::cout<<std::endl;
  for(k=0;k<=L;k++){
    Matrix A3(std::pow(2,L),std::pow(2,L));
    Matrix B3(std::pow(2,L),std::pow(2,L));
    for(i=0;i<N;i++){
      time5=TimeTransposeBlocking(A3,i,std::pow(2,k));
      sum5+=time5;
      sum52+=time5*time5;
      time6=TimeMultiplicationBlocking(A3,B3,i,std::pow(2,k));
      sum6+=time6;
      sum62+=time6*time6;
    }
    A3.Destroyer(); B3.Destroyer();
    double mean5 = sum5/N;
    double sigma5 = std::sqrt(N*std::abs(sum52/N-mean5*mean5)/(N-1));
    double mean6 = sum6/N;
    double sigma6 = std::sqrt(N*std::abs(sum62/N-mean6*mean6)/(N-1));

    std::cout<<std::pow(2,k)<<"\t"<<mean5<<"\t"<<sigma5<<"\t"<<mean6<<"\t"<<sigma6<<std::endl;
  }
  /*
  Matrix A(4,4),B(4,4),C(4,4);
  std::cout.precision(6);
  std::cout.setf(std::ios::scientific);
  A.FillingRandom(1); B.FillingRandom(3); C=A*B;
  A.Print();std::cout<<std::endl;
  B.Print();std::cout<<std::endl;
  C.Print();std::cout<<std::endl;
  A.MultiplicationBlocking(B,1,1); std::cout<<"\n"<<std::endl;
  A.MultiplicationBlocking(B,2,2); std::cout<<"\n"<<std::endl;
  A.MultiplicationBlocking(B,4,4); std::cout<<"\n"<<std::endl;
  */
  return 0;
}
