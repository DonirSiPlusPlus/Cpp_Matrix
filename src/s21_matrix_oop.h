#include <math.h>

#include <iostream>

class S21Matrix {
 private:
  int rows_, cols_;  // Rows and columns
  double** matrix_;  // Pointer to the memory where the matrix is allocated

 public:
  /* Destructor & Constructors */
  ~S21Matrix();
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);

  /* Operations */
  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  double Determinant();
  S21Matrix CalcComplements();
  S21Matrix InverseMatrix();
  S21Matrix Minor(int row, int col);

  /* Overloading */
  bool operator==(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const double num);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const double num);
  S21Matrix& operator*=(const S21Matrix& other);
  double& operator()(int i, int j);

  /* Accesors & Mutators */
  int GetRows();
  int GetCols();
  void SetRows(int rows);
  void SetCols(int cols);

 private:
  double TriangleDeter();
  bool SwapRows(S21Matrix& Triangle, int row);
};
