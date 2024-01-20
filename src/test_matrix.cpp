#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

using namespace std;

TEST(ConstructorTest, default_) {
  S21Matrix A;
  int i = A.GetRows();
  int j = A.GetCols();
  EXPECT_EQ(i, j);
}

TEST(ConstructorTest, parameters) {
  S21Matrix A(4, 5);
  int i = A.GetRows();
  int j = A.GetCols();
  EXPECT_EQ(i, 4);
  EXPECT_EQ(j, 5);
}

TEST(ConstructorTest, parameters_throw_1) {
  EXPECT_ANY_THROW({ S21Matrix A(0, 3); });
}

TEST(ConstructorTest, parameters_throw_2) {
  EXPECT_ANY_THROW({ S21Matrix A(-3, -1); });
}

TEST(ConstructorTest, parameters_throw_3) {
  EXPECT_ANY_THROW({ S21Matrix A(10, -100); });
}

TEST(ConstructorTest, copy) {
  S21Matrix A(4, 5);
  S21Matrix B(A);
  int i = B.GetRows();
  int j = B.GetCols();
  EXPECT_EQ(i, 4);
  EXPECT_EQ(j, 5);
}

TEST(ConstructorTest, move) {
  S21Matrix A(2, 2);
  S21Matrix B(move(A));
  EXPECT_EQ(2, B.GetRows());
  EXPECT_EQ(2, B.GetCols());
  EXPECT_EQ(0, A.GetRows());
  EXPECT_EQ(0, A.GetCols());
}

TEST(OperationTest, EqMatrix) {
  S21Matrix A(2, 2);
  S21Matrix B(3, 2);
  EXPECT_FALSE(A.EqMatrix(B));

  S21Matrix C(5, 6);
  S21Matrix D(5, 6);
  for (int i = 0; i < C.GetRows(); i++)
    for (int j = 0; j < C.GetCols(); j++) C(i, j) = D(i, j) = i + j;

  EXPECT_TRUE(C.EqMatrix(D));
  C(2, 3) = 100;
  EXPECT_FALSE(C.EqMatrix(D));
}

TEST(OperationTest, Sum) {
  S21Matrix a(3, 3), b(3, 3), d(3, 3);
  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = 3;
  a(1, 0) = 0, a(1, 1) = 4, a(1, 2) = 5;
  a(2, 0) = 0, a(2, 1) = 0, a(2, 2) = 6;

  b(0, 0) = 1, b(0, 1) = 0, b(0, 2) = 0;
  b(1, 0) = 2, b(1, 1) = 0, b(1, 2) = 0;
  b(2, 0) = 3, b(2, 1) = 4, b(2, 2) = 1;

  d(0, 0) = 2, d(0, 1) = 2, d(0, 2) = 3;
  d(1, 0) = 2, d(1, 1) = 4, d(1, 2) = 5;
  d(2, 0) = 3, d(2, 1) = 4, d(2, 2) = 7;

  a.SumMatrix(b);
  EXPECT_TRUE(a.EqMatrix(d));
}

TEST(OperationTest, Sum_Throw) {
  S21Matrix a(3, 3), b(10, 6);
  EXPECT_THROW(a.SumMatrix(b), logic_error);  // !!!! leak
}

TEST(OperationTest, Sub) {
  S21Matrix a(3, 3), b(3, 3), d(3, 3);
  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = 3;
  a(1, 0) = 0, a(1, 1) = 4, a(1, 2) = 5;
  a(2, 0) = 0, a(2, 1) = 0, a(2, 2) = 6;

  b(0, 0) = 1, b(0, 1) = 0, b(0, 2) = 0;
  b(1, 0) = 2, b(1, 1) = 0, b(1, 2) = 0;
  b(2, 0) = 3, b(2, 1) = 4, b(2, 2) = 1;

  d(0, 0) = 0, d(0, 1) = 2, d(0, 2) = 3;
  d(1, 0) = -2, d(1, 1) = 4, d(1, 2) = 5;
  d(2, 0) = -3, d(2, 1) = -4, d(2, 2) = 5;

  a.SubMatrix(b);
  EXPECT_TRUE(a.EqMatrix(d));
}

TEST(OperationTest, Sub_Throw) {
  S21Matrix a(3, 5), b(9, 5);
  EXPECT_THROW(a.SubMatrix(b), logic_error);
}

TEST(OperationTest, MulNumber) {
  S21Matrix a(3, 3), b(3, 3);
  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = -0.03;
  a(1, 0) = 0, a(1, 1) = 4, a(1, 2) = -5;
  a(2, 0) = 0, a(2, 1) = 0, a(2, 2) = 6;

  b(0, 0) = 10, b(0, 1) = 20, b(0, 2) = -0.3;
  b(1, 0) = 0, b(1, 1) = 40, b(1, 2) = -50;
  b(2, 0) = 0, b(2, 1) = 0, b(2, 2) = 60;

  a.MulNumber(10);
  EXPECT_TRUE(a.EqMatrix(b));
}

TEST(OperationTest, MulMatrix_1) {
  S21Matrix a(3, 2), b(2, 3), d(3, 3);
  a(0, 0) = 1, a(0, 1) = 4;
  a(1, 0) = 2, a(1, 1) = 5;
  a(2, 0) = 3, a(2, 1) = 6;

  b(0, 0) = 1, b(0, 1) = -1, b(0, 2) = 1;
  b(1, 0) = 2, b(1, 1) = 3, b(1, 2) = 4;

  d(0, 0) = 9, d(0, 1) = 11, d(0, 2) = 17;
  d(1, 0) = 12, d(1, 1) = 13, d(1, 2) = 22;
  d(2, 0) = 15, d(2, 1) = 15, d(2, 2) = 27;

  a.MulMatrix(b);
  EXPECT_TRUE(a.EqMatrix(d));
}

TEST(OperationTest, MulMatrix_2) {
  S21Matrix a(4, 2), b(2, 4), d(4, 4);
  a(0, 0) = 1, a(0, 1) = 4, a(1, 0) = 2, a(1, 1) = 5;
  a(2, 0) = 3, a(2, 1) = 6;
  a(3, 0) = 5, a(3, 1) = 3;

  b(0, 0) = 1, b(0, 1) = -1, b(0, 2) = 1, b(0, 3) = 6.045;
  b(1, 0) = 2, b(1, 1) = 3, b(1, 2) = 4, b(1, 3) = 3.5;

  d(0, 0) = 9, d(0, 1) = 11, d(0, 2) = 17, d(0, 3) = 20.045;
  d(1, 0) = 12, d(1, 1) = 13, d(1, 2) = 22, d(1, 3) = 29.59;
  d(2, 0) = 15, d(2, 1) = 15, d(2, 2) = 27, d(2, 3) = 39.135;
  d(3, 0) = 11, d(3, 1) = 4, d(3, 2) = 17, d(3, 3) = 40.725;

  a.MulMatrix(b);
  EXPECT_TRUE(a.EqMatrix(d));
}

TEST(OperationTest, Mul_Throw_1) {
  S21Matrix a(4, 5), b(4, 5);
  EXPECT_THROW(a.MulMatrix(b), logic_error);
}

TEST(OperationTest, Mul_Throw_2) {
  S21Matrix a(4, 5);
  EXPECT_THROW(a.MulNumber(NAN), logic_error);
}

TEST(OperationTest, Transpose_1) {
  S21Matrix a(3, 2), b(2, 3);

  a(0, 0) = 1, a(0, 1) = 2;
  a(1, 0) = 3, a(1, 1) = 4;
  a(2, 0) = 5, a(2, 1) = 6;

  b(0, 0) = 1, b(0, 1) = 3, b(0, 2) = 5;
  b(1, 0) = 2, b(1, 1) = 4, b(1, 2) = 6;

  S21Matrix result = a.Transpose();
  EXPECT_TRUE(result.EqMatrix(b));
}

TEST(OperationTest, Transpose_2) {
  S21Matrix a(2, 2), b(2, 2);

  a(0, 0) = 1, a(0, 1) = 2;
  a(1, 0) = 4, a(1, 1) = 6;

  b(0, 0) = 1, b(0, 1) = 4;
  b(1, 0) = 2, b(1, 1) = 6;

  S21Matrix result = a.Transpose();
  EXPECT_TRUE(result.EqMatrix(b));
}

TEST(OperationTest, Deter_1) {
  S21Matrix a(1, 1);

  a(0, 0) = 1;

  double expected_result = 1;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_2) {
  S21Matrix a(1, 1);

  a(0, 0) = -1234.9;

  double expected_result = -1234.9;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_3) {
  S21Matrix a(1, 1);

  a(0, 0) = 0;

  double expected_result = 0;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_4) {
  S21Matrix a(2, 2);

  a(0, 0) = 312, a(0, 1) = 20;
  a(1, 0) = 654, a(1, 1) = 12;

  double expected_result = -9336;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_5) {
  S21Matrix a(2, 2);

  a(0, 0) = 3, a(0, 1) = 2;
  a(1, 0) = 6, a(1, 1) = 1;

  double expected_result = -9;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_6) {
  S21Matrix a(2, 2);

  a(0, 0) = 3, a(0, 1) = 0;
  a(1, 0) = 6, a(1, 1) = 1;

  double expected_result = 3;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_7) {
  S21Matrix a(3, 3);

  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = 3;
  a(1, 0) = 4, a(1, 1) = 0, a(1, 2) = 6;
  a(2, 0) = 7, a(2, 1) = 8, a(2, 2) = 9;

  double expected_result = 60;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_8) {
  S21Matrix a(4, 4);

  a(0, 0) = 0, a(0, 1) = 2, a(0, 2) = 3, a(0, 3) = 3;
  a(1, 0) = 4, a(1, 1) = 0, a(1, 2) = 6, a(1, 3) = 3;
  a(2, 0) = 7, a(2, 1) = 8, a(2, 2) = 0, a(2, 3) = 3;
  a(3, 0) = 7, a(3, 1) = 8, a(3, 2) = 32, a(3, 3) = 3;

  double expected_result = -3648;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_9) {
  S21Matrix a(4, 4);

  a(0, 0) = 0, a(0, 1) = 2, a(0, 2) = 3, a(0, 3) = 3;
  a(1, 0) = 4, a(1, 1) = 0, a(1, 2) = 6, a(1, 3) = 3;
  a(2, 0) = 7, a(2, 1) = 8, a(2, 2) = 0, a(2, 3) = 3;
  a(3, 0) = 7, a(3, 1) = 0, a(3, 2) = 32, a(3, 3) = 3;

  double expected_result = -2856;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_10) {
  S21Matrix a(4, 4);

  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = 3, a(0, 3) = 5;
  a(1, 0) = 3, a(1, 1) = 4, a(1, 2) = 0, a(1, 3) = 4;
  a(2, 0) = 8, a(2, 1) = 6, a(2, 2) = 0, a(2, 3) = 4;
  a(3, 0) = 6, a(3, 1) = 8, a(3, 2) = 9, a(3, 3) = 4;

  double expected_result = 510;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_11) {
  S21Matrix a(4, 4);

  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = 3, a(0, 3) = 4;
  a(1, 0) = 2, a(1, 1) = 3, a(1, 2) = 4, a(1, 3) = 1;
  a(2, 0) = 3, a(2, 1) = 4, a(2, 2) = 1, a(2, 3) = 2;
  a(3, 0) = 4, a(3, 1) = 1, a(3, 2) = 2, a(3, 3) = 3;

  double expected_result = 160;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_12) {
  S21Matrix a(4, 4);

  a(0, 0) = -2, a(0, 1) = 1, a(0, 2) = 3, a(0, 3) = 2;
  a(1, 0) = 3, a(1, 1) = 0, a(1, 2) = -1, a(1, 3) = 2;
  a(2, 0) = -5, a(2, 1) = 2, a(2, 2) = 3, a(2, 3) = 0;
  a(3, 0) = 4, a(3, 1) = -1, a(3, 2) = 2, a(3, 3) = -3;

  double expected_result = -80;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_13) {
  S21Matrix a(3, 3);

  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = 3;
  a(1, 0) = 4, a(1, 1) = 5, a(1, 2) = 6;
  a(2, 0) = 7, a(2, 1) = 8, a(2, 2) = 9;

  double expected_result = 0;
  double actual_result = a.Determinant();

  EXPECT_TRUE(fabs(expected_result - actual_result) < 1e-7);
}

TEST(OperationTest, Deter_14) {
  S21Matrix a(5, 5);

  a(0, 0) = 4.0, a(0, 1) = 1.0, a(0, 2) = 1.0, a(0, 3) = 2.0, a(0, 4) = 1.0;
  a(1, 0) = 1.0, a(1, 1) = 2.0, a(1, 2) = -1.0, a(1, 3) = 1.0, a(1, 4) = 1.0;
  a(2, 0) = 3.0, a(2, 1) = 1.0, a(2, 2) = 1.0, a(2, 3) = 1.0, a(2, 4) = 1.0;
  a(3, 0) = 2.0, a(3, 1) = 1.0, a(3, 2) = 1.0, a(3, 3) = 4.0, a(3, 4) = 1.0;
  a(4, 0) = 2.0, a(4, 1) = -1.0, a(4, 2) = 1.0, a(4, 3) = 1.0, a(4, 4) = 5.0;

  double result = a.Determinant();
  EXPECT_DOUBLE_EQ(result, 64.0);
}

TEST(OperationTest, Deter_Throw) {
  S21Matrix a(4, 5);
  EXPECT_THROW(a.Determinant(), logic_error);
}

TEST(OperationTest, Calc_Compl_1) {
  S21Matrix a(3, 3);

  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = 3;
  a(1, 0) = 0, a(1, 1) = 4, a(1, 2) = 2;
  a(2, 0) = 5, a(2, 1) = 2, a(2, 2) = 1;

  S21Matrix result = a.CalcComplements();

  EXPECT_DOUBLE_EQ(result(0, 0), 0);
  EXPECT_DOUBLE_EQ(result(0, 1), 10);
  EXPECT_DOUBLE_EQ(result(0, 2), -20);
  EXPECT_DOUBLE_EQ(result(1, 0), 4);
  EXPECT_DOUBLE_EQ(result(1, 1), -14);
  EXPECT_DOUBLE_EQ(result(1, 2), 8);
  EXPECT_DOUBLE_EQ(result(2, 0), -8);
  EXPECT_DOUBLE_EQ(result(2, 1), -2);
  EXPECT_DOUBLE_EQ(result(2, 2), 4);
}

TEST(OperationTest, Calc_Compl_2) {
  double m[3][3] = {{1, 2, 3}, {0, 4, 2}, {5, 2, 1}};
  double r[3][3] = {{0, 10, -20}, {4, -14, 8}, {-8, -2, 4}};

  S21Matrix a(3, 3), res(3, 3);

  for (int i = 0; i < a.GetRows(); i++)
    for (int j = 0; j < a.GetCols(); j++) a(i, j) = m[i][j];

  for (int i = 0; i < res.GetRows(); i++)
    for (int j = 0; j < res.GetCols(); j++) res(i, j) = r[i][j];

  S21Matrix my_res = a.CalcComplements();

  EXPECT_TRUE(my_res.EqMatrix(res));
}

TEST(OperationTest, CCompl_Throw_1) {
  S21Matrix a(3, 10);
  EXPECT_THROW(a.CalcComplements(), logic_error);
}

TEST(OperationTest, CCompl_Throw_2) {
  S21Matrix a(1, 1);
  EXPECT_THROW(a.CalcComplements(), logic_error);
}

TEST(OperationTest, inverse_1) {
  S21Matrix A(3, 3), res(3, 3);

  A(0, 0) = 2.0, A(0, 1) = 5.0, A(0, 2) = 7.0;
  A(1, 0) = 6.0, A(1, 1) = 3.0, A(1, 2) = 4.0;
  A(2, 0) = 5.0, A(2, 1) = -2.0, A(2, 2) = -3.0;

  res(0, 0) = 1.0, res(0, 1) = -1.0, res(0, 2) = 1.0;
  res(1, 0) = -38.0, res(1, 1) = 41.0, res(1, 2) = -34.0;
  res(2, 0) = 27.0, res(2, 1) = -29.0, res(2, 2) = 24.0;

  S21Matrix my_res = A.InverseMatrix();
  EXPECT_TRUE(my_res.EqMatrix(res));
}

TEST(OperationTest, inverse_2) {
  S21Matrix a(1, 1);

  a(0, 0) = 5;

  S21Matrix my_res = a.InverseMatrix();
  EXPECT_TRUE(fabs(my_res(0, 0) - 0.2 < 1e-7));
}

TEST(OperationTest, inverse_3) {
  S21Matrix A(4, 4), res(4, 4);

  A(0, 0) = 1, A(0, 1) = 2, A(0, 2) = 3, A(0, 3) = 4;
  A(1, 0) = 2, A(1, 1) = 3, A(1, 2) = 4, A(1, 3) = 1;
  A(2, 0) = 3, A(2, 1) = 4, A(2, 2) = 1, A(2, 3) = 2;
  A(3, 0) = 4, A(3, 1) = 1, A(3, 2) = 2, A(3, 3) = 3;

  res(0, 0) = -0.225, res(0, 1) = 0.025, res(0, 2) = 0.025, res(0, 3) = 0.275;
  res(1, 0) = 0.025, res(1, 1) = 0.025, res(1, 2) = 0.275, res(1, 3) = -0.225;
  res(2, 0) = 0.025, res(2, 1) = 0.275, res(2, 2) = -0.225, res(2, 3) = 0.025;
  res(3, 0) = 0.275, res(3, 1) = -0.225, res(3, 2) = 0.025, res(3, 3) = 0.025;

  S21Matrix my_res = A.InverseMatrix();

  EXPECT_TRUE(my_res.EqMatrix(res));
}

TEST(OperationTest, inverse_4) {
  S21Matrix A(4, 4), res(4, 4);

  A(0, 0) = -2, A(0, 1) = 1, A(0, 2) = 3, A(0, 3) = 2;
  A(1, 0) = 3, A(1, 1) = 0, A(1, 2) = -1, A(1, 3) = 2;
  A(2, 0) = -5, A(2, 1) = 2, A(2, 2) = 3, A(2, 3) = 0;
  A(3, 0) = 4, A(3, 1) = -1, A(3, 2) = 2, A(3, 3) = -3;

  res(0, 0) = -0.1, res(0, 1) = 23. / 80, res(0, 2) = 9. / 80,
         res(0, 3) = 1. / 8;
  res(1, 0) = -0.7, res(1, 1) = 71. / 80, res(1, 2) = 73. / 80,
         res(1, 3) = 1. / 8;
  res(2, 0) = 0.3, res(2, 1) = -9. / 80, res(2, 2) = -7. / 80,
         res(2, 3) = 1. / 8;
  res(3, 0) = 0.3, res(3, 1) = 1. / 80, res(3, 2) = -17. / 80,
         res(3, 3) = -1. / 8;

  S21Matrix my_res = A.InverseMatrix();
  EXPECT_TRUE(my_res.EqMatrix(res));
}

TEST(OperationTest, Inverse_Throw_1) {
  S21Matrix a(5, 5);
  EXPECT_THROW(a.InverseMatrix(), logic_error);
  a.SetCols(6);
  EXPECT_THROW(a.InverseMatrix(), logic_error);
}

TEST(OperationTest, Inverse_Throw_2) {
  S21Matrix a(1, 1);
  EXPECT_THROW(a.InverseMatrix(), logic_error);
}

TEST(OverloadTest, Sum) {
  S21Matrix a(3, 3), b(3, 3), d(3, 3);
  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = 3;
  a(1, 0) = 0, a(1, 1) = 4, a(1, 2) = 5;
  a(2, 0) = 0, a(2, 1) = 0, a(2, 2) = 6;

  b(0, 0) = 1, b(0, 1) = 0, b(0, 2) = 0;
  b(1, 0) = 2, b(1, 1) = 0, b(1, 2) = 0;
  b(2, 0) = 3, b(2, 1) = 4, b(2, 2) = 1;

  d(0, 0) = 2, d(0, 1) = 2, d(0, 2) = 3;
  d(1, 0) = 2, d(1, 1) = 4, d(1, 2) = 5;
  d(2, 0) = 3, d(2, 1) = 4, d(2, 2) = 7;

  a = a + b;
  EXPECT_TRUE(a.EqMatrix(d));
}

TEST(OverloadTest, Sum_Throw) {
  S21Matrix a(4, 5), b(4, 6);
  EXPECT_THROW(a + b, logic_error);
}

TEST(OverloadTest, Sub) {
  S21Matrix a(3, 3), b(3, 3), d(3, 3);
  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = 3;
  a(1, 0) = 0, a(1, 1) = 4, a(1, 2) = 5;
  a(2, 0) = 0, a(2, 1) = 0, a(2, 2) = 6;

  b(0, 0) = 1, b(0, 1) = 0, b(0, 2) = 0;
  b(1, 0) = 2, b(1, 1) = 0, b(1, 2) = 0;
  b(2, 0) = 3, b(2, 1) = 4, b(2, 2) = 1;

  d(0, 0) = 0, d(0, 1) = 2, d(0, 2) = 3;
  d(1, 0) = -2, d(1, 1) = 4, d(1, 2) = 5;
  d(2, 0) = -3, d(2, 1) = -4, d(2, 2) = 5;

  a = a - b;
  EXPECT_TRUE(a.EqMatrix(d));
}

TEST(OverloadTest, Sub_Throw) {
  S21Matrix a(4, 5), b(4, 6);
  EXPECT_THROW(a - b, logic_error);
}

TEST(OverloadTest, MulMatrix_1) {
  S21Matrix a(3, 3), b(3, 3);
  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = -0.03;
  a(1, 0) = 0, a(1, 1) = 4, a(1, 2) = -5;
  a(2, 0) = 0, a(2, 1) = 0, a(2, 2) = 6;

  b(0, 0) = 10, b(0, 1) = 20, b(0, 2) = -0.3;
  b(1, 0) = 0, b(1, 1) = 40, b(1, 2) = -50;
  b(2, 0) = 0, b(2, 1) = 0, b(2, 2) = 60;

  a = a * 10;
  EXPECT_TRUE(a.EqMatrix(b));
}

TEST(OverloadTest, MulMatrix_2) {
  S21Matrix a(3, 2), b(2, 3), d(3, 3);
  a(0, 0) = 1, a(0, 1) = 4;
  a(1, 0) = 2, a(1, 1) = 5;
  a(2, 0) = 3, a(2, 1) = 6;

  b(0, 0) = 1, b(0, 1) = -1, b(0, 2) = 1;
  b(1, 0) = 2, b(1, 1) = 3, b(1, 2) = 4;

  d(0, 0) = 9, d(0, 1) = 11, d(0, 2) = 17;
  d(1, 0) = 12, d(1, 1) = 13, d(1, 2) = 22;
  d(2, 0) = 15, d(2, 1) = 15, d(2, 2) = 27;

  a = a * b;
  EXPECT_TRUE(a.EqMatrix(d));
}

TEST(OverloadTest, MulNumber) {
  S21Matrix a(3, 2), b(3, 2);
  a(0, 0) = 1, a(0, 1) = 4;
  a(1, 0) = 2, a(1, 1) = 5;
  a(2, 0) = 3, a(2, 1) = 6;

  b(0, 0) = 10, b(0, 1) = 40;
  b(1, 0) = 20, b(1, 1) = 50;
  b(2, 0) = 30, b(2, 1) = 60;

  a = a * 10;
  EXPECT_TRUE(a.EqMatrix(b));
}

TEST(OverloadTest, Mul_Throw) {
  S21Matrix a(4, 5), b(4, 4);
  EXPECT_THROW(a * b, logic_error);
}

TEST(OverloadTest, Compare) {
  S21Matrix A(2, 2);
  S21Matrix B(3, 2);
  EXPECT_FALSE(A == B);

  S21Matrix C(5, 6);
  S21Matrix D(5, 6);
  for (int i = 0; i < C.GetRows(); i++)
    for (int j = 0; j < C.GetCols(); j++) C(i, j) = D(i, j) = i + j;

  EXPECT_TRUE(C == D);
}

TEST(OverloadTest, Appropriation) {
  S21Matrix A(2, 2);
  S21Matrix B(10, 7);
  A = B;
  EXPECT_TRUE(A == B);

  S21Matrix C(5, 6);
  S21Matrix D;
  for (int i = 0; i < C.GetRows(); i++)
    for (int j = 0; j < C.GetCols(); j++) C(i, j) = i + j;

  D = C;

  EXPECT_TRUE(C == D);
}

TEST(OverloadTest, Add_Assign) {
  S21Matrix a(3, 3), b(3, 3), d(3, 3);
  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = 3;
  a(1, 0) = 0, a(1, 1) = 4, a(1, 2) = 5;
  a(2, 0) = 0, a(2, 1) = 0, a(2, 2) = 6;

  b(0, 0) = 1, b(0, 1) = 0, b(0, 2) = 0;
  b(1, 0) = 2, b(1, 1) = 0, b(1, 2) = 0;
  b(2, 0) = 3, b(2, 1) = 4, b(2, 2) = 1;

  d(0, 0) = 2, d(0, 1) = 2, d(0, 2) = 3;
  d(1, 0) = 2, d(1, 1) = 4, d(1, 2) = 5;
  d(2, 0) = 3, d(2, 1) = 4, d(2, 2) = 7;

  a += b;
  EXPECT_TRUE(a.EqMatrix(d));
}

TEST(OverloadTest, AddA_Throw) {
  S21Matrix a(4, 5), b(4, 6);
  EXPECT_THROW(a += b, logic_error);
}

TEST(OverloadTest, Diff_Assign) {
  S21Matrix a(3, 3), b(3, 3), d(3, 3);
  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = 3;
  a(1, 0) = 0, a(1, 1) = 4, a(1, 2) = 5;
  a(2, 0) = 0, a(2, 1) = 0, a(2, 2) = 6;

  b(0, 0) = 1, b(0, 1) = 0, b(0, 2) = 0;
  b(1, 0) = 2, b(1, 1) = 0, b(1, 2) = 0;
  b(2, 0) = 3, b(2, 1) = 4, b(2, 2) = 1;

  d(0, 0) = 0, d(0, 1) = 2, d(0, 2) = 3;
  d(1, 0) = -2, d(1, 1) = 4, d(1, 2) = 5;
  d(2, 0) = -3, d(2, 1) = -4, d(2, 2) = 5;

  a -= b;
  EXPECT_TRUE(a.EqMatrix(d));
}

TEST(OverloadTest, DiffA_Throw) {
  S21Matrix a(4, 5), b(4, 3);
  EXPECT_THROW(a -= b, logic_error);
}

TEST(OverloadTest, MulNumber_Assign) {
  S21Matrix a(3, 2), b(3, 2);
  a(0, 0) = 1, a(0, 1) = 4;
  a(1, 0) = 2, a(1, 1) = 5;
  a(2, 0) = 3, a(2, 1) = 6;

  b(0, 0) = 10, b(0, 1) = 40;
  b(1, 0) = 20, b(1, 1) = 50;
  b(2, 0) = 30, b(2, 1) = 60;

  a *= 10;
  EXPECT_TRUE(a.EqMatrix(b));
}

TEST(OverloadTest, Mul_Assign_1) {
  S21Matrix a(3, 3), b(3, 3);
  a(0, 0) = 1, a(0, 1) = 2, a(0, 2) = -0.03;
  a(1, 0) = 0, a(1, 1) = 4, a(1, 2) = -5;
  a(2, 0) = 0, a(2, 1) = 0, a(2, 2) = 6;

  b(0, 0) = 10, b(0, 1) = 20, b(0, 2) = -0.3;
  b(1, 0) = 0, b(1, 1) = 40, b(1, 2) = -50;
  b(2, 0) = 0, b(2, 1) = 0, b(2, 2) = 60;

  a *= 10;
  EXPECT_TRUE(a.EqMatrix(b));
}

TEST(OverloadTest, Mul_Assign_2) {
  S21Matrix a(3, 2), b(2, 3), d(3, 3);
  a(0, 0) = 1, a(0, 1) = 4;
  a(1, 0) = 2, a(1, 1) = 5;
  a(2, 0) = 3, a(2, 1) = 6;

  b(0, 0) = 1, b(0, 1) = -1, b(0, 2) = 1;
  b(1, 0) = 2, b(1, 1) = 3, b(1, 2) = 4;

  d(0, 0) = 9, d(0, 1) = 11, d(0, 2) = 17;
  d(1, 0) = 12, d(1, 1) = 13, d(1, 2) = 22;
  d(2, 0) = 15, d(2, 1) = 15, d(2, 2) = 27;

  a *= b;
  EXPECT_TRUE(a.EqMatrix(d));
}

TEST(OverloadTest, MulA_Throw) {
  S21Matrix a(4, 5), b(4, 6);
  EXPECT_THROW(a *= b, logic_error);
}

TEST(OverloadTest, Index_Throw) {
  S21Matrix a(10, 2);
  EXPECT_THROW(a(1, 3), out_of_range);
  EXPECT_THROW(a(11, 1), out_of_range);
  a.SetRows(5);
  EXPECT_THROW(a(6, 2), out_of_range);
  a.SetCols(1);
  EXPECT_THROW(a(4, 2), out_of_range);
}

TEST(SetterTest, SetRows) {
  S21Matrix d(4, 4);
  d(0, 0) = 99, d(0, 1) = 11, d(0, 2) = 17, d(0, 3) = 20.045;
  d(1, 0) = 12, d(1, 1) = 13, d(1, 2) = 22, d(1, 3) = 29.59;
  d(2, 0) = 15, d(2, 1) = 15, d(2, 2) = 27, d(2, 3) = 39.135;
  d(3, 0) = 11, d(3, 1) = 4, d(3, 2) = 17, d(3, 3) = 40.725;

  d.SetRows(9);
  EXPECT_EQ(9, d.GetRows());
  for (int i = 4; i < d.GetRows(); i++)
    for (int j = 0; i < d.GetCols(); j++) EXPECT_EQ(0, d(i, j));

  d.SetRows(1);
  EXPECT_EQ(1, d.GetRows());
}

TEST(SetterTest, SetCols) {
  S21Matrix d(4, 4);
  d(0, 0) = 99, d(0, 1) = 11, d(0, 2) = 17, d(0, 3) = 20.045;
  d(1, 0) = 12, d(1, 1) = 13, d(1, 2) = 22, d(1, 3) = 29.59;
  d(2, 0) = 15, d(2, 1) = 15, d(2, 2) = 27, d(2, 3) = 39.135;
  d(3, 0) = 11, d(3, 1) = 4, d(3, 2) = 17, d(3, 3) = 40.725;

  d.SetCols(9);
  EXPECT_EQ(9, d.GetCols());

  for (int i = 0; i < d.GetRows(); i++)
    for (int j = 4; j < d.GetCols(); j++) EXPECT_EQ(0, d(i, j));

  d.SetCols(1);
  EXPECT_EQ(1, d.GetCols());
}

TEST(SetterTest, SetRows_Throw) {
  S21Matrix d(4, 4);
  EXPECT_ANY_THROW({ d.SetRows(-5); });
}

TEST(SetterTest, SetCols_Throw) {
  S21Matrix d(4, 4);
  EXPECT_ANY_THROW({ d.SetCols(0); });
}

int main() {
  testing::InitGoogleTest();

  return RUN_ALL_TESTS();
}