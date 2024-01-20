#include "s21_matrix_oop.h"

/* Деструктор */
S21Matrix::~S21Matrix() {
  if (matrix_) {
    for (int i = 0; i < this->rows_; i++) {
      if (matrix_[i]) delete[] matrix_[i];
    }
    delete[] matrix_;
  }
  rows_ = cols_ = 0;
}

/* Базовый конструктор */
S21Matrix::S21Matrix() {
  rows_ = cols_ = 1;
  this->matrix_ = new double*();
  this->matrix_[0] = new double();
}

/* Конструктор с параметрами строк и столбцов */
S21Matrix::S21Matrix(int rows, int cols) {
  if (rows > 0 && cols > 0) {
    this->rows_ = rows;
    this->cols_ = cols;
    this->matrix_ = new double*[rows]();
    for (int i = 0; i < rows; i++) this->matrix_[i] = new double[cols]();
  } else {
    throw std::underflow_error(
        "Число строк и столбцов должны быть больше нуля");
  }
}

/* Конструктор копирования */
S21Matrix::S21Matrix(const S21Matrix& other) {
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->matrix_ = new double*[rows_];
  for (int i = 0; i < this->rows_; i++)
    this->matrix_[i] = new double[this->cols_]();

  for (int i = 0; i < this->rows_; i++)
    for (int j = 0; j < this->cols_; j++)
      this->matrix_[i][j] = other.matrix_[i][j];
}

/* Конструктор переноса */
S21Matrix::S21Matrix(S21Matrix&& other) {
  this->rows_ = other.rows_;
  this->cols_ = other.rows_;
  this->matrix_ = other.matrix_;
  other.matrix_ = nullptr;
  other.rows_ = other.cols_ = 0;
}

/* Сравнение матриц */
bool S21Matrix::EqMatrix(const S21Matrix& other) {
  bool Equal_Status(true);
  if (this->rows_ == other.rows_ && this->cols_ == other.cols_) {
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        if (fabs(this->matrix_[i][j] - other.matrix_[i][j]) > 1e-7) {
          Equal_Status = false;
          break;
        }
      }
    }
  } else {
    Equal_Status = false;
  }
  return Equal_Status;
}

/* Сложение Матриц */
void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (this->rows_ == other.rows_ && this->cols_ == other.cols_) {
    for (int i = 0; i < this->rows_; i++)
      for (int j = 0; j < this->cols_; j++)
        this->matrix_[i][j] += other.matrix_[i][j];

  } else {
    throw std::logic_error("Матрицы не совпадают");
  }
}

/* Вычитание матриц */
void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (this->rows_ == other.rows_ && this->cols_ == other.cols_) {
    for (int i = 0; i < this->rows_; i++)
      for (int j = 0; j < this->cols_; j++)
        this->matrix_[i][j] -= other.matrix_[i][j];

  } else {
    throw std::logic_error("Матрицы не совпадают");
  }
}

/* Умножение матрицы на число */
void S21Matrix::MulNumber(const double num) {
  if (num == num) {
    for (int i = 0; i < this->rows_; i++)
      for (int j = 0; j < this->cols_; j++) this->matrix_[i][j] *= num;

  } else {
    throw std::invalid_argument("Num is NAN");
  }
}

/* Перемножение матриц */
void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (this->cols_ == other.rows_) {
    S21Matrix result(this->rows_, other.cols_);
    for (int i = 0; i < result.rows_; i++)
      for (int j = 0; j < result.cols_; j++)
        for (int c = 0; c < this->cols_; c++)
          result.matrix_[i][j] += this->matrix_[i][c] * other.matrix_[c][j];

    *this = result;
  } else {
    throw std::logic_error("Вычисление невозможно");
  }
}

/* Транспонирование */
S21Matrix S21Matrix::Transpose() {
  S21Matrix result(this->cols_, this->rows_);
  for (int i = 0; i < result.rows_; i++)
    for (int j = 0; j < result.cols_; j++)
      result.matrix_[i][j] = this->matrix_[j][i];

  return result;
}

/* Определитель */
double S21Matrix::Determinant() {
  double result(0);
  if (rows_ == cols_) {
    if (rows_ == 1) {
      result = matrix_[0][0];
    } else if (rows_ == 2) {
      result += matrix_[0][0] * matrix_[1][1];
      result -= matrix_[0][1] * matrix_[1][0];
    } else if (rows_ == 3) {
      result += matrix_[0][0] * matrix_[1][1] * matrix_[2][2];
      result += matrix_[0][1] * matrix_[1][2] * matrix_[2][0];
      result += matrix_[0][2] * matrix_[1][0] * matrix_[2][1];
      result -= matrix_[0][2] * matrix_[1][1] * matrix_[2][0];
      result -= matrix_[0][0] * matrix_[1][2] * matrix_[2][1];
      result -= matrix_[0][1] * matrix_[1][0] * matrix_[2][2];
    } else {
      result = TriangleDeter();
    }
  } else {
    throw std::logic_error("Матрица должна быть квадратной");
  }
  if (result != result) throw std::invalid_argument("Determinant is NAN");
  return result;
}

/* Матрица алгебраических дополнений */
S21Matrix S21Matrix::CalcComplements() {
  S21Matrix result;
  if (this->rows_ == this->cols_) {
    if (this->rows_ == 1) {
      throw std::logic_error("Вычисление невозможно");
    } else {
      result.SetRows(this->rows_);
      result.SetCols(this->cols_);
      for (int i = 0; i < this->rows_; i++) {
        for (int j = 0; j < this->cols_; j++) {
          S21Matrix minor = this->Minor(i, j);
          int sign = (i + j + 2) % 2 > 0 ? -1 : 1;
          result.matrix_[i][j] = minor.Determinant() * sign;
        }
      }
    }
  } else {
    throw std::logic_error("Матрица должна быть квадратной");
  }
  return result;
}

/* Обратная матрица */
S21Matrix S21Matrix::InverseMatrix() {
  S21Matrix result;
  if (this->rows_ != this->cols_) {
    throw std::logic_error("Матрица должна быть квадратной");
  } else if (this->rows_ == 1) {
    if (fabs(this->matrix_[0][0]) > 1e-6) {
      S21Matrix Inverse_Matrix(1, 1);
      Inverse_Matrix.matrix_[0][0] = 1. / this->matrix_[0][0];
      result = Inverse_Matrix;
    } else {
      throw std::logic_error("Вычисление невозможно");
    }
  } else {
    double Determinant_Of_This = this->Determinant();
    if (fabs(Determinant_Of_This) > 1e-6) {
      Determinant_Of_This = 1 / Determinant_Of_This;
      S21Matrix Inverse_Matrix(this->CalcComplements().Transpose());
      Inverse_Matrix.MulNumber(Determinant_Of_This);
      result = Inverse_Matrix;
    } else {
      throw std::logic_error("Вычисление невозможно");
    }
  }
  return result;
}

/* Перегрузка = (Присвоение) */
S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    this->~S21Matrix();
    this->matrix_ = new double*[other.rows_];
    for (int i = 0; i < other.rows_; i++) {
      this->matrix_[i] = new double[other.cols_];
      for (int j = 0; j < other.cols_; j++)
        this->matrix_[i][j] = other.matrix_[i][j];
    }
    this->rows_ = other.rows_;
    this->cols_ = other.cols_;
  }
  return *this;
}

/* Перегрузка + (Сложение) */
S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

/* Перегрузка - (вычитание) */
S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

/* Перегрузка * (умножение на число) */
S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix result(*this);
  result.MulNumber(num);
  return result;
}

/* Перегрузка * (перемножение матриц) */
S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

/* Перегрузка == (сравнение) */
bool S21Matrix::operator==(const S21Matrix& other) { return EqMatrix(other); }

/* Перегрузка += (Сложение) */
S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

/* Перегрузка -= (вычитание) */
S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

/* Перегрузка *= (умножение на число) */
S21Matrix& S21Matrix::operator*=(const double num) {
  this->MulNumber(num);
  return *this;
}

/* Перегрузка *= (перемножение матриц) */
S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

/* Перегрузка () (Индексирование) */
double& S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || j >= cols_)
    throw std::out_of_range("index is out of range");
  return matrix_[i][j];
}

/* Сеттер строк */
void S21Matrix::SetRows(int rows) {
  if (rows > 0) {
    S21Matrix result(rows, this->cols_);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < this->cols_; j++) {
        if (i >= this->rows_) {
          result.matrix_[i][j] = 0;
        } else {
          result.matrix_[i][j] = this->matrix_[i][j];
        }
      }
    }
    *this = result;
  } else {
    throw std::underflow_error("Число строк должно быть больше нуля");
  }
}

/* Сeттер столбцов */
void S21Matrix::SetCols(int cols) {
  if (cols > 0) {
    S21Matrix result(this->rows_, cols);
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < cols; j++) {
        if (j >= this->cols_) {
          result.matrix_[i][j] = 0;
        } else {
          result.matrix_[i][j] = this->matrix_[i][j];
        }
      }
    }
    *this = result;
  } else {
    throw std::underflow_error("Число столбцов должно быть больше нуля");
  }
}

int S21Matrix::GetRows() { return rows_; }
int S21Matrix::GetCols() { return cols_; }

/* Построение минора матрицы */
S21Matrix S21Matrix::Minor(int row, int col) {
  S21Matrix minor(this->rows_ - 1, this->cols_ - 1);
  int row_minor(0), col_minor(0);  // счетчики для прохода по минору
  for (int i = 0; i < this->rows_; i++) {
    if (i == row) continue;
    col_minor = 0;
    for (int j = 0; j < this->cols_; j++) {
      if (j == col) continue;
      minor.matrix_[row_minor][col_minor++] = this->matrix_[i][j];
    }
    row_minor++;
  }
  return minor;
}

/* Вычисление определителя через треугольный вид */
double S21Matrix::TriangleDeter() {
  double Determinant(1);
  int Swaps_Count(0);
  S21Matrix Triangle(*this);
  for (int i = 0; i < Triangle.rows_; i++) {
    if (fabs(matrix_[i][i]) < 1e-7) {
      bool was_swap = SwapRows(Triangle, i);
      if (!was_swap) {
        // если элемент на гл. диагонали остался 0
        Determinant = 0;
        break;
      }
      Swaps_Count++;
    }
    // зануление всего под гл. дигональю
    for (int j = i + 1; j < Triangle.rows_; j++) {
      double tmp = Triangle.matrix_[j][i];
      for (int c = i; c < Triangle.cols_; c++) {
        Triangle.matrix_[j][c] -=
            tmp * (Triangle.matrix_[i][c] / Triangle.matrix_[i][i]);
      }
    }
  }

  for (int i = 0; i < Triangle.rows_; i++)
    Determinant *= Triangle.matrix_[i][i];

  if (Swaps_Count % 2) Determinant *= -1;
  return Determinant;
}

/* Перестановка строк в матрице */
bool S21Matrix::SwapRows(S21Matrix& Triangle, int row) {
  bool check_swap(false);
  for (int i = Triangle.rows_ - 1; i > row; i--) {
    if (fabs(Triangle.matrix_[i][row]) > 1e-7) {
      for (int j = 0; j < Triangle.cols_; j++) {
        std::swap(Triangle.matrix_[row][j], Triangle.matrix_[i][j]);
      }
      check_swap = true;
      break;
    }
  }
  return check_swap;
}
