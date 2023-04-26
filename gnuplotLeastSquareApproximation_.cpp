#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

class ColumnVector;

class Matrix {
public:
    int rows, columns;
    double **values;

    Matrix(int rows, int cols);

    Matrix(const Matrix &);

    ~Matrix();

    friend std::ostream &operator<<(std::ostream &out, const Matrix &m);

    friend std::istream &operator>>(std::istream &in, Matrix &m);

    Matrix &operator=(const Matrix &obj);

    Matrix operator+(const Matrix &obj);

    Matrix operator-(const Matrix &obj);

    Matrix operator*(const Matrix &obj) const;

    ColumnVector operator*(const ColumnVector &obj);

    Matrix transpose() const;
};

void printTwoMatrices(Matrix &a, Matrix &b);

class ColumnVector : public Matrix {
public:
    explicit ColumnVector(int n) : Matrix(n, 1) {};

    ColumnVector operator+(const ColumnVector &obj);

    ColumnVector operator-(const ColumnVector &obj);

    ColumnVector operator*(const ColumnVector &obj);
};

class SquareMatrix : public Matrix {
public:
    explicit SquareMatrix(int n) : Matrix(n, n) {};

    SquareMatrix operator+(const SquareMatrix &obj);

    SquareMatrix operator-(const SquareMatrix &obj);

    SquareMatrix operator*(const SquareMatrix &obj);

    SquareMatrix transpose();

    SquareMatrix gaussianElimination();

    SquareMatrix inverse();

    ColumnVector solveLinear(ColumnVector &columnVector);

    double determinant();
};

class IdentityMatrix : public SquareMatrix {
public:
    explicit IdentityMatrix(int n);
};

class EliminationMatrix : public SquareMatrix {
public:
    explicit EliminationMatrix(int row, int column, SquareMatrix &a);
};

class PermutationMatrix : public SquareMatrix {
public:
    explicit PermutationMatrix(int n, int row1, int row2);
};

Matrix::Matrix(int rows, int cols) {
    this->rows = rows;
    this->columns = cols;
    values = new double *[rows];
    for (int i = 0; i < rows; i++) {
        values[i] = new double[this->columns];
    }
}

Matrix::Matrix(const Matrix &m) {
    rows = m.rows;
    columns = m.columns;
    values = new double *[rows];
    for (int i = 0; i < rows; ++i) {
        values[i] = new double[columns];
        for (int j = 0; j < columns; ++j) {
            values[i][j] = m.values[i][j];
        }
    }
}

Matrix::~Matrix() {
    for (int i = 0; i < rows; ++i) {
        delete[] values[i];
    }
    delete[] values;
}

std::ostream &operator<<(std::ostream &out, const Matrix &m) {
    for (int i = 0; i < m.rows; ++i) {
        for (int j = 0; j < m.columns; ++j) {
            if (abs(m.values[i][j]) < 0.00005) {
                cout << "0.0000 ";
                continue;
            }
            cout << fixed << setprecision(4) << m.values[i][j] << " ";
        }
        cout << "\n";
    }
    return out;
}

std::istream &operator>>(std::istream &in, Matrix &m) {
    for (int i = 0; i < m.rows; ++i) {
        for (int j = 0; j < m.columns; ++j) {
            cin >> m.values[i][j];
        }
    }
    return in;
}

Matrix &Matrix::operator=(const Matrix &obj) {
    if (this == &obj)
        return *this;

    for (int i = 0; i < rows; ++i) {
        delete[] values[i];
    }
    delete[] values;
    this->rows = obj.rows;
    this->columns = obj.columns;
    values = new double *[rows];
    for (int i = 0; i < rows; ++i) {
        values[i] = new double[columns];
        for (int j = 0; j < columns; ++j) {
            this->values[i][j] = obj.values[i][j];
        }
    }
    return *this;
}

Matrix Matrix::operator+(const Matrix &obj) {
    Matrix m = *this;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            m.values[i][j] = this->values[i][j] + obj.values[i][j];
        }
    }
    return m;
}

Matrix Matrix::operator-(const Matrix &obj) {
    Matrix m = *this;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            m.values[i][j] = this->values[i][j] - obj.values[i][j];
        }
    }
    return m;
}

Matrix Matrix::operator*(const Matrix &obj) const {
    Matrix m = Matrix(this->rows, obj.columns);
    for (int i = 0; i < this->rows; ++i) {
        for (int j = 0; j < obj.columns; ++j) {
            m.values[i][j] = 0;
            for (int k = 0; k < this->columns; k++) {
                m.values[i][j] += this->values[i][k] * obj.values[k][j];
            }
        }
    }
    return m;
}


ColumnVector Matrix::operator*(const ColumnVector &obj) {
    Matrix b = *(Matrix *) &obj;
    Matrix a = *(Matrix *) this;
    Matrix res = a * b;
    ColumnVector vector = *(ColumnVector *) &res;
    return vector;
}


Matrix Matrix::transpose() const {
    Matrix m = Matrix(columns, rows);
    for (int i = 0; i < this->rows; i++)
        for (int j = 0; j < this->columns; j++)
            m.values[j][i] = this->values[i][j];
    return m;
}

ColumnVector ColumnVector::operator+(const ColumnVector &obj) {
    Matrix a = *(ColumnVector *) this;
    Matrix b = *(ColumnVector *) &obj;
    Matrix d = a + b;
    return *(ColumnVector *) &d;
}

ColumnVector ColumnVector::operator-(const ColumnVector &obj) {
    Matrix b = *(Matrix *) this;
    Matrix a = *(Matrix *) &obj;
    Matrix e = b - a;
    return *(ColumnVector *) &e;
}

ColumnVector ColumnVector::operator*(const ColumnVector &obj) {
    Matrix c = *(ColumnVector *) this;
    Matrix a = *(ColumnVector *) &obj;
    Matrix f = c * a;
    return *(ColumnVector *) &f;
}

SquareMatrix SquareMatrix::operator+(const SquareMatrix &obj) {
    Matrix a = *(Matrix *) this;
    Matrix b = *(Matrix *) &obj;
    Matrix d = a + b;
    return *(SquareMatrix *) &d;
}

SquareMatrix SquareMatrix::operator-(const SquareMatrix &obj) {
    Matrix b = *(Matrix *) this;
    Matrix a = *(Matrix *) &obj;
    Matrix e = b - a;
    return *(SquareMatrix *) &e;
}

SquareMatrix SquareMatrix::operator*(const SquareMatrix &obj) {
    Matrix c = *(Matrix *) this;
    Matrix a = *(Matrix *) &obj;
    Matrix f = c * a;
    return *(SquareMatrix *) &f;
}

SquareMatrix SquareMatrix::transpose() {
    Matrix a = *(Matrix *) this;
    Matrix g = a.transpose();
    return *(SquareMatrix *) &g;
}

double SquareMatrix::determinant() {
    Matrix eliminated = gaussianElimination();
    double det = 1;
    for (int i = 0; i < rows; ++i) {
        det *= eliminated.values[i][i];
    }
    return det;
}

IdentityMatrix::IdentityMatrix(int n) : SquareMatrix(n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            this->values[i][j] = i == j ? 1 : 0;
        }
    }
}

EliminationMatrix::EliminationMatrix(int row, int column, SquareMatrix &a) : SquareMatrix(a.rows) {
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.rows; j++) {
            this->values[i][j] = i == j ? 1 : 0;
        }
    }
    this->values[row][column] = -a.values[row][column] / a.values[column][column];
}

PermutationMatrix::PermutationMatrix(int n, int row1, int row2) : SquareMatrix(n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            this->values[i][j] = i == j ? 1 : 0;
        }
    }
    for (int i = 0; i < n; ++i) {
        double tmp = this->values[row1][i];
        this->values[row1][i] = this->values[row2][i];
        this->values[row2][i] = tmp;
    }
}

SquareMatrix SquareMatrix::gaussianElimination() {
    SquareMatrix result = *this;
    int step = 1;
    for (int i = 0; i < rows; i++) {
        int maxRow = i;
        double maxV = abs(result.values[i][i]);
        for (int j = i + 1; j < rows; j++) {
            if (abs(result.values[j][i]) > maxV) {
                maxV = abs(result.values[j][i]);
                maxRow = j;
            }
        }

        if (maxRow != i) {
            auto permutationMatrix = PermutationMatrix(rows, i, maxRow);
            result = permutationMatrix * result;
            cout << "step #" << step++ << ": permutation" << endl;
            cout << result;
        }

        for (int j = i + 1; j < rows; j++) {
            if (result.values[j][i] == 0) continue;
            auto eliminationMatrix = EliminationMatrix(j, i, result);
            result = eliminationMatrix * result;
            cout << "step #" << step++ << ": elimination" << endl;
            cout << result;
        }
    }
    return result;
}

SquareMatrix SquareMatrix::inverse() {
    SquareMatrix original = *this;
    SquareMatrix augmented = IdentityMatrix(rows);

    for (int i = 0; i < rows; i++) {
        int maxRow = i;
        double maxV = abs(original.values[i][i]);
        for (int j = i + 1; j < rows; j++) {
            if (abs(original.values[j][i]) > maxV) {
                maxV = abs(original.values[j][i]);
                maxRow = j;
            }
        }

        if (maxRow != i) {
            auto permutationMatrix = PermutationMatrix(rows, i, maxRow);
            original = permutationMatrix * original;
            augmented = permutationMatrix * augmented;
        }

        for (int j = i + 1; j < rows; j++) {
            if (original.values[j][i] == 0) continue;
            auto eliminationMatrix = EliminationMatrix(j, i, original);
            original = eliminationMatrix * original;
            augmented = eliminationMatrix * augmented;
        }
    }

    for (int i = rows - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            if (original.values[j][i] == 0) continue;
            auto eliminationMatrix = EliminationMatrix(j, i, original);
            original = eliminationMatrix * original;
            augmented = eliminationMatrix * augmented;
        }
    }

    for (int i = 0; i < rows; i++) {
        double factor = original.values[i][i];
        original.values[i][i] = 1;
        for (int j = 0; j < rows; j++) {
            augmented.values[i][j] /= factor;
        }
    }

    return augmented;
}

void printAnswer(Matrix a, ColumnVector b, FILE* pipe){
    cout << "A:" << endl;
    cout << a;
    auto A_T_mul_A_raw = a.transpose() * a;
    auto A_T_mul_A = *((SquareMatrix *)&A_T_mul_A_raw);
    cout << "A_T*A:" << endl;
    cout << A_T_mul_A;
    auto A_T_mul_A_degree_min_1 = A_T_mul_A.inverse();
    cout << "(A_T*A)^-1:" << endl;
    cout << A_T_mul_A_degree_min_1;
    auto A_T_mul_b = a.transpose() * b;
    cout << "A_T*b:" << endl;
    cout << A_T_mul_b;
    cout << "x~:" << endl;
    auto coefficients = (Matrix) A_T_mul_A_degree_min_1 * A_T_mul_b;

    if (pipe != NULL) {

        if (coefficients.rows == 2){
            fprintf(pipe, "plot %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n",
                    coefficients.values[0][0], coefficients.values[1][0]);
            for (int j = 0; j < 100; ++j) {
                double x = j;
                double y;
                y = coefficients.values[0][0] * x + coefficients.values[1][0];

                fprintf(pipe, "%f\t%f\n", x, y);
            }
        }

        if (coefficients.rows == 3){
            fprintf(pipe, "plot %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n",
                    coefficients.values[0][0], coefficients.values[1][0], coefficients.values[2][0]);
            for (int j = 0; j < 100; ++j) {
                double x = j;
                double y;
                y = coefficients.values[0][0] * x * x + coefficients.values[1][0] * x + coefficients.values[2][0];

                fprintf(pipe, "%f\t%f\n", x, y);
            }
        }

        if (coefficients.rows == 4){
            fprintf(pipe, "plot %lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n",
                    coefficients.values[0][0], coefficients.values[1][0], coefficients.values[2][0], coefficients.values[3][0]);
            for (int j = 0; j < 100; ++j) {
                double x = j;
                double y;
                y = coefficients.values[0][0] * x * x * x + coefficients.values[1][0] * x * x
                    + coefficients.values[2][0] * x + coefficients.values[3][0];

                fprintf(pipe, "%f\t%f\n", x, y);
            }
        }

        fprintf(pipe, "e\n");
        fflush(pipe);

    }
}

int main() {

    FILE* pipe = _popen(GNUPLOT_NAME, "w");

    int n;
    cin >> n;

    Matrix a(n, 2);

    ColumnVector b(n);
    ColumnVector t(n);

    for (int i = 0; i < n; ++i) {
        double t1, b1;
        cin >> t1 >> b1;
        t.values[i][0] = t1;
        b.values[i][0] = b1;

        a.values[i][0] = 1;
        a.values[i][1] = t1;
    }

    int degree;
    cin >> degree;

    if (degree>1){
        Matrix a1(n, degree+1);
        for (int i = 0; i < n; ++i) {
            double tmp = a.values[i][1];
            for (int j = 0; j < degree+1; ++j) {
                if (j==0){
                    a1.values[i][0] = 1;
                } else if (j==1){
                    a1.values[i][j] = tmp;
                } else {
                    tmp *= tmp;
                    a1.values[i][j] = tmp;
                }
            }
        }
        printAnswer(a1, b, pipe);
    } else{
        printAnswer(a, b, pipe);
    }

    _pclose(pipe);

    return 0;
}