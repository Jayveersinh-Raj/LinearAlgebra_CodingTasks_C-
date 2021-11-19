#include <bits/stdc++.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
double **createMatrix(int,int);
void print(double **,int,int);
class Matrix
{
private:
    int rows;
    int columns;
    double **mat;
public:
    Matrix() {}

    Matrix(int rows,int columns)
    {
        this->rows = rows;
        this->columns = columns;
        mat = new double* [rows];
        for(int i=0; i<rows; i++)
        {
            mat[i] = new double [columns];
            memset(mat[i],0,columns);
        }
    };

    void setRows(int r)
    {
        this->rows = r;
    }

    int getRows()
    {
        return this->rows;
    }

    void setColumns(int c)
    {
        this->columns = c;
    }

    int getColumns()
    {
        return this->columns;
    }

    void setMatrix(double ** matrix)
    {
        this->mat = matrix;
    }

    double ** getMatrix()
    {
        return mat;
    }

    void addEntry(int row,int col,double value)
    {
        mat[row][col] = value;
    }

    Matrix operator +(Matrix  &obj)
    {
        Matrix ans;
        double **res = createMatrix(rows,columns);
        res = addMatrices(obj.getMatrix(),obj.rows,obj.columns);
        ans.setMatrix(res);
        ans.setRows(obj.rows);
        ans.setColumns(obj.columns);
        return ans;
    }

    double ** addMatrices(double **M,int rows,int columns)
    {

        if(rows!=this->rows || columns!=this->columns)
        {
            cout<<"Error: the dimensional problem occurred"<<endl;
            return nullptr;
        }
        else
        {
            double ** result = 0;
            result = createMatrix(rows,columns);
            for(int i=0; i<rows; i++)
            {
                for(int j=0; j<columns; j++)
                {
                    result[i][j] = this->mat[i][j] + M[i][j];
                }
            }
            return result;
        }
    }

    Matrix operator -(Matrix  &obj)
    {
        Matrix ans;
        double **res = createMatrix(rows,columns);
        res = subtractMatrices(obj.getMatrix(),obj.rows,obj.columns);
        ans.setMatrix(res);
        ans.setRows(obj.rows);
        ans.setColumns(obj.columns);
        return ans;
    }

    double ** subtractMatrices(double **M,int rows,int columns)
    {

        if(rows!=this->rows || columns!=this->columns)
        {
            cout<<"Error: the dimensional problem occurred"<<endl;
            return 0;
        }
        else
        {
            double ** result = 0;
            result = createMatrix(rows,columns);
            for(int i=0; i<rows; i++)
            {
                for(int j=0; j<columns; j++)
                {
                    result[i][j] = this->mat[i][j] - M[i][j];
                }
            }
            return result;
        }
    }

    Matrix operator *(Matrix  &obj)
    {
        Matrix ans;
        double **res = createMatrix(rows,columns);
        res = multiplyMatrices(obj.getMatrix(),obj.rows,obj.columns);
        ans.setMatrix(res);
        ans.setRows(rows);
        ans.setColumns(obj.columns);
        return ans;
    }

    double ** multiplyMatrices(double **M,int rows,int columns)
    {
        if(this->columns!=rows)
        {
            cout<<"Error: the dimensional problem occurred"<<endl;
            return 0;
        }
        else
        {
            double ** result = 0;
            result = createMatrix(this->rows,columns);
            for(int i=0; i<this->rows; i++)
            {
                for(int j=0; j<columns; j++)
                {
                    double res = 0;
                    for(int k=0; k<rows; k++)
                    {
                        res+=this->mat[i][k]*M[k][j];
                    }
                    result[i][j] = res;
                }

            }
            return result;
        }
    }

    Matrix transpose()
    {
        Matrix ans;
        double **res = createMatrix(columns,rows);
        for(int i=0; i<rows; i++)
            for(int j=0; j<columns; j++)
            {
                res[j][i] = this->mat[i][j];
            }
        ans.setMatrix(res);
        ans.setRows(columns);
        ans.setColumns(rows);
        return ans;
    }

    void print()
    {
        if(mat!=0)
            for(int i=0; i<rows; i++)
            {
                for(int j=0; j<columns; j++)
                {
                    if(fabs(mat[i][j]) > 0.0000001)
                        cout<<setprecision(4)<<fixed<<this->mat[i][j]<<' ';
                    else cout<<"0.0000"<<' ';
                }
                cout<<endl;
            }
    }
    friend void operator << (ostream &out,  Matrix &s);
    friend void operator >> (istream &in,  Matrix &s);
};

class ColumnVector : public Matrix
{
private:
public:
    ColumnVector() {}
    ColumnVector(int dim) : Matrix(dim,1) {}
    void multiply(double coef)
    {
        for(int i=0; i<getRows(); i++)
        {
            addEntry(i,0,getMatrix()[i][0]*coef);
        }
    }

    double diffLength(ColumnVector V)
    {
        double res = 0;
        for(int i=0; i<V.getRows(); i++)
        {
            res += (V.getMatrix()[i][0] - getMatrix()[i][0]) * (V.getMatrix()[i][0] - getMatrix()[i][0]);
        }
        return sqrt(res);
    }
    void operator = (Matrix  obj)
    {
        setMatrix(obj.getMatrix());
        setRows(obj.getRows());
        setColumns(obj.getColumns());
    }

    Matrix operator * (Matrix & obj)
    {
        Matrix firstOperand(this->getRows(),this->getColumns());
        firstOperand.setMatrix(getMatrix());
        Matrix secondOperand (obj.getRows(),obj.getColumns());
        secondOperand.setMatrix(obj.getMatrix());
        Matrix result = firstOperand * secondOperand;
        ColumnVector res (firstOperand.getRows());
        res.setMatrix(result.getMatrix());
        res = result;
        return res;
    }
};

class SquareMatrix : public Matrix
{
public:
    SquareMatrix() {}

    SquareMatrix(int dim) : Matrix(dim,dim) {}

    void operator = (Matrix  obj)
    {
        setMatrix(obj.getMatrix());
        setRows(obj.getRows());
        setColumns(obj.getColumns());
    }

    SquareMatrix operator + (SquareMatrix & obj)
    {
        Matrix firstOperand(this->getRows(),this->getColumns());
        firstOperand.setMatrix(getMatrix());
        Matrix secondOperand (obj.getRows(),obj.getColumns());
        secondOperand.setMatrix(obj.getMatrix());
        Matrix result = firstOperand + secondOperand;
        SquareMatrix res;
        res = result;
        return res;
    }

    SquareMatrix operator - (SquareMatrix & obj)
    {
        Matrix firstOperand(this->getRows(),this->getColumns());
        firstOperand.setMatrix(getMatrix());
        Matrix secondOperand (obj.getRows(),obj.getColumns());
        secondOperand.setMatrix(obj.getMatrix());
        Matrix result = firstOperand - secondOperand;
        SquareMatrix res;
        res = result;
        return res;
    }

    SquareMatrix operator * (SquareMatrix & obj)
    {
        Matrix firstOperand(this->getRows(),this->getColumns());
        firstOperand.setMatrix(getMatrix());
        Matrix secondOperand (obj.getRows(),obj.getColumns());
        secondOperand.setMatrix(obj.getMatrix());
        Matrix result = firstOperand * secondOperand;
        SquareMatrix res;
        res = result;
        return res;
    }

    ColumnVector operator * (ColumnVector & obj)
    {
        Matrix firstOperand(this->getRows(),this->getColumns());
        firstOperand.setMatrix(getMatrix());
        Matrix secondOperand (obj.getRows(),obj.getColumns());
        secondOperand.setMatrix(obj.getMatrix());
        Matrix result = firstOperand * secondOperand;
        ColumnVector res;
        res = result;
        return res;
    }
    double det();
    double** invert();
};

void  operator << (ostream &out,  Matrix &s)
{
    s.print();
}

void operator >> (istream &in,  Matrix &s)
{
    double value;
    for(int i=0; i<s.getRows(); i++)
        for(int j=0; j<s.getColumns(); j++)
        {
            cin>>value;
            s.addEntry(i,j,value);
        }
}

class IdentityMatrix : public SquareMatrix
{
public:
    IdentityMatrix() {}

    IdentityMatrix(int dim)
    {
        setRows(dim);
        setColumns(dim);
        setMatrix(createMatrix(dim,dim));
        for(int i=0; i<dim; i++)
            for(int j=0; j<dim; j++)
            {
                if(i==j)
                {
                    addEntry(i,i,1);
                }

                else
                {
                    addEntry(i,j,0);
                }
            }
    }
};

class EliminationMatrix : public SquareMatrix
{
public:
    EliminationMatrix() {}
    EliminationMatrix(int dim)
    {
        setRows(dim);
        setColumns(dim);
        setMatrix(createMatrix(dim,dim));
        for(int i=0; i<dim; i++)
            for(int j=0; j<dim; j++)
            {
                if(i==j)
                {
                    addEntry(i,i,1);
                }

                else
                {
                    addEntry(i,j,0);
                }
            }
    }

    void nullify (Matrix A, int i, int j, int pivRow, int pivCol)
    {

        double pivot = A.getMatrix()[pivRow][pivCol];
        addEntry(i,j, -1 * (A.getMatrix()[i][j] / pivot));
    }

    Matrix operator * (Matrix & obj)
    {
        Matrix firstOperand(this->getRows(),this->getColumns());
        firstOperand.setMatrix(getMatrix());
        Matrix secondOperand (obj.getRows(),obj.getColumns());
        secondOperand.setMatrix(obj.getMatrix());
        Matrix result = firstOperand * secondOperand;
        SquareMatrix res;
        res = result;
        return res;
    }
};

class PermutationMatrix : public SquareMatrix
{
public:
    PermutationMatrix(int dim, int i1, int i2)
    {
        setRows(dim);
        setColumns(dim);
        setMatrix(createMatrix(dim,dim));
        for(int i=0; i<dim; i++)
            for(int j=0; j<dim; j++)
            {
                if(i==j&&i!=i1&&i!=i2)
                {
                    addEntry(i,i,1);
                }

                else
                {
                    addEntry(i,j,0);
                }
            }
        addEntry(i1,i2,1);
        addEntry(i2,i1,1);
    }

    Matrix operator * (Matrix & obj)
    {
        Matrix firstOperand(this->getRows(),this->getColumns());
        firstOperand.setMatrix(getMatrix());
        Matrix secondOperand (obj.getRows(),obj.getColumns());
        secondOperand.setMatrix(obj.getMatrix());
        Matrix result = firstOperand * secondOperand;
        SquareMatrix res;
        res = result;
        return res;
    }
};

double SquareMatrix::det ()
{
    double** grid = getMatrix();
    SquareMatrix matrix(getRows());
    matrix.setMatrix(grid);
    double result = 1;
    int i = 0;
    int j = 0;
    int steps = 1;
    while (i<getRows())
    {
        double currentMax = fabs(matrix.getMatrix()[i][j]);
        int currentMaxRow = i;
        for(int row = i + 1; row < getRows(); row++)
        {
            if(fabs(matrix.getMatrix() [row][j]) > currentMax)
            {
                currentMax = fabs(matrix.getMatrix()[row][j]);
                currentMaxRow = row;
            }
        }

        if(currentMax == 0 )
        {

            return 0;
        }

        if(currentMaxRow!=i)
        {
            PermutationMatrix permutation(getRows(),i,currentMaxRow);
            matrix = permutation * matrix;
            steps++;
        }

        double flag = 0;
        for(int row = i + 1; row<getRows(); row++)
        {
            flag += fabs(matrix.getMatrix()[row][j]);
        }

        if(flag!=0)
            for(int row = i + 1; row < getRows(); row++)
            {
                EliminationMatrix elimination(getRows());
                elimination.nullify(matrix,row,j,i,j);
                matrix = elimination * matrix;
                steps++;
            }
        i++;
        j++;
    }
    grid = matrix.getMatrix();
    for(int i=0; i<getColumns(); i++)
    {
        result *= grid[i][i];
    }
    return result;
}

double** SquareMatrix::invert()
{
    double** grid = getMatrix();
    SquareMatrix matrix(getRows());
    matrix.setMatrix(grid);
    IdentityMatrix I (getRows());
    Matrix augmented (getRows(), 2 * getColumns());
    double** augmentedMatrix = createMatrix(getRows(), 2 * getColumns());
    for(int i=0; i<getRows(); i++)
        for(int j=0; j<getColumns(); j++)
        {
            augmentedMatrix[i][j] = grid[i][j];
        }

    for(int i=0; i<getRows(); i++)
        for(int j=getColumns(); j<2 * getColumns(); j++)
        {
            augmentedMatrix[i][j] = I.getMatrix()[i][j - getColumns()];
        }
    augmented.setMatrix(augmentedMatrix);
    int i = 0;
    int j = 0;
    int steps = 0;
    steps++;
    while (i<getRows())
    {
        double currentMax = fabs(augmented.getMatrix()[i][j]);
        int currentMaxRow = i;
        for(int row = i + 1; row < getRows(); row++)
        {
            if(fabs(augmented.getMatrix() [row][j]) > currentMax)
            {
                currentMax = fabs(augmented.getMatrix()[row][j]);
                currentMaxRow = row;
            }
        }

        if(currentMax == 0 )
        {

            return nullptr;
        }
        if(currentMaxRow!=i)
        {
            PermutationMatrix permutation(getRows(),i,currentMaxRow);
            augmented = permutation * augmented;
            steps++;
        }

        double flag = 0;
        for(int row = i + 1; row<getRows(); row++)
        {
            flag += fabs(augmented.getMatrix()[row][j]);
        }

        if(flag!=0)
            for(int row = i + 1; row < getRows(); row++)
            {
                EliminationMatrix elimination(getRows());
                elimination.nullify(augmented,row,j,i,j);
                augmented = elimination * augmented;
                steps++;

            }
        i++;
        j++;
    }

    i--;
    j--;
    while (i >= 0)
    {
        double flag = 0;
        for(int row = i - 1; row >= 0; row--)
        {
            flag += fabs(augmented.getMatrix()[row][j]);
        }
        if(flag!=0)
            for(int row = i - 1; row >= 0; row--)
            {
                EliminationMatrix elimination(getRows());
                elimination.nullify(augmented,row,j,i,j);
                augmented = elimination * augmented;
                steps++;
            }
        i--;
        j--;
    }
    augmentedMatrix = augmented.getMatrix();
    for(int i=0; i<getRows(); i++)
    {
        double den = augmentedMatrix[i][i];
        augmentedMatrix[i][i] = 1;
        for(int j=getColumns(); j<2 * getColumns(); j++)
        {
            augmentedMatrix[i][j] /= den;
        }
    }
    augmented.setMatrix(augmentedMatrix);
    double** res = createMatrix(getRows(),getColumns());
    for(int i=0; i<getRows(); i++)
    {
        for(int j=getColumns(); j<2 * getColumns(); j++)
        {
            res[i][j-getColumns()] = augmentedMatrix[i][j];
        }
    }
    return res;
}

int main()
{
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
    int m = 15,n = 4;
    double x[] = {-6.5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,6.5};
    double y[] = {2.3,5.321,5.43,2.56,1.276,0.635,5.432,-2.334,-2.432,3.4,4.3,1.54,6.45,7.34,1.57,2.657};
    double** arr = createMatrix(m,1);
    for(int i=0; i<m; i++)
    {
        arr[i][0] = y[i];
    }
    double** grid = createMatrix(m,n+1);
    Matrix A (m,n+1);
    SquareMatrix B(A.getColumns());
    SquareMatrix B_inverse(A.getColumns());
    ColumnVector b (m);
    ColumnVector beta(A.getColumns());
    ColumnVector answer(A.getColumns());
    b.setMatrix(arr);
    for(int i=0; i<m; i++)
    {
        grid[i][0] = 1;
    }
    for(int j=1; j<=n; j++)
    {
        for(int i=0; i<m; i++)
        {
            grid[i][j] = x[i] * grid[i][j-1];
        }
    }
    A.setMatrix(grid);
    Matrix A_t(n+1,m);
    A_t = A.transpose();
    B = A_t * A;
    double** inv = createMatrix(A.getColumns(),A.getColumns());
    inv =  B.invert();
    B_inverse.setMatrix(inv);
    beta = A_t * b;
    answer = B_inverse * beta;
    if(pipe != nullptr){
        fprintf(pipe, "%s\n", "plot '-' using 1:2 title 'plot' with lines");
        for(int i=0;i<m;i++){
            double xpoint = x[i];
            double ypoint = 0;
            double temp = xpoint;
            ypoint += answer.getMatrix()[0][0];
            for(int j=1;j<answer.getRows();j++){
                ypoint += answer.getMatrix()[j][0] * temp;
                temp *= xpoint;
            }
            fprintf(pipe, "%f\t%f\n",xpoint,ypoint);
        }
        fprintf(pipe, "%s\n","e");
        fflush(pipe);
        _pclose(pipe);
    }
    return 0;
}
double** createMatrix(int rows,int columns)
{
    double **res = new double* [rows];
    for(int i=0; i<rows; i++)
    {
        res[i] = new double [columns];
    }
    return res;
}

