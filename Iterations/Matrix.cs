using System;
using System.Linq;
using System.Text;

namespace MatrixLibrary
{
    public class Matrix : ICloneable
    {
        public double[][] Source { get; set; }

        public int Width => Source[0].Length;
        public int Height => Source.Length;

        public bool IsSquare => Width == Height;

        public Matrix(double[][] matrix)
        {
            Source = matrix;
        }

        public Matrix(double[] matrix)
        {
            Source = new double[matrix.Length][];
            for (var i = 0; i < matrix.Length; i++)
            {
                Source[i] = new double[1];
                Source[i][0] = matrix[i];
            }
        }

        public Matrix(int height, int width)
        {
            var matrix = new double[height][];
            for (var i = 0; i < height; i++)
                matrix[i] = new double[width];
            Source = matrix;
        }

        public Matrix(string[] matrix)
        {
            Source =
                matrix
                .Select(
                    s => s
                    .Split()
                    .Where(x => !string.IsNullOrEmpty(x))
                    .Select(x => double.Parse(x.Replace('.', ',')))
                    .ToArray())
                .ToArray();
        }

        public double this[int i, int j]
        {
            get { return Source[i][j]; }
            set { Source[i][j] = value; }
        }

        public Matrix Transpose()
        {
            var matrix = new Matrix(Width, Height);

            for (var i = 0; i < Width; i++)
                for (var j = 0; j < Height; j++)
                    matrix[i, j] = this[j, i];

            return matrix;
        }

        /// <summary>
        /// Меняет строки местами
        /// </summary>
        /// <param name="row1"></param>
        /// <param name="row2"></param>
        public void SwapRows(int row1, int row2)
        {
            var temp = Source[row1];
            Source[row1] = Source[row2];
            Source[row2] = temp;
            
        }

        /// <summary>
        /// Вычисляет максимальный элемент в столбце
        /// </summary>
        /// <param name="column"></param>
        /// <param name="startRow"></param>
        /// <returns></returns>
        public int GetMaxRowElement(int column, int startRow = 0)
        {
            var max = startRow;
            for (var i = startRow; i < Height; i++)
                if (Math.Abs(this[i, column]) > Math.Abs(this[max, column])) max = i;
            return max;
        }

        /// <summary>
        /// Вычисляет квадратную норму матрицы
        /// </summary>
        /// <returns></returns>
        public double GetCubicNorm()
        {
            var max = 0.0;
            for (var i = 0; i < Height; i++)
            {
                var sum = 0.0;
                for (var j = 0; j < Width; j++)
                    sum += Math.Abs(this[i, j]);
                if (sum > max) max = sum;
            }

            return max;
        }

        public double GetNorm()
        {
            return Source.Max(x => Math.Abs(x.Max()));
        }

        /// <summary>
        /// Вычисляет октаэдрическую норму матрицы
        /// </summary>
        /// <returns></returns>
        public double GetOctahedralNorm()
        {
            var max = 0.0;
            for (var j = 0; j < Width; j++)
            {
                var sum = 0.0;
                for (var i = 0; i < Height; i++)
                    sum += Math.Abs(this[i, j]);
                if (sum > max) max = sum;
            }

            return max;
        }

        /// <summary>
        /// Вычисляет единичную матрицу заданной размерности
        /// </summary>
        /// <param name="dimension"></param>
        /// <returns></returns>
        public static Matrix GetIdentityMatrix(int dimension)
        {
            var matrix = new Matrix(dimension, dimension);
            for (var i = 0; i < dimension; i++)
                matrix[i, i] = 1;
            return matrix;
        }

        /// <summary>
        /// Оператор умножения матриц
        /// </summary>
        /// <param name="multiplicand"></param>
        /// <param name="multiplier"></param>
        /// <returns></returns>
        public static Matrix operator *(Matrix multiplicand, Matrix multiplier)
        {
            if (multiplicand.Width != multiplier.Height)
                throw new ArgumentException();

            var matrix = new Matrix(multiplicand.Height, multiplier.Width);

            for (var i = 0; i < multiplicand.Height; i++)
                for (var j = 0; j < multiplier.Width; j++)
                {
                    for (var k = 0; k < multiplier.Height; k++)
                        matrix[i, j] += multiplicand[i, k] * multiplier[k, j];
                }

            return matrix;
        }

        /// <summary>
        /// Оператор умножение числа на матрицу
        /// </summary>
        /// <param name="multiplicand"></param>
        /// <param name="multiplier"></param>
        /// <returns></returns>
        public static Matrix operator *(double multiplicand, Matrix multiplier)
        {
            var matrix = new Matrix(multiplier.Height, multiplier.Width);

            for (var i = 0; i < multiplier.Height; i++)
                for (var j = 0; j < multiplier.Width; j++)
                    matrix[i, j] = multiplier[i, j] * multiplicand;

            return matrix;
        }

        /// <summary>
        /// Оператор умножение матрицы на число
        /// </summary>
        /// <param name="multiplicand"></param>
        /// <param name="multiplier"></param>
        /// <returns></returns>
        public static Matrix operator *(Matrix multiplicand, double multiplier)
        {
            return multiplier * multiplicand;
        }

        /// <summary>
        /// Оператор вычитания матриц
        /// </summary>
        /// <param name="subtrahend"></param>
        /// <param name="subtractor"></param>
        /// <returns></returns>
        public static Matrix operator -(Matrix subtrahend, Matrix subtractor)
        {
            if (subtrahend.Width != subtractor.Width ||
                subtrahend.Height != subtractor.Height)
                throw new ArgumentException();

            var result = new Matrix(subtractor.Height, subtractor.Width);
            for (var i = 0; i < subtrahend.Height; i++)
                for (var j = 0; j < subtractor.Width; j++)
                    result[i, j] = subtrahend[i, j] - subtractor[i, j];

            return result;
        }

        /// <summary>
        /// Оператор сложения матриц
        /// </summary>
        /// <param name="term1"></param>
        /// <param name="term2"></param>
        /// <returns></returns>
        public static Matrix operator +(Matrix term1, Matrix term2)
        {
            if (term1.Width != term2.Width ||
                term1.Height != term2.Height)
                throw new ArgumentException();

            var result = new Matrix(term2.Height, term2.Width);
            for (var i = 0; i < term1.Height; i++)
                for (var j = 0; j < term2.Width; j++)
                    result[i, j] = term1[i, j] + term2[i, j];

            return result;
        }

        public Vector ToVector()
        {
            return new Vector(Source.Select(x => x[0]).ToArray());
        }

        public override string ToString()
        {
            var builder = new StringBuilder();
            for (var i = 0; i < Height; i++)
            {
                for (var j = 0; j < Width; j++)
                {
                    builder.Append($"{this[i, j], 8:0.####}");
                    builder.Append(" ");
                }
                builder.Append("\r\n");
            }
            return builder.ToString();
        }
        
        public object Clone()
        {
            var array = new double[Height][];
            for (var i = 0; i < Height; i++)
                array[i] = new double[Width];
            for (var i = 0; i < Height; i++)
                for (var j = 0; j < Width; j++)
                    array[i][j] = Source[i][j];

            return new Matrix(array);
        }
    }
}
