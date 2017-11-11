using System;
using System.Windows.Forms;
using System.IO;
using System.Collections.Generic;

namespace SLAU
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            SolveButton.Enabled = true;
        }
        public const int  n=5; //Размерности
        public double[,] A = new double[n+1 , n + 2]; //Исходная, единичная
        public double[] b = new double[2];  //Вектор b
        public double[] solution; //Для решения
        public int[] kol;
        const double E = 0.0001; //Точность
        public double chis;
        List<string> strList = new List<string>(); 

        //Вывод вектора
        public void Write(double[] vec)
        {
            for (int i = 0; i < n; i++)
            {
                string str1 = "";
                str1 = String.Format("{0,-18}", vec[i]);
                strList.Add(str1);
            }
            strList.Add("");
        }

        //Выводим матрицу
        public void WriteMas(double[,] Mas)
        {
            if (Mas == A)
                strList.Add("Матрица А:");
            for (int i = 0; i <= n; i++)
            {
                string str1 = "";
                for (int j = 0; j < n+2; j++)
                    str1 += String.Format("{0,-15} ", Mas[i, j]);
                strList.Add(str1);
            }
            strList.Add("");
        }

        //Сохраняем решение в файл
        public void SaveFile()
        {
            SaveFileDialog saveFile = new SaveFileDialog();
            saveFile.Filter = "Текстовый(*.txt)|*.txt";
            saveFile.DefaultExt = "txt";
            saveFile.Title = "Сохранение решения";

            if (saveFile.ShowDialog() == DialogResult.OK)
            {
                FileStream fs = new FileStream(saveFile.FileName, FileMode.Create);
                StreamWriter st = new StreamWriter(fs);
                foreach (string str in strList)
                {
                        st.WriteLine(str);
                }
                st.Close();
            }
            SolveBox.Items.Add("Решение сохранено в " + saveFile.FileName);
            strList.Clear();
        }       

        //Нахождение транспонированной матрицы
        public double[,] Tran(double[,] Mas)
        {
            double[,] C = Mas;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < i; j++)
                    if (i != j)
                    {
                        double c = C[i, j];
                        C[i, j] = C[j, i];
                        C[j, i] = c;
                    }
            return C;
        }

        //Нахождение обратной матрицы
        public double[,] Reverse(double[,] Mas)
        {
            double o = Determ(Mas);
            double[,] C = new double[n, n];

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    C[i, j] = Mas[i, j];

            C = Tran(AlgAdd(C));
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    C[i, j] = C[i, j] / o;
            return C;
        }

        //Нахождение матрицы алгебраических дополнений
        public double[,] AlgAdd(double[,] Mas)
        {
            double[,] Alg = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    Alg[i, j] = (i + j) % 2 == 0 ? Determ(GetMinor(Mas, i, j)) : (-1) * Determ(GetMinor(Mas, i, j));
            return Alg;
        }

        //Определитель квадратной матрицы
        public double Determ(double[,] matrix)
        {
            if (matrix.GetLength(0) != matrix.GetLength(1)) throw new Exception(" Число строк в матрице не совпадает с числом столбцов");
            double det = 0;
            int Rank = matrix.GetLength(0);
            if (Rank == 1) det = matrix[0, 0];
            if (Rank == 2) det = matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];
            if (Rank > 2)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    det += Math.Pow(-1, 0 + j) * matrix[0, j] * Determ(GetMinor(matrix, 0, j));
                }
            }
            return det;
        }

        public double [] Raund(double[] vec)
        {
            for (int i = 0; i < n; i++)
                vec[i] = Math.Round(vec[i], 4);
            return vec;
        }
        //Найти минор
        public double[,] GetMinor(double[,] matrix, int row, int column)
        {
            if (matrix.GetLength(0) != matrix.GetLength(1)) throw new Exception(" Число строк в матрице не совпадает с числом столбцов");
            double[,] buf = new double[matrix.GetLength(0) - 1, matrix.GetLength(0) - 1];
            for (int i = 0; i < matrix.GetLength(0); i++)
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    if ((i != row) || (j != column))
                    {
                        if (i > row && j < column) buf[i - 1, j] = matrix[i, j];
                        if (i < row && j > column) buf[i, j - 1] = matrix[i, j];
                        if (i > row && j > column) buf[i - 1, j - 1] = matrix[i, j];
                        if (i < row && j < column) buf[i, j] = matrix[i, j];
                    }
                }
            return buf;
        }

        //Разность векторов
        public double[] Subtraction(double[] vec1, double[] vec2)
        {
            double[] C = new double[n];
            for (int i = 0; i < n; i++)
                C[i] = vec1[i] - vec2[i];
            return C;
        }

        //Скалярное произведение векторов
        public double MultiScal(double[] A, double[] B)
        {
            double sum = A[0] * B[0];
            for (int i = 1; i < n; i++)
                sum += A[i] * B[i];
            return sum;
        }

        //Для умножения
        public double Sum(int i, int j, double[,] A, double[,] B)
        {
            double sum = 0;
            for (int k = 0; k < n; k++)
                sum += A[i, k] * B[k, j];
            return sum;
        }

        public double Sum(int i, double[,] A, double[] B)
        {
            double sum = 0;
            for (int k = 0; k < n; k++)
                sum += A[i, k] * B[k];
            return sum;
        }

        public double Sum(int i, double[] A, double[,] B)
        {
            double sum = 0;
            for (int k = 0; k < n; k++)
                sum += A[k] * B[k, i];
            return sum;
        }
        //умножение вектора на число
        public double[] Multiplication(double A, double[] B)
        {
            if (B != null)
            {
                double[] C = new double[n];
                for (int i = 0; i < n; i++)                    
                        C[i] = B[i]*A;
                return C;
            }
            return null;
        }

        //Умножение матрицы на вектор
        public double[] Multiplication(double[,] A, double[] B)
        {
            if (A != null && B != null)
            {
                double[] C = new double[n];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        C[i] = Sum(i, A, B);
                return C;
            }
            return null;
        }

        //Умножение матриц
        public double[,] Multiplication(double[,] A, double[,] B)
        {
            if (A != null && B != null)
            {
                double[,] C = new double[n, n];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        C[i, j] = Sum(i, j, A, B);
                return C;
            }
            return null;
        }

        //Нахождение нормы матрицы (max)
        public double Norma(double[,] mas)
        {
            double max = 0, temp = 0;
            for (int i = 0; i < n; i++)
            {
                temp = 0;
                for (int j = 0; j < n; j++)
                    temp += Math.Abs(mas[i, j]);
                if (temp > max)
                    max = temp;
            }
            return max;
        }

        //Нахождение нормы вектора (max)
        public double Norma(double[] vec)
        {
            double max = 0;
            for (int i = 0; i < n; i++)
                if (max < Math.Abs(vec[i]))
                    max = Math.Abs(vec[i]);
            return max;
        }

        //Для функции
        public double f_13(double x)
        {
            return Math.Pow(3, x) + 2 * x - 5;
        }

        public double f_derivative_13(double x)
        {
            return (Math.Pow(3, x) * Math.Log(3) + 2);
        }

       
        public double f_22(double x)
        {            
            return Math.Atan(x) + 2 * x - 1;
        }

        public double f_derivative_22(double x)
        {
            return (1 / (x * x + 1) + 2);
        }

       
        public double[] f_13_reverse(double[] x_)
        {
            var temp = new double[n];
            temp[0] = 1 - Math.Sin(x_[1]) / 2;
            temp[1] = 0.7 - Math.Cos(x_[0] - 1);
            return temp;
        }

        public double[] f_22_reverse(double[] x_)
        {
            var temp = new double[n];
            //temp[0] = Math.Acos(0.8 - x_[1]) + 1;
            //temp[1] = Math.Acos(x_[0] - 2);
            temp[0] = 2 + Math.Cos(x_[1]);
            temp[1] = 0.8 - Math.Cos(x_[0] - 1);
            return temp;
        }    
        public void create_table(int type)
        {
            int a = 1, b = 2;
            double h = (double)(b - a) / n;
            for(int i=0;i<=n;i++)
            {
                A[i, 0] = a + i * h;
            }          
            if (type == 13)
            {
                for (int i = 0; i <= n; i++)
                {
                    A[i, 1] = f_13(A[i, 0]);
                }
            }
            else
            {
                for (int i = 0; i <= n; i++)
                {
                    A[i, 1] = f_22(A[i, 0]);
                }
            }
            for (int i = 2; i <= n + 1; i++)//отвечает за столбцы
            { 
                //проход по элементам в столбце
                for (int j = 0; j <= n - i + 1; j ++)
                {
                    A[j, i] = (A[j + 1, i - 1] - A[j, i - 1]) / (A[j+i - 1, 0] - A[j , 0]);
                }

            }
            
            WriteMas(A);


        }
        
        //Метод Ньютона
        public void Niuton()
        {          
           
            strList.Add("");
            strList.Add("Метод Ньютона:");           
            if (rb_13.Checked)
            {
               
            }
            else
            {
                
            }              
        }

        
        //Решение уравнения
        private void SolveButton_Click(object sender, EventArgs e)
        {
            //Niuton();   
                     
            create_table(13);
            SaveFile();
        }

       
    }
}
