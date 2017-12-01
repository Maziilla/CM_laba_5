using System;
using System.Windows.Forms;
using System.IO;
using System.Collections.Generic;
using Rychusoft.NumericalLibraries.Integral;
using MatrixLibrary;
using System.Text;

namespace SLAU
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            SolveButton.Enabled = true;
        }
        public const int a = 1, b = 2;
        public const int  n=5; //Размерности
        public double[,] A = new double[n+1 , n + 2]; //Исходная, единичная       
        public double[] solution,RightPart; //Для решения
        public int[] kol;
        public double[,] MatrForM = new double[n - 1, n - 1], MateForSpline = new double[n + 1, n + 1];
        const double E = 0.0001; //Точность
        public double chis;
        List<string> strList = new List<string>();
        Dictionary<double, double> pointsDict;
        Vector coefficient;
        double[,] matrix;

        //Вывод вектора
        public void Write(double[] vec)
        {
            for (int i = 0; i < n-1; i++)
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
            strList.Add("");
            for (int i = 0; i <= n; i++)
            {
                string str1 = "";
                for (int j = 0; j < n+2; j++)
                    str1 += String.Format("{0,-15} ", Mas[i, j]);
                strList.Add(str1);
            }
            strList.Add("");
        }

        public void WriteMas_(double[,] Mas,int a,int b)
        {
            strList.Add("");
            for (int i = 0; i < a; i++)
            {
                string str1 = "";
                for (int j = 0; j < b ; j++)
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

        public double Sum(int i, double[] A, double[,] B)
        {
            double sum = 0;
            for (int k = 0; k < n; k++)
                sum += A[k] * B[k, i];
            return sum;
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

        //Для функции
        public double f(double x)
        {
            if(rb_13.Checked)
                return Math.Pow(3, x) + 2 * x - 5;            
            return Math.Atan(x) + 2 * x - 1;
        }

        public double f_derivative(double x)
        {
            if (rb_13.Checked)
                return (Math.Pow(3, x) * Math.Log(3) + 2);            
            return (1 / (x * x + 1) + 2);
        }
        //d^4/dx^4 arccotx+2x-1 для вольфрама
        public double f_derivative_4(double x)
        {
            if (rb_13.Checked)
                return (Math.Pow(3, x) * Math.Pow(Math.Log(3), 4));            
            return (24 * x * (x * x - 1) / Math.Pow((x * x + 1), 4));
        }
        public double f_derivative_5(double x)
        {
            if (rb_13.Checked)
                return (Math.Pow(3, x) * Math.Pow(Math.Log(3), 5));
            return (-24 * (5 * x * x * x * x - 10 * x * x + 1) / Math.Pow((x * x + 1), 5));
        }
        public double f_derivative_6(double x)
        {
            if (rb_13.Checked)
                return (Math.Pow(3, x) * Math.Pow(Math.Log(3),6));
            return (240 * x * (3 * x * x * x * x - 10 * x * x + 3) / Math.Pow((x * x + 1), 6));
        } 
        public long Factorial(int n)
        {
            long temp = 1;
            for (int i = 1; i <= n; i++)
                temp *= i;
            return temp;
        }        
        public double[] ProgonkaSLAU(double [,] Matr,double[] Right)
        {
            double[] res = new double[n-1], a = new double[n - 1], B = new double[n - 1];
            int N1 = n - 2;
            double y;
            y = Matr[0,0];
            a[0] = -Matr[0,1] / y;
            B[0] = Right[0] / y;
            for (int i = 1; i < N1; i++)
            {
                y = Matr[i,i] + Matr[i,i -1] * a[i - 1];
                a[i] = -Matr[i,i + 1] / y;
                B[i] = (Right[i] - Matr[i,i-1] * B[i - 1]) / y;
            }
            B[N1] = (Right[N1] - Matr[N1, N1 - 1] * B[N1 - 1]) / (Matr[N1, N1] + Matr[N1, N1 - 1] * a[N1 - 1]);
            res[N1] = (Right[N1] - Matr[N1,N1 - 1] * B[N1 - 1]) / (Matr[N1,N1] + Matr[N1,N1 - 1] * a[N1 - 1]);
            for (int i = N1 - 1; i >= 0; i--)
            {
                res[i] = a[i] * res[i + 1] + B[i];
            }
            return res;
        }

        //Метод Ньютона
        public void Niuton()
        {
            strList.Add("Метод Ньютона");

            double h = (double)(b - a) / n;
            for (int i = 0; i <= n; i++)
            {
                A[i, 0] = a + i * h;
                A[i, 1] = f(A[i, 0]);
            }           
            for (int i = 2; i <= n + 1; i++)//отвечает за столбцы
            {
                //проход по элементам в столбце
                for (int j = 0; j <= n - i + 1; j++)
                {
                    A[j, i] = (A[j + 1, i - 1] - A[j, i - 1]) / (A[j + i - 1, 0] - A[j, 0]);
                }
            }
            strList.Add("Таблица разностей");
            WriteMas(A);
            double M6, temp;            
            M6 = f_derivative_6(a);
            for (double i = a+h; i <= b; i += h)
            {
                temp = f_derivative_6(i);
                if (temp > M6)
                    M6 = temp;
            }
            strList.Add(" M6 =  " + M6);            
            //Оценки
            double[,] A1 = new double[5,5];
            h = (double)(b - a) / n;
            for (int i = 0; i < n; i++)
            {                              
                A1[i, 0] = a + (i + 0.5) * h;                
                double w_temp = 1;
                //P
                A1[i, 2] = A[0, 1]; 
                for (int j=1;j<=n;j++)
                {
                    w_temp = 1;
                    for (int p = 0; p < j; p++)
                    {
                        w_temp *= (A1[i, 0] - A[p, 0]);
                    }
                    A1[i, 2] += A[0, j + 1] * w_temp;
                }
                //
                double w = 1;
                for (int j = 0; j <= n; j++)
                {
                    w *= (A1[i, 0] - A[j, 0]);
                }
                w = Math.Abs(w);                
                A1[i, 1] = f(A1[i, 0]); 
                A1[i, 4] = (M6 * w) / Factorial(n + 1);
                A1[i, 3] = Math.Abs(A1[i, 1] - A1[i, 2]);
            }
            strList.Add("Оценка погрешности");
            strList.Add("       х       |      f(x)      |     Pn(x)     |   Реальная   |       Оценка");
            WriteMas_(A1,n,n);

        }

        //Метод кубических сплайнов 1го дефекта
        public void CubSplain()
        {
            strList.Add("Метод кубических сплайнов 1го дефекта");
            double[,] splains = new double[n, n];
            double h = (double)(b - a) / n;
            double nu = h / (2 * h), lambda = nu;
            for (int i = 0; i <= n; i++)
                MateForSpline[i, 0] = a + (i * h);
            for (int i = 0; i <= n; i++)               
                MateForSpline[i, 1] = f(MateForSpline[i, 0]);                
            for (int i = 0; i < n - 1; i++)
            {
                MatrForM[i, i] = 4;
                if (i - 1 >= 0) MatrForM[i, i-1] = 1;
                if (i + 1 < n - 1) MatrForM[i, i+1] = 1;
            }
            strList.Add("Трёхдиагональная матрица");
            WriteMas_(MatrForM, 4, 4);
            RightPart = new double[n - 1];
            for (int i = 1; i < n; i++)
            {
                RightPart[i - 1] = 3 * (MateForSpline[i + 1, 1] - MateForSpline[i - 1, 1]) / h;
            }            
            RightPart[0] -= f_derivative(MateForSpline[0, 0]);
            RightPart[3] -= f_derivative(MateForSpline[n, 0]);            
            strList.Add("Правая часть:");
            Write(RightPart);
            var answer = ProgonkaSLAU(MatrForM, RightPart);
            //Write(answer);
            for (int i = 1; i < n; i++)
                MateForSpline[i,2] = answer[i - 1];
            double M5,M4;
            double temp;           
            MateForSpline[0, 2] = f_derivative(MateForSpline[0, 0]);
            MateForSpline[n, 2] = f_derivative(MateForSpline[n, 0]);
            M5 = Math.Abs(f_derivative_5(a));
            M4 = Math.Abs(f_derivative_4(a));
            for (double i = a + h; i <= b; i += h)
            {
                temp = Math.Abs(f_derivative_5(i));
                if (temp > M5)
                    M5 = temp;
                temp = Math.Abs(f_derivative_4(i));
                if (temp > M4)
                    M4 = temp;
            }
            for (int i = 0; i <= n; i++)
                    MateForSpline[i, 1] = f_derivative(MateForSpline[i, 0]);
            strList.Add(" M5 =  " + M5);
            strList.Add(" M4 =  " + M4);
            temp = M5 / 60 * h * h * h * h;
            for(int i=0;i<=n;i++)
            {
                MateForSpline[i, 3] = Math.Abs(MateForSpline[i, 1] - MateForSpline[i, 2]);
                MateForSpline[i, 4] = temp;
            }
            strList.Add("Оценки погрешностей в точках интерполяции");
            strList.Add("       х       |      f(x)      |      m[i]     |   Реальная   |       Оценка");
            WriteMas_(MateForSpline, n+1, n);

            for (int i = 0; i < n; i++)
            {
                splains[i, 0] = (a + 0.1) + h * i;
                splains[i, 4] = (M4 / 384 + M5 * h / 240) * h * h * h * h;
            }
            Func<double, double> z0 = x => (1 + 2 * x) * (1 - x) * (1 - x);
            Func<double, double> z1 = x => x * (1 - x) * (1 - x);           
            for (int i = 0; i < n; i++)
            {
                splains[i, 1] = f(splains[i, 0]);
                double t = (splains[i, 0] - MateForSpline[i, 0]) / h;
                splains[i, 2] = f(MateForSpline[i, 0]) * z0(t) + f(MateForSpline[i + 1, 0]) * z0(1 - t) + h * (z1(t) * MateForSpline[i, 3] + z1(1 - t) * MateForSpline[i + 1, 3]);
                splains[i, 3] = Math.Abs(splains[i, 1] - splains[i, 2]);
            }
            strList.Add("Оценки погрешностей");
            strList.Add("       х       |      f(x)      |    S31(f;x)   |   Реальная   |       Оценка");
            WriteMas_(splains, n, n);
        }

        //Метод обратной интерполяции
        public void obrNiuton()

        {
            strList.Add("Обратный метод Ньютона");

            double h = (double)(b - a) / n;
            for (int i = 0; i <= n; i++)
            { 
                A[i, 0] = a + i * h;  
                A[i, 1] = f(A[i, 0]);
                double temp_revers = A[i, 1];
                A[i, 1] = A[i, 0];
                A[i, 0] = temp_revers;
                
            }   
            for (int i = 2; i <= n + 1; i++)//отвечает за столбцы
            {
                //проход по элементам в столбце
                for (int j = 0; j <= n - i + 1; j++)
                {
                    A[j, i] = (A[j + 1, i - 1] - A[j, i - 1]) / (A[j + i - 1, 0] - A[j, 0]);
                }
            }
            strList.Add("Таблица разностей");
            WriteMas(A);         
            //P
            double w_temp, answer = A[0, 1], c;
            if (rb_13.Checked)
                c = 5.1;
            else
                c = 3.1;
            for (int j = 1; j <= n; j++)
            {
                w_temp = 1;
                for (int p = 0; p < j; p++)
                {
                    w_temp *= (c - A[p, 0]);
                }
            answer += A[0, j + 1] * w_temp;
            }              
            strList.Add(" х  = "+ answer);
            w_temp = Math.Abs(c - f(answer));
            strList.Add(" Невязка равна " + w_temp);            
        }

        //Строим таблицу
        private Dictionary<double,double> GetPoints(double a, double b, int n, double d = 0)
        {
            var Dict = new Dictionary<double, double>();
            var h = (b - a) / n;
            for (var i = 0; i <= n; i++)
            {
                var x = a + (i + d) * h;
                if (x > b) break;
                Dict.Add(x, f(x));
            }
            strList.Add($"Отрезок аппроксимации [{a}, {b}], h = {h}\r\n");
            strList.Add($"Таблица: \r\n");
            foreach (var point in Dict)
                strList.Add($"{point.Key,5}|{point.Value,10}\r\n");
            return Dict;
        }

        //Таблица скалярных произведений-Для скалярного случая
        public double[,] GetDiscreteMatrix()
        {
            var functionList = new List<Func<double, double>>
            {
                x => 1,
                x => x,
                x => x * x,
                f
            }; //Список функций
            var matrix = new double[functionList.Count, functionList.Count];
            for (var i = 0; i < functionList.Count; i++)
                for (var j = 0; j < functionList.Count; j++)
                    foreach (var point in pointsDict)
                        matrix[i, j] += functionList[i](point.Key) * functionList[j](point.Key);

            return matrix;
        }

        //Функция в строке
        private string functionString()
        {
            if (rb_13.Checked)
                return "3^x+2*x-5";
            else
                return "actg(x)+2*x-1";
        }

        //Таблица скалярных произведений-Для дискретного случая
        public double[,] GetContinuousMatrix()
        {
            var functionList = new List<Func<double, double>>
            {
                x => 1,
                x => x,
                x => x * x,
                f
            }; //Список функций
            var matrix = new double[functionList.Count, functionList.Count];
            for (var i = 0; i < functionList.Count; i++)
                for (var j = i; j < functionList.Count; j++)
                {
                    var f = i == 0 ? "1" : i == functionList.Count - 1 ? functionString() : $"x^{i}";
                    var g = j == 0 ? "1" : j == functionList.Count - 1 ? functionString() : $"x^{j}";
                    var result = new Integral($"({f})*({g})", a, b).ComputeIntegral();
                    matrix[i, j] = result;
                    matrix[j, i] = result;
                }
            return matrix;
        }

        //Коэф
        private Vector GetCoefficient()
        {
            var n = matrix.GetLength(0) - 1;
            var a = new double[n][];
            for (var i = 0; i < n; i++)
                a[i] = new double[n];
            for (var i = 0; i < n; i++)
                for (var j = 0; j < n; j++)
                    a[i][j] = matrix[i, j];
            var b = new double[n];
            for (var i = 0; i < n; i++)
                b[i] = matrix[i, n];

            var sys = new ConjugateGradient(new Matrix(a), new Matrix(b));
            return sys.Iteration();
        }

        public double Solve(double x)
        {
            var functionList = new List<Func<double, double>>
            {
                y => 1,
                y => y,
                y => y * y,
                f
            }; //Список функций
            var g = 0.0;
            for (var i = 0; i < coefficient.Length; i++)
                g += coefficient[i] * functionList[i](x);

            return g;
        }

        public double DiscreteValueError()
        {
            var functionList = new List<Func<double, double>>
            {
                x => 1,
                x => x,
                x => x * x,
                f
            }; //Список функций
            var g = 0.0;
            var ff = 0.0;
            foreach (var point in pointsDict)
            {
                var x = point.Key;
                g += Solve(x) * Solve(x);
                ff += f(x) * f(x);
            }
            return ff - g;
        }

        public void DiscreteErrorOut()
        {
            strList.Add($"Оценка погрешности: {DiscreteValueError()}\r\n");
        }

        public void ContinuousErrorOut()
        {
            strList.Add($"Оценка погрешности: {ContinuousValueError()}\r\n");
        }

        public double ContinuousValueError()
        {
            var f = new Integral($"({functionString()})*({functionString()})", a, b).ComputeIntegral();
            var gString = $"{coefficient[0]}";
            for (var i = 1; i < coefficient.Length; i++)
                gString += $"{(coefficient[i] >= 0 ? $"+{coefficient[i]}" : $"{coefficient[i]}")} * x^{i}";
            var g = new Integral($"({gString})*({gString})", a, b).ComputeIntegral();
            return f - g;
        }

        public void CoefficientOut()
        {
            strList.Add($"Таблица коэффициентов: \r\n");
            var builder = new StringBuilder();
            for (var i = 0; i < coefficient.Length; i++)
                builder.Append($"C{i}|{coefficient[i],11:N8}\r\n");
            strList.Add(builder.ToString());
        }

        private void MatrixOut()
        {
            var lenght = matrix.GetLength(0);
            strList.Add($"Таблица скалярных произведений: \r\n");
            var builder = new StringBuilder();
            builder.Append("  ");
            for (var i = 0; i < lenght - 1; i++)
                builder.Append($"|{$"g{i}",11:N8}");
            builder.Append($"|{$"f",11:N8}\r\n");

            for (var i = 0; i < lenght; i++)
            {
                builder.Append(i == lenght - 1 ? " f" : $"g{i}");
                for (var j = 0; j < lenght; j++)
                    builder.Append($"|{matrix[i, j],11:N8}");
                builder.Append("\r\n");
            }
            strList.Add(builder.ToString());
        }

        //Метод среднеквадратичных приближений
        public void AverageSquare()
        {
            pointsDict= GetPoints(a, b, n);
            strList.Add("Дискретный случай:\r\n");
            matrix = GetDiscreteMatrix();
            coefficient = GetCoefficient();
            MatrixOut();
            CoefficientOut();
            DiscreteErrorOut();
            strList.Add("Непрерывный случай:\r\n");
            matrix = GetContinuousMatrix();
            coefficient = GetCoefficient();
            MatrixOut();
            CoefficientOut();
            ContinuousErrorOut();
        }

        //Решение уравнения
        private void SolveButton_Click(object sender, EventArgs e)
        {
            Niuton();
            CubSplain();
            AverageSquare();
            obrNiuton();
            SaveFile();
        }
    }
}
