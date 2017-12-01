using System.Collections.Generic;
using System.Text;

namespace MatrixLibrary
{
    public class ConjugateGradient 
    {
        private Matrix A;
        private Vector B;
        private double e;


        public ConjugateGradient(Matrix a, Matrix b)
        {
            A = a;
            B = b.ToVector();
            e = 0.0001;
        }

        private double ValueNorm(Queue<Vector> queue)
        {
            var x = queue.ToArray();
            var q = (x[2] - x[1]).GetCubicNorm() /
                    (x[1] - x[0]).GetCubicNorm();
            return double.IsInfinity(q) ? double.NaN : q;
        }

        private double ValueError(Queue<Vector> queue)
        {
            var x = queue.ToArray();
            var q = ValueNorm(queue);
            return q / (1 - q) * (x[2] - x[1]).GetCubicNorm();
        }

        public Vector Iteration()
        {
            var X = (Vector)B.Clone();
            var queue = new Queue<Vector>();
            queue.Enqueue(X);
            queue.Enqueue(X);
            var a = 1.0;
            var r = A * X - B;
            var T = (r * r) / ((A * r) * r);
            var step = 1;
            queue.Enqueue(X = a * X - T * a * r);
            while (r.GetCubicNorm() > e)
            {
                queue.Dequeue();
                var newR = A * X - B;
                var newT = (newR * newR) / ((A * newR) * newR);
                a = 1 / (1 - (newT * (newR * newR)) / (T * a * (r * r)));
                var array = queue.ToArray();
                queue.Enqueue(X = a * array[1] + (1 - a) * array[0] - newT * a * newR);
                r = newR;
                T = newT;
            }

            return X;
        }
    }
}
