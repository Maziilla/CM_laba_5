using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatrixLibrary
{
    public class Vector : ICloneable
    {
        private double[] source;
        public int Length => source.Length;

        public Vector(int length)
        {
            source = new double[length];
        }
        
        public Vector(double[] array)
        {
            source = array;
        }
        
        public double this[int i]
        {
            get { return source[i]; }
            set { source[i] = value; }
        }

        public double GetCubicNorm()
        {
            return source.Max(x => Math.Abs(x));
        }

        public double GetNorm()
        {
            return Math.Sqrt(source.Sum(x => x * x));
        }

        public static double operator *(Vector multiplicand, Vector multiplier)
        {
            if (multiplicand.Length != multiplier.Length)
                throw new ArgumentException();

            var result = 0.0;
            for (var i = 0; i < multiplicand.Length; i++)
                result += multiplicand[i] * multiplier[i];
            return result;
        }
        
        public static Vector operator *(Matrix multiplicand, Vector multiplier)
        {
            if (multiplicand.Width != multiplier.Length)
                throw new ArgumentException();

            var vector = new Vector(multiplier.Length);
            for (var i = 0; i < multiplier.Length; i++)
                for (var j = 0; j < multiplier.Length; j++) 
                    vector[i] += multiplicand[i, j] * multiplier[j];
            return vector;
        }
        
        public static Vector operator -(Vector subtrahend, Vector subtractor)
        {
            if (subtrahend.Length != subtractor.Length)
                throw new ArgumentException();

            var vector = new Vector(subtractor.Length);
            for (var i = 0; i < vector.Length; i++)
                vector[i] = subtrahend[i] - subtractor[i];
            return vector;
        }
        
        public static Vector operator +(Vector term1, Vector term2)
        {
            if (term1.Length != term2.Length)
                throw new ArgumentException();

            var vector = new Vector(term2.Length);
            for (var i = 0; i < vector.Length; i++)
                vector[i] = term1[i] + term2[i];
            return vector;
        }
        
        public static Vector operator *(double multiplicand, Vector multiplier)
        {
            var vector = new Vector(multiplier.Length);
            for (var i = 0; i < vector.Length; i++)
                vector[i] = multiplicand * multiplier[i];
            return vector;
        }
        
        public object Clone()
        {
            return new Vector((double[])source.Clone());
        }

        public override string ToString()
        {
            var builder = new StringBuilder();
            builder.Append("(");
            for (var i = 0; i < Length; i++)
                builder.Append($"{this[i]}, ");
            builder.Remove(builder.Length - 2, 2);
            builder.Append(")");
            return builder.ToString();
        }
    }
}
