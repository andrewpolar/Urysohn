using System;
using System.Collections.Generic;
using System.Text;

namespace Urysohn
{
    static class Helper
    {
        public static void ShowMatrix(double[][] M)
        {
            for (int i = 0; i < M.GetLength(0); i++)
            {
                for (int j = 0; j < M[i].Length; j++)
                {
                    Console.Write("{0:0.0000} ", M[i][j]);
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        public static double[][] GetProduct(double [][] Left, double[][] Right)
        {
            int N = Left.GetLength(0);
            double[][] P = new double[N][];
            for (int k = 0; k < N; ++k)
            {
                P[k] = new double[N];
            }

            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < N; ++j)
                {
                    P[i][j] = 0.0;
                    for (int k = 0; k < N; ++k)
                    {
                        P[i][j] += Left[i][k] * Right[k][j];
                    }
                }
            }
            return P;
        }

        public static void ShowBasisValues(Basis basis)
        {
            int N = 10;
            for (int i = 0; i < basis.splines.Count - 1; i++)
            {
                for (int j = 0; j < N; ++j)
                {
                    double dist = j * 1.0 / N;
                    double point = basis.GetValue(i, dist);

                    Console.WriteLine("{0:0.0000} ", point);
                }
            }
        }
    }
}
