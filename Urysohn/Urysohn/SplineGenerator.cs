using System;
using System.Collections.Generic;
using System.Text;

namespace Urysohn
{
    internal class SplineGenerator
    {
        //Code from https://visualstudiomagazine.com/articles/2024/01/03/matrix-inverse.aspx
        public double[][] MatInverseQR(double[][] M)
        {
            // inverse using QR decomp (Householder algorithm)
            double[][] Q;
            double[][] R;
            MatDecomposeQR(M, out Q, out R);

            // TODO: check if determinant is zero (no inverse)

            double[][] Rinv = MatInverseUpperTri(R); // std algo
            double[][] Qtrans = MatTranspose(Q);  // is inv(Q)
            return MatProduct(Rinv, Qtrans);

            // ----------------------------------------------------
            // 10 helpers: MatDecomposeQR, MatMake, MatTranspose,
            // MatInverseUpperTri, MatIdentity, MatCopy, VecNorm,
            // VecDot, VecToMat, MatProduct
            // ----------------------------------------------------

            static void MatDecomposeQR(double[][] mat,
              out double[][] q, out double[][] r)
            {
                // QR decomposition, Householder algorithm.
                int m = mat.Length;  // assumes mat is nxn
                int n = mat[0].Length;  // check m == n

                double[][] Q = MatIdentity(m);
                double[][] R = MatCopy(mat);
                int end = n - 1;
                // if (m == n) end = n - 1; else end = n;

                for (int i = 0; i < end; ++i)
                {
                    double[][] H = MatIdentity(m);
                    double[] a = new double[n - i];
                    int k = 0;
                    for (int ii = i; ii < n; ++ii)
                        a[k++] = R[ii][i];

                    double normA = VecNorm(a);
                    if (a[0] < 0.0) { normA = -normA; }
                    double[] v = new double[a.Length];
                    for (int j = 0; j < v.Length; ++j)
                        v[j] = a[j] / (a[0] + normA);
                    v[0] = 1.0;

                    double[][] h = MatIdentity(a.Length);
                    double vvDot = VecDot(v, v);
                    double[][] alpha = VecToMat(v, v.Length, 1);
                    double[][] beta = VecToMat(v, 1, v.Length);
                    double[][] aMultB = MatProduct(alpha, beta);

                    for (int ii = 0; ii < h.Length; ++ii)
                        for (int jj = 0; jj < h[0].Length; ++jj)
                            h[ii][jj] -= (2.0 / vvDot) * aMultB[ii][jj];

                    // copy h into lower right of H
                    int d = n - h.Length;
                    for (int ii = 0; ii < h.Length; ++ii)
                        for (int jj = 0; jj < h[0].Length; ++jj)
                            H[ii + d][jj + d] = h[ii][jj];

                    Q = MatProduct(Q, H);
                    R = MatProduct(H, R);
                } // i

                q = Q;
                r = R;
            } // QR decomposition 

            static double[][] MatInverseUpperTri(double[][] U)
            {
                int n = U.Length;  // must be square matrix
                double[][] result = MatIdentity(n);

                for (int k = 0; k < n; ++k)
                {
                    for (int j = 0; j < n; ++j)
                    {
                        for (int i = 0; i < k; ++i)
                        {
                            result[j][k] -= result[j][i] * U[i][k];
                        }
                        result[j][k] /= U[k][k];
                    }
                }
                return result;
            }

            static double[][] MatTranspose(double[][] m)
            {
                int nr = m.Length;
                int nc = m[0].Length;
                double[][] result = MatMake(nc, nr);  // note
                for (int i = 0; i < nr; ++i)
                    for (int j = 0; j < nc; ++j)
                        result[j][i] = m[i][j];
                return result;
            }

            static double[][] MatMake(int nRows, int nCols)
            {
                double[][] result = new double[nRows][];
                for (int i = 0; i < nRows; ++i)
                    result[i] = new double[nCols];
                return result;
            }

            static double[][] MatIdentity(int n)
            {
                double[][] result = MatMake(n, n);
                for (int i = 0; i < n; ++i)
                    result[i][i] = 1.0;
                return result;
            }

            static double[][] MatCopy(double[][] m)
            {
                int nRows = m.Length; int nCols = m[0].Length;
                double[][] result = MatMake(nRows, nCols);
                for (int i = 0; i < nRows; ++i)
                    for (int j = 0; j < nCols; ++j)
                        result[i][j] = m[i][j];
                return result;
            }

            static double[][] MatProduct(double[][] matA,
              double[][] matB)
            {
                int aRows = matA.Length;
                int aCols = matA[0].Length;
                int bRows = matB.Length;
                int bCols = matB[0].Length;
                if (aCols != bRows)
                    throw new Exception("Non-conformable matrices");

                double[][] result = MatMake(aRows, bCols);

                for (int i = 0; i < aRows; ++i) // each row of A
                    for (int j = 0; j < bCols; ++j) // each col of B
                        for (int k = 0; k < aCols; ++k)
                            result[i][j] += matA[i][k] * matB[k][j];

                return result;
            }

            static double VecDot(double[] v1, double[] v2)
            {
                double result = 0.0;
                int n = v1.Length;
                for (int i = 0; i < n; ++i)
                    result += v1[i] * v2[i];
                return result;
            }

            static double VecNorm(double[] vec)
            {
                int n = vec.Length;
                double sum = 0.0;
                for (int i = 0; i < n; ++i)
                    sum += vec[i] * vec[i];
                return Math.Sqrt(sum);
            }

            static double[][] VecToMat(double[] vec,
              int nRows, int nCols)
            {
                double[][] result = MatMake(nRows, nCols);
                int k = 0;
                for (int i = 0; i < nRows; ++i)
                    for (int j = 0; j < nCols; ++j)
                        result[i][j] = vec[k++];
                return result;
            }

        } // MatInverseQR

        //https://people.clas.ufl.edu/kees/files/CubicSplines.pdf
        public double[][] GenerateTriDiagonal(int N, double[] h)
        {
            double[][] M = new double[N][];
            for (int i = 0; i < N; ++i)
            {
                M[i] = new double[N];
            }

            M[0][0] = 1.0;
            for (int j = 1; j < N; ++j)
            {
                M[0][j] = 0.0;
            }
            for (int i = 1; i < N - 1; ++i)
            {
                for (int j = 0; j < N; ++j)
                {
                    if (i == j) M[i][j] = 2.0 * (h[i - 1] + h[i]);
                    else if (1 == i - j) M[i][j] = h[i - 1];
                    else if (1 == j - i) M[i][j] = h[j - 1];
                    else M[i][j] = 0.0;
                }
            }
            for (int j = 0; j < N - 1; ++j)
            {
                M[N - 1][j] = 0.0;
            }
            M[N - 1][N - 1] = 1.0;

            return M;
        }

        //https://people.clas.ufl.edu/kees/files/CubicSplines.pdf
        public (double[] a, double[] b, double[] c, double[] d) MakeSplines(double[][] A, double[] y, double[] h)
        {
            int N = y.Length;
 
            double[] z = new double[N];
            z[0] = 0.0;
            for (int i = 1; i < N - 1; ++i)
            {
                z[i] = 3.0 * (y[i + 1] - y[i]) / h[i] - 3.0 * (y[i] - y[i - 1]) / h[i - 1];
            }
            z[N - 1] = 0.0;

            double[] v = new double[N];
            for (int i = 0; i < N; ++i)
            {
                v[i] = 0.0;
                for (int j = 0; j < N; ++j)
                {
                    v[i] += A[i][j] * z[j];
                }
            }

            double[] a = new double[N - 1];
            for (int i = 0; i < a.Length; ++i)
            {
                a[i] = y[i];
            }

            double[] b = new double[N - 1];
            for (int i = 0; i < N - 1; ++i)
            {
                b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * v[i] + v[i + 1]) / 3.0;
            }

            double[] c = new double[N - 1];
            for (int i = 0; i < c.Length; ++i)
            {
                c[i] = v[i];
            }

            double[] d = new double[N - 1];
            for (int i = 0; i < N - 1; ++i)
            {
                d[i] = (v[i + 1] - v[i]) / 3.0 / h[i];
            }

            return (a, b, c, d);
        }

        public void SelfTest()
        {
            double[] x = new double[] { 0.0, 1.0, 2.0, 2.5 };
            double[] y = new double[] { 0.0, 1.0, 8.0, 9.0 };
            double[] h = new double[x.Length - 1];
            for (int i = 0; i < h.Length; ++i)
            {
                h[i] = x[i + 1] - x[i];
            }
            double[][] M = GenerateTriDiagonal(x.Length, h);
            double[][] R = MatInverseQR(M);
            (double[] a, double[] b, double[] c, double[] d) = MakeSplines(R, y, h);
            for (int j = 0; j < a.Length; ++j)
            {
                Console.WriteLine("{0:0.000} {1:0.000} {2:0.000} {3:0.000}", a[j], b[j], c[j], d[j]);
            }
            Console.WriteLine();

            //expected result
            //0.000 - 1.091 -0.000  2.091
            //1.000   5.182  6.273 -4.455
            //8.000   4.364 -7.091  4.727
        }
    }
}
