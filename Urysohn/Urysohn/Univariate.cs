using System;
using System.Collections.Generic;
using System.Text;

namespace Urysohn
{
    internal class Spline
    {
        public double a { get; set; }
        public double b { get; set; }
        public double c { get; set; }
        public double d { get; set; }

        public Spline(double A, double B, double C, double D)
        {
            this.a = A; this.b = B; this.c = C; this.d = D; 
        }
    }

    internal class Basis
    {
        public List<Spline> splines = new List<Spline>();

        public void AddSpline(double A, double B, double C, double D)
        {
            splines.Add(new Spline(A, B, C, D));
        }

        public double GetValue(int spline, double relativeDistance)
        {
            return splines[spline].a + splines[spline].b * relativeDistance +
                splines[spline].c * relativeDistance * relativeDistance +
                splines[spline].d * relativeDistance * relativeDistance * relativeDistance;
        }
    }

    internal class Univariate
    {
        private int _points;
        List<Basis> _basisList = new List<Basis>();
        double[] _coefficients = null;
        double _xmin;
        double _xmax;
        double _ymin;
        double _ymax;
        double _deltax;
        Random _rnd = new Random();

        public Univariate(double xmin, double xmax, double ymin, double ymax, int points)
        {
            _points = points;
            _xmin = xmin;
            _xmax = xmax;
            _xmin -= 0.01 * (_xmax - _xmin);
            _xmax += 0.01 * (_xmax - _xmin);
            _deltax = (_xmax - _xmin) / (_points - 1);
            _ymin = ymin;
            _ymax = ymax;
            Initialize();
        }

        private void PopulateBasisFunctions(SplineGenerator sg, double[][] R, double[] h)
        {
            for (int i = 0; i < _points; ++i)
            {
                double[] e = new double[_points];
                for (int j = 0; j < _points; ++j)
                {
                    e[j] = 0.0;
                }
                e[i] = 1.0;

                (double[] a, double[] b, double[] c, double[] d) = sg.MakeSplines(R, e, h);

                Basis basis = new Basis();
                for (int j = 0; j < a.Length; ++j)
                {
                    basis.AddSpline(a[j], b[j], c[j], d[j]);
                }
                _basisList.Add(basis);
            }
        }

        private void InitializeCoefficients()
        {
            _coefficients = new double[_basisList.Count];
            for (int i = 0; i < _coefficients.Length; ++i)
            {
                _coefficients[i] = _rnd.Next(10, 1000) / 1000.0 * (_ymax - _ymin) + _ymax;
                _coefficients[i] /= _coefficients.Length;
            }
        }

        private void FitDefinition(double x)
        {
            if (x < _xmin)
            {
                x = x - 0.01 * (_xmax - _xmin);
                _deltax = (_xmax - x) / (_points - 1);
                _xmin = x;
            }
            else if (x > _xmax)
            {
                x = x + 0.01 * (_xmax - _xmin);
                _deltax = (x - _xmin) / (_points - 1);
                _xmax = _xmin + (_points - 1) * _deltax;
            }
        }

        public (int k, double relative) GetSplineAndRelative(double x)
        {
            int k = (int)((x - _xmin) / _deltax);
            if (k > _points - 2) k = _points - 2;
            double relative = (x - (_xmin + _deltax * k)) / _deltax;
            return (k, relative);
        }

        public double GetFunctionValue(double x)
        {
            FitDefinition(x);

            (int k, double relative) = GetSplineAndRelative(x);  

            double v = 0.0;
            for (int i = 0; i < _basisList.Count; i++)
            {
                v += _basisList[i].GetValue(k, relative) * _coefficients[i];
            }
            return v;
        }

        public void Update(double x, double delta, double mu)
        {
            FitDefinition(x);

            (int k, double relative) = GetSplineAndRelative(x);

            double[] vectorX = new double[_basisList.Count];
            for (int i = 0; i < vectorX.Length; i++)
            {
                vectorX[i] = _basisList[i].GetValue(k, relative);
            }

            for (int i = 0; i < vectorX.Length; i++)
            {
                _coefficients[i] += delta * mu * vectorX[i];
            }
        }

        private void Initialize()
        { 
            SplineGenerator sg = new SplineGenerator();
            double[] h = new double[_points - 1];
            for (int i = 0; i < h.Length; i++)
            {
                h[i] = 1.0;
            }
            double[][] M = sg.GenerateTriDiagonal(_points, h);
            double[][] R = sg.MatInverseQR(M);
            PopulateBasisFunctions(sg, R, h);
            InitializeCoefficients();
        }
    }
}
