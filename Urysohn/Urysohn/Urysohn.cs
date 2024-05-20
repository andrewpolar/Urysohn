using System;
using System.Collections.Generic;
using System.Text;

namespace Urysohn
{
    internal class Urysohn
    {
        private List<Univariate> _univariateList = new List<Univariate>();

        public Urysohn(double[] xmin, double[] xmax, double targetMin, double targetMax, int[] points)
        {
            double ymin = targetMin / points.Length;
            double ymax = targetMax / points.Length;
            for (int i = 0; i < points.Length; ++i)
            {
                Univariate univariate = new Univariate(xmin[i], xmax[i], ymin, ymax, points[i]);
                _univariateList.Add(univariate);
            }
        }

        public void Clear()
        {
            _univariateList.Clear();
        }

        public void Update(double delta, double[] inputs, double mu)
        {
            int i = 0;
            foreach (Univariate uni in _univariateList)
            {
                uni.Update(inputs[i++], delta, mu);
            }
        }

        public double GetU(double[] inputs)
        {
            double f = 0.0;
            int i = 0;
            foreach (Univariate uni in _univariateList)
            {
                f += uni.GetFunctionValue(inputs[i++]);
            }
            return f;
        }
    }
}
