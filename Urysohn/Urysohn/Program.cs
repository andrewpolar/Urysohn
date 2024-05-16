//This is iterative identification of the discrete Urysohn operator by input/output.
//This method is applicable for identification of Kolmogorov-Arnold representation.
//Developed by Andrew Polar and Mike Poluektov.
//Published:
//https://www.sciencedirect.com/science/article/abs/pii/S0016003220301149
//https://www.sciencedirect.com/science/article/abs/pii/S0952197620303742
//https://arxiv.org/abs/2305.08194

//Usual execution print out
//epoch 0, relative error 0.307858
//epoch 10, relative error 0.026019
//epoch 20, relative error 0.011861
//epoch 30, relative error 0.005931
//epoch 40, relative error 0.003509
//epoch 50, relative error 0.002689
//epoch 60, relative error 0.002464
//epoch 70, relative error 0.002408
//epoch 80, relative error 0.002395
//epoch 90, relative error 0.002391
//Time for identification .17 seconds
//relative error for validation data 0.0024

using System;
using System.Collections.Generic;

namespace Urysohn
{
    internal class Program
    {
        static double Function(double x, double y, double z)
        {
            return Math.Exp(Math.Sin(x)) + Math.Exp(Math.Cos(y)) + Math.Sin(z) / (1.0 + z);
        }

        static (List<double[]> inputs, List<double> target) GenerateData()
        {
            double xmin = 2.0;
            double xmax = 9.0;
            double ymin = 1.0;
            double ymax = 8.0;
            double zmin = 4.0;
            double zmax = 9.0;
            int N = 1000;
            Random rand = new Random();

            List<double[]> inputs = new List<double[]>();
            List<double> target = new List<double>();

            for (int i = 0; i < N; ++i)
            {
                double arg1 = rand.Next(10, 1000) / 1000.0 * (xmax - xmin) + xmin;
                double arg2 = rand.Next(10, 1000) / 1000.0 * (ymax - ymin) + ymin;
                double arg3 = rand.Next(10, 1000) / 1000.0 * (zmax - zmin) + zmin;
                inputs.Add(new double[] { arg1, arg2, arg3 });
                target.Add(Function(arg1, arg2, arg3));
            }

            return (inputs, target);
        }

        static (double[] xmin, double[] xmax, double targetMin, double targetMax) FindMinMax(List<double[]> inputs, List<double> target)
        {
            int size = inputs[0].Length;
            double[] xmin = new double[size];
            double[] xmax = new double[size];

            for (int i = 0; i < size; ++i)
            {
                xmin[i] = double.MaxValue;
                xmax[i] = double.MinValue;
            }

            for (int i = 0; i < inputs.Count; ++i)
            {
                for (int j = 0; j < inputs[i].Length; ++j)
                {
                    if (inputs[i][j] < xmin[j]) xmin[j] = inputs[i][j];
                    if (inputs[i][j] > xmax[j]) xmax[j] = inputs[i][j];
                }
            }

            double targetMin = double.MaxValue;
            double targetMax = double.MinValue;
            for (int j = 0; j < target.Count; ++j)
            {
                if (target[j] < targetMin) targetMin = target[j];
                if (target[j] > targetMax) targetMax = target[j];
            }

            return (xmin, xmax, targetMin, targetMax);
        }

        static void Main(string[] args)
        {
            //Generation data
            (List<double[]> inputs, List<double> target) = GenerateData();
            (double[] xmin, double[] xmax, double targetMin, double targetMax) = FindMinMax(inputs, target);

            DateTime start = DateTime.Now;
            //Identification
            Urysohn urysohn = new Urysohn(xmin, xmax, targetMin, targetMax, new int[] { 10, 10, 10 });
            int epochs = 100;
            for (int epoch = 0; epoch < epochs; ++epoch)
            {
                double error = 0.0;
                int cnt = 0;
                for (int i = 0; i < inputs.Count; ++i)
                {
                    double m = urysohn.GetU(inputs[i]);
                    double delta = target[i] - m;
                    urysohn.Update(delta, inputs[i], 0.01);
                    error += delta * delta;
                    ++cnt;
                }
                error /= cnt;
                error = Math.Sqrt(error);
                error /= (targetMax - targetMin);
                if (0 == epoch % 10)
                {
                    Console.WriteLine("epoch {0}, relative error {1:0.000000}", epoch, error);
                }
            }
            DateTime end = DateTime.Now;
            TimeSpan duration = end - start;
            double time = duration.Minutes * 60.0 + duration.Seconds + duration.Milliseconds / 1000.0;
            Console.WriteLine("Time for identification {0:####.00} seconds", time);

            //Validation
            (List<double[]> inputs_test, List<double> target_test) = GenerateData();

            double error_test = 0.0;
            int cnt_test = 0;
            for (int i = 0; i < inputs_test.Count; ++i)
            {
                double m = urysohn.GetU(inputs_test[i]);
                double delta = target_test[i] - m;
                error_test += delta * delta;
                ++cnt_test;
            }
            error_test /= cnt_test;
            error_test = Math.Sqrt(error_test);
            error_test /= (targetMax - targetMin);
            Console.WriteLine("relative error for validation data {0:0.0000}", error_test);
        }
    }
}
