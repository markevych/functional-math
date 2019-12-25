using System;
using System.Collections.Generic;

namespace Binomial
{
    class Program
    {
        static void Main(string[] args)
        {
            double spotPrice = 100.0; 
            double exercisePrice = 100.0; 
            double interestRate = 0.025;
            double upMovement = 1.05;
            double downMovement = 1 / upMovement;

            int noPeriods = 2;

            var onePeriodEuropeanCall = OptionPriceCallEuropeanBinomialSinglePeriod(spotPrice, exercisePrice, interestRate, upMovement, downMovement);
            var twoPeriodEuropeanCall = OptionPriceCallEuropeanBinomialMultiPeriodGivenud(spotPrice, exercisePrice, interestRate, upMovement, downMovement, noPeriods);

            Console.WriteLine($"One period European Call: {Math.Round(onePeriodEuropeanCall, 5)}");
            Console.WriteLine($"Two period European Call: {Math.Round(twoPeriodEuropeanCall, 5)}");

            Console.ReadKey();
        }

        static double OptionPriceCallEuropeanBinomialSinglePeriod(double spotPrice, double exercisePrice, double interestRate, double upMovement, double downMovement)
        {
            var putPriceUp = (Math.Exp(interestRate) - downMovement) / (upMovement - downMovement);
            var putPriceDown = 1.0 - putPriceUp;
            var callPriceUp = Math.Max(0.0, (upMovement * spotPrice - exercisePrice));
            var callPriceDown = Math.Max(0.0, (downMovement * spotPrice - exercisePrice));

            return Math.Exp(-interestRate) * (putPriceUp * callPriceUp + putPriceDown * callPriceDown);
        }

        static double OptionPriceCallEuropeanBinomialMultiPeriodGivenud(double spotPrice, double exercisePrice, double interestRate, double upMovement, double downMovement, int noPeriods)
        {
            var inverseOfInterestRates = Math.Exp(-interestRate);
            var upMovementSqaure = upMovement * upMovement;
            var putPriceUp = (Math.Exp(interestRate) - downMovement) / (upMovement - downMovement);
            var putPriceDown = 1.0 - putPriceUp;
            var prices = new List<double>(new double[noPeriods + 1]);

            prices[0] = spotPrice * Math.Pow(downMovement, noPeriods);

            for (int i = 1; i <= noPeriods; i++)
            {
                prices[i] = upMovementSqaure * prices[i - 1];
            }

            var callValues = new List<double>(new double[noPeriods + 1]);

            for (int i = 0; i <= noPeriods; i++)
            {
                callValues[i] = Math.Max(0.0, (prices[i] - exercisePrice));
            }

            for (int step = noPeriods - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; i++)
                {
                    callValues[i] = (putPriceUp * callValues[i + 1] + putPriceDown * callValues[i]) * inverseOfInterestRates;
                }
            }

            return callValues[0];
        }

        List<List<double>> BuildBinomialTree(double spotPrice, double upMovement, double downMovement, int noSteps)
        {
            var tree = new List<List<double>>();

            for (int i = 0; i < noSteps; i++)
            {
                var treeNode = new List<double>(new double[i]);

                for (int j = 0; j < i; j++)
                {
                    treeNode[j] = spotPrice * Math.Pow(upMovement, j) * Math.Pow(downMovement, i - j - 1);
                }

                tree.Add(treeNode);
            }

            return tree;
        }
    }
}
