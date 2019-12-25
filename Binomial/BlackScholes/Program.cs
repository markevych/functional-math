using System;

namespace BlackScholes
{
    class Program
    {
        static void Main(string[] args)
        {
            double S = 50; double K = 50; double r = 0.10;
            double sigma = 0.30; double time = 0.50, C = 2.5;
            Console.Write(" Black Scholes call price = ");
            Console.WriteLine(option_price_call_black_scholes(S, K, r, sigma, time));
            Console.WriteLine("----------------------------------------------");


            Console.WriteLine(" Black Scholes call partial derivatives ");
            double Delta, Gamma, Theta, Vega, Rho;
            option_price_partials_call_black_scholes(S, K, r, sigma, time, out Delta, out Gamma, out Theta, out Vega, out Rho);
            Console.WriteLine(" Delta = {0}", Delta);
            Console.WriteLine(" Gamma = {0}", Gamma);
            Console.WriteLine(" Theta = {0}", Theta);
            Console.WriteLine(" Vega = {0}", Vega);
            Console.WriteLine(" Rho = {0}", Rho);
            Console.WriteLine("----------------------------------------------");

            Console.Write(" Black Scholes implied volatility using Newton search = ");
            Console.WriteLine(option_price_implied_volatility_call_black_scholes_newton(S, K, r, time, C));
            Console.Write(" Black Scholes implied volatility using bisections = ");
            Console.WriteLine(option_price_implied_volatility_call_black_scholes_bisection(S, K, r, time, C));
            Console.ReadLine();
        }

        public static double n(double z)
        {
            return (1.0 / Math.Sqrt(2.0 * Math.PI)) * Math.Exp(-0.5 * z * z);
        }

        public static double N(double z)
        {
            double b1 = 0.31938153;
            double b2 = -0.356563782;
            double b3 = 1.781477937;
            double b4 = -1.821255978;
            double b5 = 1.330274429;
            double p = 0.2316419;
            double c2 = 0.3989423;

            if (z > 6.0) { return 1.0; }; // this guards against overflow 
            if (z < -6.0) { return 0.0; };
            double a = Math.Abs(z);
            double t = 1.0 / (1.0 + a * p);
            double b = c2 * Math.Exp((-z) * (z / 2.0));
            double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
            n = 1.0 - b * n;
            if (z < 0.0) n = 1.0 - n;
            return n;
        }

        public static double option_price_call_black_scholes(double S, // spot (underlying) price
                                                             double K, // strike (exercise) price,
                                                             double r, // interest rate
                                                             double sigma, // volatility
                                                             double time)
        { // time to maturity
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
            double d2 = d1 - (sigma * time_sqrt);
            return S * N(d1) - K * Math.Exp(-r * time) * N(d2);
        }

        public static void option_price_partials_call_black_scholes(double S, // spot price
                                                                    double K, // Strike (exercise) price,
                                                                    double r, // interest rate
                                                                    double sigma, // volatility
                                                                    double time, // time to maturity
                                                                    out double Delta, // partial wrt S
                                                                    out double Gamma, // second prt wrt S
                                                                    out double Theta, // partial wrt time
                                                                    out double Vega, // partial wrt sigma
                                                                    out double Rho) // partial wrt r
        {
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
            double d2 = d1 - (sigma * time_sqrt);
            Delta = N(d1);
            Gamma = n(d1) / (S * sigma * time_sqrt);
            Theta = -(S * sigma * n(d1)) / (2 * time_sqrt) - r * K * Math.Exp(-r * time) * N(d2);
            Vega = S * time_sqrt * n(d1);
            Rho = K * time * Math.Exp(-r * time) * N(d2);
        }


        public static double option_price_implied_volatility_call_black_scholes_bisection(double S,
                                                                                          double K,
                                                                                          double r,
                                                                                          double time,
                                                                                          double option_price)
        {
            if (option_price < 0.99 * (S - K * Math.Exp(-time * r)))
            { // check for arbitrage violations.
                return 0.0; // Option price is too low if this happens
            }

            // simple binomial search for the implied volatility.
            // relies on the value of the option increasing in volatility
            double ACCURACY = 1.0e-5; // make this smaller for higher accuracy
            int MAX_ITERATIONS = 100;
            double HIGH_VALUE = 1e10;
            double ERROR = -1e40;

            // want to bracket sigma. first find a maximum sigma by finding a sigma
            // with a estimated price higher than the actual price.
            double sigma_low = 1e-5;
            double sigma_high = 0.3;

            double price = option_price_call_black_scholes(S, K, r, sigma_high, time);

            while (price < option_price)
            {
                sigma_high = 2.0 * sigma_high; // keep doubling.

                price = option_price_call_black_scholes(S, K, r, sigma_high, time);

                if (sigma_high > HIGH_VALUE) return ERROR; // panic, something wrong
            }

            for (int i = 0; i < MAX_ITERATIONS; i++)
            {
                double sigma = (sigma_low + sigma_high) * 0.5;

                price = option_price_call_black_scholes(S, K, r, sigma, time);

                double test = (price - option_price);

                if (Math.Abs(test) < ACCURACY)
                {
                    return sigma;
                };
                if (test < 0.0)
                {
                    sigma_low = sigma;
                }
                else
                {
                    sigma_high = sigma;
                }
            }

            return ERROR;
        }

        public static double option_price_implied_volatility_call_black_scholes_newton(double S,
                                                                                       double K,
                                                                                       double r,
                                                                                       double time,
                                                                                       double option_price)
        {
            if (option_price < 0.99 * (S - K * Math.Exp(-time * r))) // check for arbitrage violations. Option price is too low if this happens
            {
                return 0.0;
            }

            int MAX_ITERATIONS = 100;
            double ACCURACY = 1.0e-5;
            double t_sqrt = Math.Sqrt(time);

            double sigma = (option_price / S) / (0.398 * t_sqrt); // find initia

            for (int i = 0; i < MAX_ITERATIONS; i++)
            {
                double price = option_price_call_black_scholes(S, K, r, sigma, time);
                double diff = option_price - price;

                if (Math.Abs(diff) < ACCURACY)
                {
                    return sigma;
                }

                double d1 = (Math.Log(S / K) + r * time) / (sigma * t_sqrt) + 0.5 * sigma * t_sqrt;
                double vega = S * t_sqrt * n(d1);
                sigma = sigma + diff / vega;
            }

            return -99e10; // something screwy happened, should throw 
        }
    }
}
