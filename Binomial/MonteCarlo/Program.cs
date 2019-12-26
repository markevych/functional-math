using System;

namespace MonteCarlo
{
    class Program
    {
        static void Main(string[] args)
        {
            double S = 100.0; double K = 100.0; double r = 0.1; double sigma = 0.25;
            double time = 1.0; int no_sims = 5000;
            Console.WriteLine(" call: black scholes price = {0}", option_price_call_black_scholes(S, K, r, sigma, time));
            Console.WriteLine(" simulated price = {0}", option_price_call_european_simulated(S, K, r, sigma, time, no_sims));
            //Console.WriteLine(" put: black scholes price = {0}", option_price_put_black_scholes(S, K, r, sigma, time));
            //Console.WriteLine(" simulated price = {0}", option_price_put_european_simulated(S, K, r, sigma, time, no_sims));
            Console.WriteLine("----------------------------------------------");

            Console.WriteLine(" call: bs delta = {0}", option_price_delta_call_black_scholes(S, K, r, sigma, time));
            Console.WriteLine(" sim delta = {0}", option_price_delta_call_european_simulated(S, K, r, sigma, time, no_sims));
            //Console.WriteLine(" put: bs delta = {0}", option_price_delta_put_black_scholes(S, K, r, sigma, time)); 
            //Console.WriteLine(" sim delta = {0}", option_price_delta_put_european_simulated(S, K, r, sigma, time, no_sims));
            Console.WriteLine("----------------------------------------------");

            Console.WriteLine("Black Scholes call option price = {0}", option_price_call_black_scholes(S, K, r, sigma, time));
            Console.WriteLine("Simulated call option price = {0}", derivative_price_simulate_european_option_generic(S, K, r, sigma, time, payoff_call, no_sims));
            //Console.WriteLine("Black Scholes put option price = {0}", option_price_put_black_scholes(S, K, r, sigma, time));
            //Console.WriteLine("Simulated put option price = {0}", derivative_price_simulate_european_option_generic(S, K, r, sigma, time, payoff_put, no_sims));

        }


        public static double payoff_call(double S, double K)
        {
            return Math.Max(0.0, S - K);
        }

        public static double payoff_put(double S, double K)
        {
            return Math.Max(0.0, K - S);
        }

        public static double derivative_price_simulate_european_option_generic(double S,
                                                                               double K,
                                                                               double r,
                                                                               double sigma,
                                                                               double time,
                                                                               Func<double, double, double> payoff,
                                                                               int no_sims) 
        {
            double sum_payoffs = 0;
            for (int n = 0; n < no_sims; n++) {
                double S_T = simulate_lognormal_random_variable(S, r, sigma, time);
                sum_payoffs += payoff(S_T, K);
            }
            return Math.Exp(r* time) * (sum_payoffs / no_sims);
        }

        public static double random_normal()
        {
            double U1, U2, V1 = 0, V2;
            double S = 2;
            var random = new Random();
            while (S >= 1)
            {
                U1 = random.NextDouble();
                U2 = random.NextDouble();
                V1 = 2.0 * U1 - 1.0;
                V2 = 2.0 * U2 - 1.0;
                S = Math.Pow(V1, 2) + Math.Pow(V2, 2);
            };
            double X1 = V1 * Math.Sqrt((-2.0 * Math.Log(S)) / S);
            return X1;
        }

        public static double option_price_delta_call_black_scholes(double S, // spot price
                                                                   double K, // Strike (exercise) price,
                                                                   double r, // interest rate
                                                                   double sigma, // volatility
                                                                   double time)
        { // time to maturity
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
            double delta = N(d1);
            return delta;
        }

        public static double simulate_lognormal_random_variable(double S, // current value of variable
                                                                double r, // interest rate
                                                                double sigma, // volatitily
                                                                double time) // time to final date
        {
            double R = (r - 0.5 * Math.Pow(sigma, 2)) * time;
            double SD = sigma * Math.Sqrt(time);
            return S * Math.Exp(R + SD * random_normal());
        }

        public static double option_price_call_european_simulated(double S,
                                                                   double K,
                                                                   double r,
                                                                   double sigma,
                                                                   double time,
                                                                   int no_sims)
        {
            double R = (r - 0.5 * Math.Pow(sigma, 2)) * time;
            double SD = sigma * Math.Sqrt(time);
            double sum_payoffs = 0.0;
            for (int n = 1; n <= no_sims; n++)
            {
                double S_T = S * Math.Exp(R + SD * random_normal());
                sum_payoffs += Math.Max(0.0, S_T - K);
            };
            return Math.Exp(-r * time) * (sum_payoffs / no_sims);
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

        public static double option_price_delta_call_european_simulated(double S,
                                                          double K,
                                                          double r,
                                                          double sigma,
                                                          double time,
                                                          int no_sims)
        {
            double R = (r - 0.5 * Math.Pow(sigma, 2)) * time;
            double SD = sigma * Math.Sqrt(time);
            double sum_payoffs = 0.0;
            double sum_payoffs_q = 0.0;
            double q = S * 0.01;

            for (int n=1; n <= no_sims; n++) {
                double Z = random_normal();
                double S_T = S * Math.Exp(R + SD * Z);
                sum_payoffs += Math.Max(0.0, S_T - K);
                double S_T_q = (S + q) * Math.Exp(R + SD * Z);
                sum_payoffs_q += Math.Max(0.0, S_T_q - K);
             }

            double c = Math.Exp(r * time) * (sum_payoffs / no_sims);
            double c_q = Math.Exp(r* time) * (sum_payoffs_q / no_sims);

            return (c_q - c) / q;
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
    }
}
