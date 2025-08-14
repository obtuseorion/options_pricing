#include <iostream>
#include <cmath>
#include<iomanip>

class BlackScholesModel {
    private:

        double cumulativeNormalDist(double x) const {
            return 0.5*(1.0+ erf(x/sqrt(2.0)));
        }

        double normalPDF(double x) const {
            return (1.0/ sqrt(2.0* M_PI)) * (exp(-0.5*x*x));
        }

        double calculateD1(double S, double K, double r, double sigma, double T) const {
            return (log(S/K) +(r + 0.5 *sigma*sigma)*T)/(sigma *(sqrt(T)));
        }

        double calculateD2(double d1, double sigma, double T) const {
            return d1 - sigma*sqrt(T);
        }
    
    public:

        double callPrice(double S, double K, double r, double sigma, double T) const {
            if (T <= 0) return std::max(S-K, 0.0);
            
            double d1 = calculateD1(S, K, r, sigma, T);
            double d2 = calculateD2(d1, sigma, T);
            
            double N_d1 = cumulativeNormalDist(d1);
            double N_d2 = cumulativeNormalDist(d2);

            return (S* N_d1) - (K * exp(-r*T) * N_d2);
        }

        double putPrice(double S, double K, double r, double sigma, double T) const {
            if (T <= 0) return std::max(K-S, 0.0);

            double d1 = calculateD1(S, K, r, sigma, T);
            double d2 = calculateD2(d1, sigma, T);

            double N_d1 = cumulativeNormalDist(-d1);
            double N_d2 = cumulativeNormalDist(-d2);

            return K * exp(-r*T) * N_d2 - S*N_d1;
        }

        struct Greeks {
            double delta;
            double gamma;
            double theta;
            double vega;
            double rho;
        };

        Greeks calculate_call_Greeks(double S, double K, double r, double sigma, double T) const {
            Greeks greeks;

            if (T <= 0) {

                greeks.delta = (S>K) ? 1.0:0.0;
                greeks.gamma = 0.0;
                greeks.theta = 0.0;
                greeks.vega = 0.0;
                greeks.rho = 0.0;

                return {greeks};
            }

            double d1 = calculateD1(S, K, r, sigma, T);
            double d2 = calculateD2(d1, sigma, T);

            double N_d1 = cumulativeNormalDist(d1);
            double N_d2 = cumulativeNormalDist(d2);
            double n_d1 = normalPDF(d1);

            greeks.delta = N_d1;
            greeks.gamma = n_d1/(S*sigma*sqrt(T));
            greeks.theta =(-(S*n_d1*sigma)/(2*sqrt(T))- r*K*exp(-r*T) * N_d2)/365.0;
            greeks.vega = (S* n_d1 * sqrt(T))/ 100.0;
            greeks.rho = (K*T*exp(-r*T)* N_d2)/ 100.0;

            return {greeks};
        }

        Greeks calculate_put_Greeks(double S, double K, double r, double sigma, double T) const {
            Greeks greeks; 

            if (T <= 0) {

                greeks.delta = (S<K) ? -1.0:0.0;
                greeks.gamma = 0.0;
                greeks.theta = 0.0;
                greeks.vega = 0.0;
                greeks.rho = 0.0;

                return {greeks};
            }

            double d1 = calculateD1(S, K, r, sigma, T);
            double d2 = calculateD2(d1, sigma, T);

            double N_d1 = cumulativeNormalDist(-d1);
            double N_d2 = cumulativeNormalDist(-d2);
            double n_d1 = normalPDF(d1);

            greeks.delta = -N_d1;
            greeks.gamma = n_d1 / (S*sigma*sqrt(T));
            greeks.theta = (-(S*n_d1*sigma)/(2*sqrt(T))+ r*K*exp(-r*T) * N_d2)/365.0;
            greeks.vega = (S* n_d1 * sqrt(T))/ 100.0;
            greeks.rho = (-K*T*exp(-r*T)* N_d2)/ 100.0;

            return {greeks};
        }
};

