#include <iostream>
#include <iomanip>
#include <cmath>
#include "black_scholes.cpp"


int main() {
    BlackScholesModel bs;
    

    double S = 100.0;    // Current stock price
    double K = 105.0;    // Strike price
    double r = 0.05;     // Risk-free rate (5%)
    double sigma = 0.20; // Volatility (20%)
    double T = 0.25;     // Time to expiration (3 months)
    
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Black-Scholes Options Pricing Model\n";
    std::cout << "====================================\n";
    std::cout << "Stock Price (S): $" << S << "\n";
    std::cout << "Strike Price (K): $" << K << "\n";
    std::cout << "Risk-free Rate (r): " << r * 100 << "%\n";
    std::cout << "Volatility (σ): " << sigma * 100 << "%\n";
    std::cout << "Time to Expiration (T): " << T << " years\n\n";
    
    
    double callPr = bs.callPrice(S, K, r, sigma, T);
    double putPr = bs.putPrice(S, K, r, sigma, T);
    
    std::cout << "OPTION PRICES:\n";
    std::cout << "Call Price: $" << callPr << "\n";
    std::cout << "Put Price: $" << putPr << "\n\n";
    
    // Verify put-call parity: C - P = S - K*e^(-rT)
    double putCallParity = callPr - putPr;
    double theoreticalParity = S - K * exp(-r * T);
    std::cout << "Put-Call Parity Check:\n";
    std::cout << "C - P = " << putCallParity << "\n";
    std::cout << "S - Ke^(-rT) = " << theoreticalParity << "\n";
    std::cout << "Difference: " << abs(putCallParity - theoreticalParity) << "\n\n";
    
   
    auto callGreeks = bs.calculate_call_Greeks(S, K, r, sigma, T);
    auto putGreeks = bs.calculate_put_Greeks(S, K, r, sigma, T);
    
    std::cout << "GREEKS:\n";
    std::cout << "                Call        Put\n";
    std::cout << "Delta:      " << std::setw(8) << callGreeks.delta 
              << "    " << std::setw(8) << putGreeks.delta << "\n";
    std::cout << "Gamma:      " << std::setw(8) << callGreeks.gamma 
              << "    " << std::setw(8) << putGreeks.gamma << "\n";
    std::cout << "Theta:      " << std::setw(8) << callGreeks.theta 
              << "    " << std::setw(8) << putGreeks.theta << "\n";
    std::cout << "Vega:       " << std::setw(8) << callGreeks.vega 
              << "    " << std::setw(8) << putGreeks.vega << "\n";
    std::cout << "Rho:        " << std::setw(8) << callGreeks.rho 
              << "    " << std::setw(8) << putGreeks.rho << "\n";
    
    return 0;
}
