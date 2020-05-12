# RobGARCHBoot

Robust bootstrap forecast densities for GARCH models (Trucíos et al; 2017). 
This R package provides the forecast densities for returns and volatilities which are useful to obtain forecast intervals and risk measures. 
The package also provides the robust GARCH estimator of Boudt et al. (2013) using the modification proposed by Trucíos et al. (2017)

For applications of the bootstrap procedure, see:

- Trucíos, C., Hotta, L. K., & Ruiz, E. (2017). Robust bootstrap forecast densities for GARCH returns and volatilities. Journal of Statistical Computation and Simulation, 87(16), 3152-3174.
- Trucíos, C. (2019). Forecasting Bitcoin risk measures: A robust approach. International Journal of Forecasting, 35(3), 836-847.

For applications using the robust estimator, see:

- Trucíos, C. (2019). Forecasting Bitcoin risk measures: A robust approach. International Journal of Forecasting, 35(3), 836-847.
- Trucíos, C., Hotta, L. K., and Valls, P. (2019). On the robustness of the principal volatility components. Journal of Empirical Finance, 52(1), 201-219.
- Trucíos, C., Tiwari, A. K., & Alqahtani, F. (2020). Value-at-Risk and Expected Shortfall in Cryptocurrencies' Portfolio: A Vine Copula-based Approach. Applied Economics, 52(24), 2580-2593.


## Installation
RobGARCHBoot is available on CRAN, but you can install the latest version using these commands:

install.packages("devtools")


devtools::intall_github("ctruciosm/RobGARCHBoot")
