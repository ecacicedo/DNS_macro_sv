# Dynamic Nelson-Siegel (DNS) macro sv
Implementation of the DNS macro-sv yield curve model

In this version of the Dynamic Nelson-Siegel (Diebold & Li, 2006) the latent factors and macroeconomic factors follow a VAR(1) specification. The yield curve follows the Nelson-Siegel equation with three latent factors.

The yield residuals and latent factors volatility follow a log-normal model.

The 'model.txt' file contains the script for fitting the DNS macro-sv model. The Diebold & Li (2006) Dynamic Nelson-Siegel yields only (DNS) model is used as reference.

The 'bugsmodel.txt' file contains the script for the WinBUGS algorithm.

Data consists of monthly data for three macroeconomic variables (cpi, inndustrial capacity and base rate) and a term structure of interest rates for 14 maturities.
