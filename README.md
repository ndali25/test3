What drives the recent surge in inflation?
The historical decomposition roller coaster
- READ ME FILE -
Drago Bergholt∗ Fabio Canova† Francesco Furlanetto‡
Nicol`o Maffei-Faccioli§ P˚al Ulvedal¶
April 9, 2025
Overview
This document describes the replication files for What drives the recent surge in inflation? The historical decomposition roller coaster. The software used for the estimations is

This code is run via MATLAB version 2022b. The codes are written so that it should be easy for readers to
follow. The replicator should expect to let the code run for about 2 hours. A description
of the data and the codes are given below.
Data availability
Statement about Rights
• We certify that the authors of the manuscript have legitimate access to and permission to use the data used in this manuscript.
• We certify that the authors of the manuscript have documented permission to redistribute/publish the data contained within this replication package.
∗Norges Bank. P.O. Box 1179 Sentrum, 0107 Oslo, Norway. E-mail: drago.bergholt@norges-bank.no.
†BI Norwegian Business School, Nydalsveien 37, 0484 Oslo, Norway. Corresponding author. E-mail:
fabio.canova@bi.no.
‡Norges Bank. P.O. Box 1179 Sentrum, 0107 Oslo, Norway. E-mail: francesco.furlanetto@norgesbank.no.
§Norges Bank. P.O. Box 1179 Sentrum, 0107 Oslo, Norway. E-mail: nicolo.maffei-faccioli@norgesbank.no.
¶Nord University Business School P.O. Box 2501, 7729 Steinkjer, Norway. E-mail:
pal.b.ulvedal@nord.no.
1
Details on each Data Source
All the data used for the estimations are stored in the folder Data. This folder contains
four data files in .xlsx format:
1. US bivariate.xlsx : Contains the series used in the main (bivariate) VAR model for
the US.
• GDPC1 : U.S. Bureau of Economic Analysis, Real Gross Domestic Product,
retrieved from FRED, Federal Reserve Bank of St. Louis: https://fred.
stlouisfed.org/series/GDPC1, May 5, 2023.
• GDPDEF: U.S. Bureau of Economic Analysis, Gross Domestic Product: Implicit Price Deflator, retrieved from FRED, Federal Reserve Bank of St. Louis:
https://fred.stlouisfed.org/series/GDPDEF, May 5, 2023.
2. US large.xlsx : Contains the series used in the larger VAR model for the US.
• GDPC1 : U.S. Bureau of Economic Analysis, Real Gross Domestic Product,
retrieved from FRED, Federal Reserve Bank of St. Louis: https://fred.
stlouisfed.org/series/GDPC1, May 5, 2023.
• GDPDEF: U.S. Bureau of Economic Analysis, Gross Domestic Product: Implicit Price Deflator, retrieved from FRED, Federal Reserve Bank of St. Louis:
https://fred.stlouisfed.org/series/GDPDEF, May 5, 2023.
• GPDIC1 : U.S. Bureau of Economic Analysis, Real Gross Private Domestic
Investment, retrieved from FRED, Federal Reserve Bank of St. Louis: https:
//fred.stlouisfed.org/series/GPDIC1, May 5, 2023.
• AHETPI: U.S. Bureau of Labor Statistics, Average Hourly Earnings of Production and Nonsupervisory Employees, Total Private, retrieved from FRED,
Federal Reserve Bank of St. Louis: https://fred.stlouisfed.org/series/
AHETPI, May 5, 2023.
• FEDFUNDS: Board of Governors of the Federal Reserve System (US), Federal
Funds Effective Rate, retrieved from FRED, Federal Reserve Bank of St. Louis:
https://fred.stlouisfed.org/series/FEDFUNDS, May 5, 2023.
3. EuroArea FRED.xlsx : Contains the series used in the model for the euro area.
• HICP: Eurostat, Harmonized Index of Consumer Prices: All-Items HICP for
Euro Area (19 Countries) [CP0000EZ19M086NEST], retrieved from FRED,
Federal Reserve Bank of St. Louis: https://fred.stlouisfed.org/series/
CP0000EZ19M086NEST, June 14, 2023.
2
• INPR: Organization for Economic Co-operation and Development, Production:
Industry: Total Industry Excluding Construction for Euro Area (19 Countries)
[EA19PRINTO01IXOBSAM], retrieved from FRED, Federal Reserve Bank
of St. Louis: https://fred.stlouisfed.org/series/EA19PRINTO01IXOBSAM,
June 14, 2023.
4. IMF data different countries.xlsx : Contains the series used in the VAR models for
Sweden, Norway, Canada and Australia.
• SE:NGDP R SA XDC : IMF, Gross domestic product (GDP), Constant prices,
Seasonally adjusted (SA), Domestic currency, Sweden, retrieved from IMF National Economic Accounts (NEA) Dataset (Series name: SWE.B1GQ.Q.SA.XDC.Q):
https://data.imf.org/en/Data-Explorer?datasetUrn=IMF.STA:QNEA(7.0.
0), June 10, 2024.
• SE:PCPI IX : IMF, Consumer Price Index (CPI), Sweden, retrieved from IMF
Consumer Price Index (CPI) Dataset (Series name: SWE.CPI. T.IX.Q): https:
//data.imf.org/en/Data-Explorer?datasetUrn=IMF.STA:CPI(3.0.1), June
10, 2024.
• NO:NGDP R SA XDC : IMF, Gross domestic product (GDP), Constant prices,
Seasonally adjusted (SA), Domestic currency, Norway, retrieved from IMF National Economic Accounts (NEA) Dataset (Series name: NOR.B1GQ.Q.SA.XDC.Q):
https://data.imf.org/en/Data-Explorer?datasetUrn=IMF.STA:QNEA(7.0.
0), June 10, 2024.
• NO:PCPI IX : IMF, Consumer Price Index (CPI), Norway, retrieved from IMF
Consumer Price Index (CPI) Dataset (Series name: NOR.CPI. T.IX.Q): https:
//data.imf.org/en/Data-Explorer?datasetUrn=IMF.STA:CPI(3.0.1), June
10, 2024.
• CA:NGDP R SA XDC : IMF, Gross domestic product (GDP), Constant prices,
Seasonally adjusted (SA), Domestic currency, Canada, retrieved from IMF National Economic Accounts (NEA) Dataset (Series name: CAN.B1GQ.Q.SA.XDC.Q):
https://data.imf.org/en/Data-Explorer?datasetUrn=IMF.STA:QNEA(7.0.
0), June 10, 2024.
• CA:PCPI IX : IMF, Consumer Price Index (CPI), Canada, retrieved from
IMF Consumer Price Index (CPI) Dataset (Series name: CAN.CPI. T.IX.Q):
https://data.imf.org/en/Data-Explorer?datasetUrn=IMF.STA:CPI(3.0.
1), June 10, 2024.
• AU:NGDP R SA XDC : IMF, Gross domestic product (GDP), Constant prices,
Seasonally adjusted (SA), Domestic currency, Australia, retrieved from IMF
National Economic Accounts (NEA) Dataset (Series name: AUS.B1GQ.Q.SA.XDC.Q):
3
https://data.imf.org/en/Data-Explorer?datasetUrn=IMF.STA:QNEA(7.0.
0), June 10, 2024.
• AU:PCPI IX : IMF, Consumer Price Index (CPI), Australia, retrieved from
IMF Consumer Price Index (CPI) Dataset (Series name: AUS.CPI. T.IX.Q):
https://data.imf.org/en/Data-Explorer?datasetUrn=IMF.STA:CPI(3.0.
1), June 10, 2024.
Computational requirements
The code requires at least 2.5 GB of free disk space. It takes around 48 hours to run on
a laptop with the following specifications:
• Processor: 11th Gen Intel(R) Core(TM) i5-1135G7 @ 2.40 GHz (1 processor, 4 cores)
• RAM: 16 GB
• Operating system: Windows 11
Software requirements
The following software is required:
1. Matlab (code was run with Matlab Release 2022b)
Controlled randomness
The seed for the random number generator is set at the beginning of each Matlab script.
Description of the codes
The following scripts should be run to replicate the results:
1. macro VAR diffuse.m This code estimates the bivariate VAR model for U.S. output
and inflation (see Section for data series), using a diffuse prior over the sample period
from 1983:Q1 to 2022:Q4. The structure of the code is the following: It loads the
data (line 20), estimates the reduced form model (lines 37 - 52), identifies a demand
and a supply shock, using sign restrictions on impact (lines 54 - 60), computes the
historical decompositions for all draws (line 63), and finds the draws closest to the
point-wise median IRFs (line 66). Finally, it plots Figure 1, Figure 2, Figure 3 (a)
and (b), and Figure 10 (a).
2. macro VAR diffuse 83 19.m This code does the same as macro VAR diffuse.m, except that it uses a sample from 1983:Q1 to 2019:Q4. It plots Figure A-3.
4
3. macro VAR diffuse bq.m This code does the same as macro VAR diffuse.m, except
that it identifies the structural shocks using a Blanchard-Quah decomposition. It
plots Figure A-1 (a).
4. macro VAR diffuse cholesky.m This code does the same as macro VAR diffuse.m,
except that it identifies the structural shocks using a Cholesky decomposition. It
plots Figure A-1 (b)
5. macro VAR diffuse long sample.m This code does the same as macro VAR diffuse.m,
except that it uses a sample from 1949:Q1 to 2022:Q4.It plots Figure 3 (e).
6. macro VAR diffuse no constant.m This code does the same as macro VAR diffuse.m,
except that it de-means the data prior to estimation, and estimates the model without a constant. It plots Figure 9.
7. macro VAR EA diffuse.m This code estimates the bivariate VAR model for euro
area output and inflation (see Section for data series), using a diffuse prior over
the sample period from 2001:M1 to 2023:M3. The code follows the same structure
as macro VAR diffuse.m. It plots Figure F-1 and Figure F-2.
8. macro VAR EA SUR.m This code does the same as macro VAR EA diffuse.m, except that it uses the single-unit-root prior for estimation of the reduced form model.
It plots Figure 7.
9. macro VAR minnesota.m This code does the same as macro VAR diffuse.m, except
that it uses a Minnesota-like prior for estimation of the reduced form model. It plots
Figure 3 (d) and Figure 10 (c).
10. macro VAR NIW.m This code does the same as macro VAR diffuse.m, except that
it uses a Normal-Inverse-Wishart prior for estimation of the reduced form model. It
plots Figure 3 (c) and Figure 10 (b)
11. macro VAR SUR.m This code does the same as macro VAR diffuse.m, except that
it uses a single-unit-root prior for estimation of the reduced form model. It plots
Figure 6, Figure 10 (d) and Figure A-2.
12. macro VAR sur AU.m This code estimates the bivariate VAR model for Australian
output and inflation (see Section for data series), using a single-unit-root prior over
the sample period from 1993:Q1 to 2022:Q4. The structural shocks are identified
using sign restrictions on impact. It plots Figure 8 (e)
13. macro VAR sur CA.m This code estimates the bivariate VAR model for Canadian
output and inflation (see Section for data series), using a single-unit-root prior over
the sample period from 1993:Q1 to 2022:Q4. The structural shocks are identified
using sign restrictions on impact. It plots Figure 8 (c)
5
14. macro VAR sur NO.m This code estimates the bivariate VAR model for Norwegian
output and inflation (see Section for data series), using a single-unit-root prior over
the sample period from 1993:Q1 to 2022:Q4. The structural shocks are identified
using sign restrictions on impact. It plots Figure 8 (b)
15. macro VAR sur SE.m This code estimates the bivariate VAR model for Swedish
output and inflation (see Section for data series), using a single-unit-root prior over
the sample period from 1993:Q1 to 2022:Q4. The structural shocks are identified
using sign restrictions on impact. It plots Figure 8 (a)
16. macro VAR US 5variables diffuse.m This code estimates the larger VAR model for
U.S. output, inflation, investment, real wage and the Fed funds rate (see Section for
data series), using a diffuse prior over the sample period from 1983:Q1 to 2022:Q4.
The structural shocks are identified using sign restrictions on impact (see identifying
restrictions in table E-1 in the appendix of the paper). It plots Figure 3 (f) and Figre
E-1.
17. macro VAR US 5variables SUR.m This code does the same as macro VAR US 5variables diffuse.m,
except that it uses a single-unit-root prior for estimation of the reduced form model.
It plots Figure E-2 and Figure E-3
18. simulation VAR.m This code produces artificial data based on a known data generating process, which has the structure of a bivariate VAR(1) model. Two different
datasets are generated (lines 23 - 60): one where variable 1 is more persistent, and
one where variable 1 is less persistent. It estimates a VAR(1) models for different
sample lengths and with different priors (lines 76 - 341). Finally, it plots Figure 4,
Figure C-2 and Figure C-3.
19. simulation VAR different delta.m This code estimates VAR(1) models on the artificial (simulated) data, using the single-unit-root prior with different values for δ.
Produces Figure 5 (a) and (b).
20. simulation VAR diffuse ar02.m This code does the same as simulation VAR.m, except that it has a lower coefficient on variable 1’s own lag. It plots Figure C-1.
21. simulation VAR dispersion.m This code estimates VAR(1) models on the artificial
(simulated) data, using the single-unit-root prior with different values for Y¯
0 and δ,
and calculates the relative dispersion between the posterior draws. It plots Figure 5
(d).
22. simulation VAR marginal likelihood different delta.m This code estimates VAR(1)
models on the artificial (simulated) data, using the single-unit-root prior with different values for δ, and calculates the log-marginal likelihood. It plots Figure 5 (c).
6
23. simulation VAR min sur forecasting.m This code estimates VAR(1) models recursively on the artificial (simulated) data, using both the Minnesota-like prior and the
single-unit-root prior. It produces forecasts for each iteration, which are evaluated
against the actual data, by computing the root mean squared errors (RMSE). It
produces Table C-1.
Instructions to Replicators
• All relevant scripts that runs estimations are stored in the main folder.
• The functions and supplementary code that are used by these scripts are stored in
the folder Functions
• All data files are in .xlsx format and stored in the folder Data
• The output files are stored in the folder Figures
Figures and programs
Comments in the relevant scripts indicates which figure is produced by each cell. An
overview of relevant code is also provided in Table 1, below.
7
Table 1: List of tables and figures and relevant code
Figure/Table Program Line number
Figure 1 macro VAR diffuse.m 72-126
Figure 2 a macro VAR diffuse.m 128-164
Figure 2 b, left macro VAR diffuse.m 166-195
Figure 2 b, right macro VAR diffuse.m 197-222
Figure 3 a macro VAR diffuse bq.m 105-133
Figure 3 b macro VAR diffuse cholesky.m 99-127
Figure 3 c macro VAR NIW.m 52-80
Figure 3 d macro VAR minnesota.m 55-81
Figure 3 e macro VAR diffuse long sample.m 71-96
Figure 3 f macro VAR US 5variables diffuse.m 186-211
Figure 4 a simulation VAR.m 340-391
Figure 4 b simulation VAR.m 393-445
Figure 5 a simulation VAR different delta.m 165-191
Figure 5 b simulation VAR different delta.m 192-217
Figure 5 c simulation VAR marginal likelihood different delta.m 132-152
Figure 5 d simulation VAR dispersion.m 193-212
Figure 6 a macro VAR sur.m 54-86
Figure 6 b, left macro VAR sur.m 88-113
Figure 6 b, right macro VAR sur.m 115-137
Figure 7 a macro VAR EA SUR.m 54-86
Figure 7 b, left macro VAR EA SUR.m 88-114
Figure 7 b, right macro VAR EA SUR.m 116-138
Figure 8 a 1 macro VAR sur SE.m 106-137
Figure 8 a 2 macro VAR sur NO.m 106-138
Figure 8 a 3 macro VAR sur CA.m 107-138
Figure 8 a 4 macro VAR sur AU.m 107-139
Figure 8 b 1 macro VAR sur SE.m 141-163
Figure 8 b 2 macro VAR sur NO.m 142-164
Figure 8 b 3 macro VAR sur CA.m 142-164
Figure 8 b 4 macro VAR sur AU.m 143-165
Figure 9 a macro VAR diffuse no constant.m 74-110
Figure 9 b, left macro VAR diffuse no constant.m 112-140
Figure 9 b, right macro VAR diffuse no constant.m 142-167
Figure 10 a macro VAR diffuse.m 224-248
Figure 10 b macro VAR NIW.m 82-105
Figure 10 c macro VAR minnesota.m 83-107
Figure 10 d macro VAR sur.m 139-162
Figure A-1 a macro VAR diffuse bq.m 135-187
Figure A-1 b macro VAR diffuse cholesky.m 129-180
Figure A-2 macro VAR sur.m 165-178
Figure A-3 macro VAR diffuse 83 19.m 70-101
Figure C-1 simulation VAR diffuse ar02.m 156-175
Figure C-2 simulation VAR.m 446-529
Figure C-3 simulation VAR.m 531-624
Figure D-1 a simulation VAR D1a.m 250-293
Figure D-1 b simulation VAR D1b.m 247-289
Figure D-1 c simulation VAR D1c.m 252-295
Figure E-1 macro VAR US 5variables diffuse.m 214-243
Figure E-2 macro VAR US 5variables sur.m 135-170
Figure E-3 macro VAR US 5variables sur.m 172-207
Figure F-1 macro VAR EA diffuse.m 160-188
Figure F-2 macro VAR EA diffuse.m 190-212
Table C-1 simulation VAR min sur forecasting.m 78-112
