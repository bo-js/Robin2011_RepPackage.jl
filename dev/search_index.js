var documenterSearchIndex = {"docs":
[{"location":"results/#Replicated-Results","page":"Results","title":"Replicated Results","text":"","category":"section"},{"location":"results/","page":"Results","title":"Results","text":"Output Name Reference in Robin (2011) Description\nfig4a_plot.png Figure 4 (a), p. 1346 Match productivity for different ability types across states\nfig4b_plot.png Figure 4 (b), p. 1346 Distribution of types across ability and employment status\nfig6_plot.png Figure 6, p. 1348 Unemployment rate as a function of the aggregate shock\nfig10_plot.png Figure 10, p. 1352 Starting and promotion wages for various ability quantiles","category":"page"},{"location":"results/","page":"Results","title":"Results","text":"(Image: Replication of Figure 4(a), p. 1346)","category":"page"},{"location":"results/","page":"Results","title":"Results","text":"(Image: Replication of Figure 4(b)), p. 1346)","category":"page"},{"location":"results/","page":"Results","title":"Results","text":"(Image: Replication of Figure 6, p. 1348)","category":"page"},{"location":"results/","page":"Results","title":"Results","text":"(Image: Replication of Figure 10, p. 1352)","category":"page"},{"location":"README/#Overview","page":"ReadMe","title":"Overview","text":"","category":"section"},{"location":"README/","page":"ReadMe","title":"ReadMe","text":"This is the ReadMe for Robin2011_RepPackage.jl, a julia package to replicated Robin (2011), \"On the Dynamics of Unemployment and Wage Distributions\", Econometrica, Vol. 79, No 5. (September 2011), 1327-1355.","category":"page"},{"location":"README/","page":"ReadMe","title":"ReadMe","text":"All code for Robin (2011) was originally written in Matlab and Stata and can be retrieved from the website of the Econometrics Society here.","category":"page"},{"location":"README/#Data-Availability","page":"ReadMe","title":"Data Availability","text":"","category":"section"},{"location":"README/","page":"ReadMe","title":"ReadMe","text":"The data used for this replication is entirely available at the website of the Econometrics Society here. The author himself has procured the data from the US Bureau of Labor Statistics.","category":"page"},{"location":"README/#Computational-Requirements","page":"ReadMe","title":"Computational Requirements","text":"","category":"section"},{"location":"README/#Software-Requirements","page":"ReadMe","title":"Software Requirements","text":"","category":"section"},{"location":"README/","page":"ReadMe","title":"ReadMe","text":"Julia (version used for this replication exercise: 1.9)\nThe following julia packages need to be installed. The versions of all packages can be found in Manifest.toml\nCSV \nCopulas \nDataFrames\nDelimitedFiles \nDistributions \nLinearAlgebra\nOptimization \nOptimizationNLopt \nPlots \nRevise ","category":"page"},{"location":"README/#Memory-and-Run-time-Requirements","page":"ReadMe","title":"Memory and Run-time Requirements","text":"","category":"section"},{"location":"README/","page":"ReadMe","title":"ReadMe","text":"The code was last run on MacBook Pro M2, 16GB, MacOS 13.06. On this device, approximate run time is about 180 seconds.","category":"page"},{"location":"README/#Description-of-the-Programs","page":"ReadMe","title":"Description of the Programs","text":"","category":"section"},{"location":"README/","page":"ReadMe","title":"ReadMe","text":"Script Content Run-time\nestmain.jl estimates the model 110 seconds\nplots.jl runs sim_wages.jl and generates plots 70 seconds\nsim_wages.jl simulates wages 65 seconds","category":"page"},{"location":"README/#Instructions-to-Replicators","page":"ReadMe","title":"Instructions to Replicators","text":"","category":"section"},{"location":"README/","page":"ReadMe","title":"ReadMe","text":"To replicated the results generated using Robin2011RepPackage.jl, change the paths in estmain.jl and simwages.jl to where you placed the replication files on your computer. Then run estmain.jl, which generates the estimated model parameters, followed by plots.jl, which runs sim_wages.jl and generates plots from the simulated data. The content of all functions used in the abovementioned script can be found in the documentation section of our GitHub page.","category":"page"},{"location":"README/#Outputs","page":"ReadMe","title":"Outputs","text":"","category":"section"},{"location":"README/","page":"ReadMe","title":"ReadMe","text":"Output Name Reference in Robin (2011) Description\nfig4a_plot.png Figure 4 (a), p. 1346 Match productivity for different ability types across states\nfig4b_plot.png Figure 4 (b), p. 1346 Distribution of types across ability and employment status\nfig6_plot.png Figure 6, p. 1348 Unemployment rate as a function of the aggregate shock\nfig10_plot.png Figure 10, p. 1352 Starting and promotion wages for various ability quantiles","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Robin2011_RepPackage","category":"page"},{"location":"#Robin2011_RepPackage","page":"Home","title":"Robin2011_RepPackage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Robin2011_RepPackage.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Robin2011_RepPackage]","category":"page"},{"location":"#Robin2011_RepPackage.SurplusVFI-Tuple{Matrix, Matrix, Matrix}","page":"Home","title":"Robin2011_RepPackage.SurplusVFI","text":"SurplusVFI(p::Matrix, z::Matrix, Π::Matrix; β::Number = 0.9466)\n\nPerforms Value Function Iteration on the Surplus Function, returns the resulting worker x state Surplus function matrix.\n\nTakes as inputs the matrix of match productivities p, the matrix of home production z, the markov transition matrix Π, and the parameter β as a keyword argument.\n\nβ should be set to equal (1-δ)/(1+r), it's default is the value used in Robin (2011).\n\n\n\n\n\n","category":"method"},{"location":"#Robin2011_RepPackage.WageVFI-Tuple{Matrix, Matrix, Matrix}","page":"Home","title":"Robin2011_RepPackage.WageVFI","text":"WageVFI(S::Matrix, Π::Matrix, z::Matrix; λ1 = 0.119, β::Number = 0.9466)\n\nThis function calculates equilibrium wages, as well as their associated value functions.\n\nWages are assigned through Bertrand Competition, resulting either in workers being paid their reservaiton wage (the monopsony wage), when being hired from unemployment, or capturing the full surplus when being poached.\n\nThe Function returns a Dict with keys:\n\n:wmin, an i x j matrix giving the monopsony wage assigned to worker j, being hired out of unemployment in state i.\n:wmax, an i x j matrix giving the full surplus wage assigned to worker j, being poached in state i.\n:Wmin, an i x j x k array giving the value of worker k, being in state i, with the monopsony wage from state j.\n:Wmax, an i x j x k array giving the value of worker k, being in state i, with the poacher's wage from state j.\n\nThe Function takes as arguments:\n\nS, the surplus matrix, calculated using SurplusVFI.\nΠ, the markov transition matrix, calculated using grids.\nz, the home production matrix, calculated using homeprod.\n\nAs well as parameters:\n\nλ1, the probability of meeting a new employer while already employed.\nβ, the discounted probability of exogenenous separation, i.e. (1 - δ)/(1+r).\n\nThe default parameter values are those used in Robin(2011).\n\n\n\n\n\n","category":"method"},{"location":"#Robin2011_RepPackage.estCrit-Tuple{Any}","page":"Home","title":"Robin2011_RepPackage.estCrit","text":"estCrit(b; M = 500, N = 100, τ = 0.5, α = 0.64, r = 0.05/4, T = 5000, burn = 1000, draw = rand(burn+T, 1), b0 = [1, 0.0226, 0.9136, 0.0579, 0.2141, 2.5296, 0.7842])\n\nestCrit takes a set of transformed parameter values b, simulates an economy of length T+burn based on an evolution of the aggregate state given in the keyword argument draw (which should be a T+burn long series of draws from a Unif(0,1) distribution), drops the first draw periods and  calculates the simulated targeted moments, before finally calculating the distance from the equivalent moments in the data, given by keyword argument b0.\n\nb should be a vector which contains a guess of the parameters to be estimated. They should be ordered and transformed in the following way:\n\nlog(ν)\nlog(μ)\nloginv(δ)\nloginv(λ0)\natanh(ρ)\nlog(σ)\nlog(z0)\nlog(C)\n\nWhere logit = x -> 1/(1+exp(x)) and loginv = x -> log(1/x - 1).\n\nExternal parameters are given by the following keyword arguments, with defaults set to those used in Robin (2011).\n\nτ is the tie-breaking probability for poachers\nα is the elasticity of wages to productivity - this doesn't matter for estimation but must be specified.\nr is the interest rate.\n\nOther keyword arguments include\n\nM the length of the worker ability grid\nN the length of the aggregate state grid\nT the number of time periods to be used in the estimation\nburn the number of additional time periods to be simulated at the beginning of the economy and then dropped to remove reliance on starting values.\nb0 a vector of the moments calculated from the data.\n\nThe targeted moments are:\n\nAverage Productivity\nStandard Deviation of Productivity\nAutocorrelation of Productivity\nAverage Unemployment Rate\nStandard Deviation of the Unemployment Rate\nKurtosis of the Unemployment Rate\nAverage Exit Rate from Unemployment\n\n\n\n\n\n","category":"method"},{"location":"#Robin2011_RepPackage.grids-Tuple{}","page":"Home","title":"Robin2011_RepPackage.grids","text":"grids(; M = 500, N = 100, ρ = 0.913702234286476, σ = 0.0257,ν::Number = 2.019365636076711, μ::Number = 5.786082109731152 )\n\nThis function takes parameter values and constructs the grids and distributions of worker types and aggregate shocks, which it returns in a Dict, the Dict keys are:\n\n:x, the grid of worker productivities\n:y, the grid of aggregate shocks,\n:Π, the Markov transition matrix,\n:l, the pdf of worker types,\n:F, the grid of aggregate shocks before transformed to fit the estimated distribution.\n\n\n\n\n\n","category":"method"},{"location":"#Robin2011_RepPackage.homeprod-Tuple{Vector, Vector}","page":"Home","title":"Robin2011_RepPackage.homeprod","text":"homeprod(x::Vector, y::Vector; B::Number = 1, C::Number = 0.726, α::Number = 0.64, z0::Number = 0.767)`\n\nTakes the worker and aggregate state grids, as well as the parameters B, C, α, and returns the worker x state matrix of home production.\n\nDefault parameter values are those used in Robin (2011).\n\n\n\n\n\n","category":"method"},{"location":"#Robin2011_RepPackage.matchprod-Tuple{Vector, Vector}","page":"Home","title":"Robin2011_RepPackage.matchprod","text":"matchprod(x::Vector, y::Vector; B::Number = 1, C::Number = 0.725767358913686)\n\nTakes the worker and aggregate state grids, as well as the parameters B and C, and returns the worker x state matrix of match productivities.\n\nDefault parameter values are those used in Robin (2011).\n\n\n\n\n\n","category":"method"},{"location":"#Robin2011_RepPackage.params_estimated-Tuple{}","page":"Home","title":"Robin2011_RepPackage.params_estimated","text":"params_estimated(; path = \"x1.txt\")\n\nThis function outputs all necessary parameters for the model, from a set of estimated parameters.\n\nAs it's input it takes a file path, which must point to a .txt file containing the parameter values correctly ordered, line separated, and transformed as in estmain.jl. \n\nThe function then performs the appropriate inverse transformations and returns the parameters as a NamedTuple. \n\nThe default file path points to those estimated by estmain.jl, stored in x1.txt. Alternatively to use those from Robin(2011), replace this with x0.txt. \n\nThe parameters must be stored in the text file in the following order with the following transformations.\n\nlog(ν)\nlog(μ)\nloginv(δ)\nloginv(λ0)\natanh(ρ)\nlog(σ)\nlog(z0)\nlog(C)\n\nWhere logit = x -> 1/(1+exp(x)) and loginv = x -> log(1/x - 1).\n\n\n\n\n\n","category":"method"},{"location":"#Robin2011_RepPackage.unemp_path-Tuple{Matrix, Vector, Integer}","page":"Home","title":"Robin2011_RepPackage.unemp_path","text":"function unemp_path(S::Matrix, statet::Vector, T::Integer; λ0::Number = 0.994544861919718, δ::Number = 0.041563759920623)\n\nThis function calculates the unemployment path for a given evolution of the aggregate state.\n\nThe Function takes as arguments:\n\nS, the surplus matrix, calculated using SurplusVFI\nstatet, a given evolution of the aggregate state index through time\nT, the number of periods\n\nAs well as parameters:\n\nλ0, the rate at which an unemployed worker meets a firm \nδ, the exogeneous job offer separatioin rate\n\nThe default parameter values are those used in Robin(2011).\n\nThe function returns a matrix, where each row gives a time period, and the columns give the unemployment rate for each worker type.\n\n\n\n\n\n","category":"method"},{"location":"#Robin2011_RepPackage.wage_dens_path-Tuple{Matrix, Matrix, Dict, Vector, Matrix, Vector, Integer}","page":"Home","title":"Robin2011_RepPackage.wage_dens_path","text":"function wage_dens_path(S::Matrix, uxt::Matrix, wd::Dict, l::Vector, U::Matrix, statet::Vector, T::Integer; λ0::Number = 0.9945, λ1::Number = 0.1193, δ::Number = 0.0416)\n\nWages are assigned through Bertrand Competition, resulting either in workers being paid their  reservation wage (the monopsony wage), when being hired from unemployment, or capturing the full surplus when being poached.  However, since workers do not receive a poaching offer every period, you may observe multiple wages at every single point in time. Furthermore, changes in the state may make it so that either the worker or the firm have a credible threat to leave/ fire leading to renegotiation of the wage.\n\nThis function calculates the wage density path for a given evolution of the aggregate state, according to the laws of motion derived in Robin (2011), giving the density of workers earning each of the possible wages at each point in time.\n\nThe Function takes as arguments:\n\nS, the surplus matrix, calculated using SurplusVFI\nuxt, the unemployment path, calculated using unemp_path\nwd, a Dict including the minimum and the maximum wage in each period as well as their current value, calculated using WageVFI\nl, the number of workers of each type, calculated using grids\nU, the value of unemployment \nstatet, a given evolution of the aggregate state index through time\nT, the number of periods\n\nAs well as parameters:\n\nλ0, the rate at which an unemployed worker meets a firm \n\nThe default parameter values are those used in Robin(2011).\n\nThe function returns a four dimensional array, where each entry is the measure of workers earning that wage.  The first index denotes the time period, the second the state that wage is assigned, the third the worker type that wage is assigned to, and the fourth wether the wage is the monopsony wage [:,:,:1], or the full surplus wage [:, :, :, 2].\n\n\n\n\n\n","category":"method"}]
}
