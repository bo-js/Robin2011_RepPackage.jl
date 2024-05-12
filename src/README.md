### Overview

This is the ReadMe for Robin2011_RepPackage.jl, a julia package to replicated Robin (2011), "On the Dynamics of Unemployment and Wage Distributions", Econometrica, Vol. 79, No 5. (September 2011), 1327-1355.

All code for Robin (2011) was originally written in Matlab and Stata and can be retrieved from the website of the Econometrics Society [here](https://www.econometricsociety.org/publications/econometrica/2011/09/01/dynamics-unemployment-and-wage-distributions/supp/9070_data%20and%20programs_0.zip).



### Data Availability 

The data used for this replication is entirely available at the website of the Econometrics Society [here](https://www.econometricsociety.org/publications/econometrica/2011/09/01/dynamics-unemployment-and-wage-distributions/supp/9070_data%20and%20programs_0.zip). The author himself has procured the data from the US Bureau of Labor Statistics.

### Computational Requirements
#### Software Requirements
- Julia (version used for this replication exercise: 1.9)
- The following julia packages need to be installed. The versions of all packages can be found in Manifest.toml
  - CSV 
  - Copulas 
  - DataFrames
  - DelimitedFiles 
  - Distributions 
  - LinearAlgebra
  - Optimization 
  - OptimizationNLopt 
  - Plots 
  - Revise 
  
#### Memory and Run-time Requirements
The code was last run on MacBook Pro M2, 16GB, MacOS 13.06. On this device, approximate run time is about 180 seconds.

#### Description of the Programs

|Script|Content|Run-time|
|---|---|---|
|estmain.jl|estimates the model| 110 seconds |
|plots.jl|runs sim_wages.jl and generates plots| 70 seconds |
|sim_wages.jl|simulates wages| 65 seconds |

#### Instructions to Replicators
To replicated the results generated using Robin2011_RepPackage.jl, change the paths in estmain.jl and sim_wages.jl to where you placed the replication files on your computer. Then run estmain.jl, which generates the estimated model parameters, followed by plots.jl, which runs sim_wages.jl and generates plots from the simulated data. The content of all functions used in the abovementioned script can be found in the documentation section of our GitHub page.

#### Outputs

|Output Name|Reference in Robin (2011)|Description|
|---|---|---|
|fig4a_plot.png| Figure 4 (a), p. 1346| Match productivity for different ability types across states |
|fig4b_plot.png| Figure 4 (b), p. 1346| Distribution of types across ability and employment status|
|fig6_plot.png| Figure 6, p. 1348| Unemployment rate as a function of the aggregate shock|
|fig10_plot.png| Figure 10, p. 1352| Starting and promotion wages for various ability quantiles |

