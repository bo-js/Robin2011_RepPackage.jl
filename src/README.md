Replication exercise for "Computational Economics" taught by [Prof. Florian Oswald](https://floswald.github.io/), SciencesPo Master of Economics, 2023/2024. Authors: [Bo Jacobs-Strom](https://github.com/bo-js) and [Lewin Nolden](https://github.com/lewinno).

The course website is available [here](https://floswald.github.io/NumericalMethods/).
### Overview

This is the ReadMe for Robin2011_RepPackage.jl, a julia package to replicate Robin (2011), "On the Dynamics of Unemployment and Wage Distributions", Econometrica, Vol. 79, No 5. (September 2011), 1327-1355.

All code for Robin (2011) was originally written in Matlab and Stata and can be retrieved from the website of the Econometrics Society [here](https://www.econometricsociety.org/publications/econometrica/2011/09/01/dynamics-unemployment-and-wage-distributions/supp/9070_data%20and%20programs_0.zip).



### Data Availability 

The data used for this replication is entirely available at the website of the Econometrics Society [here](https://www.econometricsociety.org/publications/econometrica/2011/09/01/dynamics-unemployment-and-wage-distributions/supp/9070_data%20and%20programs_0.zip). The author himself has procured the data from the US Bureau of Labor Statistics.

To replicate Robin (2011) using this Julia package, replicators are required to download the data from the original replication package under the link above and place the `USquarterly.raw` datafile in the data folder of this repo.



### Computational Requirements
#### Software Requirements
- Julia (version used for this replication exercise: 1.9)
- The following julia packages need to be installed. The versions of all packages can be found in Manifest.toml.
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

The correct package versions can be installed by running the following code in the Julia REPL while in the root directory of the repo:
```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
  
#### Memory and Run-time Requirements
The code was last run on MacBook Pro M2, 16GB, MacOS 13.06. On this device, approximate run time is about 180 seconds.

#### Description of the Programs

|Script|Content|Run-time|
|---|---|---|
|estmain.jl|estimates the model| 120 seconds |
|plots.jl|runs sim_wages.jl and generates plots| 70 seconds |
|sim_wages.jl|simulates wages| 65 seconds |

#### Instructions to Replicators
To replicate the results generated using Robin2011_RepPackage.jl, change the paths in estmain.jl and sim_wages.jl to where you placed the replication files on your computer. Make sure the paths lead to the data files procured from the original replication package available at the Econometrics Society's website (see above). Then run estmain.jl, which generates the estimated model parameters, followed by plots.jl, which runs sim_wages.jl and generates plots from the simulated data. The content of all functions used in the above mentioned script can be found in the documentation section of our GitHub page [here](https://bo-js.github.io/Robin2011_RepPackage.jl/dev/).

#### Outputs

|Output Name|Reference in Robin (2011)|Description|
|---|---|---|
|fig4a_plot.png| Figure 4 (a), p. 1346| Match productivity for different ability types across states |
|fig4b_plot.png| Figure 4 (b), p. 1346| Distribution of types across ability and employment status|
|fig6_plot.png| Figure 6, p. 1348| Unemployment rate as a function of the aggregate shock|
|fig10_plot.png| Figure 10, p. 1352| Starting and promotion wages for various ability quantiles |


