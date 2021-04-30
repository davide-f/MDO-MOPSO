# MDO-MOPSO [![View Multiple Design Options - MOPSO (MDO-MOPSO) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://it.mathworks.com/matlabcentral/fileexchange/82174-multiple-design-options-mopso-mdo-mopso)

Improved version of Multiple Objective Particle Swarm Optimization (MOPSO), that:
1. has advanced convergence criteria,
2. tracks the entire history of the simulated points including auxiliary quantities,
3. supports parallel computing by specifying a paremter,
4. identifies Multiple Design Options: points close to the Pareto frontier for a given tolerance,
5. presents a more user-friendly method to set options is also included, and
6. supports vectorized functions.

The two advanced convergence criteria are based on:
1. the spread measure and 
2. quadratic mean of the crouding distances, which enhance the standard criterion based on maximum number of iterations.

Similarly to MATLAB optimization functions, the code implements options for the solver with default values; this version enables the user to change only the desired parameter in a simple way, as also shown in the example.
More details are reported in the referenced paper or please leave a comment below.

The proposed formulation automatically supports vectorized functions or not-vectorized function; for the non-vectorized functions, the methodology also supports parallel computing by simply enabling a parameter.

Cite as:
> Fioriti, D., Lutzemberger, G., Poli, D., Duenas-Martinez, P., & Micangeli, A. (2021). Coupling economic multi-objective optimization and multiple design options: A business-oriented approach to size an off-grid hybrid microgrid. International Journal of Electrical Power & Energy Systems, 127, 106686. Elsevier BV. Retrieved from https://doi.org/10.1016%2Fj.ijepes.2020.106686

bibtex:
> @article{Fioriti_2021, doi = {10.1016/j.ijepes.2020.106686}, url = {https://doi.org/10.1016%2Fj.ijepes.2020.106686}, year = 2021, month = {may}, publisher = {Elsevier {BV}}, volume = {127}, pages = {106686}, author = {Davide Fioriti and Giovanni Lutzemberger and Davide Poli and Pablo Duenas-Martinez and Andrea Micangeli}, title = {Coupling economic multi-objective optimization and multiple design options: A business-oriented approach to size an off-grid hybrid microgrid}, journal = {International Journal of Electrical Power {\&} Energy Systems} }

Code inspired from Multiple Objective Particle Swarm Optimization by Víctor Martínez-Cagigal https://it.mathworks.com/matlabcentral/fileexchange/62074-multi-objective-particle-swarm-optimization-mopso
