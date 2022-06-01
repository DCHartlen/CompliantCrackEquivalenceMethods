# Crack Equivalence Methods for DCB and ENF

A set of methods for evaluating the resistance curve (R-curve) behaviour of double cantilever beam (DCB) and end notch flexure (ENF) test specimens. Unlike traditional methods, these techniques do not require the crack tip to be tracked during testing, simplifying experimental setups and data processing. 

Methods for both DCB and ENF require knowledge of some material properties prior for accurate extraction of resistance behaviour. 

# Methods Used

Two techniques are implimented for DCB analysis. Both are broadly equivalent, but the method of de Gracia et al. is slightly more accurate and robust for large displacements. 
> de Moura, M. F. S. F., Morais, J. J. L., & Dourado, N. (2008). A new data reduction scheme for mode I wood fracture characterization using the double cantilever beam test. <i>Engineering Fracture Mechanics</i>, <i>75</i>(13), 3852–3865. https://doi.org/10.1016/J.ENGFRACMECH.2008.02.006

>de Gracia, J., Boyano, A., Arrese, A., & Mujika, F. (2015). A new approach for determining the R-curve in DCB tests without optical measurements. <i>Engineering Fracture Mechanics</i>, <i>135</i>, 274–285. https://doi.org/10.1016/j.engfracmech.2015.01.016

Only a single method is implimented for ENF analysis. 

>Xavier, J., Oliveira, M., Morais, J. J. L., & de Moura, M. F. S. F. (2014). Determining mode II cohesive law of Pinus pinaster by combining the end-notched flexure test with digital image correlation. <i>Construction and Building Materials</i>, <i>71</i>, 109–115. https://doi.org/10.1016/j.conbuildmat.2014.08.021

# Usage

Methods are based on object oriented programming. Example scripts can be found in the 'Examples' directory. 


