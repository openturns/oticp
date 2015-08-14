from openturns import *
from integralcompoundpoisson import *
from math import *

# First, build the atom distribution
d = IntegralUserDefined([1, 2, 4, 7], [0.1, 0.2, 0.3, 0.4])
# Second, build Poisson's theta
theta = 20.0

# Build the CompoundPoisson distribution
ICP = OTICP(d, theta)

# Compute the PDF and CDF of the CompoundPoisson distribution
x = 75
print "CDF(", x, ")=", ICP.computeCDF(x)
print "PDF(", x, ")=", ICP.computePDF(x)

# Draw the PDF
