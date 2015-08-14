from __future__ import print_function
from openturns import *
from oticp import *
from math import *

drawFlag = False
support = [1, 2, 4, 7]
weights = [0.1, 0.2, 0.3, 0.4]

# First, build the atom distribution
d = IntegralUserDefined(support, weights)
# Second, build Poisson's theta
theta = 20.0
# Build the CompoundPoisson distribution
ICP = IntegralCompoundPoisson(d, theta)
ICPSample = ICP.getSample(100000)
userDefined = IntegralUserDefinedFactory().buildAsIntegralUserDefined(ICPSample)
for i in range(1, 181):
    udPDF = userDefined.computePDF(i)
    icpPDF = ICP.computePDF(i)
    udCDF = userDefined.computeCDF(i)
    icpCDF = ICP.computeCDF(i)
    udCDFc = userDefined.computeComplementaryCDF(i)
    icpCDFc = ICP.computeComplementaryCDF(i)
    print("error rel % (PDF, CDF, CDF c)=")
    print(NumericalPoint([udPDF / icpPDF - 1.0, udCDF / icpCDF - 1.0, udCDFc / icpCDFc - 1.0]) * 100)
    print("error abs (PDF, CDF, CDF c)=")
    print(NumericalPoint([fabs(udPDF - icpPDF), fabs(udCDF - icpCDF), fabs(udCDFc - icpCDFc)]))

if drawFlag:
    xMin = -10.0
    xMax = 200.0
    g = ICP.drawPDF(xMin, xMax)
    drawable = userDefined.drawPDF(xMin, xMax, 10 * int(xMax - xMin + 1)).getDrawable(0)
    drawable.setColor("blue")
    drawable.setLineStyle("dotted")
    g.addDrawable(drawable)
    g.draw("ICPPDF", 1280, 1024, 1)
    g = ICP.drawCDF(xMin, xMax)
    drawable = userDefined.drawCDF(xMin, xMax, 10 * int(xMax - xMin + 1)).getDrawable(0)
    drawable.setColor("blue")
    drawable.setLineStyle("dotted")
    g.addDrawable(drawable)
    g.draw("ICPCDF", 1280, 1024, 1)
