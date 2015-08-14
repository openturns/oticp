from __future__ import print_function
from openturns import *
from oticp import *

support = [1, 2, 3, 4]
weights = [0.1, 0.2, 0.3, 0.8]

d = IntegralUserDefined(support, weights)
print("d=", d)
print("support=", d.getSupport())
print("weigths=", d.getWeights())
print("normalized weigths=", d.getNormalizedWeights())
print("realization=", d.getRealization())
print("sample=")
print(d.getSample(20))
print("pdf(2)=", d.computePDF(2))
print("cdf(2)=", d.computeCDF(2))
