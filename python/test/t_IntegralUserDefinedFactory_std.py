from __future__ import print_function
from openturns import *
from oticp import *

support = [1, 2, 3, 4]
weights = [0.1, 0.2, 0.3, 0.8]

d = IntegralUserDefined(support, weights)
print("d=", d)
sample = d.getSample(10000)
d2 = IntegralUserDefinedFactory().buildAsIntegralUserDefined(sample)
print("d2=", d2)
