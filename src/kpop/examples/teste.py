from kpop import *
import numpy as np
from matplotlib import pyplot as plt

p1 = Population.make_random(100, 100)
p2 = Population.make_random(100, 100)

d = np.array(p1 + p2)
(p1 + p2).plot.pca()
plt.show()
