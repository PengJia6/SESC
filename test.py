import numpy as np
import matplotlib.pyplot as plt
lens=np.random.exponential(scale=100,size=1000).round()
plt.hist(lens,bins=100)
print(lens)
plt.show()
