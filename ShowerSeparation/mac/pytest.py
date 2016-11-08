import matplotlib.pyplot as plt

import numpy as np


circle = np.array([[473.046, 2268.6],[473.007, 2269.23],[472.889, 2269.85],[472.695, 2270.44],[472.428, 2271.01],[472.091, 2271.54],[471.691, 2272.03],[471.233, 2272.46],[470.726, 2272.82],[470.175, 2273.13],[469.591, 2273.36],[468.983, 2273.51],[468.36, 2273.59],[467.732, 2273.59],[467.109, 2273.51],[466.501, 2273.36],[465.917, 2273.13],[465.367, 2272.82],[464.859, 2272.46],[464.402, 2272.03],[464.001, 2271.54],[463.665, 2271.01],[463.397, 2270.44],[463.203, 2269.85],[463.086, 2269.23],[463.046, 2268.6],[463.086, 2267.98],[463.203, 2267.36],[463.397, 2266.76],[463.665, 2266.19],[464.001, 2265.66],[464.402, 2265.18],[464.859, 2264.75],[465.367, 2264.38],[465.917, 2264.08],[466.501, 2263.85],[467.109, 2263.69],[467.732, 2263.61],[468.36, 2263.61],[468.983, 2263.69],[469.591, 2263.85],[470.175, 2264.08],[470.726, 2264.38],[471.233, 2264.75],[471.691, 2265.18],[472.091, 2265.66],[472.428, 2266.19],[472.695, 2266.76],[472.889, 2267.36],[473.007, 2267.98] ])

vertex  = np.array([ 468.046, 2268.6 ] )

plt.plot(circle[:,0],circle[:,1],'bo')
plt.plot(vertex[0],vertex[1],'gx')
plt.show()

