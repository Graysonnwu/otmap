import matplotlib.pyplot as plt
import numpy as np

dat = np.loadtxt('x8.dat')

plt.figure(figsize=(8, 6), facecolor='black')

# for i, flag in enumerate(dat[:,4]):
#     if flag == -1:
#         dat[i,4] = 1
#     else:
#         dat[i,4] = 0

plt.scatter(-dat[:,3], -dat[:,2], s=0.1, color='white', alpha=.2)
# plt.scatter(dat[:,0], dat[:,2], s=0.1, color='white', alpha=.5)
# plt.scatter(-dat[:,3]*dat[:,4], -dat[:,2]*dat[:,4], s=0.1, color='white', alpha=1)

plt.axis('equal')
plt.axis('off')

plt.show()