import numpy as np
import matplotlib.pyplot as plt

# Load the columns of Chi2Data.txt as indidual lists
a2 = np.loadtxt("Chi2Data.txt", unpack = True)[0]
Der1 = np.loadtxt("Chi2Data.txt", unpack = True)[1]
DerN1 = np.loadtxt("Chi2Data.txt", unpack = True)[2]

# Sorting the data
data = np.loadtxt("Chi2Data.txt", dtype = 'd, d, d, d')
sortdata = sorted(data, key = lambda tup: tup[0])
a2sorted = [x[0] for x in sortdata]
ratio = [x[3] for x in sortdata]

# Print the number of points generated
print("The number of points generated is:", len(a2sorted))

# Plot the analytical and the numerical gradient vs a1
fig = plt.figure()

plt.subplot(211)
plt.plot(a2, Der1, 'o', color = 'r', label = 'Analytical')
plt.plot(a2, DerN1, 'o', color = 'b', label = 'Numerical')
plt.title('Gradient Xi2')
plt.ylabel('dXi2/a1 (a1-h)')
plt.legend()

plt.subplot(212)
plt.plot(a2sorted, ratio, color = 'g', label = 'ratio')
plt.plot(a2sorted, [1 for i in range(len(a2sorted))], color = 'r')
plt.xlabel('a1-h')
plt.legend()

# plt.show()
fig.savefig('gradientX2.png')
