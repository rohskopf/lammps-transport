import numpy as np
import matplotlib.pyplot as plt

fh = open("S0St.dat", 'r')

line = fh.readline()
# loop over lines until you reach last section of S0St.dat, which is time averaged across all previous sections.
while ("250000 400" not in line):
    line = fh.readline()

ndat = int([float(x) for x in line.split()][1])
time = []
corr = []
for i in range(0,ndat):
    line = fh.readline()
    line_split = line.split()
    timestep = int(line_split[1])
    y1 = float(line_split[3])
    y2 = float(line_split[4])
    y3 = float(line_split[5])
    y_avg = (y1+y2+y3)/3
    time.append(timestep*0.0005)
    corr.append(y_avg)
time = np.array(time)
corr = np.array(corr)
# normalize the correlation
corr = corr/corr[0]

fh.close()

plt.plot(time, corr)
plt.xlim([time[0], time[-1]])
plt.xlabel("Correlation time (ps)")
plt.ylabel("Normalized stress correlation")
plt.savefig("correlation_plot.png", dpi=500)
