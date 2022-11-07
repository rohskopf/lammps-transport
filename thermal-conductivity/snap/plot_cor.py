"""
This script loops over heat flux autocorrelation files calculated from LAMMPS. For each file, we 
get the last autocorrelation time window, which is averaged over all previous time windows; this 
is our final autocorrelation that is time-integrated to calculate thermal conductivity. That 
time-integration is handled in `in.run` for us. 

Here, we simply plot the heat flux autocorrelation vs. correlation time for each file. This shows 
us whether the autocorrelation has converged or not.

NOTE: There are some things that need to be hard-coded here:

      `final_autocor_str`: (str) Output of the last autocorrelation block in J0Jt.dat 
      `dt`: (float) Timestep setting from LAMMPS input script
       
       Also change x-axis time units on matplotlib plot.

Run with:

    python plot_cor.py 

"""

import numpy as np
import matplotlib.pyplot as plt

"""
filenames = ["e_5/J0St.dat",
             "2p4kmps_158GPa_3650K_1000_preheat/e_3/J0St.dat",
             "3p5kmps_6840K_288GPa/e_3/J0St.dat",
             "3p5kmps_7200K_317GPa_600_preheat/e_3/J0St.dat"]
"""

filenames = ["J0Jt.dat"]

final_autocor_str = "100000 200"
dt = 4.0 # timestep setting from LAMMPS input script

for filename in filenames:

    fh = open(filename, 'r')

    line = fh.readline()
    # loop over lines until you reach last section of S0St.dat, which is time averaged across all previous sections.
    while (final_autocor_str not in line):
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
        time.append(timestep*dt)
        corr.append(y_avg)
    time = np.array(time)
    corr = np.array(corr)
    # normalize the correlation
    corr = corr/corr[0]

    fh.close()

    plt.plot(time, corr)
    plt.xlim([time[0], time[-1]])
plt.xlabel("Correlation time (ps)")
plt.ylabel("Normalized HFACF")
plt.savefig("correlation_plot.png", dpi=500)
