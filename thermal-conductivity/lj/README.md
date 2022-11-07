# Simple LJ example

Run this example with:

    lmp < in.run

And it'll give you the thermal conductivity value at the end of the simulation. This value comes 
from averaging many heat flux autocorrelation windows throughout the simulation. To be sure that 
we calculated thermal conductivity properly, we must look at the heat flux autocorrelation function 
(HFACF). This can be plotted for the last time-averaged window by doing

    python plot_cor.py

As you can see from the output `correlation_plot.png`, the HFACF converged to near-zero within the 
correlation time window. This means that our system is in equilibrium and we can rely on the 
time-integral of this correlation to give us thermal conductivity. That integral is alredy 
calculated in the `in.run` LAMMPS script, but this plot tells us that the calculation is properly
converged:

![alt text](https://github.com/rohskopf/lammps-transport/blob/main/thermal-conductivity/lj/correlation_plot.png?raw=true)


