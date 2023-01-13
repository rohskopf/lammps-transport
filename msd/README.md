Run the LAMMPS script with

    lmp -in in.run

to create unwrapped trajectory files. Then do a time-averaged MSD vs. time 
calculation using FFTs:

    python calc_msd.py

