import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, write
import ase

print(ase.units.fs)

# Load lammps trajectory file with velocities.
# Seems like ASE sorts by atom id automatically.
# NOTE: ASE converts to its own units internally, so A/ps velocities from LAMMPS convert to A/fs in ASE.
#       However, ASE's definition of velocity is a bit weird, so we multiply by ase.units.fs below.
#       Then multiply by 100 to get units back in A/ps.

frames = read("dump.lammpstrj", format="lammps-dump-text", index=":")

nframes = len(frames)
print(f"{nframes} frames")
natoms = len(frames[0])
print(f"{natoms} atoms")

# create time data of velocities (nframes, natoms, 3)

v = np.zeros((nframes,natoms,3))
x = np.zeros((nframes,natoms,3))
types = np.zeros((nframes, natoms)) # atom types
indx_frame = 0
for atoms in frames:
    v[indx_frame] = atoms.get_velocities()*ase.units.fs*1e3
    x[indx_frame] = atoms.get_positions()
    types[indx_frame] = atoms.numbers
    indx_frame += 1

# Now we have velocity signal with shape (nframes, natoms, 3).
# To take average power spectrum of each component, let's convert this to shape (nframes, natoms*3).

v = v.reshape((nframes, natoms*3))

# Calculate power spectrum of each velocity component signal.
# First define sampling interval; data collected every this many ps.
# E.g. 0.5 fs timestep, dumping every 5 timesteps --> sampling_interval = 0.5e-3 * 5
time_step = 0.5e-3 # timestep in fs
sampling_interval = 0.5e-3 * 5
sampling_frequency = 1.0 / sampling_interval
time_points = np.arange(0, nframes*sampling_interval, sampling_interval)

# Determine frequency space axis
tp_count = len(time_points)
values      = np.arange(int(tp_count/2))
time_period  = tp_count/sampling_frequency
frequencies = values/time_period

# Loop through velocity components and calculate power spectrums
# TODO: List comprehension could be faster.
#       This could be parallelized with an all_reduce operation for averaging.
len_ps = int(tp_count/2) # length of power spectrum is half so we exclude sampling frequency
ps_all = np.zeros((len_ps, 3*natoms))
for i in range(0,3*natoms):
  # Get the signal
  amplitude = v[:,i] # Need to index +1 beacuse first column is time.
  mean_amp = np.mean(amplitude)
  amplitude = amplitude - mean_amp
  # Frequency domain representation
  fourier_transform = np.fft.fft(amplitude)/len(amplitude)           # Normalize amplitude
  fourier_transform = fourier_transform[range(int(len(amplitude)/2))] # Exclude sampling frequency
  # Populate the collection of power spectrums
  ps_all[:,i] = abs(fourier_transform)**2 # power spectrum

ps_avg = np.mean(ps_all, axis=1)
plt.plot(frequencies, ps_avg)
plt.xlim([0,20])
plt.show()

