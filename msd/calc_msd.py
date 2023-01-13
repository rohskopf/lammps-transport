import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, write

def autocorrFFT(x):
    N=len(x)
    F = np.fft.fft(x, n=2*N)  #2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res= (res[:N]).real   #now we have the autocorrelation in convention B
    n=N*np.ones(N)-np.arange(0,N) #divide res(m) by (N-m)
    return res/n #this is the autocorrelation in convention A

def msd_fft(r, n):
    """
    Following section 4 (4.1, 4.2) in this paper:
    https://www.neutron-sciences.org/articles/sfn/abs/2011/01/sfn201112010/sfn201112010.html

    Args:
        r : Numpy array of positions (nframes, natoms, 3)
        n : Index of particular atom to calculate signal for
    """
    #N=len(r)
    N = np.shape(r)[0] # number of frames
    D=np.square(r).sum(axis=2) # sum over Cartesian coordinates
                               # now size (nframes, natoms)

    #D = np.square(r)[:,0] # x-component

    #print(np.shape(D))

    # append a zero onto each atom signal
    #D=np.append(D,0) 
    nframes = np.shape(D)[0]
    natoms = np.shape(D)[1]
    npad = ((0, 1), (0, 0))
    D = np.pad(D, pad_width=npad, mode='constant', constant_values=0)
    # D is now size (nframes+1, natoms), since we appended a zero for each atom

    # calculate for each atom n
    #n = 0
    S2=sum([autocorrFFT(r[:, n, a]) for a in range(r.shape[2])]) # size (nframes,) for a particular atom

    #print(np.shape(D))

    # sum over all times for all atoms
    #Q=2*D.sum() 
    Q=2*D.sum(axis=0) # now size (natoms,)

    S1=np.zeros(N)
    for m in range(N):
        Q[n]=Q[n]-D[m-1,n]-D[N-m,n]
        S1[m]=Q[n]/(N-m)
    return S1-2*S2

# load unwrapped xyz file
# seems like ASE sorts by atom id automatically

frames = read("dump.lammpstrj", format="lammps-dump-text", index=":")

nframes = len(frames)
print(f"{nframes} frames")
natoms = len(frames[0])
print(f"{natoms} atoms")

# create time data of positions (nframes, natoms, 3)

positions = np.zeros((nframes,natoms,3))
types = np.zeros((nframes, natoms)) # atom types
indx_frame = 0
for atoms in frames:
    positions[indx_frame] = atoms.get_positions()
    types[indx_frame] = atoms.numbers
    indx_frame += 1

#print(types[0])
#print(atoms.get_tags())
#print(atoms.numbers)
#print(atoms)
#print(positions[0])
#assert(False)
#print(positions[1,5,2])
#print(frames[1].get_positions()[5,2])
#assert(positions[1,5,2] == frames[1].get_positions()[5,2])
#print(np.shape(positions))

natoms_Li = 19 # Li ions are indexed from 0 to this index
msd_atoms = [msd_fft(positions, n) for n in range(0,natoms_Li)]
msd_atoms = np.array(msd_atoms) # size (natoms, nframes)
msd = np.mean(msd_atoms, axis=0) # size (nframes,)

dt = 2e-3 # 2 fs
nsample = 1 # sample every this many timesteps
time = np.zeros(nframes)
for t in range(0,nframes):
    time[t] += t*nsample*dt

indx_half = int(nframes/2)

tau = np.array([time[:indx_half]]).T
msd = np.array([msd[:indx_half]]).T
data = np.concatenate((tau,msd), axis=1)
np.savetxt("msd.dat", data)

plt.plot(time[:indx_half], msd[0:indx_half])
plt.xlabel(r"$\Delta \tau (ps)$")
plt.ylabel(r"MSD $(\AA^2)$")
plt.savefig("msd.png", dpi=500)
