### Si vibrational density of states

Run the MD simulation with: (WARNING: Generates ~0.1 GB file)

    lmp -in in.run

Process the velocities and calculate the power spectrum with:

    python calc_dos.py
