# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 17:11:27 2020

@author: Hayden
"""
import freud

# Make a 10x10 square box (for 2-dimensional systems).
freud.box.Box.square(10)

# Make a 10x10x10 cubic box.
freud.box.Box.cube(10)



# Create a 10x10 square box from a list of two items.
freud.box.Box.from_box([10, 10])

# Create a 10x10x10 cubic box from a list of three items.
freud.box.Box.from_box([10, 10, 10])

# Create a triclinic box from a list of six items (including tilt factors).
freud.box.Box.from_box([10, 5, 2, 0.1, 0.5, 0.7])

# Create a triclinic box from a dictionary.
freud.box.Box.from_box(dict(Lx=8, Ly=7, Lz=10, xy=0.5, xz=0.7, yz=0.2))

# Directly call the constructor.
freud.box.Box(Lx=8, Ly=7, Lz=10, xy=0.5, xz=0.7, yz=0.2, dimensions=3)


positions = ...  # Read positions from trajectory file.
op = freud.order.Hexatic(k=6)
op.compute(
    system=({'Lx': 5, 'Ly': 5, 'dimensions': 2}, positions),
    neighbors=dict(r_max=3)
)





## MD ANALYSIS library - Useful for loading gromacs data
import glob
from MDAnalysis import Universe

testdir = 'C:\\Users\\Hayden\\Documents\\Patey_Lab\\Results\\LiI\\Step_W_JC_NPT\\'
# Note topology is a gro file not a top file
topology = glob.glob(testdir + '*.gro')[0]
trajectory = glob.glob(testdir + '*.trr')[0] # this can be a trr or xtc file

# Creates a "universe" to contain your trajectory
u = Universe(topology, trajectory)

# Access the coordinates of a specific atom at a specific time as follows
time_index = 10000
atom_index = 5000
u.trajectory[time_index][atom_index]

print(u.trajectory[time_index][atom_index])

# Can also select atoms using the select_atoms() method
lithiums = u.select_atoms("name LI")


