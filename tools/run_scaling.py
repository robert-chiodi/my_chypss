#!/usr/bin/python

###################################################################

#Modules needed
import math
import sys
import string
import numpy as np
import shutil
import glob, os
import subprocess

# Base case is 64 Elements
# Elements at refined level calculated as 64*(2^(3*t)),
# with t refinement level (and 0 being no refinement)

print("Running job creation and submission script for scaling.\n")
cores_per_node = 56

strong_starting_cores = 32
strong_scaling_num_elem = 16777216
strong_scaling_growth = 8 # Cores at level is strong_starting_cores*(2**t)
print("Planned strong scaling study:")
print("NumElem".center(10)+"Cores".center(10)+"Nodes".center(10)+"Elem/Core".center(10))
for t in range(0,strong_scaling_growth):
  elem = strong_scaling_num_elem
  cores = int(strong_starting_cores*(2**t))
  nodes = int(math.floor(cores/cores_per_node))
  if cores - nodes*cores_per_node > 0:
    nodes = nodes + 1
  print(str(elem).center(10)+ \
        str(cores).center(10)+ \
        str(nodes).center(10)+ \
        str(elem/cores).center(10))
print("\n")

weak_starting_cores = 1
weak_scaling_num_elem_start = 32768
weak_scaling_growth = 5 # Cores at level is weak_starting_cores*(2**((t+3)*3))
print("Planned weak scaling study:\n")
print("NumElem".center(10)+"Cores".center(10)+"Nodes".center(10)+"Elem/Core".center(10))
for t in range(0,weak_scaling_growth):
  elem = weak_scaling_num_elem_start*(2**(3*t))
  cores = int(weak_starting_cores*(2**(t*3)))
  nodes = int(math.floor(cores/cores_per_node))
  if cores - nodes*cores_per_node > 0:
      nodes = nodes + 1
  print(str(elem).center(10)+ \
        str(cores).center(10)+ \
        str(nodes).center(10)+ \
        str(elem/cores).center(10))
print("\n")

print("Options: ")
print("Enter R to run scaling test (via submitting jobs with sbatch")
print("Enter D to generate all files but not submit jobs")
user_input = input()
if(user_input != "R" and user_input != "r" \
   and user_input != "D" and user_input != "d"):
  quit()

submit_job = False
if(user_input == "R" or user_input == "r"):
  submit_job = True



strong_scale_dir = "strong_scaling"
os.mkdir(strong_scale_dir)
os.chdir("./"+strong_scale_dir)

for t in range(2,strong_scaling_growth): # 7 is max for flex queue on Frontera, ends with 2340 cells per core
  dir = "simdir"+str(t)
  os.mkdir(dir)
  os.chdir("./"+dir)
  execute_command = "../../scale_test -m ../../inline-hex.mesh -rs
  6" # This is 16,777,216 Elements
  r = open("run","w")
  r.write("#!/bin/bash \n")
  r.write("#SBATCH -p flex \n")
  r.write("#SBATCH -t 00:10:00 \n")
  total_cores = strong_starting_cores*(2**(t))
  num_nodes = int(math.floor(total_cores/cores_per_node))
  if total_cores - num_nodes*cores_per_node > 0:
    num_nodes = num_nodes + 1
  r.write("#SBATCH -N "+str(int(num_nodes))+"\n")
  r.write("#SBATCH --ntasks-per-node "+str(int(cores_per_node))+"\n")
  r.write("#SBATCH -n "+str(int(total_cores))+"\n") # Number of cores/node on Frontera
  r.write("#SBATCH -J MFEM_scale \n")
  r.write("#SBATCH -o output \n")
  r.write("#SBATCH -e error \n")

  r.write("ibrun  "+execute_command+"\n")
  t_string = str(t)
  r.write("mv ./output ../output_"+t_string.zfill(4)+"\n")
  r.close()

  if(submit_job):
    command = "sbatch run"
    subprocess.check_call(command.split(), shell=False)
  os.chdir("../")
os.chdir("../")

weak_scale_dir = "weak_scaling"
os.mkdir(weak_scale_dir)
os.chdir("./"+weak_scale_dir)
for t in range(0,weak_scaling_growth):
  dir = "simdir"+str(t+3)
  os.mkdir(dir)
  os.chdir("./"+dir)
  execute_command = "../../scale_test -m ../../inline-hex.mesh -rs "+str(t+3) # 64*(2^(3*(t+3)))
  r = open("run","w")
  r.write("#!/bin/bash \n")
  r.write("#SBATCH -p flex \n")
  r.write("#SBATCH -t 00:10:00 \n")
  total_cores = weak_starting_cores*(2**(t*3))
  num_nodes = int(math.floor(total_cores/cores_per_node))
  if total_cores - num_nodes*cores_per_node > 0:
    num_nodes = num_nodes + 1
  r.write("#SBATCH -N "+str(int(num_nodes))+"\n")
  r.write("#SBATCH --ntasks-per-node "+str(int(cores_per_node))+"\n")
  r.write("#SBATCH -n "+str(int(total_cores))+"\n") # Number of cores/node on Frontera
  r.write("#SBATCH -J MFEM_scale \n")
  r.write("#SBATCH -o output \n")
  r.write("#SBATCH -e error \n")

  r.write("ibrun  "+execute_command+"\n")
  t_string = str(t)
  r.write("mv ./output ../output_"+t_string.zfill(4)+"\n")
  r.close()

  if(submit_job):
    command = "sbatch run"
    subprocess.check_call(command.split(), shell=False)
  os.chdir("../")
os.chdir("../")

