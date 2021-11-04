#!/usr/bin/env python3
import sys
import subprocess
import os

print("Test PlanarFrontNoP...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-2):
  mpicmd = mpicmd + " "+sys.argv[i]
exe = sys.argv[nargs-2]
inp = sys.argv[nargs-1]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
print("Run command: {}".format(command))

output = subprocess.check_output(command,shell=True)

#analyse PFiSM standard output
lines=output.split(b'\n')

time=0.
sfraction=0.
for line in lines:
  if line.count(b'step'):
    #print(line)
    words=line.split()
    time = eval(words[6].split(b',')[0])

  if line.count(b'fraction'):
    #print(line)
    words=line.split()
    sfraction=eval(words[2])

print("Final time: {}".format(time))
if time < 28.:
  print("Expected time not reached!!!")
  sys.exit(1)

print("Final solid fraction: {}".format(sfraction))
if sfraction < 0.52:
  print("Expected solid fraction not reached!!!")
  sys.exit(1)

sys.exit(0)
