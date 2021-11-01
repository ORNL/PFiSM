#!/usr/bin/env python
import sys
import subprocess
import os

print("Test PlanarFrontUniformT...")

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

previous_sfraction = -1.
previous_time = -1.

tol_percent = 0.025 #2.5%
tol_velocity = 6.5*tol_percent

ly = 400. #length of domain
for line in lines:
  if line.count(b'step'):
    print(line)
    words=line.split()
    time = eval(words[6][:-1])
    print(time)

  if line.count(b'fraction'):
    print(line)
    words=line.split()
    sfraction=eval(words[2])

    if previous_sfraction>0.:
      dfs=sfraction-previous_sfraction
      delt=time-previous_time
      velocity = ly*dfs/delt
      print("velocity={}".format(velocity))

      if time>20:
        if abs(velocity-6.5)>tol_velocity:
          print("Wrong velocity!")
          sys.exit(1)

    previous_sfraction = sfraction
    previous_time      = time


sys.exit(0)

