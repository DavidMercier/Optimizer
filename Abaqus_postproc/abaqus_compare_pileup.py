#!/usr/bin/env python3

import os,sys,damask,string
import numpy as np
import math
from optparse import OptionParser

scriptID   = str.replace('$Id: compare_pileup.py 154 2015-11-06 14:33:26Z chakra34 $','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Return relative difference between simulated and AFM measured surface topography
""", version = scriptID)

parser.add_option('--sim',
                  dest = 'sim_file',
                  type = 'string', metavar = 'string',
                  help = 'simulated/file1 to compare with AFM file')

parser.add_option('--error',  action="store_true",
                  dest="error",
                  help="prints out the abs error [False]")

parser.add_option('--multiplier',
                  dest = 'multiplier',
                  type = 'float',
                  help = 'multiplier for matching the values of the two datasets')

parser.add_option('--flipud',  action="store_true",
                   dest="flipud",
                   help="flips up/down - generally required for simulation [False]")
                   

parser.add_option('--step',
                  dest = 'step',
                  type = 'int',
                  help = 'step range for cross correlation')

parser.set_defaults(error      = False,
                    flipud     = False,
                    step       = 15,
                   )

(options,filenames) = parser.parse_args()


if options.multiplier == None :
  print("There has to be a multiplier!!!, aborting")
  sys.exit()


if not options.sim_file or not os.path.exists(options.sim_file):
  parser.error('No file selected for comparison')



#------------------- Reading the simulated tip from file -----------------

sim = damask.ASCIItable(name = options.sim_file, buffered = False,readonly=True)
sim.head_read()
sim.data_readArray()
if options.flipud == True: sim.data = np.flipud(sim.data)

#------------------- loop over input files ------------------------------------ 
if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              readonly = True,
                             )
  except: continue

  damask.util.report(scriptName,name)

  table.head_read()
  table.data_readArray()
  table.data *= options.multiplier
  
  xdim = table.data.shape[1]                 #---->x
  ydim = table.data.shape[0]
  
  relative_error = []
  for i in range(options.step):
    for j in range(options.step):
      error1 = table.data[i:,j:] - sim.data[:ydim-i,:xdim-j]
      abs_error1 = np.abs(table.data[i:,j:] - sim.data[:ydim-i,:xdim-j])
      error2 = table.data[:ydim-i,j:] - sim.data[i:,:xdim-j]
      abs_error2 = np.abs(table.data[:ydim-i,j:] - sim.data[i:,:xdim-j])
      error3 = table.data[i:,:xdim-j] - sim.data[:ydim-i,j:]
      abs_error3 = np.abs(table.data[i:,:xdim-j] - sim.data[:ydim-i,j:])
      error4 = table.data[:ydim-i,:xdim-j] - sim.data[i:,j:]
      abs_error4 = np.abs(table.data[:ydim-i,:xdim-j] - sim.data[i:,j:])
      relative_error.append(np.nansum(abs_error1)/np.nansum(table.data[i:,j:]))
      relative_error.append(np.nansum(abs_error2)/np.nansum(table.data[:ydim-i,j:]))
      relative_error.append(np.nansum(abs_error3)/np.nansum(table.data[i:,:xdim-j]))
      relative_error.append(np.nansum(abs_error4)/np.nansum(table.data[:ydim-i,:xdim-j]))
  
  # relative_error.sort()

  id = np.argmin(relative_error)//(options.step*4)
  jd = (np.argmin(relative_error)%(options.step*4))//4
  mark = (np.argmin(relative_error)%(options.step*4))%4
  
#----------------- printing the absolute error if needed ---------------------------#
  if options.error :
    if mark==0:
      error = table.data[id:,jd:] - sim.data[:ydim-id,:xdim-jd]
    elif mark==1:
      error = table.data[:ydim-id,jd:] - sim.data[id:,:xdim-jd]
    elif mark==2:
      error = table.data[id:,:xdim-jd] - sim.data[:ydim-id,jd:]
    elif mark==3:
      error = table.data[:ydim-id,:xdim-jd] - sim.data[id:,jd:]
    with open("diff_topo.txt",'a') as diff:
      diff.write("2 head \n")
      diff.write("relative error {}, argmin_id:{}, i:{}, j:{} mark: {}\n".format(min(relative_error),np.argmin(relative_error),id,jd,mark))
      diff.write("values \n")  
      np.savetxt(diff, error, delimiter=' ', newline='\n')  
           
  print(min(relative_error))
#  


  table.close
