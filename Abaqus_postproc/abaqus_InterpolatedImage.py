#!/usr/bin/env python3

import math,numpy as np
import scipy 
import damask
import os,sys,string
from subprocess import call 
from optparse import OptionParser
from scipy.interpolate import griddata


scriptID   = str.replace('$Id: abaqus_pileUpRescale.py 247 2019-01-10 z.zhao modified based on Aritra add_InterpolatedImage.py $','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

def radius_rescale(stop, final, R):
  if final > (1-0.5*np.sqrt(2))*R:
    if stop > (1-0.5*np.sqrt(2))*R:
      return (final+(np.sqrt(2)-1)*R)/(stop+(np.sqrt(2)-1)*R)
    elif stop <= (1-0.5*np.sqrt(2))*R:
      return (final+(np.sqrt(2)-1)*R)/np.sqrt(stop*(2*R-stop))
  elif final <= (1-0.5*np.sqrt(2))*R: 
    return np.sqrt(final*(2*R-final))/np.sqrt(stop*(2*R-stop))
    
def z_rescale(stop, final):
  fc2r3 = np.array([0.0134,0.0296,0.0485,0.0681,0.0889,0.1058,0.1260,0.1457,0.1655,0.1863,0.2055,0.2316,0.2582])
  fc2r4 = np.array([0.0136,0.03059,0.0485,0.0694,0.0909,0.1083,0.1233,0.1427,0.1611,0.1793,0.2009,0.2252,0.2505])
  fc3r2 = np.array([0.0124,0.0283,0.0476,0.0655,0.0857,0.1068,0.1247,0.1424,0.1605,0.1807,0.2074,0.2310,0.2497])
  xp = np.array([50,100,150,200,250,300,350,400,450,500,550,600,650])
  fp = (fc2r3+fc2r4+fc3r2)/3.0
  return np.interp(final,xp,fp)/np.interp(stop,xp,fp)    

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Rescale regular grid and gives the resulting image. 

""", version = scriptID)

parser.add_option('-c','--coords',
                  dest   = 'coords',
                  type = 'string', metavar = 'string',
                  help   = 'column label of xyz coordinates')
parser.add_option('--rescale',
                  dest = 'rescale',
                  type = 'float', metavar = 'float',
                  help = 'rescale data size xyz same way [%default]')
parser.add_option('--incomplete',  action="store_true",
                  dest="incomplete",
                  help="indicates simualtation does not converge completely [False]")                  
parser.add_option('--stop_depth',
                  dest = 'stop',
                  type = 'float', metavar = 'float',
                  help = 'depth where simualtion stopps in nm')
parser.add_option('--final_depth',
                  dest = 'final',
                  type = 'float', metavar = 'float',
                  help = 'depth where complete simualtion aims in nm')
parser.add_option('--radius',
                  dest = 'R',
                  type = 'float', metavar = 'float',
                  help = 'tip radius in micron')
parser.add_option('--tilt',
                  dest = 'tilt',
                  type = 'float', metavar = 'float',
                  help = 'angle for surface tilt in degree along y axis counterclockwise [%default]')                                                                                                               
parser.add_option('--rotation',
                  dest = 'rotation',
                  type = 'float', metavar = 'float',
                  help = 'angle for surface rotation in degree along z axis counterclockwise [%default]')                                                                                
parser.add_option('--grid',
                  dest = 'grid',
                  type = 'int', nargs = 2, metavar = 'int int',
                  help = 'interpolation grid')                 
parser.add_option('--size',
                  dest = 'size',
                  type = 'float', nargs = 2, metavar = 'float float',
                  help = 'interpolation size')
parser.add_option('--center',
                  dest = 'center',
                  type = 'float', nargs = 2, metavar = 'float float',
                  help = 'coordinates of interpolation patch center')
parser.add_option('-p','--pixelsize',
                  dest = 'pix_size',
                  type = 'string', metavar = 'string',
                  help = 'pixel size [20.0/255]')

#making the default values and let them show
parser.set_defaults(rescale = 1.0,
                    tilt = 0.0,
                    rotation = 0.0,
                    incomplete = False, 
                    )


(options,filenames) = parser.parse_args()
#----------------------------------------  sanity checks   ------------------------------------------------

if options.pix_size:
  options.pix_size = float(eval(options.pix_size))
  if options.grid: 
    options.size = tuple(options.pix_size * (x - 1) for x in options.grid)
  elif options.size:
    options.grid = tuple(round(x/options.pix_size + 1) for x in options.size)
    options.size      = tuple(options.pix_size * (x - 1) for x in options.grid)
  else:
    parser.error("Either dimension or size has to be specified if pixel size is given.")
else:
  if options.size and options.grid:
    options.pix_size = options.size/options.grid
  else:
    parser.error("Both dimension and size has to be specified if pixel size is not given.")

# --------------------------------------- loop over input files -------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  out_file = "out_"+os.path.basename(name)
  try:
    table = damask.ASCIItable(name     = name, 
                              outname  = os.path.join(os.path.dirname(name),out_file),
                              buffered = False)
  except: continue
  damask.util.report(scriptName,name)
# ------------------------------------------ read header and data ------------------------------------------

  table.head_read()
  table.data_readArray([options.coords])
  if len(table.data[0]) != 3:
    continue 

#-------------------------------------------- process and store output ---------------------------------------
  if options.tilt != 0.0:                                                                   #tilt before coordinates rotation
    ang_tilt = np.radians(options.tilt)
    tilt = np.zeros((len(table.data),2),dtype=float)
    tilt[:,0] = table.data[:,0]*np.cos(ang_tilt) + table.data[:,2]*np.sin(ang_tilt)
    tilt[:,1] = -table.data[:,0]*np.sin(ang_tilt) + table.data[:,2]*np.cos(ang_tilt)
    table.data[:,0] = tilt[:,0]  
    table.data[:,2] = tilt[:,1] 
    
  if options.incomplete:
    table.data[:,:2] *= radius_rescale(options.stop/1e3, options.final/1e3, options.R) 
    table.data[:,2] = (table.data[:,2]-3.0)*z_rescale(options.stop, options.final)
  else:  
    table.data[:,:2] *= options.rescale
    table.data[:,2] = (table.data[:,2]-3.0)*options.rescale  
      
  if not options.center:
    options.center = 0.5*(table.data[:,:2].max(axis=0)+table.data[:,:2].min(axis=0))
 
  if options.rotation != 0.0:
    ang_rot = np.radians(options.rotation)
    rot = np.zeros((len(table.data),2),dtype=float)
    rot[:,0] = table.data[:,0]*np.cos(ang_rot) - table.data[:,1]*np.sin(ang_rot)
    rot[:,1] = table.data[:,0]*np.sin(ang_rot) + table.data[:,1]*np.cos(ang_rot)
    table.data[:,:2] = rot[:,:2] 

  grid_x, grid_y = np.meshgrid(np.linspace(options.center[0] - 0.5 * options.size[0],
                                           options.center[0] + 0.5 * options.size[0], num=options.grid[0]),
                               np.linspace(options.center[1] - 0.5 * options.size[1],
                                           options.center[1] + 0.5 * options.size[1], num=options.grid[1]))
  grid = np.vstack((grid_x.flatten(),grid_y.flatten())).T

  interpolation = griddata(table.data[:,:2], table.data[:,2],grid,method='linear')
  table.data = np.vstack((grid_x.flatten().T,
                          grid_y.flatten().T,
                          interpolation.T)).T 
                        
                                                                                                                           
#--------------------------------------------------- output header info --------------------------------------

  table.labels_clear()
  table.labels_append(['{}_gridInterpolation'.format(1+i) for i in range(3)])
  table.info_clear()
  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.head_write()

  table.data_writeArray()

  table.close()
