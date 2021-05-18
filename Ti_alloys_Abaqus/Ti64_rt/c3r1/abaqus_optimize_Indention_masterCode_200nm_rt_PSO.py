#!/usr/bin/env python3

import numpy as np
import sys,os,shutil
import time
import Optimization
from optparse import OptionParser
from subprocess import run,check_output
import damask

#------------------------------------------------------------------------------------------------- #
_ratio = None
class optimize(Optimization.Optimization):

  def id(self,x):
    return str(self.map2space(x[:self.dimension])).translate(str.maketrans(' []','___')) 
      
#===========================================
  def fitness(self,x):
    time.sleep(np.random.randint(2,300)*2)
    if self.id(x) in self.locations:
      self.curr_locations.append(np.append(x,self.locations[self.id(x)]))
      return self.locations[self.id(x)]
    with open('{}/output_gen{}_{}.log'.format(options.root,self.generation+1,self.id(x)),'a') as file:
      file.write("\n Generation %i "%(self.generation+1))

      orientations = np.array([ [ 9.967, 75.654, 336.86 ],
                              ])
      error = []
      for orientation_run,orientation in enumerate(orientations):
        dir_loc = "{}/Generation_{}/Test_{}/Orientation_{}".format(self.root,str(self.generation+1),self.id(x),str(orientation_run+1))
        os.makedirs(dir_loc,exist_ok=True)
        msg = '\n'.join(['generation {} location {} orientation {}'.format(self.generation+1,self.id(x),str(orientation_run+1))])
        print(msg)
        file.write(msg+'\n')

        rep = {}                                                                                        # dictionary representing the parameters
        for i in range(len(x[:self.dimension])):
          rep["%i_coords"%(i+1)] = self.map2space(x[:self.dimension])[i]                                                                 # parameters are in MPa 
        for j in range(3):
          rep["{}_eulerangles".format(j+1)] = orientation[j]                                                                 # parameters are in MPa 

        if self.job_id != None:
          shutil.copy("%s/%s.inp"%(self.root,self.job_id),'%s/'%dir_loc)                               #  *.inp file for Abaqus
          print("%s/%s.inp"%(self.root,self.job_id))
          print("%s/"%dir_loc)
#           shutil.copy("%s/abaqus_postResults.py"%(self.root),'%s/'%dir_loc)                               #  post processing code
#           shutil.copy("%s/%s"%(self.root,self.exp_file),'%s/'%dir_loc)

        with open('{}/{}'.format(self.root,self.mat_file),'r') as file_in:
          contents = file_in.read()
        for key,value in list(rep.items()):
          contents = contents.replace(str(key),str(value) )
        with open("{}/material.config".format(dir_loc),'w') as f_out:
          f_out.write(contents)

  #--------- Running damask postProcessing and finding error against "experimental" data ----------- #
        print('Start checking server load')
        availability={}
        cpu_num = 4
        for i in [1,4,9]:
          cmd0 = 'ssh compute{:02}.egr.msu.edu -t "nc --recv-only localhost 51003"'.format(i)
          cores = check_output(cmd0, shell=True, universal_newlines=True)
          cmd0_1 = 'ssh compute{:02}.egr.msu.edu -t "nc --recv-only localhost 51001"'.format(i)
          load = check_output(cmd0_1, shell=True, universal_newlines=True)
          availability[i] = float(cores.split()[0]) - float(load.split()[0])  
        j = max(availability, key=availability.get)
        print('compute{:02} is least loaded'.format(j)) 
         
        cmd1 = 'ssh compute{:02}.egr.msu.edu -t "cd {};export DAMASK_NUM_THREADS={ncpus};/opt/soft/SIMULIA/2017/Commands/abaqus job={jobs} & sleep 8"'.format(j,dir_loc,ncpus=cpu_num,jobs=self.job_id)
        print(cmd1)
        run([cmd1],shell=True)
        time.sleep(5)
        while os.path.isfile('{}/{}.lck'.format(dir_loc,self.job_id)):
          time.sleep(100)
        if not os.path.isfile('{}/{}.sta'.format(dir_loc,self.job_id)):
          print("UMAT terminated before abaqus starts running")
        else:
          exit_status = open('{}/{}.sta'.format(dir_loc,self.job_id),'r').readlines()[-1].split()[4]
          if exit_status != "SUCCESSFULLY":
            print("Error in termination")
            print("+++++++++++++++++++++++++++++++++ failed fitness and points +++++++++++++++++++++++++++")
            file.write("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
            print("Error in termination")
            file.write('failed points that returned none {}\n'.format(x[:self.dimension]))
            error.append(np.nan)
        
          else:
            print("Normal Abaqus termination doing post results")                                         # Repeating the simulation with parameter values in ratio of last converged 
          
            cmd2 = 'abaqus python ./abaqus_postResults.py -p {} -f {}.odb '.format(dir_loc,self.job_id)
            run([cmd2],shell=True)
            print(cmd2)
            
            cmd3 = 'addCalculation -l "depth/nm","load/mN" -f "#z_displacement# * 1e3","#z_reactForce# * 1e-9" {}/Force_displacement_z.txt'.format(dir_loc)
            run([cmd3],shell=True)
            print(cmd3)
        
            cmd4 = 'filterTable --white="depth/nm","load/mN" < {}/Force_displacement_z.txt > {}/abaqus_result_load_disp.txt'.format(dir_loc,dir_loc)
            run([cmd4],shell=True)
            print(cmd4)
        
            cmd5 = 'chopTable.py -l "depth/nm" -v 50.0 {}/abaqus_result_load_disp.txt'.format(dir_loc)
            run([cmd5],shell=True)
            print(cmd5)
        
            cmd6 = 'InterpolateTables.py --expX "exp_depth/nm" --expY "exp_load/mN" --expfile {}/exp_load_disp_rt_{}.txt --simX "depth/nm" --simY "load/mN" --simfile {}/chopped_abaqus_result_load_disp.txt'\
                   .format(self.root,orientation_run+1,dir_loc)
            load_disp_error = check_output(cmd6, shell=True, universal_newlines=True)
            print(cmd6)
        
            cmd7 = 'abaqus_InterpolatedImage.py -c coordinate --rescale 1.0 --incomplete --stop_depth 50.0 --final_depth 423.0 --radius 1.0 --rotation 90.0 --grid 256 256 -p 10.0/255 {}/coord_after_retracting.txt'.format(dir_loc)
            run([cmd7],shell=True)
            print(cmd7)

            cmd8 = 'subset_image.py -d 256 256 --rmin 0 --rmax 50 -l 3_gridInterpolation --positive {}/out_coord_after_retracting.txt'\
                        .format(dir_loc)
            run([cmd8], shell=True)
            print(cmd8)

            cmd9 = 'abaqus_compare_pileup.py --multiplier 1e6 --flipud --sim {}/subset_sim.txt {}/exp_pileUp_rt_{}.txt '\
                        .format(dir_loc,self.root,orientation_run+1)
            pile_up_error = check_output(cmd9, shell=True, universal_newlines=True)
            print(cmd9)
        
            file.write("***\n load_disp error-- {} *** pile_up_error-- {} *** orientation --> {} \n*** \n".format(load_disp_error,pile_up_error,orientation))
            print("*** load_disp error--",load_disp_error," *** pile_up_error-- ",pile_up_error,"*** orientation -->",orientation,"***")
            error.append(0.5*float(load_disp_error) + 0.5*float(pile_up_error))
        
#------------------------------------------------------------------------------------------------- #
      if np.count_nonzero(np.isnan(error)) <= 0.5*len(error):
        avg_error = np.nanmean(error)
      else: 
        avg_error = np.mean(error)
      print("+++++++++++++++++++++++++++++++++ current fitness and points +++++++++++++++++++++++++++")
      print(float(avg_error))
      print(x[:self.dimension])
      print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
      file.write("\n +++++++++++++++++++++++++++++++++ current fitness and points +++++++++++++++++++++++++++\n")
      file.write("\n {}".format(float(avg_error)))
      file.write("\n {}".format(x[:self.dimension]))
      file.write("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
      self.curr_locations.append(np.append(x,float(avg_error)))
      self.locations[self.id(x)] = float(avg_error)
    
    return float(avg_error)
#    return np.linalg.norm(np.dot(x,x))


  def info_fitness(self,mat_file,job_id):
    self.mat_file  = mat_file
    self.job_id    = job_id
#    self.exp_file  = exp_file
  


#-------------------------------- main program starts here ----------------------------------------- #

parser = OptionParser()
parser.add_option(      '--mat',
                  dest = 'mat_file',
                  type = 'string', metavar = 'string',
                  help = ' material.config file (including suffix)')
parser.add_option(      '--job_id',
                  dest = 'job_id',
                  type = 'string', metavar = 'string',
                  help = ' job id for abaqus input file')
parser.add_option(      '--root',
                  dest = 'root',
                  type = 'string', metavar = 'string',
                  help = ' desired root of this process ')
parser.add_option('--restart',  action="store_true",
                  dest="restart",
                  help="restart optimization")
parser.add_option('-c','--concise',  action="store_true",
                  dest="concise",
                  help="concise outputs")
parser.add_option(      '--points',
                  dest = 'points_data',
                  type = 'string', metavar = 'string',
                  help = 'points for next generation ')
parser.add_option(      '--pBest',
                  dest = 'posPBest_data',
                  type = 'string', metavar = 'string',
                  help = 'posPBest for current generation ')
parser.add_option(      '--gBest',
                  dest = 'posGBest_data',
                  type = 'string', metavar = 'string',
                  help = 'posGBest for current generation ')
parser.add_option(      '--velo',
                  dest = 'velocities_data',
                  type = 'string', metavar = 'string',
                  help = 'velocities for current generation ')

#making the default values and let them show
parser.set_defaults(mat_file = 'cpTi_material.config',
                    job_id = 'Indentation_mesh1008',
                    concise = False,
                   )

(options,filenames) = parser.parse_args()

if not options.mat_file or not os.path.exists(options.mat_file):
  print("Suitable format material config (file containing parameters) is not supplied ")
if not os.path.exists('%s.inp'%options.job_id):
  parser.error('No job_file selected for abaqus')

options.root = os.path.dirname(os.path.realpath(__file__)) if options.root == None else options.root

tick = time.time()

if options.restart:

  table1 = damask.ASCIItable(name = options.points_data, buffered = False)
  table1.head_read()
  table1.data_readArray()
  
  table2 = damask.ASCIItable(name = options.posPBest_data, buffered = False)
  table2.head_read()
  table2.data_readArray()
  
  table3 = damask.ASCIItable(name = options.posGBest_data, buffered = False)
  table3.head_read()
  table3.data_readArray()
  
  table4 = damask.ASCIItable(name = options.velocities_data, buffered = False)
  table4.head_read()
  table4.data_readArray()
  

  theOptimizer = optimize( method = 'pso',
                           bounds = np.array([[10e6,2000e6],                     
                                              [10e6,2000e6],                    
                                              [10e6,4000e6]                   
                                             ]),
                           tolerance = 0.050,
                           root      = options.root,
                           concise_outputs = options.concise,
                           rigid     = True,
                           restart   = True,
                           points_rs         = table1.data[:,1:],
                           posPBest_rs       = table2.data[:,1:],                                           # particle history best pos
                           posGBest_rs       = table3.data,                                                 # global history best pos
                           velocities_rs     = table4.data[:,1:],
                           pOrder_rs         = table1.data[:,0],
                           localbest_rs      = np.array([0.51447839,0.58132769,0.53090445,0.53328172,0.51516947,\
                                                         0.579725,0.54207936,0.51486742,0.56395055,0.5327716,\
                                                         0.62838784,0.53119355,0.51694101,0.51506147,0.51631351,\
                                                         0.53126051,0.51738898,0.60547042,0.5555745,0.55399484]),        # array of history best fitnesses of each particle for current generation
                           gbest_rs          = 0.5144783867015663,                                                        # history global best fitness for current generation
                           )
else:
  theOptimizer = optimize( method = 'pso',
                           bounds = np.array([[10e6,2000e6],                     
                                              [10e6,2000e6],                    
                                              [10e6,4000e6]                   
                                             ]),
                           tolerance = 0.050,
                           root      = options.root,
                           concise_outputs = options.concise,
                           rigid     = True,
                           )

theOptimizer.info_fitness(options.mat_file,options.job_id)
theOptimizer.optimize(verbose = False)
tock = time.time()
print("Time for simulation",(tock - tick))
print("Cost {}".format(theOptimizer.cost()))
print("Best parameters and fitness {}".format(theOptimizer.best()))
with open("{}/output_{}.log".format(options.root,theOptimizer.method),'a') as file:
  file.write("\nTime for simulation {}".format(tock - tick))
  file.write("\nCost {}".format(theOptimizer.cost()))
  file.write("\nBest parameters and fitness {}".format(theOptimizer.best()))
