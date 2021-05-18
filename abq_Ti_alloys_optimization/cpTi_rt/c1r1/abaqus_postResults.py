# -*- coding: utf-8 -*-

from odbAccess import *
from abaqusConstants import *
import os,sys,string
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id: abaqus_postResults.py 4416 2016-01-20 19:42:05Z $','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Add columns containing Hertzian contact depth and force.

""", version = scriptID)
# indicate parameters
parser.add_option('-p','--pathname',
                  dest = 'pathname',
                  type = 'string', metavar='string',
                  help = 'absolute pathname where odb file exists')  
parser.add_option('-f','--filename',
                  dest = 'filename',
                  type = 'string', metavar='string',
                  help = 'name of odb file')
                
                  
(options,filenames) = parser.parse_args()


save_path = '{}'.format(options.pathname)
completeName = os.path.join(save_path, options.filename)         


odb = openOdb(path='{}'.format(completeName))

beginstep= "Indenting"
endstep=   "Retracting"

firstStep=1                                                                                             
lastStep=2                                                                                                 
steptime=[]

firstFrame = odb.steps[beginstep].frames[0]
lastFrame = odb.steps[endstep].frames[-1]

coordinate_first = firstFrame.fieldOutputs['COORD']
coordinate_last = lastFrame.fieldOutputs['COORD']
indent_surf = odb.rootAssembly.instances['PART-1-1'].\
    nodeSets['SURFACE']                                                                        #fixed 'PART-1-1' for instance, in captal letter
indenter_tip = odb.rootAssembly.instances['PART-1-1'].\
    nodeSets['TIP_TOE']                                                                        #nodal reference of the indenter tip in captal letter 

indentsurf_coord_first = coordinate_first.getSubset(region=indent_surf)
indentsurf_coord_last = coordinate_last.getSubset(region=indent_surf)

outputFileName1 = os.path.join(save_path,'coord_before_indenting.txt')
outputFile1 = open(outputFileName1,'w')                                                  #output file names
outputFile1.write('1' +'\t' +'head' +'\n')
outputFile1.write('Node_label' +'\t' +'1_coordinate' +'\t' +'2_coordinate' +'\t' +'3_coordinate' +'\n')  #X,Y,Z coordinate respectively 

outputFileName2 = os.path.join(save_path,'coord_after_retracting.txt')
outputFile2 = open(outputFileName2,'w')   
outputFile2.write('1' +'\t' +'head' +'\n')
outputFile2.write('Node_label' +'\t' +'1_coordinate' +'\t' +'2_coordinate' +'\t' +'3_coordinate' +'\n')
for v in indentsurf_coord_first.values:
#    print 'Position = ', v.position,'Type = ',v.type
    Node_label1 = v.nodeLabel
    Xcoord1 = v.data[0]
    Ycoord1 = v.data[1]
    Zcoord1 = v.data[2]
    #print Node_label1
    #print Xcoord1
    #print Ycoord1
    #print Zcoord1
    outputFile1.write(str(Node_label1)+'\t'+ str(Xcoord1)+ '\t'+ str(Ycoord1)+'\t'+ str(Zcoord1)+'\n')
outputFile1.close() 
    
for d in indentsurf_coord_last.values:
    Node_label2 = d.nodeLabel
    Xcoord2 = d.data[0]
    Ycoord2 = d.data[1]
    Zcoord2 = d.data[2]
    #print Node_label2
    #print Xcoord2
    #print Ycoord2
    #print Zcoord2
    outputFile2.write(str(Node_label2)+'\t'+ str(Xcoord2)+ '\t'+ str(Ycoord2)+'\t'+ str(Zcoord2)+'\n')
outputFile2.close() 

outputFileName3 = os.path.join(save_path,'Force_displacement_z.txt')
outputFile3 = open(outputFileName3,'w')   
outputFile3.write('1' + '\t' + 'head' + '\n')
outputFile3.write('z_reactForce' + '\t' + 'z_displacement' + '\n')		# to append the file 
for stepKey in range(firstStep, lastStep+1):
    currentStepName=list(odb.steps.keys())[stepKey-1]
           #print 'Step: ', currentStepName
    currentStep=odb.steps[currentStepName]
           # Get number of frames in current step
    numberFrames=len(currentStep.frames)
    print(numberFrames)                                                                               #number of output for each step
    
  
    for cframe in range(numberFrames):
            
      # Outputs saved in current frame                    
        currentFrame = currentStep.frames[cframe]
    #     stepTime=currentFrame.frameValue
    #     steptime.append(stepTime)

        displacementList  = currentFrame.fieldOutputs['U']
        #print displacementList
        reaktionForceList = currentFrame.fieldOutputs['RT']
        coordList         = currentFrame.fieldOutputs['COORD']

        displacement = displacementList.getSubset(region=indenter_tip).values
        u1=displacement[0].data[0]
        u2=displacement[0].data[1]
        u3=displacement[0].data[2]
        u3=abs(u3)
        #print u3
    
        reaktionForce = reaktionForceList.getSubset(region=indenter_tip).values
        rf1=reaktionForce[0].data[0]
        rf2=reaktionForce[0].data[1]
        rf3=reaktionForce[0].data[2]
        rf3=abs(rf3)
        #print rf3
        outputFile3.write(str(rf3)+'\t'+str(u3) + '\n')
outputFile3.close()

odb.close()
