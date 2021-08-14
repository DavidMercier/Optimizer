#!/usr/bin/env python3

import os,sys,damask
import os,sys,string
from optparse import OptionParser
import numpy as np
import math

scriptID   = str.replace('$Id: chopTable.py 153 2015-11-06 14:32:50Z chakra34 $','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
prints the size of the file (as ASCII table) .

""", version = scriptID)

parser.add_option('-l','--label',
                  dest = 'label',
                  type = 'string',
                  help = 'column label for chopping')

parser.add_option('-v','--value',
                  dest = 'value',
                  type = 'float',
                  help = 'retain only rows above this value for selected label')


(options, filenames) = parser.parse_args()
if options.label == None or options.value == None:
  parser.error('Requires both column label and value')

# --- loop over input files ------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              outname = os.path.join(os.path.dirname(name),"chopped_"+os.path.basename(name)),
                  			      buffered = False,
                             )
  except: continue

  damask.util.report(scriptName,name)
  table.head_read()
  col_index = table.label_index(options.label)
  table.head_write()

  outputAlive = True
  while outputAlive and table.data_read():
    outputAlive = table.data_write() and float(table.data[col_index]) != options.value

  table.close()
