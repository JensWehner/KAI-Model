#!/usr/bin/env python

from __tools__ import MyParser
from __geometry__ import geometry
import numpy as np
import sys
import copy
import subprocess as sp
import os
import errno
import shutil
import re
import numpy.linalg as lg
import time



parser=MyParser(description="A programm to model switching in ferroelectrics with interface to flexpde")
parser.add_argument("--option","-o",type=str,required=True,help="Optionfile")
args=parser.parse_args()

Geometry =geometry() 
Geometry.readoptions(args.option)
Geometry.setupGrainstructure()
Geometry.assignpoldirection()
Geometry.printstructuretofile()









