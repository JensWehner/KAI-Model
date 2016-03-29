#!/usr/bin/env python
import numpy as np
import sys
import copy
import subprocess as sp
import os
import errno
import shutil
import re
import argparse 
import lxml.etree as lxml
import numpy.linalg as lg
import time

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser=MyParser(description="A programm to run Flexpde and integrate up an equation of motion")
parser.add_argument("--option","-o",type=str,required=True,help="Optionfile")
args=parser.parse_args()




