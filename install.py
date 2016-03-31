#!/usr/bin/python

import lxml.etree as lxml
import argparse as ap
import os

pythonpaths=["lib"]
paths=["bin"]

parser=ap.ArgumentParser(description="Adds required paths to PYTHONPATH and PATH in the .bashrc file")

args=parser.parse_args()

path=os.getcwd()
home=os.environ['HOME']
print "Writing entries to {}/.bashrc".format(home)
with open(home+"/.bashrc","r") as f:
    lines=f.readlines()
with open(home+"/.bashrc","a") as f:

    for entry in paths:
        pathcommand="export PATH={}:$PATH\n".format(os.path.join(path,entry))
        if pathcommand not in lines:
            f.write(pathcommand)

    for entry in pythonpaths:
        pathcommand="export PYTHONPATH={}:$PYTHONPATH\n".format(os.path.join(path,entry))
        if pathcommand not in lines:
            f.write(pathcommand)


