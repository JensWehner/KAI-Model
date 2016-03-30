#!/usr/bin/env python

from __tools__ import MyParser
from __geometry__ import geometry
import lxml.etree as lxml



parser=MyParser(description="A programm to model switching in ferroelectrics with interface to flexpde")
parser.add_argument("--option","-o",type=str,required=True,help="Optionfile")
parser.add_argument("--restart","-r",type=str,default=False,help="Outputfile to restart from")
args=parser.parse_args()


run=Parser(args.option,"run")     
voltages=np.array((run.find("voltages").text).split(),dtype=float)
outputfile=run.find("file").text)
maxsteps=int(run.find("maxsteps").text))
outputfreq=int(run.find("output").text)




Geometry=geometry() 
Geometry.readoptions(args.option)
if args.restart==False:
    Geometry.createSeed()
    Geometry.setupGrainstructure()
    Geometry.assignpoldirection()
else:
    Geometry.readfromxml(args.restart)
    Geometry.setupGrainstructure()

Geometry.printstructuretofile()
for voltage in voltages:
    outputfile="{}_{}V.xml".format(outputfile,voltage)
    Integrator=integrator(geometry,material,maxsteps,voltage,t0,outputfreq,root)
    Integrator.setOutput(root,outputfreq)








