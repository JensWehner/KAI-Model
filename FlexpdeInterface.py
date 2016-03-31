#!/usr/bin/env python
from __tools__ import Parser
from __tools__ import MyParser
from __geometry__ import geometry
from __integrator__ import integrator
from __material__ import material
import lxml.etree as lxml
import numpy as np
import datetime


parser=MyParser(description="A programm to model switching in ferroelectrics with interface to flexpde")
parser.add_argument("--option","-o",type=str,required=True,help="Optionfile")
parser.add_argument("--restart","-r",type=str,default=False,help="Outputfile to restart from")
args=parser.parse_args()


run=Parser(args.option,"run")     
voltages=np.array((run.find("voltages").text).split(),dtype=float)
outfile=run.find("outfile").text
maxsteps=int(run.find("maxsteps").text)
outputfreq=int(run.find("output").text)
t0=0


Geometry=geometry() 
Geometry.readoptions(args.option)
if args.restart==False:
    Geometry.createSeed()
    Geometry.setupGrainstructure()
    Geometry.assignpoldirection()
else:
    Geometry.readfromxml(args.restart)
    Geometry.setupGrainstructure()

Material=material()
Material.readoptions(args.option)
if args.restart==False:
    Material.readfromxml(args.restart)
Material.setupGrainproperties(Geometry)

Geometry.printstructuretopng()
for i,voltage in enumerate(voltages):
    print "Starting calculation for V={}".format(voltage)
    print "Voltage {} of {}".format(i+1,len(voltages))
    start=datetime.datetime.now()
    outputfile="{}_{}V.xml".format(outfile,voltage)
    print "Writing output to \'{}\'".format(outputfile)
    runlog = lxml.Element("runlog",start=start.ctime(),V="{:1.3f}".format(voltage),steps="{:d}".format(maxsteps))
    Geometry.writetoxml(runlog)
    Material.writetoxml(runlog)
    with open(outputfile, 'w') as f:
            f.write(lxml.tostring(runlog, pretty_print=True))
    Integrator=integrator(Geometry,Material,maxsteps,voltage,t0)
    Integrator.readoptions(args.option)
    Integrator.setOutput(runlog,outputfreq,outputfile)
    Integrator.IntegrateDgl(maxsteps)
    end=datetime.datetime.now()
    runlog.set("end",end.ctime())
    with open(outputfile, 'w') as f:
            f.write(lxml.tostring(runlog, pretty_print=True))
    print "Finished calculation for V={}".format(voltage)








