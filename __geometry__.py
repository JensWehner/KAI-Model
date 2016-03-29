import numpy as np
import pyvoro
from __tools__ import Parser
import lxml.etree as lxml

class geometry:

    def __init__(self):
        self.grains=[]
        self.numberofgrains=0
        self.angle=0
        self.width=0
        self.safetyscale=0.95
        self.cells=None
        
    def readoptions(self,optionfile):
        geometry = Parser(optionfile,"geometry")     
        self.numberofgrains=int(geometry.find("numofgrains").text)
        self.angle=float(geometry.find("anglewindow").text)
        self.width=float(geometry.find("width").text)
        self.safetyscale=float(geometry.find("safetyscale").text)

    def setupGrainstructure(self):
        a0=np.random.rand(self.numberofgrains,2)
        a1=(a0*self.safetyscale)+(1.0-self.safetyscale)/2.0
        box=[[0,self.width],[0,self.width]]
        #see https://pypi.python.org/pypi/pyvoro/1.3.2
        self.cells=pyvoro.compute_2d_voronoi(a1,box,2.0)
        
        
          
        
        
