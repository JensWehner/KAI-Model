import numpy as np
import pyvoro
from __tools__ import Parser
import lxml.etree as lxml
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class geometry:

    def __init__(self):
        self.grains=[]
        self.numberofgrains=0
        self.angle=0
        self.width=0
        self.safetyscale=0.95
        self.cells=None
        self.angles=None
        
    def readoptions(self,optionfile):
        geometry = Parser(optionfile,"geometry")     
        self.numberofgrains=int(geometry.find("numofgrains").text)
        self.angle=float(geometry.find("anglewindow").text)*np.pi/180.0  #convert to radian
        self.width=float(geometry.find("width").text)
        self.safetyscale=float(geometry.find("safetyscale").text)

    def setupGrainstructure(self):
        a0=np.random.rand(self.numberofgrains,2)
        self.seeds=(a0*self.safetyscale)+(1.0-self.safetyscale)/2.0
        self.box=[[0,self.width],[0,self.width]]
        #see https://pypi.python.org/pypi/pyvoro/1.3.2
        self.cells=pyvoro.compute_2d_voronoi(self.seeds,self.box,2.0)
        self.Volumes=[]
        for cell in self.cells:
            self.Volumes.append(cell['volume'])
        self.printtoscreen()
       

    def printtoscreen(self):
        print self.cells
        print "Volume: ",np.sum(self.Volumes)

    def printstructuretofile(self):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        ax1.add_patch(patches.Rectangle((0.0, 0.0),self.width,self.width,fill=False))
        for cell in self.cells:
            coordinates=np.array(cell['vertices'])
            ax1.add_patch(patches.Polygon(coordinates,fill=False))    
        ax1.scatter(self.width*self.seeds[:,0],self.width*self.seeds[:,1])        
        scaling=self.width/np.sqrt(self.numberofgrains)
        if type(self.angles)!= None:
            for dx,dy,seed in zip(scaling*np.cos(self.angles),scaling*np.sin(self.angles),self.seeds):
               ax1.add_patch(patches.Arrow(self.width*seed[0],self.width*seed[1],dx,dy,width=0.25*scaling,color="black"))
                    
        fig1.savefig('structure.png', dpi=90, bbox_inches='tight')     

    def assignpoldirection(self):
        self.angles=(np.random.rand(self.numberofgrains,1).T)[0]*self.angle-self.angle/2.0+np.pi/2.0
        
        
          
        
        
