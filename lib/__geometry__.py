import numpy as np
import pyvoro
from __flextools__ import Parser
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
        self.volume=0
        self.boundarylayer=0.003
        
    def readoptions(self,optionfile):
        geometry = Parser(optionfile,"geometry")     
        self.numberofgrains=int(geometry.find("numofgrains").text)
        self.angle=float(geometry.find("anglewindow").text)*np.pi/180.0  #convert to radian
        self.width=float(geometry.find("width").text)
        self.safetyscale=float(geometry.find("safetyscale").text)
        self.boundarylayer=float(geometry.find("boundarylayer").text)
        self.volume=self.width**2

    def createSeed(self):
        a0=np.random.rand(self.numberofgrains,2)
        self.seeds=self.width*((a0*self.safetyscale)+(1.0-self.safetyscale)/2.0)
        self.box=np.array([[0,self.width],[self.width,0]])

    def setupGrainstructure(self):
        #see http://pypi.python.org/pypi/pyvoro/1.3.2
        self.cells=pyvoro.compute_2d_voronoi(self.seeds,self.box,2.0)
        for cell in self.cells:
            coordinates=cell['vertices']
            for coordinate in coordinates:
                #print "before",coordinate
                if coordinate[0]<self.boundarylayer*self.width:
                   
                    coordinate[0]=self.boundarylayer*self.width
                    
                elif coordinate[0]>(1-self.boundarylayer)*self.width:
                    coordinate[0]=(1-self.boundarylayer)*self.width
                #print "after",coordinate
        self.setupVolfracs()
        #self.printtoscreen()

    def setupVolfracs(self):
        Volumes=[]
        for cell in self.cells:
            Volumes.append(cell['volume'])
        Volumes=np.array(Volumes)
        self.volfrac=Volumes/np.sum(Volumes)

    def getVolumes(self):
        
        return self.volfrac*self.volume
       

    def printtoscreen(self):
        print self.cells
        print "Volume: ",np.sum(self.volfrac)*self.volume

    def printstructuretopng(self):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        ax1.add_patch(patches.Rectangle((0.0, 0.0),self.width,self.width,fill=False))
        for cell in self.cells:
            coordinates=np.array(cell['vertices'])
            ax1.add_patch(patches.Polygon(coordinates,fill=False))    
        ax1.scatter(self.width*self.seeds[:,0],self.width*self.seeds[:,1])        
        scaling=self.width/np.sqrt(self.numberofgrains)
        if type(self.angles)!= None:
            for dr,seed in zip(scaling*self.Pdirections.T,self.seeds):
               ax1.add_patch(patches.Arrow(self.width*seed[0],self.width*seed[1],dr[0],dr[1],width=0.25*scaling,color="black"))
                    
        fig1.savefig('structure.png', dpi=90, bbox_inches='tight')     

    def assignpoldirection(self):
        angles=(np.random.rand(self.numberofgrains,1).T)[0]*self.angle-self.angle/2.0+np.pi/2.0
        self.Pdirections=np.array([np.cos(angles),np.sin(angles)])

    def writetoxml(self,root):
        geo=lxml.SubElement(root,"geometry")
        box=lxml.SubElement(geo,"box")
        for vector in self.box.T:
            lxml.SubElement(box,"vector",x="{:1.5f}".format(vector[0]),y="{:1.5f}".format(vector[1]))
        grains=lxml.SubElement(geo,"grains")
        for i,(seed,pdirection) in enumerate(zip(self.seeds,self.Pdirections.T)):
            grain=lxml.SubElement(grains,"grain",id="{:d}".format(i+1),x="{:1.5f}".format(seed[0]),y="{:1.5f}".format(seed[1]),px="{:1.5f}".format(pdirection[0]),py="{:1.5f}".format(pdirection[1]))
            


    def readfromxml(self,xmlfile):
        return
    
        
        
          
        
        
