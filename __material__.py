from __geometry__ import geometry
import numpy as np
from __tools__ import Parser
import lxml.etree as lxml
class material:

    def __init__(self):
        self.epsilon=None
        self.tau=None
        self.Ea=None
        self.Ps=None
        self.beta=None
        self.e=np.identity(2)


    def readoptions(self,optionfile):
        material = Parser(optionfile,"material")     
        e=np.array((material.find("epsilon").text).split(),dtype=float)
        self.e[0,0]=e[0]
        self.e[1,1]=e[1]
        self.e[1,0]=e[2]
        self.e[0,1]=e[2]
        self.tau=float(material.find("tau").text)
        self.Ea=float(material.find("Ea").text)
        self.Ps=float(material.find("Ea").text)
        self.beta=float(material.find("beta").text)

    def readfromxml(self,xmlfile):
        return

    def writetoxml(self,root):
        material=lxml.SubElement(root,"material")
        e=lxml.SubElement(material,"epsilon")
        e.text="{} {} {}".format(self.e[0,0],self.e[1,1],self.e[1,0])
        tau=lxml.SubElement(material,"tau")
        tau.text="{}".format(self.tau)
        Ea=lxml.SubElement(material,"Ea")
        Ea.text="{}".format(self.Ea)
        Ps=lxml.SubElement(material,"Ps")
        Ps.text="{}".format(self.Ps)
        beta=lxml.SubElement(material,"beta")
        beta.text="{}".format(self.beta)
    
    def setupGrainproperties(self,geometry):
        vector1=geometry.Pdirections
        vector2=np.array([-vector1[1],vector1[0]])
        transmatrix=np.array([vector1,vector2])
        self.epsilon=np.einsum('ihl,ij,jkl->hkl',transmatrix,self.e,transmatrix)
        print self.epsilon.shape
       
