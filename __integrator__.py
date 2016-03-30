from __tools__ import Parser
from __geometry__ import geometry
import scipy.constants as const
import subprocess 
import time
class integrator():

    def __init__(self,geometry,material,maxstep,voltage,t0):
        self.geometry=geometry
        self.material=material
        self.maxstep=maxstep
        self.t0=t0
        self.voltage=voltage
        self.createmesh=False
    
    def readoptions(self,optionfile):
        integrator = Parser(optionfile,"integrator")     
        self.errlim0=float(geometry.find("errlim0").text)
        self.delta0=float(geometry.find("delta0").text)
        self.threads=int(geometry.find("threads").text)
        self.safetyfactor=float(geometry.find("safetyfactor").text)
        self.flexpdeinputfile="flex.in"
        self.flexdemeshfile="flex.mesh"
        self.flexdeoutputfile="flex.out"

    def setOutput(self,root,outputfreq)
        self.root=root
        self.outputfreq=outputfreq

    def IntegrateDgl(self):
        P0=self.material.Ps*np.ones(self.geometry.numberofgrains)
        P=P0
        t=self.t0
        deltat=100.*self.delta0
        
        for step in xrange(maxstep):
            start=time.time()
            if step==0:
                self.createmesh=True 
            else:
                self.createmesh=False         
            deltaP,deltat,E=self.rk_variabletimestep(self.voltage,P,t,deltat)
            P=P+deltaP 
            t=t+deltat  
            stop=time.time() 
            print "{} of {}, calctime={}:t={} P={}".format(step+1,maxstep,start-stop,P*self.material.Ps*self.geometry.Pdirections)
            if (step%self.outputfreq)==0:
                self.appendOutput(t,P,E)
        return

    def rk46(self,voltage,P,t,deltat):
        k1,E=self.evaluatedP(voltage,P,t)
        k1=k1*deltat
        k2,ignored=self.evaluatedP(voltage,P+0.2*k1,t+0.2*deltat)
        k2=k2*deltat
        k3,ignored=self.evaluatedP(voltage,P+(0.075*k1+0.225*k2),t+0.3*deltat)
        k3=k3*deltat
        k4,ignored=self.evaluatedP(voltage,P+(0.3*k1-0.9*k2+1.2*k3),t+0.6*deltat)
        k4=k4*deltat
        k5,ignored=self.evaluatedP(voltage,P+(-11./54*k1+2.5*k2-70./27*k3+35./27*k4),t+deltat)
        k5=k5*deltat
        k6,ignored=self.evaluatedP(voltage,P+(1631./55296*k1+175./512*k2+575./13824*k3+44275./110592*k4+253./4096*k5),t+7./8*deltat)
        k6=k6*deltat
        deltay=37./378*k1+250./621*k3+125./594*k4+512./1771*k6
		deltaystar=2825./27648*k1+18575./48384*k3+13525./55296*k4+277./14336*k5+0.25*k6
        delta1=np.amax(deltay-deltaystar)
        return deltaystar,delta1,E


    def rk_variabletimestep(self,voltage,P,t,deltat):
        deltaystar,delta1,E=self.rk46(voltage,P,t,deltat)
        while delta1>self.delta0):
            deltat=Safety*deltat*(np.absolute(self.delta0/delta1))**0.25
            deltaystar,delta1,E=self.rk46(voltage,P,t,deltat)
        else:
			deltatnew=Safety*deltat*(np.absolute(delta0/delta1))**0.2
			if deltatnew>10*deltat:
				deltat=10*deltat
			else:
				deltat=deltatnew
        return deltaP,deltat,E
        
    def evaluatedP(self,voltage,P,t):
        E=self.evaluateE(voltage,P)
        beta0=self.materials.beta
        Ps=self.geometry.Pdirections
        dP=(np.sign(E[1])-P)*(np.exp(-(1./np.absolute(np.dot(E,Ps)))))#**beta0*(beta0*t**(beta0-1.0))
        return dP,E

    def evaluateE(self,voltage,P):
        self.writeFlexpdeinput(voltage,P,self.createmesh)
        subprocess.call(["flexpde6n","-Q",self.flexpdeinputfile])
        E=self.readFlexpdeoutput()
        return E

    def writeFlexpdeinput(self,voltage,P,createmesh):
        with open(self.flexpdeinputfile,"w") as f:
            f.write("TITLE \'input{}\'\n".format(self.geometry.numberofgrains))
            f.write("SELECT\nerrlim=1e-{}\n".format(self,errlim0))
            f.write("SELECT THREADS = {}\n".format(self.threads))
            f.write("COORDINATES cartesian2\n")
            f.write("VARIABLES\nphi\n")
            f.write("DEFINITIONS\neaa={}\necc={}\neac={}\nP=VECTOR(0,0)\nE_x=-dx(phi)\nE_y=-dy(phi)\n".format(self.material.eaa,self.material.ecc,self.material.eac)
	        f.write("e0={}\n".format(const.epsilon_0))
	        for j in xrange(1,self.geometry.numberofgrains+1):
                f.write("E{0}x= VOL_INTEGRAL(E_x, \'Grain{0}\')\n".format(j))
                f.write("E{0}y= VOL_INTEGRAL(E_y, \'Grain{0}\')\n".format(j))
            if not createmesh:
	            f.write("transfermesh(\'{}\')\n".format(self.flexdemeshfile))
	        f.write("EQUATIONS\nphi: e0*(dx(eaa*dx(phi)+eac*dy(phi))+dy(ecc*dy(phi)+eac*dx(phi)))=0\n")
            averagee=np.sqrt(self.material.eaa*self.material.ecc) #no idea why I used this
            f.write("\nBOUNDARIES \nREGION 1    \'box\'\neaa={0}\necc={0}\n".format(averagee)
            f.write("START(0,0)\nPERIODIC(x+{0},y)LINE TO (0,{0})\nVALUE(phi)=0 LINE TO ({0},{0})\nnobc(phi)\nLINE TO ({0},0)\nVALUE(phi)={1} LINE TO CLOSE\n".format(self.geometry.width,voltage))
            P=self.geometry.Pdirections*self.material.Ps*P
	        for j,cell in enumerate(self.geometry.cells):
		        f.write("REGION {} \'Grain {}\'\n".format(j+2,j+1))
		        f.write("P=VECTOR({},{})\n".format(P[0,j],P[1,j]) 
		        f.write("eaa={}\necc={}\neac={}\n".format(self.material.epsilon[0,j],self.material.epsilon[1,j],self.material.epsilon[2,j])
                vertices=cell['vertices']
		        f.write("START({},{})NATURAL(phi)=normal(+ps)\n".format(coordinates[0][0],coordinates[0][1]))
		        for i,vertex in enumerate(vertices):
                    if i==0:
                        f.write("START({},{})NATURAL(phi)=normal(+ps)\n".format(vertex[0],vertex[1]))
                    else:
                        f.write("LINE TO({},{})NATURAL(phi)=normal(+ps)\n".format(vertex[0],vertex[1]))
	            f.write("LINE TO CLOSE\n")
	        f.write("PLOTS\nSUMMARY EXPORT FILE \'{}\'\nREPORT(\'avenergy\')\n".format(self.flexdeoutputfile)
	        for j in xrange(1,self.geometry.numberofgrains+1):
	            f.write("REPORT(E{0}x)\nREPORT(E{0}y)\n".format(j))
            if createmesh:
                f.write("transfer() file=\"{}\"".format(self.flexdemeshfile)))
	        f.write("END")

    def readFlexpdeoutput(self):
        with open(self.flexdeoutputfile,'r') as f:
	        content=f.readlines()
	        a=[]
	        b=[]
	        for line in content:
		        if len(line) > 2:
			        if line[1]=="E" and "x" in line:
				        a.append(float(line.split()[1]))
			        elif line[1]=="E" and "y" in line:			
					    b.append(float(line.split()[1]))		
	        E=np.array([a,b])/self.material.Ea/self.geometry.getVolumes()
	    return E


