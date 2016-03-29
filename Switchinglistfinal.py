#!/usr/bin/env python
import sys
import os
import numpy as np
from scipy.spatial import Delaunay
import subprocess 
import shutil
import time
import scipy.constants as const
print "Start"
inputdat = sys.argv[1]          

laufzeitanfang=time.time()

def inputauslesen(inputdat):
	f=open(inputdat,'r')
	content=f.readlines()
	f.close()
	print "los"
	numofpoints=(content[0].split()[0])
	steps=int(content[1].split()[0])+1
	t0=float(content[2].split()[0])
	maxSpannung=float(content[3].split()[0])
	minSpannung=float(content[4].split()[0])
	Spannungsschritte=int(content[5].split()[0])
	Voltagerange=list(np.linspace(minSpannung,maxSpannung,num=Spannungsschritte))
	#Voltagerange=np.array([16.89795918,  17.02040816,  17.14285714,  17.26530612,17.3877551 ,  17.51020408,  17.63265306,  17.75510204,17.87755102,  18.])
	delta0=float(content[6].split()[0])
	Safety=float(content[7].split()[0])
	eaa=float(content[8].split()[0])
	ecc=float(content[9].split()[0])
	Ea0=float(content[10].split()[0])
	tau0=float(content[11].split()[0])
	ps0=float(content[12].split()[0])
	beta0=float(content[13].split()[0])
	d0=float(content[14].split()[0])
	errlim0=int(content[15].split()[0])
	return numofpoints, steps , t0, Voltagerange, delta0, Safety, eaa, ecc, Ea0, tau0, ps0, beta0, d0, errlim0

def dot2(u, v):
	return u[0]*v[0] + u[1]*v[1]

def geometrieerstellen(numofpoints):
	nnumofpoints=int(numofpoints)
	a0=np.random.rand(nnumofpoints,2)
        a1=(a0*0.93)+(1.0-0.93)/2.0 
        b1=np.array([a1[:,0],-a1[:,1]]).T
        c1=np.array([-a1[:,0],a1[:,1]]).T
        d1=np.array([-a1[:,0]+2.0,a1[:,1]]).T
        e1=np.array([a1[:,0],-a1[:,1]+2.0]).T
        f1=np.vstack((a1,b1,c1,d1,e1))
	tri = Delaunay(f1)
	p = tri.points[tri.vertices]
	A = p[:,0,:].T
	B = p[:,1,:].T
	C = p[:,2,:].T
	a = A - C
	b = B - C

	def cross2(u, v, w):
                """u x (v x w)"""
    		return dot2(u, w)*v - dot2(u, v)*w

	def ncross2(u, v):
                """|| u x v ||^2"""
    	        return sq2(a)*sq2(b) - dot2(a,b)**2

	def sq2(u):
   		return dot2(u, u)

	cc = np.round(cross2(sq2(a) * b - sq2(b) * a, a, b) / (2*ncross2(a, b)) + C,decimals=12)
        
    def unique(a):
            #loescht doppelte Punkte aus den Polygonen
            order = np.lexsort(a.T)
            a = a[order]
            diff = np.diff(a, axis=0)
            ui = np.ones(len(a), 'bool')
            ui[1:] = (diff != 0).any(axis=1) 
            return a[ui]

	neighborlist=[]
	for i in range(np.shape(f1)[0]):

		a=np.where(tri.vertices==i)[0] 
		b=unique(cc[:,a].T).T
		if b.shape[1] > 2 and all((b >= 0.0).flatten())==True and all((b <= 1.0).flatten())==True:
		        neighborlist.append(b)
			
	def winkel(P,b):
		c=b-P
		theta=1.-(c[0]/np.sqrt(np.square(c[0])+np.square(c[1])))
		d=b[:,np.argsort(theta)]
		return d
	
	
	def complexhull(arrayofpoints):
		c=arrayofpoints
		a=np.argmin(c[1])
		b=np.delete(c,a,1)
		P=(np.atleast_2d(c[:,a])).T
		d=winkel(P,b)
		return np.hstack((P,d))

	sortedlist=[]
	for i in range(len(neighborlist)):
                nx=np.where(neighborlist[i][0]==0,0.003,np.where(neighborlist[i][0]==1,0.997,neighborlist[i][0]))
		n=np.vstack((nx,neighborlist[i][1]))
		a=complexhull(n)
		sortedlist.append(a)
	return sortedlist, len(sortedlist)

def materialkennwerte(anzahlkorn,eaa,ecc):
	#materialkennwerte erstellen winkel (0,2pi)
	#alpha=(np.random.rand(len(sortedlist),1).T)[0]*np.pi #rotationswinkel
	#materialkennwerte erstellen winkel (pi/4,3/4pi)
	alpha=(np.random.rand(anzahlkorn,1).T)[0]*np.pi/2+np.pi/4 #rotationswinkel
	beta=alpha-np.pi/2
	Ps=np.vstack((np.cos(alpha),-np.sin(alpha)))		#polarisation
	orthoPs=np.vstack((-Ps[1],Ps[0]))	#einheitsvektorder senkrecht auf Ps steht
	epsilon=np.vstack((eaa*np.square(np.cos(beta))+ecc*np.square(np.sin(beta)),ecc*np.square(np.cos(beta))+eaa*np.square(np.sin(beta)),(ecc-eaa)*np.cos(beta)*np.sin(beta))) #eaa,ecc,eac
	return Ps, orthoPs, epsilon

def ausleseninput(datei,d0): #methode zum auslesen einer vorherigen input.pde
        f=open(datei,'r')
	content=f.readlines()
	f.close()
        eaa=[]
        ecc=[]
        eac=[]
        px=[]
        py=[]
        pointx=[]
        pointy=[]
        sortedlist=[]
        check=False
        for line in content:
                if "REGION 2" in line:
                        check=True
                if "eaa" in line and check==True:
                        eaa.append(float(line.split("=")[1]))
                elif "ecc" in line and check==True:
                        ecc.append(float(line.split("=")[1]))
                elif "eac" in line and check==True:
                        eac.append(float(line.split("=")[1]))
                elif "Ps=VECTOR" in line and check==True:
                        px.append(float(line.split("(")[1].split(",")[0]))
                        py.append(float(line.split("(")[1].split(",")[1].split(")")[0]))
                elif "START" in line and check==True:
                        pointx.append(float(line.split("(")[1].split(",")[0]))
                        pointy.append(float(line.split("(")[1].split(",")[1].split(")")[0]))
                elif "LINE TO" in line and "TO CLOSE" not in line and check==True:
                        pointx.append(float(line.split("(")[1].split(",")[0]))
                        pointy.append(float(line.split("(")[1].split(",")[1].split(")")[0]))
                elif "LINE TO CLOSE" in line and check==True:
                        a=np.array([pointx,pointy])/d0
                        sortedlist.append(a)
                        pointx=[]
                        pointy=[]
                        
        p=np.array([px,py])
        ps=p/np.sqrt(dot2(p,p))
        orthops=np.vstack((-ps[1],ps[0]))
        epsilon=np.array([eaa,ecc,eac])
        anzahlkorn=len(sortedlist)
	print anzahlkorn
        return sortedlist,anzahlkorn,ps,epsilon

def flexpdeeingabedat(intfile,anzahlkorn,sortedlist,eaa,ecc,epsilon,Ps,Voltage,meshfile,errlim,d0): #erzeugt input.pde, ausfuehren von input.pde erzeugt die datei Ewerte
	f=open(intfile,"w")
	f.write("TITLE \'input"+str(anzahlkorn)+"\'\nSELECT\nerrlim=1e-"+str(errlim)+"\nSELECT THREADS = 4\nCOORDINATES cartesian2\nVARIABLES\nphi\nDEFINITIONS\neaa="+str(eaa)+"\necc="+str(ecc)+"\neac=0\nPs=VECTOR(0,0)\nE_x=-dx(phi)\nE_y=-dy(phi)\n")
	f.write("e0="+str(const.epsilon_0)+"\n")
	for j in range(anzahlkorn): #schreiben der Integrale
		n=str(j+1)
		f.write("E"+n+"x= VOL_INTEGRAL(E_x, \'Korn"+n+"\')\nE"+n+"y= VOL_INTEGRAL(E_y, \'Korn"+n+"\')\n")
	f.write("transfermesh(\'"+meshfile+"\')\n")
	f.write("EQUATIONS\nphi: e0*(dx(eaa*dx(phi)+eac*dy(phi))+dy(ecc*dy(phi)+eac*dx(phi)))=0\n\nBOUNDARIES \nREGION 1    \'box\'\neaa="+str(np.sqrt(eaa*ecc))+"\necc="+str(np.sqrt(eaa*ecc))+"\nSTART(0,0)\nPERIODIC(x+"+str(d0)+",y)LINE TO (0,"+str(d0)+")\nVALUE(phi)=0 LINE TO ("+str(d0)+","+str(d0)+")\nnobc(phi)\nLINE TO ("+str(d0)+",0)\nVALUE(phi)="+str(Voltage)+" LINE TO CLOSE\n")
	for j in range(anzahlkorn):
		f.write("REGION "+str(j+2)+" \'Korn"+str(j+1)+"\'\n")
		f.write("Ps=VECTOR("+str(Ps[0,j])+","+str(Ps[1,j])+")\n") 
		f.write("eaa="+str(epsilon[0,j])+"\necc="+str(epsilon[1,j])+"\neac="+str(epsilon[2,j])+"\n")
		f.write("START("+str(d0*sortedlist[j][0,0])+","+str(d0*sortedlist[j][1,0])+")NATURAL(phi)=normal(+ps)\n")
		for k in range(sortedlist[j].shape[1]-1):
			f.write("LINE TO("+str(d0*sortedlist[j][0,k+1])+","+str(d0*sortedlist[j][1,k+1])+")NATURAL(phi)=normal(+ps)\n")
		f.write("LINE TO CLOSE\n")
	f.write("PLOTS\nSUMMARY EXPORT FILE \'Ewerte\'\nREPORT(\'avenergy\')\n")
	for j in range(anzahlkorn):
		n=str(j+1)
		f.write("REPORT(E"+n+"x)\nREPORT(E"+n+"y)\n")
	f.write("END")
	f.close()
	
def auslesenE(datei,Ea0):
	f=open(datei,'r')
	content=f.readlines()
	f.close()
	a=[]
	b=[]
	for line in content:
		if len(line) > 2:
			if line[1]=="E" and "x" in line:
				a.append(float(line.split()[1]))
			else:
				if line[1]=="E" and "y" in line:			
					b.append(float(line.split()[1]))		
	E=np.array([a,b])/Ea0
	return E

def auslesenV(datei):
	f=open(datei,'r')
	content=f.readlines()
	f.close()
	a=[]
	for line in content:
		if len(line) > 2:
			if line[1]=="V":
				a.append(float(line.split()[1]))			
	V=np.array(a)
	return V

def volumenauslesen(volfile,anzahlkorn,sortedlist,meshfile,eaa,ecc): #Auslesen der Volumen der Koerner sowie die erzeugung der Mesh.dat mit dem fertigen Netz
#Volumenauslesen
	f=open(volfile,"w")
	f.write("TITLE \'vol"+str(anzahlkorn)+"\'\nSELECT\nerrlim=1e-5\nSELECT THREADS = 4\nCOORDINATES cartesian2\nVARIABLES\nphi\nDEFINITIONS\neaa="+str(eaa)+"\necc="+str(ecc)+"\neac=0\nPs=VECTOR(0,0)")
	f.write("e0="+str(const.epsilon_0)+"\n")
	for j in range(anzahlkorn): #schreiben der Integrale
		n=str(j+1)
		f.write("V"+n+"= VOL_INTEGRAL(1, \'Korn"+n+"\')\n")
	f.write("EQUATIONS\nphi: e0*(dx(eaa*dx(phi)+eac*dy(phi))+dy(ecc*dy(phi)+eac*dx(phi)))=0\n\nBOUNDARIES \nREGION 1    \'box\'\neaa="+str(np.sqrt(eaa*ecc))+"\necc="+str(np.sqrt(eaa*ecc))+"\nSTART(0,0)\nPERIODIC(x+"+str(d0)+",y)LINE TO (0,"+str(d0)+")\nVALUE(phi)=0 LINE TO ("+str(d0)+","+str(d0)+")\nnobc(phi)\nLINE TO ("+str(d0)+",0)\nVALUE(phi)=0 LINE TO CLOSE\n")
	for j in range(anzahlkorn):
		f.write("REGION "+str(j+2)+" \'Korn"+str(j+1)+"\'\n")
		f.write("START("+str(d0*sortedlist[j][0,0])+","+str(d0*sortedlist[j][1,0])+")NATURAL(phi)=normal(+ps)\n")
		for k in range(sortedlist[j].shape[1]-1):
			f.write("LINE TO("+str(d0*sortedlist[j][0,k+1])+","+str(d0*sortedlist[j][1,k+1])+")NATURAL(phi)=normal(+ps)\n")
		f.write("LINE TO CLOSE\n")
	f.write("PLOTS\n SUMMARY EXPORT FILE \'Vol\'\nREPORT(\'listofvolumes\')\n")
	for j in range(len(sortedlist)):
		n=str(j+1)
		f.write("REPORT(V"+n+")\n")
	f.write("transfer() file=\""+meshfile+"\"")
	f.write("END")
	f.close()
	subprocess.call(["flexpde6n","-Q",volfile])
	V=auslesenV("Vol")
	Vges=V.sum()
	return V,Vges

def verteilungsfunktionEx0(Eraw,Vrel,anzahlkorn):
	Ex=Eraw[0]
	A,binedges=np.histogram(Ex,bins=int(np.sqrt(anzahlkorn)),weights=Vrel)
	C=A/A.sum()
	bincenterpoints=np.delete(binedges,-1)+np.diff(binedges)/2.
	Zex=np.vstack((bincenterpoints,C)).T
	return Zex,binedges

def verteilungsfunktionEy0(Eraw,Vrel,anzahlkorn):
	Ey=Eraw[1]
	A,binedges=np.histogram(Ey,bins=int(np.sqrt(anzahlkorn)),weights=Vrel)
	C=A/A.sum()
	bincenterpoints=np.delete(binedges,-1)+np.diff(binedges)/2.
	Zey=np.vstack((bincenterpoints,C)).T
	return Zey,binedges

def verteilungsfunktionEmag0(Eraw,Vrel,anzahlkorn):
	Emag=np.sqrt(dot2(Eraw,Eraw))
	A,binedges=np.histogram(Emag,bins=int(np.sqrt(anzahlkorn)),weights=Vrel)
	C=A/A.sum()
	bincenterpoints=np.delete(binedges,-1)+np.diff(binedges)/2.
	Zemag=np.vstack((bincenterpoints,C)).T
	return Zemag,binedges

def verteilungsfunktionEx(Eraw,Vrel,binedges):
	Ex=Eraw[0]
	A,bla=np.histogram(Ex,bins=binedges,weights=Vrel)
	C=A/A.sum()
	bincenterpoints=np.delete(binedges,-1)+np.diff(binedges)/2.
	Zex=np.vstack((bincenterpoints,C)).T
	return Zex

def verteilungsfunktionEy(Eraw,Vrel,binedges):
	Ey=Eraw[1]
	A,bla=np.histogram(Ey,bins=binedges,weights=Vrel)
	C=A/A.sum()
	bincenterpoints=np.delete(binedges,-1)+np.diff(binedges)/2.
	Zey=np.vstack((bincenterpoints,C)).T
	return Zey

def verteilungsfunktionEmag(Eraw,Vrel,binedges):
	Emag=np.sqrt(dot2(Eraw,Eraw))
	A,bla=np.histogram(Emag,bins=binedges,weights=Vrel)
	C=A/A.sum()
	bincenterpoints=np.delete(binedges,-1)+np.diff(binedges)/2.
	Zemag=np.vstack((bincenterpoints,C)).T
	return Zemag

def f1(pscalar,t,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim):
	flexpdeeingabedat(workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,-pscalar*Ps*ps0,Voltage,meshfile,errlim,d0)
	subprocess.call(["flexpde6n","-Q",workfile])
	Eraw=auslesenE("Ewerte",Ea0)
	#pnorm=p/np.sqrt(dot2(p,p))
	E=Eraw/V
	f=(np.sign(E[1])-pscalar)*(np.exp(-(1./np.absolute(dot2(E,Ps)))))#**beta0*(beta0*t**(beta0-1.0))
	#f=(pnorm*np.sign(dot2(E,p))-p)*(np.exp(-(1/np.absolute(dot2(E,p)))))**beta0*(beta0*t**(beta0-1.0))
	return f,E

def f2(pscalar,t,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim):
	flexpdeeingabedat(workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,-pscalar*Ps*ps0,Voltage,meshfile,errlim,d0)
	subprocess.call(["flexpde6n","-Q",workfile])
	Eraw=auslesenE("Ewerte",Ea0)
	#pnorm=p/np.sqrt(dot2(p,p))
	E=Eraw/V
	f=(np.sign(E[1])-pscalar)*(np.exp(-(1./np.absolute(dot2(E,Ps)))))#**beta0*(beta0*t**(beta0-1.0))
	#f=(pnorm*np.sign(dot2(E,p))-p)*(np.exp(-(1/np.absolute(dot2(E,p)))))**beta0*(beta0*t**(beta0-1.0))
	return f	


#Beginn des eigentlichen programms

numofpointsstring,steps,t0,Voltagerange,delta0,Safety,eaa,ecc,Ea0,tau0,ps0,beta0,d0,errlim0=inputauslesen(inputdat)
print numofpointsstring,steps

if ".pde" in numofpointsstring:
        Kornecken,anzahlkorn,Ps,epsilon=ausleseninput(numofpointsstring,d0)
	p0=Ps
  
else: 
	numofpoints=int(numofpointsstring)
        Kornecken,anzahlkorn=geometrieerstellen(numofpoints)
        Ps,orthoPs,epsilon=materialkennwerte(anzahlkorn,eaa,ecc)
        p0=Ps
	
print "av. grainsize", d0/np.sqrt(float(anzahlkorn)), "m"	
        
volfile="vol."+str(anzahlkorn)+".pde"
meshfile="mesh"+str(anzahlkorn)+".dat"
print "Start"
V,Vges=volumenauslesen(volfile,anzahlkorn,Kornecken,meshfile,eaa,ecc)
Vrel=V/Vges

print "Spannungen",Voltagerange
print "Anzahl Koerner",anzahlkorn

for Voltage in Voltagerange:
        
        print "Spannung",Voltage
        Kennung=str(anzahlkorn)+"K"+str(Voltage)+"V"
        intfile="input"+Kennung+".pde"
	errlim=errlim0
        flexpdeeingabedat(intfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,Ps*ps0,Voltage,meshfile,errlim,d0)
        workfile="temp"+Kennung+".pde"
        shutil.copyfile(intfile, workfile)
        t=t0
        #deltat=1.0
        logfile="log"+Kennung+".txt"
        zemag="zemag"+Kennung+".txt"
        zex="zex"+Kennung+".txt"
        zey="zey"+Kennung+".txt"
        ts=np.array([t])
        p=p0
	pscalar=np.sqrt(dot2(p,p))*np.sign(p[1])
	print pscalar
	
	deltat=100.*delta0

        for i in range(steps):
		errlim=errlim0
		k1stern,E=f1(pscalar,t,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim)
		k1=k1stern*deltat
		k2=f2(pscalar+0.2*k1,t+0.2*deltat,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim)*deltat
                k3=f2(pscalar+(0.075*k1+0.225*k2),t+0.3*deltat,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim)*deltat
		k4=f2(pscalar+(0.3*k1-0.9*k2+1.2*k3),t+0.6*deltat,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim)*deltat
		k5=f2(pscalar+(-11./54*k1+2.5*k2-70./27*k3+35./27*k4),t+deltat,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim)*deltat
		k6=f2(pscalar+(1631./55296*k1+175./512*k2+575./13824*k3+44275./110592*k4+253./4096*k5),t+7./8*deltat,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim)*deltat
		deltay=37./378*k1+250./621*k3+125./594*k4+512./1771*k6
		deltaystar=2825./27648*k1+18575./48384*k3+13525./55296*k4+277./14336*k5+0.25*k6
		delta1=np.amax(deltay-deltaystar)
		errlimcount=False
		print delta1
		while delta1>delta0:
			if errlimcount==True:
				errlim=errlim #errlim+1 
			deltat=Safety*deltat*(np.absolute(delta0/delta1))**0.25
			k1=f2(pscalar,t,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim)*deltat
			k2=f2(pscalar+0.2*k1,t+0.2*deltat,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim)*deltat
			k3=f2(pscalar+(0.075*k1+0.225*k2),t+0.3*deltat,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim)*deltat
			k4=f2(pscalar+(0.3*k1-0.9*k2+1.2*k3),t+0.6*deltat,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim)*deltat
			k5=f2(pscalar+(-11./54*k1+2.5*k2-70./27*k3+35./27*k4),t+deltat,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim)*deltat
			k6=f2(pscalar+(1631./55296*k1+175./512*k2+575./13824*k3+44275./110592*k4+253./4096*k5),t+7./8*deltat,workfile,anzahlkorn,Kornecken,eaa,ecc,epsilon,ps0,Ps,d0,Voltage,meshfile,beta0,V,errlim)*deltat
			deltay=37./378*k1+250./621*k3+125./594*k4+512./1771*k6
			deltaystar=2825./27648*k1+18575./48384*k3+13525./55296*k4+277./14336*k5+0.25*k6
			delta1=np.amax(deltay-deltaystar)
			errlimcount=True
			print delta1
		else:
			deltatnew=Safety*deltat*(np.absolute(delta0/delta1))**0.2
			if deltatnew>10*deltat:
				deltat=10*deltat
			else:
				deltat=deltatnew
		deltap=deltaystar
		
                Edepol=E #vorsicht Edepol kann nicht berechnet werden, da Eext nicht bekannt ist.
                depolmax=np.amax(np.sqrt(dot2(Edepol,Edepol)))
		p=-Ps*pscalar
		pxges=(sum(p[0]*V))/Vges
		pyges=(sum(p[1]*V))/Vges

                if t==t0:
	                Zex,binedgesx=verteilungsfunktionEx0(E,Vrel,anzahlkorn)
	                Zey,binedgesy=verteilungsfunktionEy0(E,Vrel,anzahlkorn)
	                Zemag,binedgesmag=verteilungsfunktionEmag0(E,Vrel,anzahlkorn)
	                log=np.array([t,pxges,pyges,depolmax])
	                logzex=np.hstack((t*np.ones((np.shape(Zex)[0],1)),Zex))
	                logzey=np.hstack((t*np.ones((np.shape(Zey)[0],1)),Zey))
	                logzemag=np.hstack((t*np.ones((np.shape(Zemag)[0],1)),Zemag))
                        a=0       
                else:
                        D=np.array([t,pxges,pyges,depolmax])
	                log=np.vstack((log,D))
                        if a==5: #(t-a)>=100: #
	                        Zex=verteilungsfunktionEx(E,Vrel,binedgesx)
	                        Zemag=verteilungsfunktionEmag(E,Vrel,binedgesmag)
	                        Zey=verteilungsfunktionEy(E,Vrel,binedgesy)
	                        logzex=np.vstack((logzex,np.hstack((t*np.ones((np.shape(Zex)[0],1)),Zex))))
	                        logzey=np.vstack((logzey,np.hstack((t*np.ones((np.shape(Zey)[0],1)),Zey))))
	                        logzemag=np.vstack((logzemag,np.hstack((t*np.ones((np.shape(Zemag)[0],1)),Zemag))))
				np.savetxt('log.txt',log,'%0.5f',' ','&')
				np.savetxt('logzinpol.txt',logzex,'%0.5f',' ','&')
				np.savetxt('logzsenkpol.txt',logzey,'%0.5f',' ','&')
				np.savetxt('logzemag.txt',logzemag,'%0.5f',' ','&')
                                a=0
                        else:
                                a=a+1
	

                print "Zeit",t,"Step "+str(i)+" of "+str(steps-1),"Voltage "+str(Voltagerange.index(Voltage)+1)+" of "+str(len(Voltagerange))
		pscalar=pscalar+deltap
		print pscalar
                t=t+deltat
		if t>1.e12 or ((p+Ps)>-0.01).all()==True:
			break
	
        #Erstellung der Log Dateien
        np.savetxt('log.txt',log,'%0.5f',' ','&')
	    np.savetxt('logzex.txt',logzex,'%0.5f',' ','&')
        np.savetxt('logzey.txt',logzey,'%0.5f',' ','&')
	    np.savetxt('logzemag.txt',logzemag,'%0.5f',' ','&')
	    f=open('log.txt','r')
        content=f.readlines()
        f.close()

        contents=content[0].split('&')
        del contents[-1]
        f=open(logfile,'w')
        for line in contents:
                f.write(line+'\n')
        f.close() 

        f=open('logzex.txt','r')
        content=f.readlines()
        f.close()

        contents=content[0].split('&')
        del contents[-1]
        f=open(zex,'w')
        a=float(contents[0].split()[0])
        for line in contents:
                if a==float(line.split()[0]):
                        f.write(line+'\n')
                else:
                        f.write('\n'+line+'\n')
                a=float(line.split()[0])        
        f.close() 

        f=open('logzey.txt','r')
        content=f.readlines()
        f.close()

        contents=content[0].split('&')
        del contents[-1]
        f=open(zey,'w')
        a=float(contents[0].split()[0])
        for line in contents:
                if a==float(line.split()[0]):
                        f.write(line+'\n')
                else:
                        f.write('\n'+line+'\n')
                a=float(line.split()[0])
        f.close() 

        f=open('logzemag.txt','r')
        content=f.readlines()
        f.close()

        contents=content[0].split('&')
        del contents[-1]
        f=open(zemag,'w')
        a=float(contents[0].split()[0])
        for line in contents:
                if a==float(line.split()[0]):
                        f.write(line+'\n')
                else:
                        f.write('\n'+line+'\n')
                a=float(line.split()[0])
        f.close() 

        os.remove("log.txt")
        os.remove('logzex.txt')
        os.remove('logzey.txt')
        os.remove('logzemag.txt')
        os.remove("Ewerte")

laufzeitende=time.time()
laufzeit=laufzeitende-laufzeitanfang
print "Runtime ",laufzeit,"s"

