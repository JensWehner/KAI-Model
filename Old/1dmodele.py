#!/usr/bin/env python
import numpy as np
import sys
import scipy.constants as const

ssteps =100000 #sys.argv[1]    #Zeit zu der hochintegriert werden soll in tau0
steps=int(ssteps)+1
Ea0=100.0e6 #V/m
tau0=1.e-10 #s
ps0=0.3 #c/m^2
beta0=1.0
d0=200.e-9 #m
d1=2.e-9
t=0.0
ef=600.
ei=60.0
deltatrange=[0.5]#,0.05,0.01,0.005,0.001]
Voltage=np.linspace(0.0,10.0,num=150)#[4.,3.,2.]
Eext=Voltage/(d0+d1*ef/ei)
#Eext=Voltage/d0
#tauext=np.exp((1./np.absolute(Eext/Ea0)))

def f(t,p,ps0,d1,d0,ef,ei,Ea0,beta0):
	Etot=(Eext-d1*p*ps0/(const.epsilon_0*ei*d0*(1.+d1/d0*ef/ei)))/Ea0
	a=(np.sign(Etot)-p)*(np.exp(-(1./np.absolute(Etot))))**beta0*(beta0*t**(beta0-1.0))
	return a
	
for faktordeltat in deltatrange:
	t=0.0
	p=-1.0*np.ones(np.shape(Voltage))
	logfile="logtextrk3"+str(faktordeltat)+"t.txt"
	for n in range(steps):
		Etot=(Eext-d1*p*ps0/(const.epsilon_0*ei*d0*(1.+d1/d0*ef/ei)))/Ea0
		#Etot=(Eext-d1*p*ps0/(const.epsilon_0*ei*d0))/Ea0
		tau=np.exp(1./np.absolute(np.max(Etot)))
		deltat=tau*faktordeltat
		#q=-ps0+2.*ps0*(1.-np.exp(-(t/tauext)**beta0))
		tarray=t*np.ones(np.shape(Voltage))
		if t==0.0:
			log=np.vstack((tarray,Eext,p)).T	
		else:
			D=np.vstack((tarray,Eext,p)).T
			log=np.vstack((log,D))
		#print deltap
		k1=f(t,p,ps0,d1,d0,ef,ei,Ea0,beta0)
		k2=f(t+0.5*deltat,p+0.5*deltat*k1,ps0,d1,d0,ef,ei,Ea0,beta0)
		k3=f(t+deltat,p-deltat*k1+2.0*deltat*k2,ps0,d1,d0,ef,ei,Ea0,beta0)
		deltap=deltat*(1./6*k1+2./3*k2+1./6*k3)
		p=p+deltap
		t=t+deltat
		#print t
		print n,"steps of",(steps-1)
		if t>10.e13:
			break

	np.savetxt('log.txt',log,'%0.5f',' ','&')

	f=open('log.txt','r')
	content=f.readlines()

	f.close()
	contents=content[0].split('&')
	del contents[-1]
	f=open(logfile,'w')
	for line in contents:
		f.write(line+'\n')
	f.close() 

        
        
        
        
