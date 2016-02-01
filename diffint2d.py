#!/usr/bin/env python
import numpy as np
import sys
from scipy import interpolate
minSpannung=12.
maxSpannung=15.0612244898
Spannungsschritte=26
numofpoints=100
t0=0.1
tmax=1.e12
d0=1.e-5 #m
Voltagerange=list(np.linspace(minSpannung,maxSpannung,num=Spannungsschritte))

for Voltage in Voltagerange:

	filename="log100K"+str(Voltage)+"V.txt"
	a=np.loadtxt(filename)
	E=Voltage/d0
	n0=np.log10(t0)
	nmax=np.log10(tmax)
	times=np.logspace(n0,nmax,numofpoints)
	texp=np.amax(a[:,0])
	times2=np.delete(times,np.where(times>texp))
	interpolationfkt=interpolate.UnivariateSpline(a[:,0],a[:,2],k=3,s=0) 
	intervalues=interpolationfkt.__call__(times2,nu=0)
	Elist=E*np.ones_like(intervalues)
	result=np.vstack((times2,Elist,intervalues)).T
	outfile="int"+str(np.around(Voltage,decimals=2))+"V.txt"
	np.savetxt('ausgabe.txt',result,'%0.20f',' ','&')
	f=open('ausgabe.txt','r')
	content=f.readlines()
	f.close()
	contents=content[0].split('&')
	del contents[-1]
	f=open(outfile,'w')
	for line in contents:
		f.write(line+'\n')
	f.close() 
	if Voltage==Voltagerange[0]:
		complet=result
	else:
		complet=np.vstack((complet,result))
	print "fertig"


print "geladen"
timeslices=np.logspace(2,9,8)
for time in timeslices:
	timeused=complet[np.argmin(np.absolute(complet[:,0]-time)),0]
	print "timeused",timeused
	b=complet[np.where(complet[:,0]==timeused)]
	print b
	Emax=np.amax(complet[:,1])
	Emin=np.amin(complet[:,1])
	#print "los"
	c=(b[:,1].flatten())[1:]
	d=(b[:,2].flatten())[1:]
	
	interpolationfkt=interpolate.UnivariateSpline(c,d,k=3,s=0.01) #0=t 1=E 2=p
	dPdE=interpolationfkt.__call__(c,nu=1)
	P=interpolationfkt.__call__(c,nu=0)
	#print dPdE
	#print "check"
	dpdEE=dPdE*c
	Emax=c[np.argmax(dpdEE)]
	#print dpdEE
	log=(np.vstack((c,c/Emax,d,P,dpdEE))).T

	np.savetxt('ausgabe.txt',log,'%0.20f',' ','&')
	logfile=str(timeused)+".txt"
	f=open('ausgabe.txt','r')
	content=f.readlines()

	f.close()
	contents=content[0].split('&')
	del contents[-1]
	f=open(logfile,'w')
	for line in contents:
		f.write(line+'\n')
	f.close() 
	
	


        
        
        
        
