#!/usr/bin/env python
import numpy as np
import sys
from scipy import interpolate
filename =sys.argv[1]   
a=np.loadtxt(filename)
print "geladen"
timeslices=[1.e2,5.e2,1.e3,5.e3,1.e4,5.e4,1.e5,5.e5,1.e6]
for time in timeslices:
	timeused=a[np.argmin(np.absolute(a[:,0]-time)),0]
	print "timeused",timeused
	b=a[np.where(a[:,0]==timeused)]
	#print b
	Emax=np.amax(a[:,1])
	Emin=np.amin(a[:,1])
	ZahlE=150
	print "los"
	c=(b[:,1].flatten())[1:]
	d=(b[:,2].flatten())[1:]
	
	interpolationfkt=interpolate.InterpolatedUnivariateSpline(c,d,k=3) #0=t 1=E 2=p
	dPdE=interpolationfkt.__call__(c,nu=1)
	print dPdE
	print "check"
	dpdEE=dPdE*c
	Emax=c[np.argmax(dpdEE)]
	print dpdEE
	log=(np.vstack((c,c/Emax,dpdEE))).T

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
	
	


        
        
        
        
