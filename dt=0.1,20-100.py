import numpy as np
import math
import matplotlib.pyplot as plt

##########################
#theoretical solution#######
##########################
k = 1000 #kg/m
m = 150 #kg
u0 = 0.03 #m
v0 = 0
wd = 2.58
wn = 2.58
w = 0.5
#print(wn,w)
ki = 0.02
c = 2*.02*m*wn
p0 = 750
Z = 1 - (w/wn)**2
Y = 2*ki*(w/wn)
#print(Z,Y)
D = (p0/k)*(Z) / ((Z**2)+(Y**2))
C = (p0/k)*(Y) / ((Z**2)+(Y**2))
#print(C,D)
A = u0 - D
B = ((A*ki*wn) - C*w) / wd
#print ("a is ",A,B)
u = []
v= []
uf = [] #free vibration
vf =[] 
t = np.linspace(0,20,201)

#print(t)

for i in t :
    ut = (np.exp(-ki*wn*t))*(A*np.cos(wd*t) + B*np.sin(wd*t))+ C*np.sin(w*t) + D*np.cos(w*t)
    vt = (-ki*wn*(np.exp(-ki*wn*t))*(A*np.cos(wd*t) + B*np.sin(wd*t))) +(np.exp(-ki*wn*t) * ((B*wd*np.cos(wd*t) - A*wd*np.sin(wd*t)))) + (C*w*np.cos(w*t) - D*w*np.sin(w*t))
u.append(ut)
v.append(vt)
unn = [j for sub in u for j in sub]
#####################################################################
u0_20 = (np.exp(-ki*wn*20))*(A*np.cos(wd*20) + B*np.sin(wd*20))+ C*np.sin(w*20) + D*np.cos(w*20)
v0_20 = (-ki*wn*(np.exp(-ki*wn*20))*(A*np.cos(wd*20) + B*np.sin(wd*20))) +(np.exp(-ki*wn*20) * ((B*wd*np.cos(wd*20) - A*wd*np.sin(wd*20)))) + (C*w*np.cos(w*20) - D*w*np.sin(w*20))
#########################################################################
tt = np.linspace(0,80,801)#free vibration


for i in tt :
     uft = (np.exp(-ki*wn*tt))*(u0_20*np.cos(wd*tt) + (((v0_20)+ki*wn*u0_20)/wd)*np.sin(wd*tt))
     vft = (((-ki*wn)*np.exp(-ki*wn*tt))*((u0_20*np.cos(wd*tt) + (((v0_20)+ki*wn*u0_20)/wd)*np.sin(wd*tt))))
     + ((-u0_20*np.sin(wd*tt)+(v0_20+(ki*wn*u0_20*np.cos(wd*tt))-wd*u0_20*np.sin(wd*tt))*np.exp(-ki*wn*tt)))
print(tt)

uf.append(uft)
uff =[j for sub in uf for j in sub]
print(uf)
#plt.subplot(211)
#plt.plot(t,unn,'r')
#plt.grid()
#plt.subplot(212)
#plt.plot(tt,uff,'c')
#plt.grid()
#plt.xlabel("time")
#plt.ylabel("displacment")

#plt.show()
###############################################################################################################
######## numerical method ###################
########## CENTRAL DIFFERENCE METHOD  #################
#####################################################################
dt = 0.1
p = [0 for x in range(801)]
#tn= np.linspace(0,20,201)
#print(tn)
#for i in tn :
 #   p0 = 750*(np.cos(0.5*i))
  #  p.append(p0)
#print(p)
Khat = (m/(dt*dt))+(c/(2*dt))
print(Khat)
a = (m/(dt*dt))-(c/(2*dt)) #01 = delta T
b = (k-((2*m)/(dt*dt)))
udotdot0 = (750 - v0*c - k*u0)/m #shetab dar lahze sefr
#u__1 = u0 - (deltat*v0) + ((((deltat)**2)/2)* udotdot0) #jabejayi -1
print(a,b,udotdot0)
##########################################
#u20 = -0.6794 u19=-0.7673
u1 = [0 for x in range(799)]
u1.insert(0,-0.6794)
u1.insert(1,-.5473417)
#print(len(u1))
#print(u1)
#print(p)
phat = [0 for x in range(800)]
phat.insert(0,-8252.49)

for i in range(1,801):
    phat[i] = (p[i]) - (a*u1[i-1]) - (b*u1[i])
    for n in range(0,800):
        u1[n+1] = (phat[n])/Khat
rounded_list = [ round(elem, 4) for elem in u1 ]
#print('phat',phat)
print(rounded_list)


#######################################
plt.plot(tt,uff,tt,rounded_list)
plt.grid()
plt.legend(["theoretical", "numerical"],loc='center left', bbox_to_anchor=(0.8, 1.05))
plt.xlabel("time(s)")
plt.ylabel("displacement (m)")
plt.title("dt = 0.1s-CENTRAL DIFFERENCE METHOD",loc = 'left')
#plt.savefig('dt=0.1,20-100.pdf')
plt.show()
#################################################
###################################################
################################################
################################################################
##########################################################################33






































