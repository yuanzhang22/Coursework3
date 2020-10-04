# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 00:41:44 2020

@author: JoanneZ
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 18:31:17 2019

@author: Yuan Zhang
"""

import numpy as np
import scipy as sp
import scipy.io as sio
import matplotlib.pyplot as plt
import scipy.signal as sig
from sympy import diff
import scipy.stats as scs
import math
#load data
data=sio.loadmat('cw1_data2.mat')
dFonF=data['dFonF']
data1=sio.loadmat('cw1_data1.mat')
mouse_theta = data1 ['mouse_theta']
fs = data1['fs'].flatten()[0]    # frequency is the same as coursework 1
row =np.size(dFonF,0)#75
col =np.size(dFonF,1)#7420
t = np.linspace(0,(col-1)*1/fs,col) #build time serious matrix

## question 1

def findpeak(y):
    dy=np.diff(y)
    thres=3
    length_data=len(y)
    peak_width=65   #多次实验得到160可以去除
    l = []
    for i in range(1,length_data-1):
            if y[i-1] <= y[i] and y[i] >= y[i+1] and y[i]>thres: #找到极值且大于阈值，阈值的设置要解释
                l.append(i);
       
    s = []
    len_l=len(l)
    if len_l > 2 :  #if find two peak position is below the width threshold, I think they belongs to the same peaks and only keep the larger value peak's position.
       for j in range(0,len_l-1):   #当一个峰值内有两个极值点,去除在同一个峰内的极值点    
          if l[j+1]-l[j]>peak_width:
               if y[j] >= y[j+1]:
                  s.append(l[j])
               else:
                  s.append(l[j+1]);         
    else:
        s=l.copy()
    
    #find the beginning point
    min1=s.copy()     #use copy in case of change the value of s
    for i in range(0,len(s)):
        lag=s[i]    #从极大值点像向前面倒退，第一个做差值是非正数的点就是起始点
        while lag >= 0:
               lag =lag-1; 
               if dy[lag] < 0:
                  min1[i]=lag; #起点就是从后往前第一个小于零的数字
                  break
              
    #plt.plot(min1,y[min1],'y*')
          
    return s,min1
#s means the peak location
#min1 means the beginning of the peak location


#yo=dFonF[3]
#y=sig.savgol_filter(yo,299,6,mode='wrap')#Too much noise may affect the determination of peak position.
#    #Because Savitzky-Golay filter has the advantage of preserving the area, position and width of peaks, so I choose this filter and find peak position on the filtered curve.
#    #Apply a Savitzky-Golay filter to an array in order to smooth data and remove noise.     
#[s,min1]=findpeak(y)
#
##yo1=dFonF[3]
##y1=sig.savgol_filter(yo1,299,6,mode='wrap')#Too much noise may affect the determination of peak position.
##    #Because Savitzky-Golay filter has the advantage of preserving the area, position and width of peaks, so I choose this filter and find peak position on the filtered curve.
##    #Apply a Savitzky-Golay filter to an array in order to smooth data and remove noise.     
##[s1,min1]=findpeak(y1)
##
##yo2=dFonF[35]
##y2=sig.savgol_filter(yo2,299,6,mode='wrap') 
##[s2,min2]=findpeak(y2)
##
##yo3=dFonF[73]
##y3=sig.savgol_filter(yo3,299,6,mode='wrap')    
##[s3,min3]=findpeak(y3)
##
##plt.figure(2)
##plt.plot(yo1,label='Cell 3')
##plt.plot(yo2,label='Cell 35')
##plt.plot(yo3,label='Cell 73')
##plt.plot(s1,yo1[s1],'r*',label='The peak of a detected calcium transient')
##plt.plot(min1,yo1[min1],'y*',label='The beginning of a detected calcium transient')
##plt.plot(s2,yo2[s2],'r*')
##plt.plot(min2,yo2[min2],'y*')
##plt.plot(s3,yo3[s3],'r*')
##plt.plot(min3,yo3[min3],'y*')
##plt.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)#add label
##plt.title('Three selected cells with detected calcuim transient(original signal)')
##plt.xlabel('time/s (t interval is 1/f =0.0323 s) ')
##plt.ylabel('Calcuim transient value')
##
##plt.figure(1)
##plt.plot(y1,label='Cell 3')
##plt.plot(y2,label='Cell 35')
##plt.plot(y3,label='Cell 73')
##plt.plot(s1,y1[s1],'r*',label='The peak of a detected calcium transient(Filtered signal)')
##plt.plot(min1,y1[min1],'y*',label='The beginning of a detected calcium transient(Filtered signal)')
##plt.plot(s2,y2[s2],'r*')
##plt.plot(min2,y2[min2],'y*')
##plt.plot(s3,y3[s3],'r*')
##plt.plot(min3,y3[min3],'y*')
##plt.legend(loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)#add label
##plt.title('Three selected cells with detected calcuim transient(Filtered signal)')
##plt.xlabel('time/s (t interval is 1/f =0.0323 s) ')
##plt.ylabel('Calcuim transient value')
#
#
#
#
#
#
##plt.figure(2)
##plt.plot(y)
##plt.plot(s,y[s],'r*')
##plt.plot(min1,y[min1],'y*')
##plt.title('beginning of a spike')
##plt.xlabel('******')
##plt.ylabel('******')
#
##plt.figure(2)
##plt.plot(yo)
##plt.plot(s,yo[s],'r*')
##plt.plot(min1,yo[min1],'y*')
##plt.title('')
##plt.xlabel('******')
##plt.ylabel('******')
##
#### Considering the last value should not be defined as an spike.Becaues the following activity of the neuron is not showing.It may be a spike or the rising period of a spike.This remains unknown.
##
#### question 2
##
##
##turn = [] #turn means how many turns/loops of the circle
##for i in range(1,len(mouse_theta)):
##    delta=mouse_theta[i]-mouse_theta[i-1]
##    if (delta < -200):          # a loop exists the adjacent angles diff by 200 degree(not 360 because there are some inevitabel errors in angular value)
##     turn.append(i+1)
##     
##plt.figure()
##plt.plot(t,mouse_theta,label='angular position')
##plt.plot(t[turn],mouse_theta[turn],'r*')
##plt.title('all loops is found and show the end of loop with red asterisk')
##     
###calculate which loop
##loop = np.zeros(len(s));
##for j in range(0,len(s)):
##    a=0
##    for i in range(0,len(turn)):
##      if (s[j] > turn[i]):
##       loop[j] = loop[j]+1;   #loop
##
##plt.figure()
##plt.plot(mouse_theta[s],loop,'b.')
##plt.xlim(0, 360);
#
###写def
#def response_peak(mouse_theta,s):
##trial indicates which loop around the track the mouse is on
##loop position indicates the binned spatial angle of the mouse’s location, and response is 0 or 1 as above.
#    turn = [] #turn means how many turns/loops of the circle
#    for i in range(1,len(mouse_theta)):
#        delta=mouse_theta[i]-mouse_theta[i-1]
#        if (delta < -200):          # a loop exists the adjacent angles diff by 200 degree(not 360 because there are some inevitabel errors in angular value)
#         turn.append(i+1)
#    #calculate which loop
#    loop = np.zeros(len(s));
#    for j in range(0,len(s)):
#        for i in range(0,len(turn)):
#          if (s[j] > turn[i]):
#           loop[j] = loop[j]+1;   #loop
#    peak_angle=mouse_theta[s]
#    
#    pos=np.zeros(len(mouse_theta))#position
#    for i in range(len(mouse_theta)):
#        pos[s]=1
#        
#    return loop,peak_angle,pos
#
#[loop,peak_angle,pos]=response_peak(mouse_theta,s)
#plt.figure()
#plt.plot(peak_angle,loop,'b.')
#plt.xlim(0, 360);
#
#
##
##
### Question 3
#
#res=np.zeros(len(yo)) 
#res[s]=1 #res is the binary response variable which means the spike location is 1 and others is 0
#bin_number=20 #the number of bin is 20 here
#ang=np.round(mouse_theta/bin_number) #angle in 20 bins and rounding
#n=len(ang)
###Mutual information can be equivalently expressed as I(x,y)=Hx+Hy-H(x,y)
##Hx=scs.entropy(res)
##Hy=scs.entropy(angle_bin)
#
#px=np.array([1-len(s)/len(yo),len(s)/len(yo)])#px is probabilty of event x = 0 or 1.  0 means no spike,1 means spike.
#py=np.zeros(bin_number)
#pxy=np.zeros((2,bin_number))
#p=1/n
#for i in range(n):
#    py[int(ang[i])] += p  # py is probabilty of event y that which angle value occurs.
#    pxy[int(res[i]),int(ang[i])] +=p   #pxy means the joint probability of x and y event happen together.
## the pxy has 40 probability.with or without spikes and the corresponding angle value.
#Ixy=0 #Ixy is mutual information
#for i in range(2):
#    for j in range(bin_number):
#        if pxy[i,j] != 0:
#            #Ixy=pxy[i,j]*math.log(pxy[i,j]/(px[i]*p[j])[,10])
#            I = math.log(pxy[i,j]/(px[i]*py[j]),2)  # If the log base 2 is used, the units of mutual information are bits.
#            Ixy += I*pxy[i,j]
#print(Ixy)
#
#
####define Ixy
#def mutual_information(y,s,mouse_theta):
#    res=np.zeros(len(y)) 
#    res[s]=1 #res is the binary response variable which means the spike location is 1 and others is 0
#    bin_number=20 #the number of bin is 20 here
#    ang=np.round(mouse_theta/bin_number) #angle in 20 bins and rounding
#    n=len(ang)
#    ##Mutual information 
#    px=np.array([1-len(s)/len(yo),len(s)/len(yo)])#px is probabilty of event x = 0 or 1.  0 means no spike,1 means spike.
#    py=np.zeros(bin_number)
#    pxy=np.zeros((2,bin_number))
#    p=1/n
#    for i in range(n):
#        py[int(ang[i])] += p  # py is probabilty of event y that which angle value occurs.
#        pxy[int(res[i]),int(ang[i])] +=p   #pxy means the joint probability of x and y event happen together.
#    # the pxy has 40 probability.with or without spikes and the corresponding angle value.
#    Ixy=0 #Ixy is mutual information
#    for i in range(2):
#        for j in range(bin_number):
#            if pxy[i,j] != 0:
#                #Ixy=pxy[i,j]*math.log(pxy[i,j]/(px[i]*p[j])[,10])
#                I = math.log(pxy[i,j]/(px[i]*py[j]),2)  # If the log base 2 is used, the units of mutual information are bits.
#                Ixy += I*pxy[i,j]
#    return Ixy
#
#Ixy1=mutual_information(yo,s,mouse_theta)
#print(Ixy1)
###question 4
# entropy function
def entropy(x):
    idx = x>0
    H = -np.sum( x[idx] * np.log2(x[idx]) )
    return H


#####
H = [] #H is Uncertainty (entropy) of the spatial location variable
for i in range(0,75):
    yo=dFonF[i,:]
    y=sig.savgol_filter(yo,299,6,mode='wrap')
    [s,min1]=findpeak(y)
    res=np.zeros(len(yo)) 
    res[s]=1 #res is the binary response variable which means the spike location is 1 and others is 0
    bin_number=20 #the number of bin is 20 here
    ang=np.round(mouse_theta/bin_number) #angle in 20 bins and rounding
    n=len(ang)
    px=np.array([1-len(s)/len(yo),len(s)/len(yo)])#px is probabilty of event x = 0 or 1.  0 means no spike,1 means spike.
    py=np.zeros(bin_number)
    pxy=np.zeros((2,bin_number))
    p=1/n
    for i in range(n):
        py[int(ang[i])] += p  # py is probabilty of event y that which angle value occurs.
        pxy[int(res[i]),int(ang[i])] +=p   #pxy means the joint probability of x and y event happen together.
    
    #H is Uncertainty (entropy) of the spatial location variable
    Hres=entropy(px)
    H.append(Hres); 
Hmax=np.max(H)
print('Max entropy is ',Hmax)
#
#
### question 5
#N = 10
###cs,bins = np.histogram(s,N)
##plt.hist(flights['arr_delay'], color = 'blue', edgecolor = 'black',bins = int(180/5))
#
#
##If the log base 2 is used, the units of mutual information are bits.
