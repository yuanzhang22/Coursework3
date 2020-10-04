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
yo=dFonF[3]

y=sig.savgol_filter(yo,299,6,mode='wrap')#Too much noise may affect the determination of peak position.
#Because Savitzky-Golay filter has the advantage of preserving the area, position and width of peaks, so I choose this filter and find peak position on the filtered curve.
#Apply a Savitzky-Golay filter to an array in order to smooth data and remove noise.
dy=np.diff(y)
max=np.max(y)
min=np.min(y)

#if dy[i]>0,dy[i+1]<0,y[i]>thres  
thres=3
length_data=len(y)
peak_width=65   #多次实验得到160可以去除
#s=find_peak(y,length_data,thres,peak_width)
###
#def find_peak(y,length_data,thres,peak_width):#定义函数找到极值
   # l = []
l = []
for i in range(1,length_data-1):
        if y[i-1] <= y[i] and y[i] >= y[i+1] and y[i]>thres: #找到极值且大于阈值，阈值的设置要解释
            l.append(i);
#        elif y[i]==y[i-1] and y[i]>thres:
#            l.append(i);      
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

#q=np.maximum(l[j],l[j-1]);    #补上最后一个点#不用补上当peak_width=60
#s.append(q);            #补上最后一个点
#   return s
plt.figure(1)
plt.plot(y)
plt.plot(s,y[s],'r*')
#find the beginning point
min1=s.copy()     #use copy in case of change the value of s
for i in range(0,len(s)):
 lag=s[i]    #从极大值点像向前面倒退，第一个做差值是非正数的点就是起始点
 while lag >= 0:
      lag =lag-1; 
      if dy[lag] < 0:
          min1[i]=lag; #起点就是从后往前第一个小于零的数字
          break

plt.plot(min1,y[min1],'y*')
plt.title('beginning of a spike')
plt.xlabel('******')
plt.ylabel('******')

plt.figure(2)
plt.plot(yo)
plt.plot(s,yo[s],'r*')
plt.plot(min1,yo[min1],'y*')
plt.title('')
plt.xlabel('******')
plt.ylabel('******')

## Considering the last value should not be defined as an spike.Becaues the following activity of the neuron is not showing.It may be a spike or the rising period of a spike.This remains unknown.

## question 2


turn = [] #turn means how many turns/loops of the circle
for i in range(1,len(mouse_theta)):
    delta=mouse_theta[i]-mouse_theta[i-1]
    if (delta < -200):          # a loop exists the adjacent angles diff by 200 degree(not 360 because there are some inevitabel errors in angular value)
     turn.append(i+1)
     
plt.figure()
plt.plot(t,mouse_theta,label='angular position')
plt.plot(t[turn],mouse_theta[turn],'r*')
plt.title('all loops is found and show with red asterisk')


## example单个数据计算第几个loop
#a=s[1]
#loop=1
#for i in range(0,len(turn)):
#    if (a > turn[i]):
#     loop = loop+1;
     
#calculate which loop
loop = np.zeros(len(s));
for j in range(0,len(s)):
    a=0
    for i in range(0,len(turn)):
      if (s[j] > turn[i]):
       loop[j] = loop[j]+1;   #loop

plt.figure()
plt.plot(mouse_theta[s],loop,'b.')
plt.xlim(0, 360);


## Question 3

res=np.zeros(len(yo)) 
res[s]=1 #res is the binary response variable which means the spike location is 1 and others is 0
bin_number=20 #the number of bin is 20 here
ang=np.round(mouse_theta/bin_number) #angle in 20 bins and rounding
n=len(ang)
##Mutual information can be equivalently expressed as I(x,y)=Hx+Hy-H(x,y)
#Hx=scs.entropy(res)
#Hy=scs.entropy(angle_bin)

px=np.array([1-len(s)/len(yo),len(s)/len(yo)])#px is probabilty of event x = 0 or 1.  0 means no spike,1 means spike.
py=np.zeros(bin_number)
pxy=np.zeros((2,bin_number))
p=1/n
for i in range(n):
    py[int(ang[i])] += p  # py is probabilty of event y that which angle value occurs.
    pxy[int(res[i]),int(ang[i])] +=p   #pxy means the joint probability of x and y event happen together.
# the pxy has 40 probability.with or without spikes and the corresponding angle value.
Ixy=0 #Ixy is mutual information
for i in range(2):
    for j in range(bin_number):
        if pxy[i,j] != 0:
            #Ixy=pxy[i,j]*math.log(pxy[i,j]/(px[i]*p[j])[,10])
            I = math.log(pxy[i,j]/(px[i]*py[j]),2)  # If the log base 2 is used, the units of mutual information are bits.
            Ixy += I*pxy[i,j]
print(Ixy)


##question 4

# define the findpeak function to find the location of peak and the beginning of the peak
def findpeak(y):
    # y is input signal
    dy=np.diff(y) #calculate the difference 
    thres=3  # thres means threshold of peak valve.The choice of numerical value is obtained after many attempts.
    length_data=len(y)
    peak_width=65   #peak_width is the threshold distance between peaks.The choice of numerical value is obtained after many attempts.
    l = []
    for i in range(1,length_data-1):
            if y[i-1] <= y[i] and y[i] >= y[i+1] and y[i]>thres: #find the maximum value and meet the conditions greater than the threshold
                l.append(i);
       
    s = []
    len_l=len(l)
    if len_l > 2 :  #if find two peak position is below the width threshold, I think they belongs to the same peaks and only keep the larger value peak's position.
       for j in range(0,len_l-1):  #whether the two maximum points are the same peak.
          if l[j+1]-l[j]>peak_width:   
               if y[j] >= y[j+1]:
                  s.append(l[j])
               else:
                  s.append(l[j+1]);         #if the distance between the peaks less than the peak_width threshold ,we think the two peaks belong to on calcium transient and only keep the larger peak location.
    else:
        s=l.copy()
    
    #find the beginning point
    min1=s.copy()     #use copy in case of change the value of s
    for i in range(0,len(s)):
        lag=s[i]    #Backward from the maximum point image, the first point where the difference is non-positive is the starting point
        while lag >= 0:
               lag =lag-1; 
               if dy[lag] < 0:
                  min1[i]=lag; ##The starting point is the first number from dy less than zero from back to front
                  break
              
    #plt.plot(min1,y[min1],'y*')
          
    return s,min1
# entropy function
def entropy(x):
    idx = x>0
    H = -np.sum( x[idx] * np.log2(x[idx]) )
    return H
#Uncertainty (entropy) of the spatial location variable
#Hres=entropy(px)
#print('Entropy of spatial location is ',Hres)
row =np.size(dFonF,0)#75
col =np.size(dFonF,1)#7420
Hres=[]
ang=np.round(mouse_theta/bin_number) #angle in 20 bins and rounding
p=1/col
#for i in range(0,75):
i=8
y=dFonF[i,:]
[s,min]=findpeak(y)
px=np.array([1-len(s)/len(yo),len(s)/len(yo)])#px is probabilty of event x = 0 or 1.  0 means no spike,1 means spike.
py=np.zeros(20)
pxy=np.zeros((2,20))
for i in range(col):
    py[int(ang[i])] += p  # py is probabilty of event y that which angle value occurs.
H=entropy(py) 
Hres.append(H);
print(Hres)
    
    
## question 5
N = 10
##cs,bins = np.histogram(s,N)
#plt.hist(flights['arr_delay'], color = 'blue', edgecolor = 'black',bins = int(180/5))



