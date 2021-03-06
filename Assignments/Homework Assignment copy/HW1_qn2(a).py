#!/usr/bin/env python
# coding: utf-8

# ## Qn 2(a)
# ## Name: James Morillo
# ## Matric: U1740375H

# In[ ]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 11:10:42 2020

@author: jameselijah
"""
import numpy as np
import matplotlib.pyplot as plt
import math

"""Initialisation"""
# N = 10 particles inside a unit box
N = 10
# b = 0.1
b = 0.1
# random initial positions for particles
x = np.random.random(N)
y = np.random.random(N)
tmin = 10
dt = []

vx = [2.0*np.random.random() - 1.0 for i in range (N)]
vy = [2.0*np.random.random() - 1.0 for i in range (N)]
def exec(x, vx,vy, y, b, Lx = 1, Ly = 1):
    l_overlap = []
    # next initialize an empty array for the overlapping pairs
    for i in range(N-1):
        for j in range(i+1,N):
            #FIRST CHECK OVERLAP
            if (x[i]-x[j])**2 + (y[i]-y[j])**2 < 4*b**2:
                l_overlap.append([i,j])
                print(l_overlap)
                
                A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                D = B*B - 4.0*A*C
                temp = (-B - np.sqrt(D))/2*A
                '''
                ================================
                Check if length of overlap is more than 0: we then 
                move the particles apart
                ================================
                '''
                if len(l_overlap) > 0: 
                    f = 1.001
                    for p in l_overlap:
                        i = p[0]
                        j = p[1]
                        # move particles i and j
                        # find unit vector along line joining two centres
                        xnij = x[i] - x[j]
                        ynij = y[i] - y[j]
                        normnij = np.sqrt(xnij**2 + ynij**2)
                        xnij = xnij/normnij
                        ynij = ynij/normnij
                        # find centre of overlap
                        xc = 0.5*(x[i] + x[j])
                        yc = 0.5*(y[i] + y[j])
                        # move i and j
                        x[i] = xc + f*b*xnij
                        y[i] = yc + f*b*ynij
                        x[j] = xc - f*b*xnij
                        y[j] = yc - f*b*ynij
                        
                    '''
                    ================================
                    We then apply the collision conditions 
                    to obtain the time and append to dt
                    ================================
                    '''
                
                        
                    A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                    B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                    C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                    D = B*B - 4.0*A*C
                    temp = (-B - np.sqrt(D))/2*A
                    '''
                    ================================
                    Apply B.C
                    ================================
                    '''
                    if D < 0.0 and D == 0:
                        if x[i] <= 0.0:
                            x[i] = x[i] + Lx
                            A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                            B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                            C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                            D = B*B - 4.0*A*C
                            temp = (-B - np.sqrt(D))/2*A
                            if temp <= 0:
                                dt.append('imaginary')
                            else:
                                dt.append(temp)
                            
                        elif x[i] >= Lx:
                            x[i] = x[i] - Lx
                            A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                            B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                            C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                            D = B*B - 4.0*A*C
                            
                            temp = (-B - np.sqrt(D))/2*A
                            if temp <= 0:
                                dt.append('imaginary')
                            else:
                                dt.append(temp)
                            
                            
                        elif y[i] <= 0.0:
                            y[i] = y[i] + Ly
                            A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                            B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                            C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                            D = B*B - 4.0*A*C
                            temp = (-B - np.sqrt(D))/2*A
                            if temp <= 0:
                                dt.append('imaginary')
                            else:
                                dt.append(temp)
                        else:
                            y[i] = y[i] - Ly
                            A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                            B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                            C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                            D = B*B - 4.0*A*C
                            temp = (-B - np.sqrt(D))/2*A
                            if temp <= 0:
                                dt.append('imaginary')
                            else:
                                dt.append(temp)
                    if temp <= 0:
                        dt.append('imaginary')
                    else:
                        dt.append(temp)
                        

                # this elif is followed thru when 
                # the length of overlap is == 0
                # so then we have the BC properties
                elif D < 0.0 and D == 0:
                    if x[i] <= 0.0:
                        x[i] = x[i] + Lx
                        A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                        B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                        C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                        D = B*B - 4.0*A*C

                        temp = (-B - np.sqrt(D))/2*A
                        if temp <= 0:
                            dt.append('imaginary')
                        else:
                            dt.append(temp)
                        
                    elif x[i] >= Lx:
                        x[i] = x[i] - Lx
                        A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                        B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                        C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                        D = B*B - 4.0*A*C
                        
                        temp = (-B - np.sqrt(D))/2*A
                        if temp <= 0:
                            dt.append('imaginary')
                        else:
                            dt.append(temp)
                        
                        
                    elif y[i] <= 0.0:
                        y[i] = y[i] + Ly
                        A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                        B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                        C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                        D = B*B - 4.0*A*C
                        temp = (-B - np.sqrt(D))/2*A
                        if temp <= 0:
                            dt.append('imaginary')
                        else:
                            dt.append(temp)
                    else:
                        y[i] = y[i] - Ly
                        A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                        B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                        C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                        D = B*B - 4.0*A*C
                        temp = (-B - np.sqrt(D))/2*A
                        if temp <= 0:
                            dt.append('imaginary')
                        else:
                            dt.append(temp)
                    
                else:
                    dt.append(temp)
                
                
                '''
                To check the overlaps
                '''
                if len(l_overlap) == 0:
                    A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                    B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                    C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                    D = B*B - 4.0*A*C
                    
                    #boundary conditions
                    if D < 0.0 and D == 0:
                        if x[i] <= 0.0:
                            x[i] = x[i] + Lx
                            A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                            B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                            C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                            D = B*B - 4.0*A*C
                            '''
                            There should be a while loop somewhere 
                            but to maintan a consistent motion of the particles
                            in the box. 
                            
                            for e.g. while D <= recursively transport back the 
                            particle to the opposite wall. However, I cannot seem
                            to get the syntax.
                            '''
                            temp = (-B - np.sqrt(D))/2*A
                            if temp <= 0:
                                dt.append('imaginary')
                            else:
                                dt.append(temp)
                            
                        elif x[i] >= Lx:
                            x[i] = x[i] - Lx
                            A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                            B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                            C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                            D = B*B - 4.0*A*C
                            
                            temp = (-B - np.sqrt(D))/2*A
                            if temp <= 0:
                                dt.append('imaginary')
                            else:
                                dt.append(temp)
                            
                            
                        elif y[i] <= 0.0:
                            y[i] = y[i] + Ly
                            A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                            B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                            C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                            D = B*B - 4.0*A*C
                            temp = (-B - np.sqrt(D))/2*A
                            if temp <= 0:
                                dt.append('imaginary')
                            else:
                                dt.append(temp)
                        else:
                            y[i] = y[i] - Ly
                            A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                            B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                            C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                            D = B*B - 4.0*A*C
                            temp = (-B - np.sqrt(D))/2*A
                            if temp <= 0:
                                dt.append('imaginary')
                            else:
                                dt.append(temp)
                    else:
                        temp = (-B - np.sqrt(D))/2*A
                        if temp <= 0:
                            dt.append('imaginary')
                        else:
                            dt.append(temp)
            else:
                A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                D = B*B - 4.0*A*C
                
                if D < 0.0 and D == 0:
                    if x[i] <= 0.0:
                        x[i] = x[i] + Lx
                        A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                        B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                        C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                        D = B*B - 4.0*A*C

                        if temp < 0:
                            dt.append('imaginary')
                        else:
                            dt.append(temp)
                            
                    elif x[i] >= Lx:
                        x[i] = x[i] - Lx
                        A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                        B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                        C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                        D = B*B - 4.0*A*C
                        
                        if temp <= 0:
                            dt.append('imaginary')
                        else:
                            dt.append(temp)
                            
                        
                    elif y[i] <= 0.0:
                        y[i] = y[i] + Ly
                        A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                        B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                        C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                        D = B*B - 4.0*A*C
                        if temp <= 0:
                            dt.append('imaginary')
                        else:
                            dt.append(temp)
                            
                    else:
                        y[i] = y[i] - Ly
                        A = (vx[i] - vx[j])**2 + (vy[i] - vy[j])**2
                        B = 2.0*((x[i] - x[j])*(vx[i] - vx[j])+(y[i] - y[j])*(vy[i] - vy[j]))
                        C = (x[i] - x[j])**2 + (y[i] - y[j])**2 - 4.0*b*b
                        D = B*B - 4.0*A*C
                        if temp <= 0:
                            dt.append('imaginary')
                        else:
                            dt.append(temp)
                            
                else:
                    temp = (-B - np.sqrt(D))/2*A
                    if temp <= 0:
                        dt.append('imaginary')
                    else:
                        dt.append(temp)
                            
        
    dt_update = [x for x in dt if str(x) != 'nan'] #to remove nan
    print (dt_update)


while len(dt) < 5:
    exec(x, vx,vy, y, b, Lx = 1, Ly = 1)


# ## The imaginary D results would mean that the two 
# ## particles just grazed each other
