"""
The Poisson Process is the model we use for describing randomly occurring events and by itself, isn’t that useful. 

Using the Poisson Distribution, we can find the probability of a number of events in a time period or finding the probability of waiting some time until the next event.

The Poisson Distribution probability mass function gives the probability of observing k events in a time period given the length of the period and the average events per time:

\begin{equation}
P (k \ events \ in \ interval) = exp ( -\lambda )\frac{\lambda^k}{k!}
\end{equation}

Lambda is 
"""

import random
import numpy as np
import matplotlib.pyplot as plt
import math
import time


n = random.randint(1,10000) #total time
event = random.randint(1,100) #randomize an event but make sure event < n
time_period = random.randint(0,n )
_lamda = (event/n)*time_period
k_events_in_time_period = random.randint(0,event)


x = np.linspace(0,20,n-1)
r = np.random.uniform(0,1,n)

print ("The randomized Lambda = %f, n is randomized = %d, randomized number of event = %d" % (_lamda,n,event ))

time.sleep (1.5)

x_sum = []

#THE POISSON in the iterated loop
#IS iterated with i as time period changes as we run thorugh the simulation

for i in range(1,n):
    poisson = np.exp(-1*((event/n)*i))*((((event/n)*i)**k_events_in_time_period))/math.factorial(k_events_in_time_period)
    x_sum.append(poisson)
   
plt.plot(x ,x_sum, color = 'r')
plt.show()
