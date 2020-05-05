import sys
import traceback
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from random import randint
import math
from tkinter import *
import matplotlib

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib import style

matplotlib.use("TkAgg")

"""
class window():
    
    def __init__(self):
        self.master = Tk()
        self.master.grid()
        
        self.master.title("Lorenz Simulation")
        self.master.geometry = ("500x500")
        self.master.resizable = (False, False)

        self.master.graphframe = LabelFrame(master = self.master, text = "Phase Space")
        self.master.graphframe.grid()
        
        self.master.label1 = Label(self.master.graphframe, text="Graph", height = 3, width = 15).grid(row=5, column=0, columnspan=10)
        
        self.slider = Scale(self.master, from_=0, to=42, orient=HORIZONTAL)
        self.slider.grid(row = 0, column = 0)
        
        self.master.label2 = Label(self.master, text="r value", height = 1, width = 12).grid(row=10, column=10, columnspan=1)
        
        self.master.mainloop()
"""      

class LorenzAttractorRungeKutta(tk.Frame):
    DT            = 1e-3     # Differential interval
    STEP          = 100000   # Time step count
    X_0, Y_0, Z_0 = 0,0,0
    #print (X_0, Y_0, Z_0) #Show randomised initial conditions

    
    def __init__(self):
        global r_entry
        global sigma_entry
        global b_entry
        
        self.master = Tk()
        self.master.grid()
        
        self.master.title("Lorenz Simulation")
        self.master.geometry = ("500x500")
        self.master.resizable = (False, False)
        self.res = [[], [], []]

        self.master.graphframe = LabelFrame(master = self.master, text = "Phase Space")
        self.master.label1 = Label(self.master.graphframe, text="Graph", height = 3, width = 15).grid(row=5, column=0, columnspan=10)
        
        r_entry = Entry(self.master).grid(row = 1, column = 0)
        self.master.r_label = Label(self.master, text="r value", height = 1, width = 12).grid(row=2, column=0, columnspan=1)
                
        sigma_entry = Entry(self.master).grid(row = 1, column = 1)
        self.master.sigma_label = Label(self.master, text="sigma value", height = 1, width = 12).grid(row=2, column=1, columnspan=1)
        
        b_entry = Entry(self.master).grid(row = 1, column = 2)
        self.master.b_label = Label(self.master, text="b value", height = 1, width = 12).grid(row=2, column=2, columnspan=1)
        
        self.plot_button = Button (self.master, command = self.exec(), height = 4, width = 2).grid(row = 3, column = 1)
        
        self.master.mainloop()

    def exec(self):
        """ Loranz attractor (Runge-Kutta method) execution """
        try:
            xyz = [self.X_0, self.Y_0, self.Z_0]
            for _ in range(self.STEP):
                k_0 = self.__lorenz(xyz)
                k_1 = self.__lorenz([x + k * self.DT / 2 for x, k in zip(xyz, k_0)])
                k_2 = self.__lorenz([x + k * self.DT / 2 for x, k in zip(xyz, k_1)])
                k_3 = self.__lorenz([x + k * self.DT for x, k in zip(xyz, k_2)])
                for i in range(3):
                    xyz[i] += (k_0[i] + 2 * k_1[i] + 2 * k_2[i] + k_3[i]) \
                            * self.DT / 6.0
                    self.res[i].append(xyz[i])
                    
            f = Figure(figsize=(5,5), dpi=100)
            a = f.add_subplot(111)
            a.plot(self.res[0], self.res[1])
    
    
            canvas = FigureCanvasTkAgg(f, self)
            canvas.show()
            canvas.get_tk_widget().grid(master = self.master.graphframe, row = 0, column = 0)

        except Exception as e:
            raise

    def __lorenz(self, xyz, p=10, r=24.4+math.sin((0.01)*DT)
                 , b=8/3.0):
        """ Lorenz equation
        :param  list xyz
        :param  float  p
        :param  float  r
        :param  float  b
        :return list xyz
        """
        try:
            return [
                -p * xyz[0] + p * xyz[1],
                -xyz[0] * xyz[2] + r * xyz[0] - xyz[1],
                xyz[0] * xyz[1] - b * xyz[2]
            ]
        except Exception as e:
            raise

        

if __name__ == '__main__':
    try:
        LorenzAttractorRungeKutta()
        
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)


        
    