import sys
import traceback
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from random import randint
import math
from tkinter import *
import tkinter as tk
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
    
    #print (X_0, Y_0, Z_0) #Show randomised initial conditions

    
    def __init__(self):
        global r_entry
        global sigma_entry
        global b_entry
        
        super(LorenzAttractorRungeKutta, self).__init__()
        #self.widgets={}
        #self.grid(column=0,row=0)
        
        self.X_0, self.Y_0, self.Z_0 = 0,0,0
        
        self.root = Tk()
        self.root.grid()
        self.root.title("Lorenz Simulation")
        self.root.geometry = ("500x500")
        self.root.resizable = (False, False)
        
        self.res = [[], [], []]

        self.root.graphframe = LabelFrame(master = self.root, text = "Phase Space")
        self.root.label1 = Label(self.root.graphframe, text="Graph", height = 3, width = 15).grid(row=5, column=0, columnspan=10)
        
        #self.root.graphframe = LabelFrame(master = self.root, text = "Phase Space")
        #self.root.label1 = Label(self.root.graphframe, text="Graph", height = 3, width = 15).grid(row=5, column=0, columnspan=10)
        
        r_entry = Entry(self.root).grid(row = 1, column = 0)
        self.root.r_label = Label(self.root, text="r value", height = 1, width = 12).grid(row=2, column=0, columnspan=1)
                
        sigma_entry = Entry(self.root).grid(row = 1, column = 1)
        self.root.sigma_label = Label(self.root, text="sigma value", height = 1, width = 12).grid(row=2, column=1, columnspan=1)
        
        b_entry = Entry(self.root).grid(row = 1, column = 2)
        self.root.b_label = Label(self.root, text="b value", height = 1, width = 12).grid(row=2, column=2, columnspan=1)
        
        self.plot_button = Button (self.root, command = self.plot, height = 4, width = 2).grid(row = 3, column = 1)
        
        self.root.mainloop()

    def plot(self):
        
        tk.Frame.__init__(self, self.root)
        
        try:
            xyz = [self.X_0, self.Y_0, self.Z_0]
            print (xyz)
            for _ in range(self.STEP):
                k_0 = self.__lorenz(xyz)
                k_1 = self.__lorenz([x + k * self.DT / 2 for x, k in zip(xyz, k_0)])
                k_2 = self.__lorenz([x + k * self.DT / 2 for x, k in zip(xyz, k_1)])
                k_3 = self.__lorenz([x + k * self.DT for x, k in zip(xyz, k_2)])
                for i in range(3):
                    xyz[i] += (k_0[i] + 2 * k_1[i] + 2 * k_2[i] + k_3[i]) \
                            * self.DT / 6.0
                    self.res[i].append(xyz[i])
                    
            print (self.res[1])
            
            f = Figure(figsize=(10,10), dpi=100)
            a = f.add_subplot(111)
            a.plot(self.res[0], self.res[1])
    
    
            canvas = FigureCanvasTkAgg(f, self)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        except Exception as e:
            raise

    def __lorenz(self, xyz, p=10, r=24.4
                 , b=8/3.0):
        global r_entry
        global sigma_entry
        global b_entry
        
        
        print (sigma_entry,r_entry,b_entry)
        p = sigma_entry
        r = r_entry
        b = b_entry

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


        
    