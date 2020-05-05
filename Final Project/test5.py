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

class LorenzAttractorRungeKutta(tk.Frame):
    DT            = 1e-3     # Differential interval
    STEP          = 100000   # Time step count
    

    
    def __init__(self,*args, **kwargs):
        
        global user_r_entry
        global user_sigma_entry
        global user_b_entry
        global user_r_entry
        global user_sigma_entry
        global user_b_entry
        
        #super(LorenzAttractorRungeKutta, self).__init__()
        
        self.X_0 = 0
        self.Y_0 = 0
        self.Z_0 = 0
        
        self.root = Tk()
        self.root.grid()
        self.root.title("Lorenz Simulation")
        self.root.geometry = ("500x500")
        self.root.resizable = (False, False)
        
        self.res = [[], [], []]

        self.root.graphframe = LabelFrame(master = self.root, text = "Phase Space")
        self.root.label1 = Label(self.root.graphframe, text="Graph", height = 3, width = 15)
        self.root.label1.grid(row=0, column=0, columnspan=10)
        
        #self.root.graphframe = LabelFrame(master = self.root, text = "Phase Space")
        #self.root.label1 = Label(self.root.graphframe, text="Graph", height = 3, width = 15).grid(row=5, column=0, columnspan=10)
        user_r_entry = DoubleVar()
        r_entry = Entry(self.root, textvariable = user_r_entry).grid(row = 1, column = 0)
        self.root.r_label = Label(self.root, text="r value", height = 1, width = 12).grid(row=2, column=0, columnspan=1)
                
        user_sigma_entry = DoubleVar()
        sigma_entry = Entry(self.root, textvariable = user_sigma_entry).grid(row = 1, column = 1)
        self.root.sigma_label = Label(self.root, text="sigma value", height = 1, width = 12).grid(row=2, column=1, columnspan=1)
        
        user_b_entry = DoubleVar()
        b_entry = Entry(self.root, textvariable = user_b_entry).grid(row = 1, column = 2)
        self.root.b_label = Label(self.root, text="b value", height = 1, width = 12).grid(row=2, column=2, columnspan=1)
        
        self.plot_button1 = Button (self.root, command = self.click1, height = 2, width = 8, text = "Run").grid(row = 3, column = 1)
        self.plot_button2 = Button (self.root, command = self.root.destroy, height = 2, width = 8, text = "Quit").grid(row = 4, column = 1)
        self.root.mainloop()
        
        
    def click1(self):
        global r_info
        global sigma_info
        global b_info
        global user_r_entry
        global user_sigma_entry
        global user_b_entry
        
        r_info = user_r_entry.get()
        sigma_info = user_sigma_entry.get()
        b_info = user_b_entry.get()
        
        print (r_info,sigma_info,b_info)
        
        tk.Frame.__init__(self, self.root)
        print (self.STEP)
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
            
            
            f = Figure(figsize=(10,10), dpi=100)
            a = f.add_subplot(111)
            a.plot(self.res[0], self.res[1])
    
    
            canvas = FigureCanvasTkAgg(f, self)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        except Exception as e:
            raise
            

    def __lorenz(self, xyz):
        global r_info
        global sigma_info
        global b_info
        
        p = sigma_info
        r = r_info
        b = b_info

        try:
            return [
                -p * xyz[0] + p * xyz[1],
                -xyz[0] * xyz[2] + r * xyz[0] - xyz[1],
                xyz[0] * xyz[1] - b * xyz[2]
            ]
        except Exception as e:
            raise
    def Quit(self):
        self.root.destroy()
        


        

if __name__ == '__main__':
    try:
        LorenzAttractorRungeKutta()
        
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)


        
    