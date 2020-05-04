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

class window():
    DT            = 1e-3     # Differential interval
    STEP          = 100000   # Time step count
    X_0, Y_0, Z_0 = 0,0,0
    
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
        
        self.plot_button = Button (self.master, command = plotter, height = 4, width = 2, state=DISABLED).grid(row = 3, column = 1)
        
        self.master.mainloop()

    def plotter(self, event = None):
        f = Figure(figsize=(5,5), dpi=100)
        a = f.add_subplot(111)
        a.plot(self.res[0], self.res[1])


        canvas = FigureCanvasTkAgg(f, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
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
            #plt.axes()
            #plt.xlim([0, 50])
            #plt.ylim([-20, 20])
            plt.xlabel("t")
            plt.ylabel("Z(t)")
            plt.title('Z(t) vs t')
            
            dimentioner_res_0 = [i for i in range(len(self.res[1]))]
            plt.plot(dimentioner_res_0, self.res[1]) #draw Z-t, X-t, Y-t to analyse if stable manifold etc.
            #plt.plot(self.res[2], self.res[0]) #self.res[0] is values corresponding to x, [1] to y, [2] to z
            
            self.__plot() #calls out function to plot in 3D

            
        except Exception as e:
            raise

    def __lorenz(self, xyz):
        global r_entry
        global sigma_entry
        global b_entry
        
        p= sigma_entry 
        r= r_entry 
        b= b_entry
        try:
            return [
                -p * xyz[0] + p * xyz[1],
                -xyz[0] * xyz[2] + r * xyz[0] - xyz[1],
                xyz[0] * xyz[1] - b * xyz[2]
            ]
        except Exception as e:
            raise

    def __plot(self, t):
        
        try:
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
            #ax.title ("XX")
            ax.plot(self.res[0], self.res[1], self.res[2], color="red", lw=1)     
            
        except Exception as e:
            raise
            


if __name__ == '__main__':
    try:
        obj = window()
        #obj.exec()
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)


        
    