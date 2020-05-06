"""
Author: James Morillo
Matric: U1740375H
SPMS
Computational Physics Final Project 
Project title: "Lorenz Equations"
"""

"""
==================
Initialise Imports
=================
"""
import sys
import traceback
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from random import randint
import random
from tkinter import *
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk #allows importing of a interactive graph
from matplotlib.figure import Figure
from matplotlib import style
from mpl_toolkits.mplot3d import Axes3D #imports the 3D 
import time
import matplotlib 
matplotlib.use("TkAgg") #Backend of Matplotlib and we pull out TkAgg



class LorenzAttractorRungeKutta(tk.Frame):
    DT            = 1e-3     # Differential interval
    STEP          = 100000   # Time step count
    X_0, Y_0, Z_0 = random.random(),random.random(),random.random()
    
    
    def __init__(self,*args, **kwargs):
        
        """
        ====================
        Initiliasing global
        variables to be used
        in other functions
        ====================
        """
        
        global user_r_entry
        global user_sigma_entry
        global user_b_entry
        global user_r_entry
        global user_sigma_entry
        global user_b_entry
        
        """
        =====================
        Initialising the GUI
        window, adding title
        and initialising into
        a grid for easy grid
        management
        =====================
        """
        self.root = Tk()
        self.root.grid()
        self.root.title("Lorenz Simulation")
        self.root.geometry = ("700x700")
        self.root.resizable = (False, False)
        
        self.res = [[], [], []]
        
        user_r_entry = DoubleVar()
        r_entry = Entry(self.root, textvariable = user_r_entry).grid(row = 1, column = 0)
        self.root.r_label = Label(self.root, text="r value", height = 1, width = 12).grid(row=2, column=0, columnspan=1)
                
        user_sigma_entry = DoubleVar()
        sigma_entry = Entry(self.root, textvariable = user_sigma_entry).grid(row = 1, column = 1)
        self.root.sigma_label = Label(self.root, text="sigma value", height = 1, width = 12).grid(row=2, column=1, columnspan=1)
        
        user_b_entry = DoubleVar()
        b_entry = Entry(self.root, textvariable = user_b_entry).grid(row = 1, column = 2)
        self.root.b_label = Label(self.root, text="b value", height = 1, width = 12).grid(row=2, column=2, columnspan=1)
        
        self.plot_button1 = Button (self.root, 
                                    command = self.click1, 
                                    height = 2, 
                                    width = 8, 
                                    text = "Run").grid(row = 3, column = 1)
        
        self.plot_button2 = Button (self.root, 
                                    command = self.Quit, 
                                    height = 2, 
                                    width = 8, 
                                    text = "Quit").grid(row = 4, column = 1)
        
        
        self.root.mainloop()
        
        
    def click1(self):
        
        """
        =================================================
        This function is the main driver when the button
        "Run" is clicked, where the main execution of the
        Lorenz estimation using RK4 method.
        =================================================
        """

        
        global r_info
        global sigma_info
        global b_info
        global user_r_entry
        global user_sigma_entry
        global user_b_entry
        
        
        canvas.delete("all")
        
        try:
            timeinit = time.process_time()
            r_info = user_r_entry.get() #This obtains the user input for r
            sigma_info = user_sigma_entry.get() #This obtains the user input for sigma
            b_info = user_b_entry.get()  #This obtains the user input for b
            
            tk.Frame.__init__(self, master = self.root) # this initialises the frame we used earlier
            
            print ("User entered r = %f, sigma = %f and b = %f." % (r_info, sigma_info, b_info))
            print ("The randomised initial values are \\ X0 = %f \\ Y0 = %f \\ Z0 = %f." %(self.X_0, self.Y_0, self.Z_0))
            
            xyz = [self.X_0, self.Y_0, self.Z_0] #initialises xyz in a list using the initial values
            
            for _ in range(self.STEP): #iterates up till the STEP size then applies RK4
                k_0 = self.__lorenz(xyz)
                k_1 = self.__lorenz([x + k * self.DT / 2 for x, k in zip(xyz, k_0)])
                k_2 = self.__lorenz([x + k * self.DT / 2 for x, k in zip(xyz, k_1)])
                k_3 = self.__lorenz([x + k * self.DT for x, k in zip(xyz, k_2)])
                for i in range(3):
                    xyz[i] += (k_0[i] + 2 * k_1[i] + 2 * k_2[i] + k_3[i]) \
                            * self.DT / 6.0
                    self.res[i].append(xyz[i])
                    
                    
            fig = Figure() #Produces a 3D graph
            ax = Axes3D(fig)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
            
            ax.plot(self.res[0], self.res[1], self.res[2], color="red", lw=1) 
    
            canvas = FigureCanvasTkAgg(fig, master = self.root) #this embeds the graph into the graphical UI 
            canvas.draw()
            canvas.get_tk_widget().grid(row = 0, column = 1)
            
            timeend = time.process_time()
            timer =  timeend - timeinit
            print ("Running time: %f. seconds " % timer)


        except Exception as e:
            raise
            

    def __lorenz(self, xyz):
        global r_info
        global sigma_info
        global b_info
        
        p = sigma_info
        r = r_info
        b = b_info
        
        return [
                -p * xyz[0] + p * xyz[1], 
                -xyz[0] * xyz[2] + r * xyz[0] - xyz[1], 
                xyz[0] * xyz[1] - b * xyz[2]
                ]

            
    def Quit(self):
        print ("Program Quitting..")
        self.root.quit()


if __name__ == '__main__':
    try:
        LorenzAttractorRungeKutta()
        
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)


        
    