
import sys
import traceback
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from random import randint
import math


class LorenzAttractorRungeKutta:
    DT            = 1e-3     # Differential interval
    STEP          = 100000   # Time step count
    X_0, Y_0, Z_0 = 1,1,1  # Initial values of x, y, z
    print (X_0, Y_0, Z_0) #Show randomised initial conditions



    def __init__(self):
        self.res = [[], [], []]

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

            plt.xlabel("t")
            plt.ylabel("Z(t)")
            plt.title('Z(t) vs t')
            
            print (self.res[1])
            
            dimentioner_res_0 = [i for i in range(len(self.res[1]))]
            plt.plot(self.res[1], self.res[2]) #draw Z-t, X-t, Y-t to analyse if stable manifold etc.
            plt.show()
            
            self.__plot() #calls out function to plot in 3D

            
        except Exception as e:
            raise

    def __lorenz(self, xyz):
        p=10 
        r=24.4+math.sin((0.01)*self.DT)
        b=8/3.0
        try:
            return [
                -p * xyz[0] + p * xyz[1],
                -xyz[0] * xyz[2] + r * xyz[0] - xyz[1],
                xyz[0] * xyz[1] - b * xyz[2]
            ]
        except Exception as e:
            raise

    def __plot(self):
        """ Protting """
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
            
    def plot(self):
        plt.plot(self.res[0], self.res[1], colour = "blue", lw=1)
        plt.show()
        

if __name__ == '__main__':
    try:
        obj = LorenzAttractorRungeKutta()
        obj.exec()
        
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)
