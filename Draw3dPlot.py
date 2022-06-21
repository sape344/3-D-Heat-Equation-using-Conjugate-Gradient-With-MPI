from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pylab import *


class Draw3d:
    def __init__(self,length=int):
        self.lenght=length
        self.x_index = np.zeros(self.lenght ** 3)
        self.y_index = np.zeros(self.lenght ** 3)
        self.z_index = np.zeros(self.lenght ** 3)
        self.data=np.zeros(self.lenght ** 3)
        i=0
        for z in range( self.lenght):
            for y in range( self.lenght):
                for x in range( self.lenght):
                    self.x_index[i] = x
                    self.y_index[i] = y
                    self.z_index[i] = z
                    i+=1

      #  self.color_map.set_clim(0.0,100)


    def SetData(self,data):
        self.fig = plt.figure(figsize=(10, 10))
        self.ax = Axes3D(self.fig)
        # configuring colorbar
        self.color_map = cm.ScalarMappable(cmap='seismic')
        self.ax.set_title("3D Heatmap")
        self.ax.set_xlabel('X-axis')
        self.ax.set_ylabel('Y-axis')
        self.ax.set_zlabel('Z-axis')

        self.data=data.copy().reshape((self.lenght**3))
        self.color_map.set_array(self.data)


    def ExamplePlot(self):
        data = np.zeros((self.lenght, self.lenght, self.lenght))
        for z in range(self.lenght):
            for y in range(self.lenght):
                for x in range(self.lenght):
                    data[z, y, x] = z
        self.SetData(data)
        self.Plot()


    def Plot(self):
        img = self.ax.scatter(self.x_index, self.y_index, self.z_index, marker='s',
                         s=100, c=self.data, cmap='seismic')
        plt.colorbar(self.color_map,aspect=40,shrink=0.75)
        plt.show()



