import sys
import random
import matplotlib
import h5py as h5
import numpy as np
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.animation as animation

from matplotlib.ticker import MaxNLocator
from matplotlib import cm

from PyQt5 import QtCore
from PyQt5.QtWidgets import QSizePolicy, QWidget, QVBoxLayout, QHBoxLayout, QSlider, QLabel, QTextEdit, QPushButton


class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, dset=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.subplots_adjust(top=0.85, bottom=0.2, left=0.2)
        self.axes1 = self.fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes1.hold(False)

        self.dset = dset
        self.time = 0
        val = self.dset[...]
        self.compute_initial_figure(self.time, self.dset)

        #
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self, time, dset=None):
        pass

class MyStaticMplCanvas1D(MyMplCanvas):
    """Simple canvas with a sine plot."""
    def compute_initial_figure(self, time, dset=None):

        #> ==================================================================
        if dset != None:
            self.dset = dset
            self.time = time

        val = self.dset[...]
        self.ntime = val[:, 0, 0, 0].size

        val_1d = np.transpose(val[time, 0, 0, :])
        dx = 1.0
        nx = val_1d.size
        x = np.linspace(0,100.0,nx)

        self.ymin = val.min()
        self.ymax = 1.1 * val.max()

        title = self.dset.name.rsplit('/')[-1]

        self.line1, = self.axes1.plot(x, val_1d)
        self.axes1.axis([x.min(), x.max(), self.ymin, self.ymax])
        self.axes1.set_xlabel('x(mm)')
        self.axes1.set_ylabel(title)
        self.axes1.set_title(title)

    def update_figure(self, time):

    	#self.time = time
    	#if self.time == self.ntime:
    	#	self.time = 0

        val = self.dset[...]
        val_1d = np.transpose(val[time, 0, 0, :])

        dx = 1.0
        nx = val_1d.size
        x = np.linspace(0,100.0,nx)

        line1, = self.axes1.plot(x, val_1d)
        self.axes1.axis([x.min(),x.max(),self.ymin,self.ymax])

        #self.line1.set_ydata(val_1d)
        return [line1]

    def save_animation(self):
        self.ani = animation.FuncAnimation(self.fig, self.update_figure, self.ntime, blit=False, interval=500, repeat=False)
        self.ani.save('movie.gif', writer='imagemagick')



class MyStaticMplCanvas2D(MyMplCanvas):
    """Simple canvas with a sine plot."""
    def compute_initial_figure(self, time, dset=None):
        self.fig.clear()
        self.fig.subplots_adjust(top=0.85, bottom=0.2, left=0.2)
        self.axes1 = self.fig.add_subplot(111)
        #> ==================================================================
        if dset != None:
            self.dset = dset
            self.time = time

        val = self.dset[...]
        val_2d = np.transpose(val[time, 0, :, :])
        #print "aaa", val.shape, dset.shape

        nx = val_2d.shape[0]
        ny = val_2d.shape[1]

        dx=1.0
        dy=1.0

        self.y,self.x=np.mgrid[slice(dx,dx*(nx+1),dx), slice(dy,dy*(ny+0.5),dy)]


        levels=MaxNLocator(nbins=100).tick_values(val_2d.min(),val_2d.max())

        if(val_2d.min() == val_2d.max()):
        	ticks_val=np.linspace(val_2d.min(),val_2d.max()+1.0,5)
        else:
        	ticks_val=np.linspace(val_2d.min(),val_2d.max(),5)

        self.cf_temp1 = self.axes1.contourf(self.x,self.y,val_2d,cmap=cm.get_cmap('jet'),levels=levels)

        self.fig.colorbar(self.cf_temp1,ticks=ticks_val)



        self.axes1.axis([self.x.min(),self.x.max(),self.y.min(),self.y.max()])
        self.axes1.set_yticks(np.arange(0,self.y.max(),100))
        self.axes1.set_title('rho')
        self.axes1.set_xlabel('x(mm)')
        self.axes1.set_ylabel('y(mm)')

    def update_figure(self, time):
        self.fig.clear()
        self.fig.subplots_adjust(top=0.85, bottom=0.2, left=0.2)
        self.axes1 = self.fig.add_subplot(111)

        #> ==================================================================
        self.time = time

        val = self.dset[...]
        val_2d = np.transpose(val[time, 0, :, :])



        levels=MaxNLocator(nbins=100).tick_values(val_2d.min(),val_2d.max())

        if(val_2d.min() == val_2d.max()):
        	ticks_val=np.linspace(val_2d.min(),val_2d.max()+1.0,5)
        else:
        	ticks_val=np.linspace(val_2d.min(),val_2d.max(),5)

        self.cf_temp1 = self.axes1.contourf(self.x,self.y,val_2d,cmap=cm.get_cmap('jet'),levels=levels)
        self.fig.colorbar(self.cf_temp1,ticks=ticks_val)



        self.axes1.axis([self.x.min(),self.x.max(),self.y.min(),self.y.max()])
        self.axes1.set_yticks(np.arange(0,self.y.max(),100))
        self.axes1.set_title('rho')
        self.axes1.set_xlabel('x(mm)')
        self.axes1.set_ylabel('y(mm)')

        return self.im

    def save_animation(self):
        self.ani = animation.FuncAnimation(self.fig, self.update_figure, frames=self.ntime, blit=False, interval=500, repeat=False)
        self.ani.save('movie.gif', writer='imagemagick')



class PlotWidget(QWidget):
    """docstring for PlotWidget"""
    def __init__(self, parent = None, dataset = None):
        super(PlotWidget, self).__init__()

        self.parent = parent

        sizePolicy = QSizePolicy();
        sizePolicy.setHorizontalPolicy(QSizePolicy.Expanding);
        sizePolicy.setVerticalPolicy(QSizePolicy.Expanding);
        self.setSizePolicy(sizePolicy);

        self.plotVboxlayout = QVBoxLayout(self)

        #> the up part, including sc and slider
        self.sc = None
        if (dataset.shape[2] == 1):
            self.sc = MyStaticMplCanvas1D(self, dataset, width=10, height=4, dpi=100)
        elif (dataset.shape[1] == 1):
            self.sc = MyStaticMplCanvas2D(self, dataset, width=5, height=4, dpi=100)

        #> The slider
        self.sp_widget = QWidget(self)
        self.sp_layout = QHBoxLayout(self.sp_widget)

        self.label = QLabel("time")
        self.plainTextEdit = QTextEdit("0")
        self.plainTextEdit.setMaximumHeight(20)
        self.plainTextEdit.setMaximumWidth(100)
        self.sp = QSlider(QtCore.Qt.Horizontal)
        self.sp.setMinimum(0)
        self.sp.setMaximum(dataset[...].shape[0]-1)
        self.sp.setTickPosition(QSlider.TicksBelow)
        self.sp.setTickInterval(1)
        self.sp.valueChanged.connect(self.timeChange)
        self.sp_layout.addWidget(self.label)
        self.sp_layout.addWidget(self.plainTextEdit)
        self.sp_layout.addWidget(self.sp)

        self.button_saveAnimation = QPushButton("Save Animation")
        self.button_saveAnimation.clicked.connect(self.save_animation)


        self.plotVboxlayout.addWidget(self.sc)
        self.plotVboxlayout.addWidget(self.sp_widget)
        self.plotVboxlayout.addWidget(self.button_saveAnimation)

    def redraw(self, dset):
        t = 0
        self.sc.compute_initial_figure(t, dset)
        self.sc.draw()


    def timeChange(self):
        t = self.sp.value()
        self.sc.compute_initial_figure(t)
        self.sc.draw()
        self.plainTextEdit.setText(str(self.sp.value()))


    def save_animation(self):
        self.sc.save_animation()
        pass
