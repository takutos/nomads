## "pymas.py"
## Expanded Realtime Simulation GUI for PyCX

## The following two lines should be placed at the beginning of your simulator code:
##
## import matplotlib
## matplotlib.use('TkAgg')

from Tkinter import *
import pycxsimulator
import pylab as PL

class GUI(pycxsimulator.GUI):

    def __init__(self,title='PyCX Simulator',interval=0, conVars=[], conSetting=[]):
        self.titleText = title
        self.timeInterval = interval
        self.control = conVars
        self.initGUI(conSetting)

    def initGUI(self, conSetting=[]):
        #create root window
        self.rootWindow = Tk()
        self.rootWindow.wm_title(self.titleText)
        self.rootWindow.protocol('WM_DELETE_WINDOW',self.quitGUI)
        self.rootWindow.geometry('+30+30')

        #buttons
        self.frameButton = Frame(self.rootWindow)
        self.frameButton.grid(row=0,column=0,padx=5,pady=5,sticky=EW)
        self.runPauseString = StringVar()
        self.runPauseString.set("Run")
        self.buttonRun = Button(self.frameButton,width=11,height=1,textvariable=self.runPauseString,command=self.runEvent)
        self.buttonRun.pack(side='left')
        self.buttonRun = Button(self.frameButton,width=11,height=1,text='Step Once',command=self.stepOnce)
        self.buttonRun.pack(side='left')
        self.buttonRun = Button(self.frameButton,width=11,height=1,text='Reset',command=self.resetModel)
        self.buttonRun.pack(side='left')

        #add control panel below the buttons
        if len(self.control) > 0:
            if len(self.control) == len(conSetting):
                self.controlVars = []
                self.controlPanel = Frame(self.rootWindow)
                self.controlPanel.grid(row=1,column=0,padx=5,pady=5,sticky=EW)
                self.controlItems = []
                for i, setting in enumerate(conSetting):
                    if setting[2] == 'int':
                        self.controlVars.append(IntVar())
                    else:
                        self.controlVars.append(DoubleVar())
                    self.controlVars[i].set(setting[0])
                    self.controlItems.append(Scale(self.controlPanel,variable=self.controlVars[i],label=setting[3],command=self.updateVars,from_=setting[4],to=setting[5],resolution=setting[6],tickinterval=(setting[5]-setting[4]),orient=HORIZONTAL,relief=SOLID))
                    self.controlItems[i].pack(padx=10,pady=5,fill='both')
        self.updateVars()

    def drawModel(self):
        if self.modelFigure == None or self.modelFigure.canvas.manager.window == None:
            self.modelFigure = PL.figure()
            PL.get_current_fig_manager().window.geometry("+420+30")
            self.modelFigure.canvas.set_window_title(self.titleText)
            PL.ion()
        self.modelDrawFunc()
        self.modelFigure.canvas.manager.window.update()

    def updateVars(self, *ignore):
        if len(self.control) > 0:
            for i, x in enumerate(self.controlVars):
                self.control[i] = x.get()

    #extended start function
    def start(self,func=[]):
        if len(func)==3:
            self.modelInitFunc = func[0]
            self.modelDrawFunc = func[1]
            self.modelStepFunc = func[2]
            self.modelInitFunc()
            self.drawModel()
        self.rootWindow.mainloop()
