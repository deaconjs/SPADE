# python imports
import string
import math
import os.path
import pickle
import time
import random
import scipy.stats
import fpformat
# dependency imports
from Tkinter import *
sys.path.append(os.path.abspath('./Dependencies'))
import Pmw
from tkFileDialog import *
# internal imports
import parms
import sys
sys.path.append(os.path.abspath('./Tools/Math'))
import Distribution
sys.path.append(os.path.abspath('./Applications'))
import MolecularViewer
import MolecularSystem

class MiniViewer(Frame):
    def __init__(self, parent, width, ht):
        Frame.__init__(self, parent, borderwidth=0, width=width, height=ht)
        self.width = width
        self.height = ht
        self.parent = parent
        # create and pack a canvas, and draw the axes.
        print "building the viewer"
        self.canvas = Canvas(self, bg='white', relief=RAISED, bd=2, height=self.height, width=self.width)
        self.canvas.pack(expand=NO)
        # bind canvas clicks to the scroll function
        self.canvas.bind('<Button-1>', self.mini_scroll)
        # get some info for the panel representations
        self.panel_height = (self.height*1.0)/len(self.parent.msw.msp)
        self.panel_length = self.width
        counter = 1
        # draw the axes
        for panel in self.parent.msw.msp:
            print "creating %d"%(counter)
            self.canvas.create_line(10,self.panel_height*counter,self.panel_length-5, self.panel_height*counter, fill='black')
            print "drawing %s %s"%(self.panel_height*counter, self.panel_height*(counter-1)+5.0)
            self.canvas.create_line(10,self.panel_height*counter, 10, (self.panel_height*(counter-1))+5.0)
            counter = counter + 1
        # update the peak recognition, which applies the update of the miniviewer
        self.filled = 0
        self.parent.sumw.sumpanel.update_msp()
        # activate clickable scrolling of the MSPanels
    def mini_scroll(self, event):
        percent_move = (event.x - 10) / (self.width-10)
        percent_move = percent_move + (10.0/self.width)
        self.parent.msw.scrollAll("moveto", percent_move)
    def update(self):
        # delete any previous peaks drawn
        self.canvas.delete('peaks')
        # draw the peaks, colored by warnings
        counter = 1
        maxy = 0
        for panel in self.parent.msw.msp:
            for peak_index in panel.peaks:
                yval = panel.y[peak_index]
                if yval > maxy:
                    maxy = yval
        for panel in self.parent.msw.msp:
            for peak_index in panel.peaks:
                # get the scaled x value
                # find the ratio covered in the large spectra
                x = (panel.x[peak_index]/panel.x[len(panel.x)-1])
                # and scale it up to that of the small one
                x = x * (self.panel_length-15)
                # and the scaled y values
                offset = panel.y[peak_index]/maxy
                offset = offset * (self.panel_height-5.0)
                y = (self.panel_height*counter)-offset
                # draw the line
                self.canvas.create_line(x,self.panel_height*counter,x,y, fill='black', tags='peaks')
            counter = counter + 1

class MSWindow(Frame):
    def __init__(self, parent, frame, width, height):
        """ MSWindow acts as a holder for multiple MSPanel objects and coordinates their activity.
        """
        Frame.__init__(self, frame, borderwidth=0)
        hScaleFrame = parent.hScaleFrame
        self.type = parent.type
        self.reaction_type = parent.reaction_type
        self.parent=parent
        self.project=parent.project
        self.width=width-300
        self.height=height
        self.config(height=self.height, bg='white')
        labelfont = ('times', 9, 'bold')
        # get the file list
        if self.type == 'Analysis':
            self.filenameList = parent.project.analysis_files
        # heightScaler scales the heights of all ms plots
        Label(parent.hScaleFrame, text='scale\nyaxis', width=5, font=labelfont).pack(side=TOP)
        self.heightScaler = Scale(parent.hScaleFrame, from_=30.0, to=1.0, resolution=1.0, orient='vertical', width=15)
        self.heightScaler.set(10.0)
        self.heightScaler.pack(side=TOP,fill=Y, expand=1)
        # mess with the callback handler
        self.heightScaler.config(command=self.scaleAll)
        Label(parent.hScaleFrame, text='             ').pack(side=LEFT)
        # this is the count of files opened, the number to show
        if self.type == 'Analysis':
            filecount = 0
            for x in self.filenameList:
                filecount = filecount + len(x)
            self.mspCount = filecount
        self.msp = []               # this holds the mspanels
        self.msp_frame = Frame(self, height=self.height)
        self.msp_frame.pack(expand=YES, fill=BOTH)
        # add an x-scrollbar that handles all of the MSPanels
        self.xscrollframe = Frame(self.msp_frame)
        self.xscrollframe.pack(side=BOTTOM, expand=0, fill=X)
        self.xsbar=Scrollbar(self.xscrollframe)
        self.xsbar.config(orient='horizontal', command=self.scrollAll)
        self.xsbar.pack(side=BOTTOM, fill=X, expand=1, anchor=S)
        # and a text to put in it
        self.msp_text = Text(self.msp_frame, height=4)
        self.msp_text.height = self.height
        self.msp_text.width = self.width
        self.msp_text.parent = self
        # add a y-axis scroll bar
        self.msp_text.ysbar=Scrollbar(self.msp_frame)
        self.msp_text.ysbar.pack(side=RIGHT, fill=Y)
        self.msp_text.pack(expand=YES, fill=BOTH)
        self.msp_text.ysbar.config(command=self.msp_text.yview, cursor='arrow')
        # and add the mass spec panels
        insert_count = 0
        if self.type == 'Analysis':
            filecount = 0
            # first get a count of the files
            for x in self.filenameList:
                for y in x:
                    filecount = filecount + 1
            ms_list = []
            for x in self.filenameList:
                for y in x:
                    ms_list.append(MSPanel(self.msp_text, y, self.mspCount, self.type))
            #if len(self.project.standards_selected) > 0:
            #    print 'standardizing,' 
            #    self.standardize_MSPanels(ms_list)
            print 'normalizing'
            self.interpolate_MSPanels(ms_list)
            print 'and drawing panels'
            self.draw_MSPanels(ms_list)
            for ms in ms_list:
                ms.pack(expand=YES, fill=BOTH, side=TOP)
                self.msp_text.window_create('end', window=ms)
                insert_count = insert_count + 1
                if (insert_count < filecount-1):
                    self.msp_text.insert('end', '\n')
                self.msp.append(ms)
        self.msp_text.config(yscrollcommand=self.msp_text.ysbar.set, state=DISABLED)
        self.msp[0].plot.config(xscrollcommand=self.xsbar.set)
        
    def interpolate_MSPanels(self, mspanels):
        # first interpolate the later panels to the first
        # use a simple triangles method for now, but the vtk code can interpolate Bezier curves
        # first find the latest start and earliest end, in amus
        latest_start = -10000.0
        earliest_end = 1000000.0
        for p in mspanels:
            if p.x[0] > latest_start:
                latest_start = p.x[0]
            if p.x[-1] < earliest_end:
                earliest_end = p.x[-1]
        # now delimit the first mspanel to these bounds
        p1 = mspanels[0]
        start_ind = 0
        end_ind = len(p1.x)
        for x_ind in range(len(p1.x)):
            if p1.x[x_ind] >= latest_start:
                start_ind = x_ind
                break
        for x_ind in range(len(p1.x)-1, -1, -1):
            if p1.x[x_ind] <= earliest_end:
                end_ind = x_ind
                break
        p1.x = p1.x[start_ind:end_ind+1]
        p1.y = p1.y[start_ind:end_ind+1]
        for p2 in mspanels[1:]:
            new_xs = []
            new_ys = []
            last_2ind = 1
            for x1 in p1.x:
                # find the nearest pair of coordinates in p2
                for x2_ind in range(last_2ind-1, len(p2.x)):
                    if p2.x[x2_ind] >= x1:
                        # use triangles to interpolate
                        ratio = (x1-p2.x[x2_ind-1]) / (p2.x[x2_ind]-p2.x[x2_ind-1])
                        new_xs.append(p2.x[x2_ind-1] + ratio * (p2.x[x2_ind]-p2.x[x2_ind-1]))
                        new_ys.append(p2.y[x2_ind-1] + ratio * (p2.y[x2_ind]-p2.y[x2_ind-1]))
                        last_2ind = x2_ind
                        break
            p2.x = new_xs
            p2.y = new_ys
    
    """
    def standardize_MSPanels(self, mspanels):
        # first get any standards information available
        standards = []
        for std_ind in self.project.standards_selected:
            standards.append(string.atof(self.project.standards_weights[std_ind]))
        standards.sort()
        print 'standards', standards
        # so far, coordinates have only been generated, not drawn
        searchrange = .002    # percent weight for tolerance of x- shift
        # right now, this code assumes only one standard is selected
        
        if len(standards) > 1:
            print 'only working with one standard, deleting the others'
            standards = standards[0:0]
        standards_heights = []
        standards_indices = []
        standards_xs      = []
        for mspanel in mspanels:
            standards_heights.append([])       # each panel gets a list of heights for each standard
            standards_indices.append([])       # each panel gets a list of heights for each standard
            standards_xs.append([])
            for stdwt in standards:
                # first locate the x-index of the weight
                midind = -1
                midx   = 0.0
                xs = mspanel.x           # dereference speedup
                for x_ind in range(len(xs)):
                    if xs[x_ind] > stdwt:
                        if abs(xs[x_ind]-stdwt) < abs(xs[x_ind-1]-stdwt):
                            midind = x_ind
                            midx   = xs[x_ind]
                        else:
                            midind = x_ind-1
                            midx   = xs[x_ind-1]
                        break
                startind = -1
                startx   = 0.0
                for x_ind in range(len(xs)):
                    if midx-xs[x_ind] <= stdwt * searchrange:
                        if abs(xs[x_ind]-stdwt) < abs(xs[x_ind-1]-stdwt):
                            startind = x_ind
                            startx   = xs[x_ind]
                        else:
                            startind = x_ind-1
                            startx   = xs[x_ind-1]
                        break
                endind = -1
                endx   = 0.0
                for x_ind in range(len(xs)):
                    if xs[x_ind]-midx > stdwt * searchrange:
                        if abs(xs[x_ind]-stdwt) < abs(xs[x_ind-1]-stdwt):
                            endind = x_ind
                            endx   = xs[x_ind]
                        else:
                            endind = x_ind-1
                            endx   = xs[x_ind-1]
                        break
                maxind = 0
                maxht  = 0.0
                ys = mspanel.y
                for x_ind in range(startind, endind):
                    if ys[x_ind] > maxht:
                        maxht = ys[x_ind]
                        maxind = x_ind
                standards_heights[-1].append(string.atof(maxht))
                standards_indices[-1].append(maxind)
                standards_xs[-1].append(mspanel.x[maxind])
        # now find the lowest
        low_height = 100000.0
        for std in standards_heights:
            if std[0] < low_height:
                low_height = std[0]
        # normalize the y-axis
        panel_index = 0
        for mspanel in mspanels:
            ys = mspanel.y
            for y_ind in range(len(ys)):
                ys[y_ind] = ys[y_ind] * (low_height/standards_heights[panel_index][0])
            panel_index += 1
        # and normalize the x-axis
        panel_index = 0
        for mspanel in mspanels:
            dif = mspanel.x[standards_indices[panel_index][0]] - standards[0]            # find the offset of the panel
            xstd = mspanel.x[standards_indices[panel_index][0]]
            for x_ind in range(len(mspanel.x)):
                #mspanel.x[x_ind] = mspanel.x[x_ind] - (mspanel.x[x_ind]/xstd)*dif
                mspanel.x[x_ind] = mspanel.x[x_ind] - (mspanel.x[x_ind]/xstd)*dif
            panel_index += 1
    """
    def draw_MSPanels(self, mspanels):
        for mspanel in mspanels:
            mspanel.draw_coordinates()

    def scaleAll(self, value):
        val = self.heightScaler.get()
        for x in self.msp:
            x.rescale(x, val)
            x.rescale_y(x, val)
        if self.type == 'Analysis':
            self.parent.sumw.sumpanel.rescale(self.parent.sumw.sumpanel, val)
            self.parent.sumw.sumpanel.rescale_y(self.parent.sumw.sumpanel, val)
    def zoomUp(self, event):
        val = self.heightScaler.get()
        for x in self.msp:
            x.rescale(x, val+5)
            x.rescale_y(x, val+5)
        if self.type == 'Analysis':
            self.parent.sumw.sumpanel.rescale(x, val+5)
            self.parent.sumw.sumpanel.rescale_y(x, val+5)
        self.heightScaler.set(val+5)
    def zoomDown(self, event):
        val = self.heightScaler.get()
        for x in self.msp:
            x.rescale(x, val-5)
            x.rescale_y(x, val-5)
        if self.type == 'Analysis':
            self.parent.sumw.sumpanel.rescale(x, val-5)
            self.parent.sumw.sumpanel.rescale_y(x, val-5)
        self.heightScaler.set(val-5)
    def scrollAll(self, *args):
        if self.type == 'Analysis':
            apply(self.parent.sumw.sumpanel.plot.xview, args)
        for x in self.msp:
            apply(x.plot.xview, args)
    def scrollRight(self, event):
        self.scrollAll("scroll", 1, "units")
    def scrollLeft(self, event):
        self.scrollAll("scroll", -1, "units")
    def scrollMSW(self, *args):
        apply(self.msp_text.yview, args)
    def scrollUp(self, event):
        self.scrollMSW("scroll", -1, "units")
    def scrollDown(self, event):
        self.scrollMSW("scroll", 1, "units")

class SumWindow(Frame):
    def __init__(self, parent, frame, MSWindowList, width, ht, type):
        Frame.__init__(self, frame, borderwidth=2, height=10)
        self.width=width
        self.height=ht
        self.parent=parent
        self.type = type                            # type = Standards or PeakRec
        self.mswindowList = MSWindowList            # this is a list of filenames plotted
        self.mspCount = len(self.mswindowList)
        self.sumpanel = MSPanel(self, "", self.mspCount, "PeakSum")
        self.sumpanel.pack()
        self.sumpanel.plot.config(xscrollcommand=self.parent.xsbar.set)
        
    def scaleAll(self, value):
        #val = self.sumHeightScaler.get()
        val = self.parent.heightScaler.get()
        self.sumpanel.rescale(self.sumpanel, val)
        self.sumpanel.rescale_y(self.sumpanel, val)
    def zoomUp(self, event):
        #val = self.sumHeightScaler.get()
        val = self.parent.heightScaler.get()
        self.sumpanel.rescale(self.sumpanel, val+2)
        self.sumpanel.rescale_y(self.sumpanel, val+2)
        self.sumHeightScaler.set(val+2)
    def zoomDown(self, event):
        #val = self.sumHeightScaler.get()
        val = self.parent.heightScaler.get()
        self.sumpanel.rescale(self.sumpanel, val-2)
        self.sumpanel.rescale_y(self.sumpanel, val-2)
        self.sumHeightScaler.set(val-2)
    def scrollAll(self, *args):
        for x in self.parent.msp:
            apply(x.plot.xview, args)
        apply(self.sumpanel.plot.xview, args)
    def scrollRight(self, event):
        self.scrollAll("scroll", 1, "units")
    def scrollLeft(self, event):
        self.scrollAll("scroll", -1, "units")

class MSPanel(Frame):
    def __init__(self, parent, filename, brothers, type):
        # init
        Frame.__init__(self, parent, bg='white')
        self.all_my_labels = []       # this is where molecular weight labels get kept track of
        self.label_labels = []        #     and their x-coordinates or visa versa
        self.parent=parent
        self.reaction_type = parent.parent.reaction_type
        try:
            self.modification_system = self.parent.parent.parent.modification_system
        except AttributeError:
            print 'no modification_system found'
            self.modification_system = None
        self.filename=filename
        self.parent=parent
        self.brotherCount = brothers
        self.width=parent.width
        self.type = type            # type = Standards or Anaylsis or PeakRec
        # get some globals
        self.smoothing_gap = self.parent.parent.parent.smoothing_gap
        self.global_x_scale = self.parent.parent.parent.global_x_scale
        self.min_weight = self.parent.parent.parent.min_weight
        self.height_thresh = self.parent.parent.parent.height_thresh
        self.min_peak_height = self.parent.parent.parent.min_peak_height
        self.x_axis_mark_every = parms.get('x_axis_mark_every')
        # define height by type of window
        if (self.type == 'Analysis'):
            self.height = (self.parent.height/5.0)
        elif (self.type == 'PeakSum'):
            self.height = self.parent.height/5.0
        #self.config(width=self.width, height=self.height)
        # create a canvas to draw on
        self.plot = Canvas(self, bg='white', relief=SUNKEN, bd=0, cursor='crosshair', height=self.height, width=self.width)
        self.plot.focus_set()
        self.plot.config(highlightthickness=1)
        self.yaxis_width = 25
        # create a canvas for the first y-axis
        self.preplot = Canvas(self, relief=SUNKEN, height=self.height, width=self.yaxis_width, bd=0)
        # and one for the second y-axis
        self.postplot = Canvas(self, relief=SUNKEN, height=self.height, width=self.yaxis_width, bd=0)
        # and draw the plots and axes
        self.use_blank = 0
        if (type == 'PeakSum'):
            self.draw_sum_plot()
        else:
            self.create_coordinates()
        self.draw_x_axes()
        self.draw_y_axes()
        self.plot.config(scrollregion=(0,self.x[0], self.x[len(self.x)-1], self.height))
        # finally, pack everything up
        self.preplot.pack(expand=NO, fill=Y, side=LEFT)
        self.postplot.pack(expand=NO, fill=Y, side=RIGHT)
        self.plot.pack(expand=YES, fill=X, side=TOP)
        # used to hold heights for the sum of quanitification distributions
        self.dist_y = []
    def takefocus(self, event):
        self.plot.focus_set()
    def print_message(self, event):
        print 'got here'
    def draw_sum_plot(self):
        """ assumes the panels are normalized and standardized """
        self.x = []
        self.y = []
        for ind in range(len(self.parent.parent.msp[0].x)):
            avg_x = 0.0
            avg_y = 0.0
            for p in self.parent.parent.msp:
                avg_x += p.x[ind]/len(self.parent.parent.msp)
                avg_y += p.y[ind]/len(self.parent.parent.msp)
            self.x.append(avg_x)
            self.y.append(avg_y)
            
        i = 0
        for zz in self.y[:-1]:
            self.plot.create_line(self.x[i], self.y[i], self.x[i+1], self.y[i+1], width=1, fill='black', tags='plot')
            i = i+1
        # scale the image
        self.plot.scale("plot", 0, 0, 1, -1.5)
        self.plot.move("plot", 0, (self.height * 0.80))
        # last_mag is used to reset the magnification to zero for changing to a new magnification
        self.plot.last_mag = 1.0
        self.plot.last_scr = self.height
        scaleval = self.parent.parent.heightScaler.get()
        self.created_peaks = []
        self.created_lin_list = []
        self.recognize_peaks()
    
    def update_msp(self):
        height_thresh = parms.get('height_thresh')
        svsl_peak_thresh = parms.get('svsl_peak_thresh')
        # given a new peak selection in the sum window, the peaks in other msps are updated here
        update_range = 9                         # how far each direction will be searched 
        for x in self.parent.parent.msp:                             # for local maxima
            x.peaks = []
            x.border_test = []
            x.parent_peak = []                      # the index (i, below) of the parent peak
            x.other_peak_test = []
            x.plot.delete('largepeaks')
            x.plot.delete('smallpeaks')
        max_place = 0                             # store the height of the current peak
        max_in_range = 0                          # store the index of the parent peak
        # for each peak in the sum window, look for the corresponding peak from each
        # individual msp
        for i in range(len(self.peaks)):
            for panel in self.parent.parent.msp:
                max_in_range = panel.y[self.peaks[i]]
                max_place    = self.peaks[i]      
                test1=0
                test2=0
                # look right one, left one, right two, left two,...
                for z in range(1,update_range):
                    # if the right side hasn't already bottomed out
                    if (test1==0):
                        if (panel.y[self.peaks[i]+z] >= max_in_range):
                            # store the maximum so far encountered
                            max_in_range = panel.y[self.peaks[i]+z]
                            max_place = self.peaks[i]+z
                        elif (panel.y[self.peaks[i]+z] < (panel.y[max_place] - height_thresh)):
                            # don't look any further in this direction if a trough reaches a y-val
                            # lower than the highest point found so far minus the threshold for peak
                            # recognition
                            test1 = 1
                    # now look one more in the left direction
                    if (test2==0):
                        if (panel.y[self.peaks[i]-z] >= max_in_range):
                            max_in_range = panel.y[self.peaks[i]-z]
                            max_place = self.peaks[i]-z
                        elif (panel.y[self.peaks[i]-z] < (panel.y[max_place] - height_thresh)):
                            test2 = 1
                # if the border is the maximum there might be a problem. store a reference, to
                # later color it green as a warning
                if (max_place == self.peaks[i]-update_range or max_place == self.peaks[i]+update_range):
                    panel.border_test.append(max_place)
                if (len(panel.peaks) > 0):
                    # see if the peak is already there. if so, store a reference, to later color
                    # it red as a warning
                    for qq in panel.peaks:
                        if (max_place == qq):
                            panel.other_peak_test.append(max_place)
                            break
                # for each peak, the parent peak is stored in a corresponding index
                panel.parent_peak.append(i)
                panel.peaks.append(max_place)

        # now draw out the peaks that were recognized
        for panel in self.parent.parent.msp:
            panel.lin_list = []
            panel.lin_index = []
            for i in range(len(panel.peaks)):
                if (panel.y[panel.peaks[i]] < svsl_peak_thresh):
                    for xx in panel.border_test:
                        if (xx == panel.peaks[i]):
                            lin = panel.plot.create_line(panel.x[panel.peaks[i]],0,
                                                    panel.x[panel.peaks[i]],panel.height*0.8,
                                                    fill='green', tags='smallpeaks')
                            break
                    else:
                        for yy in panel.other_peak_test:
                            if (yy == panel.peaks[i]):
                                lin = panel.plot.create_line(panel.x[panel.peaks[i]],0,
                                                panel.x[panel.peaks[i]],panel.height*0.8,
                                                fill='red', tags='smallpeaks')
                                break
                        else:
                            lin = panel.plot.create_line(panel.x[panel.peaks[i]],0,
                                            panel.x[panel.peaks[i]],panel.height*0.8,
                                                fill='cyan', tags='smallpeaks')
                    panel.lin_index.append(i)
                    panel.lin_list.append(lin)
                else:
                    for xx in panel.border_test:
                        if (xx == panel.peaks[i]):
                            lin = panel.plot.create_line(panel.x[panel.peaks[i]],0,
                                                panel.x[panel.peaks[i]],panel.height*0.8,
                                                fill='green', tags='largepeaks')
                            break
                    else:  
                        for yy in panel.other_peak_test:
                            if (yy == panel.peaks[i]):
                                lin = panel.plot.create_line(panel.x[panel.peaks[i]],0,
                                            panel.x[panel.peaks[i]],panel.height*0.8,
                                            fill='red', tags='largepeaks')
                                break
                        else:
                            lin = panel.plot.create_line(panel.x[panel.peaks[i]],0,
                                            panel.x[panel.peaks[i]],panel.height*0.8,
                                            fill='blue', tags='largepeaks')
                    panel.lin_index.append(i)
                    panel.lin_list.append(lin)
            for id in panel.lin_list:
                panel.plot.tag_bind(id, '<Button-1>', panel.movepeak)
                panel.plot.tag_bind(id, '<ButtonRelease-1>', panel.releasepeak)

        # update the miniviewer
        if self.parent.parent.parent.minibox_filled:
            self.parent.parent.parent.minibox.update()
        if self.parent.parent.parent.sumbox != None:
            # update the sumbox
            sb = self.parent.parent.parent.sumbox
            sb.initialize_peak_fragment_associations()
            sb.updateRateWindow(0, 1)
            sb.sequence_viewer.reset_fragments(sb.get_active_fragments())
            sb.sequence_viewer.draw()

    def movepeak(self, event):
        # first just make the line dashed, and add a tag
        find_scale = 2
        id = event.widget.find_overlapping(self.plot.canvasx(event.x)-find_scale, self.plot.canvasy(event.y)-find_scale,
                                    self.plot.canvasx(event.x)+find_scale, self.plot.canvasy(event.y)+find_scale)
        lin_count = 0
        last_peak = 0
        next_peak = 0
        for i in self.lin_list:
            if (id[0] == i):
                self.plot.itemconfigure(i, stipple='gray50')
                self.plot.addtag_withtag('motion_active', i)
                
    def releasepeak(self, event):
        # first figure out which index of lin_list is tagged
        lin_id = self.plot.find_withtag('motion_active')
        # now get the index of it
        count = 0
        for i in range(len(self.lin_list)):
            count = i
            if (lin_id[0] == self.lin_list[count]):
                break
        # get the indices of the peaks before and after
        this_peak_x = self.peaks[count]
        x_before = 0                    # the self.x index of the peak before
        bef_dist = 10000
        one_before = 0
        x_after  = 0                    # the index of the peak after
        aft_dist = 10000
        one_after = 0
        for i in range(len(self.peaks)):
            if (self.peaks[i] < this_peak_x):
                if (abs(this_peak_x - self.peaks[i]) < bef_dist):
                    one_before = 1
                    bef_dist = abs(this_peak_x - self.peaks[i])
                    x_before = self.peaks[i]
            elif (self.peaks[i] > this_peak_x):
                if (abs(this_peak_x - self.peaks[i]) < aft_dist):
                    one_after = 1
                    aft_dist = abs(this_peak_x - self.peaks[i])
                    x_after  = self.peaks[i]
        # see if the mouse left the range between the previous and next peaks
        canvas_x = self.plot.canvasx(event.x)
        placement = 0
        place_idx = 0
        closest = 1000
        if (canvas_x < self.x[x_before]+2):
            placement = self.x[x_before+2]
            place_idx = x_before+2
        elif (canvas_x > self.x[x_after]-2):
            placement = self.x[x_after-2]
            place_idx = x_after-2
        else:
            # search for the closest x bin to the place the mouse was released
            for i in range(x_before+2, x_after-2):
                if (abs(canvas_x-self.x[i]) < closest):
                    placement = self.x[i]
                    place_idx = i
                    closest = abs(canvas_x-self.x[i])
        lin = self.plot.create_line(placement,0,
                                    placement,self.height*0.8,
                                    fill='blue', tags='largepeaks')
        self.plot.tag_bind(lin, '<Button-1>', self.movepeak)
        self.plot.tag_bind(lin, '<ButtonRelease-1>', self.releasepeak)
        self.plot.delete(self.lin_list[count])
        self.lin_list[count] = lin
        self.peaks[count] = place_idx

    def create_coordinates(self):
        self.plot.delete('plot')
        smoothing_gap = self.smoothing_gap = parms.get('smoothing_gap')
        global_x_scale = self.global_x_scale = parms.get('global_x_scale')
        lines = open(self.filename).readlines()
        i = 0
        x_temp = []
        y_temp = []
        # write the lines to an array for easy access drawing them
        # convert scientific notation to standard
        for zz in lines:
            if (string.find(zz, "Intensity") > 0):
                continue
            if len(zz) == 0:
                break
            coors = string.split(zz, ',')
            if string.find(coors[0], 'e') > 0:
                x_temp.append(float('0.000'))
            else:
                x_temp.append(global_x_scale*float(coors[0]))
            if string.find(coors[1], 'e') > 0:
                y_temp.append(float('0.000'))
            else:
                y_temp.append(float(coors[1]))
        if (self.use_blank):
            blank_lines = open(self.blank_file).readlines()
            j = 0
            x_blank_temp = []
            y_blank_temp = []
            for zz in blank_lines:
                if (string.find(zz, "Intensity")):
                    continue
                if len(zz) == 0:
                    break
                coors = string.split(zz, ',')
                if string.find(coors[1], 'e') > 0:
                    x = 1
                else:
                    y_temp[j] = y_temp[j]-float(coors[1])
                j = j+1
        # now that the lines are read in, run the smoothing algorithm and move the info
        # over to self.x and self.y
        self.x = []
        self.y = []
        tmpy = 0
        # iterate through the indices that will be filled, starting at the point at
        # which each coordinate can be represented by the full smoothing range
        for zz in range(smoothing_gap, len(x_temp)-smoothing_gap-1):
            self.x.append(x_temp[zz])                                 # store the original x-coor
            tmpy = 0.00
            if (smoothing_gap > 0):                                   # skip this loop if no smoothing
                for yy in range(zz-smoothing_gap, zz+smoothing_gap):    # over the range to smooth
                    tmpy = tmpy + y_temp[yy]                # collect the values
                tmpy = tmpy /((smoothing_gap*2)-1)                      # and divide by the width
            else:                                                     # if no smoothing
                tmpy = y_temp[zz]                                       # just copy the y-val directly
            self.y.append(tmpy)
        self.plot.last_mag = 1.0
        self.plot.last_scr = self.height
    
    def draw_coordinates(self):
        # and draw the lines
        i = 0
        for zz in self.y:
            if i == len(self.y)-5:
                break;
            if i == len(self.x)-5:
                break;
            self.plot.create_line(self.x[i], self.y[i], self.x[i+1], self.y[i+1], width=1, fill='black', tags='plot')
            i = i+1
        # scale the image
        self.plot.scale("plot", 0, 0, 1, -1.5)
        self.plot.move("plot", 0, (self.height * 0.80))
        # last_mag is used to reset the magnification to zero for changing to a new magnification
        scaleval = self.parent.parent.heightScaler.get()
        for x in self.parent.parent.msp:
            x.rescale(x, scaleval)
    
    def rescale(MSWindow, self, val):
        for type in ("plot", "dist_plot", "peak_plot"):
            self.plot.move(type, 0, -(self.height * 0.80))
            self.plot.scale(type, 0, 0, 1, (1.0/self.plot.last_mag))
            self.plot.scale(type, 0, 0, 1, val)
            self.plot.move(type, 0, (self.height * 0.80))
        self.plot.last_mag = val;

    def recognize_modification_peaks(self):
        """ recognizing modified peaks is significantly different from recognizing normal proteolyzed fragments
            and the recognize_peaks function seems stable enough, so I'm forking for now. Remember that some 
            modifications should affect both.
        """
        self.plot.delete('largepeaks')
        self.plot.delete('smallpeaks')
        svsl_peak_thresh = parms.get('svsl_peak_thresh')
        direction = 0
        last_peak = 0
        last_trough = 0
        self.peaks = []
        self.troughs = []
        self.peak_count = 0
        self.trough_count = 0
        smoothing_gap = self.parent.parent.parent.smoothing_gap
        global_x_scale = self.parent.parent.parent.global_x_scale
        min_weight = self.parent.parent.parent.min_weight
        height_thresh = self.parent.parent.parent.height_thresh
        min_peak_height = self.parent.parent.parent.min_peak_height
        # first recognize peaks in the zero timepoint, to confirm fragment recognition
        self.parent.parent.msp[0].recognize_peaks(0) # dont draw
        confirmation_peaks = self.parent.parent.msp[0].peaks
        # this code performs the recognition of peaks in the sum window
        for i in range(0,(len(self.x)-2*smoothing_gap)):
            if self.x[i] > global_x_scale * min_weight:
                if direction == 0:                                                      # if no direction yet
                    if self.y[last_peak] >= self.y[i] + height_thresh:                  #   if last peak >> current
                        direction = -1                                                  #     then decreasing
                    elif self.y[i] >= self.y[last_trough] + height_thresh:              #   else if current >> last trough
                        direction = 1                                                   #     then increasing
                    if self.y[last_peak] < self.y[i]:                                   #   if last peak < current 
                        last_peak = i                                                   #     last peak = current
                    elif self.y[i] < self.y[last_trough]:                               #   else if current < last trough
                        last_trough = i                                                 #     last trough = current
                elif direction == 1:                                                    # else if increasing
                    if self.y[last_peak] < self.y[i]:                                   #   if last peak < current
                        last_peak = i                                                   #     last peak = current
                    elif self.y[last_peak] >= self.y[i] + height_thresh:                #   else if last peak >> current
                        direction = -1                                                  #     direction decreasing
                        last_trough = i                                                 #     last trough = current
                        if self.y[i] > min_peak_height:                                 #     if current > min peak height
                            print '%s > %s'%(self.y[i], min_peak_height)
                            if self.modification_system == None:
                                print self.modification_system
                                self.peaks.append(last_peak)                            #       record this peak
                                self.peak_count = self.peak_count+1                     #
                            else:
                                # for changing the number of missed sites, modify the fill_proteolysis_fragments call, not this
                                fraglist = self.modification_system.get_proteolysis_fragments_within(self.x[last_peak], self.x[last_peak]*0.002)
                                if len(fraglist) > 0:
                                    for confirmation_peak in confirmation_peaks:        # peaks should be within two ranges of error of each other
                                        first = self.x[confirmation_peak]
                                        second = self.x[last_peak]
                                        if abs(first-second) < (first+second)/10000:      # ... < 2 * (((first+second)/2) / 500)
                                            self.peaks.append(last_peak)                #       record this peak
                                            self.peak_count = self.peak_count+1         #
                elif direction == -1:                                                   # else if decreasing
                    if self.y[last_trough] > self.y[i]:                                 #   if last trough > current
                        last_trough = i                                                 #     last trough = current
                    elif self.y[i] >= self.y[last_trough] + height_thresh:              #   else if current >> last trough
                        direction = 1                                                   #     direction increasing
                        last_peak = i                                                   #     last peak = current
                        self.troughs.append(last_trough)                                #     record this trough
        # make sure none have been stored in the 'manually created' list, if so, delete them
        # from the new set
        for ii in self.created_peaks:
            count = 0
            for jj in self.peaks:
                if (ii == jj):
                    del self.peaks[count]
                    break
                count = count+1
        # create lines for the peaks that were recognized in the sum window
        self.lin_list = []
        for i in range(len(self.peaks)):
            if (self.y[self.peaks[i]] < svsl_peak_thresh):
                lin = self.plot.create_line(self.x[self.peaks[i]],0,
                                        self.x[self.peaks[i]],self.height*0.8,
                                        fill='cyan', tags='smallpeaks')
                test = 0
                self.lin_list.append(lin)
            else:
                lin = self.plot.create_line(self.x[self.peaks[i]],0,
                                        self.x[self.peaks[i]],self.height*0.8,
                                        fill='blue', tags='largepeaks')
                self.lin_list.append(lin)
            
        self.plot.scale("largepeaks", 0, 0, 1, -1.5)
        self.plot.move("largepeaks", 0, (self.height * 0.80))
        self.plot.scale("smallpeaks", 0, 0, 1, -1.5)
        self.plot.move("smallpeaks", 0, (self.height * 0.80))
        for id in self.lin_list:
            self.plot.tag_bind(id, '<ButtonPress-1>', self.changePeakLabel);
        self.plot.bind('<ButtonPress-3>', self.removePeak);
        # now append created peaks to the peak list, if it isn't already there.
        for ii in self.created_peaks:
            for jj in self.peaks:
                if (ii == jj):
                    break;
            else:
                self.peaks.append(ii)
        for ii in self.created_lin_list:
            for jj in self.lin_list:
                if (ii == jj):
                    break
            else:
                self.lin_list.append(ii)
                
    def recognize_peaks(self, draw=1):
        print 'called recognize peaks on %s %s.'%(self.reaction_type, self.type)
        if self.reaction_type == 'NHS' and self.type == 'PeakSum':
            print 'calling recognize_modification peaks instead'
            self.recognize_modification_peaks()
            return
        print 'going ahead to recognize peaks'
        if draw:
            self.plot.delete('largepeaks')
            self.plot.delete('smallpeaks')
        svsl_peak_thresh = parms.get('svsl_peak_thresh')
        direction = 0
        last_peak = 0
        last_trough = 0
        self.peaks = []
        self.troughs = []
        self.peak_count = 0
        self.trough_count = 0
        smoothing_gap = self.parent.parent.parent.smoothing_gap
        global_x_scale = self.parent.parent.parent.global_x_scale
        min_weight = self.parent.parent.parent.min_weight
        height_thresh = self.parent.parent.parent.height_thresh
        min_peak_height = self.parent.parent.parent.min_peak_height
        # this code performs the recognition of peaks in the sum window
        for i in range(0,(len(self.x)-2*smoothing_gap)):
            if self.x[i] > global_x_scale * min_weight:
                if direction == 0:                                                      # if no direction yet
                    if self.y[last_peak] >= self.y[i] + height_thresh:                  #   if last peak >> current
                        direction = -1                                                  #     then decreasing
                    elif self.y[i] >= self.y[last_trough] + height_thresh:              #   else if current >> last trough
                        direction = 1                                                   #     then increasing
                    if self.y[last_peak] < self.y[i]:                                   #   if last peak < current 
                        last_peak = i                                                   #     last peak = current
                    elif self.y[i] < self.y[last_trough]:                               #   else if current < last trough
                        last_trough = i                                                 #     last trough = current
                elif direction == 1:                                                    # else if increasing
                    if self.y[last_peak] < self.y[i]:                                   #   if last peak < current
                        last_peak = i                                                   #     last peak = current
                    elif self.y[last_peak] >= self.y[i] + height_thresh:                #   else if last peak >> current
                        direction = -1                                                  #     direction decreasing
                        last_trough = i                                                 #     last trough = current
                        if self.y[i] > min_peak_height:                                 #     if current > min peak height
                            if self.modification_system == None:
                                self.peaks.append(last_peak)                            #       record this peak
                                self.peak_count = self.peak_count+1                     #
                            else:
                                fraglist = self.modification_system.get_proteolysis_fragments_within(self.x[last_peak], self.x[last_peak]*0.002)
                                if len(fraglist) > 0:
                                    self.peaks.append(last_peak)                        #       record this peak
                                    self.peak_count = self.peak_count+1                 #
                elif direction == -1:                                                   # else if decreasing
                    if self.y[last_trough] > self.y[i]:                                 #   if last trough > current
                        last_trough = i                                                 #     last trough = current
                    elif self.y[i] >= self.y[last_trough] + height_thresh:              #   else if current >> last trough
                        direction = 1                                                   #     direction increasing
                        last_peak = i                                                   #     last peak = current
                        self.troughs.append(last_trough)                                #     record this trough
        if draw:
            # make sure none have been stored in the 'manually created' list, if so, delete them
            # from the new set
            for ii in self.created_peaks:
                count = 0
                for jj in self.peaks:
                    if (ii == jj):
                        del self.peaks[count]
                        break
                    count = count+1
            # create lines for the peaks that were recognized in the sum window
            self.lin_list = []
            for i in range(len(self.peaks)):
                if (self.y[self.peaks[i]] < svsl_peak_thresh):
                    lin = self.plot.create_line(self.x[self.peaks[i]],0,
                                            self.x[self.peaks[i]],self.height*0.8,
                                            fill='cyan', tags='smallpeaks')
                    test = 0
                    self.lin_list.append(lin)
                else:
                    lin = self.plot.create_line(self.x[self.peaks[i]],0,
                                            self.x[self.peaks[i]],self.height*0.8,
                                            fill='blue', tags='largepeaks')
                    self.lin_list.append(lin)
                
            self.plot.scale("largepeaks", 0, 0, 1, -1.5)
            self.plot.move("largepeaks", 0, (self.height * 0.80))
            self.plot.scale("smallpeaks", 0, 0, 1, -1.5)
            self.plot.move("smallpeaks", 0, (self.height * 0.80))
            for id in self.lin_list:
                self.plot.tag_bind(id, '<ButtonPress-1>', self.changePeakLabel);
            self.plot.bind('<ButtonPress-3>', self.removePeak);
            # now append created peaks to the peak list, if it isn't already there.
            for ii in self.created_peaks:
                for jj in self.peaks:
                    if (ii == jj):
                        break;
                else:
                    self.peaks.append(ii)
            for ii in self.created_lin_list:
                for jj in self.lin_list:
                    if (ii == jj):
                        break
                else:
                    self.lin_list.append(ii)
        
    def addPeak(self, event_x):
        # this code responds to button clicks that manually add peaks.
        last_y = 0
        count = 0
        index = 0
        # find the index of the contact point
        for tx in self.x:
            if (self.plot.canvasx(event_x) < tx):
                break
            last_y = self.y[count]
            index = count
            count = count+1
        # see if its already in there
        test = 0
        for xx in self.peaks:
            if (xx == index):
                test = 1
                return
        self.peaks.append(index)           # add this one to the list of peaks
        self.created_peaks.append(index)   # and to created_peaks, so it isn't deleted with updates
        scaleval = self.parent.parent.heightScaler.get()
        lin = self.plot.create_line(self.x[index], 0,
                                    self.x[index], self.height*0.8,
                                    fill='purple', tags='created_peaks')
        self.plot.tag_bind(lin, '<ButtonPress-1>', self.changePeakLabel);
        #self.plot.tag_bind(lin, '<Double-1>', self.removePeak);
        # append renewable and permenant lists
        test = 0
        for xx in self.lin_list:
            if (lin == xx):
                test = 1
        if (test == 0):
            self.lin_list.append(lin)
            self.created_lin_list.append(lin)
        # fix the plot    
        self.plot.scale(lin, 0, 0, 1, -1.5)
        self.plot.move(lin, 0, (self.height * 0.80))
        # and update the rest
        self.update_msp()
            
    def removePeak(self, event):
        id = event.widget.find_overlapping( self.plot.canvasx(event.x)-1,self.plot.canvasy(event.y)+1,self.plot.canvasx(event.x)+1,self.plot.canvasy(event.y)-1)
        count = 0  # remove the peak from this list
        peak = 0
        test = 0
        if (len(id) > 0):
            for i in self.lin_list:
                if (id[0] == i):                        # if this is the line
                    my_coords = self.plot.coords(i)
                    # first get rid of any labels
                    label_count = 0
                    for ii in self.label_labels:
                        if (ii == my_coords[0]):
                            self.plot.delete(self.all_my_labels[label_count])
                            del self.label_labels[label_count]
                            del self.all_my_labels[label_count]
                        for x in self.parent.parent.msp:
                            x.plot.delete(x.all_my_labels[label_count])
                            del x.label_labels[label_count]
                            del x.all_my_labels[label_count]
                        break
                        label_count = label_count + 1
                    # remove it from created_peaks, if it is there
                    peak = self.peaks[count]
                    c_count = 0
                    for ii in self.created_peaks:
                        if (ii == peak):
                            del self.created_peaks[c_count]
                            del self.created_lin_list[c_count]
                            break
                        c_count = c_count + 1
                    # now remove the line
                    self.plot.delete(self.lin_list[count])
                    del self.peaks[count]
                    del self.lin_list[count]
                    test = 1
                    break
                count = count + 1
        if (test == 0):
            self.addPeak(event.x)
        self.update_msp()
    
    def changePeakLabel(self, event):
        find_scale = 2
        id = event.widget.find_overlapping(self.plot.canvasx(event.x)-find_scale, self.plot.canvasy(event.y)-find_scale,
                                        self.plot.canvasx(event.x)+find_scale, self.plot.canvasy(event.y)+find_scale)
        lin_count = 0
        msp_coords = []
        msp_height = []
        for i in self.lin_list:
            if (id[0] == i):
                my_coords = self.plot.coords(i)
                for peak_index in self.peaks:
                    if (round(self.x[peak_index],2) == round(my_coords[0],2)):
                        my_height = self.y[peak_index]
                        break
                # save the children, to label them as well
                for x in self.parent.parent.msp:
                    x_lin = x.lin_list[lin_count]
                    xcor = x.plot.coords(x_lin)
                    msp_coords.append(xcor[0])
                    msp_count = 0
                    for peak_index in x.peaks:
                        if (round(x.x[peak_index],2) == round(xcor[0],2)):
                            msp_height.append(x.y[peak_index])
                            break
                        msp_count = msp_count+1
                my_count = 0
                # delete a label, if found for this coordinate, and those of children
                for ii in self.label_labels:
                    if (ii == my_coords[0]):
                        self.plot.delete(self.all_my_labels[my_count])
                        del self.label_labels[my_count]
                        del self.all_my_labels[my_count]
                        for x in self.parent.parent.msp:
                            x.plot.delete(x.all_my_labels[my_count])
                            del x.label_labels[my_count]
                            del x.all_my_labels[my_count]
                        break
                    my_count = my_count + 1
                else:
                    labelFont = ('Helvetica', parms.get('label_font_size'))
                    # label the molecular weight just below the line recognizing the peak
                    my_label = str(int(round((my_coords[0]/self.global_x_scale), 0)))+'\n'+str(round(my_height,2))
                    id = self.plot.create_text(my_coords[0], -8, fill='black', text=my_label, tags='mw_label', font=labelFont)
                    self.plot.scale(id, 0, 0, 1, -1.5)
                    self.plot.move(id, -0.5, (self.height * 0.80))
                    self.all_my_labels.append(id)
                    self.label_labels.append(my_coords[0])
                    self.plot.tkraise(id)
                    # label the children
                    msp_count = 0
                    for x in self.parent.parent.msp:
                        my_label = str(int(round((msp_coords[msp_count]/self.global_x_scale), 0)))+'\n'+str(round(msp_height[msp_count],2))
                        id = x.plot.create_text(msp_coords[msp_count], -10, fill='black', text=my_label, tags='mw_label', font=labelFont)
                        x.plot.scale(id, 0, 0, 1, -1.5)
                        x.plot.move(id, -0.5, (x.height * 0.80))
                        x.all_my_labels.append(id)
                        x.label_labels.append(msp_coords[msp_count])
                        x.plot.tkraise(id)
                        msp_count = msp_count + 1
                break
            lin_count = lin_count + 1
    def draw_y_axes(self):
        self.preplot.delete('y_axis')
        self.postplot.delete('y_axis')
        labelFont = ('Helvetica', parms.get('label_font_size'))
        for zz in range(0, int(50.0 * self.height), int(round(100 * parms.get('y_axis_mark_every'),0))):
            zz = zz/100
            # pre plot first
            widget = Label(self.preplot, text=zz, fg='red', font=labelFont)
            widget.pack()
            self.preplot.create_window(10, zz, window=widget, tags='y_axis')
            self.preplot.create_line(self.yaxis_width-5, zz, self.yaxis_width, zz, fill='red', tags='y_axis')
            # post plot
            widget = Label(self.postplot, text=zz, fg='red', font=labelFont)
            widget.pack()
            self.postplot.create_window(self.yaxis_width-10, zz, window=widget, tags='y_axis')
            self.postplot.create_line(0, zz, 5, zz, fill='red', tags='y_axis')
        self.preplot.scale("y_axis", 0, 0, 1, -1.5)
        self.preplot.move("y_axis", 0, (self.height * 0.80))
        self.postplot.scale("y_axis", 0, 0, 1, -1.5)
        self.postplot.move("y_axis", 0, (self.height * 0.80))
        self.preplot.last_mag = 1.0
        self.preplot.last_scr = self.height
        scaleval = self.parent.parent.heightScaler.get()
        for x in self.parent.parent.msp:
            x.rescale_y(x, scaleval)
    def rescale_y(MSWindow, self, val):
        for plot in [self.preplot, self.postplot]:
            plot.move('y_axis', 0, -(self.height * 0.80))
            plot.scale('y_axis', 0, 0, 1, (1.0/self.preplot.last_mag))
            plot.scale('y_axis', 0, 0, 1, val)
            plot.move('y_axis', 0, (self.height * 0.80))
        self.preplot.last_mag = val;
    def delete_y_axes(self):
        self.plot.delete('y_axis')
    def draw_x_axes(self):
        self.plot.delete('x_axis')
        # draw the lines now
        self.plot.create_line(0,0,self.global_x_scale*self.x[-1],0,width=0.5, fill='red', tags='x_axis')
        labelFont = ('Helvetica', parms.get('label_font_size'))     
        for zz in range(0,int(self.global_x_scale*self.x[-1]), int(self.global_x_scale*self.x_axis_mark_every)):
            self.plot.create_line(zz, -2, zz, 2, width=0.5, fill='red', tags='x_axis')
            my_label = str(round(zz/self.global_x_scale))
            x = self.plot.create_text(zz, -10, fill='red', text = my_label, tags = 'x_axis', font=labelFont)
            self.plot.lower(x)
        self.plot.scale("x_axis", 0, 0, 1, -1.5)
        self.plot.move("x_axis", 0, (self.height * 0.80))
    def delete_x_axes(self):
        self.plot.delete('x_axis')

