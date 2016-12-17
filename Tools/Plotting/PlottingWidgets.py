# python imports
import ImageTk
import string
import math
import os.path
import time
import random
import scipy.stats
import fpformat
import pickle
import sys
import copy
# dependency imports
from Tkinter import *
from tkFileDialog import *
sys.path.append(os.path.abspath('./Dependencies'))
import Pmw
# internal imports
import parms
import MolecularSystem
import MolecularViewer
sys.path.append(os.path.abspath('./Tools/Math'))
import Distribution
sys.path.append(os.path.abspath('./Tools/Plotting'))
import PlottingTools
sys.path.append(os.path.abspath('./Tools/MassSpectrometry'))
import PeakRationalizer
import FragmentPeakTablet

class FixedMenu(Menu):
    """ This fix for the Menu delete command was produced by Sverker Nilsson, and found at
        http://sourceforge.net/tracker/index.php?func=detail&aid=1342811&group_id=5470&atid=105470
    """
    # A fix for the .delete() method in Menu.
    # To delete commands defined in the menu items deleted.
    # Also changed the comment: INDEX2 is actually INCLUDED.
    def delete_special(self, index1, index2=None):
        """Delete menu items between INDEX1 and INDEX2 (included)."""
        if index2 is None:
            index2 = index1
        # First find out what entries have defined commands.
        cmds = []
        for i in range(Menu.index(self, index1), Menu.index(self, index2)+1):
            c = str(Menu.entrycget(self, i, 'command'))
            if c in self._tclCommands:
                # I don't want to delete the command already, since it
                # seems mystical to do that while the entry is not yet deleted.
                cmds.append(c)
        # Delete the menu entries.
        top = Menu.winfo_toplevel(self)
        top.tk.call(top._w, 'delete', index1, index2)
        # Now that the menu entries have been deleted, we can delete their commands.
        for c in cmds:
            Menu.deletecommand(self, c)


class PeakProfile:
    def __init__(self, weight, heights, areas, indices):
        self.weight          = weight
        self.heights         = heights
        self.trough1_indices = indices[0]
        self.peak_indices    = indices[1]
        self.trough2_indices = indices[2]
        self.tag_ids         = {}
        self.avg_heights     = []
        self.areas           = areas
        self.avg_areas       = []

        keys = heights.keys()
        keys.sort()

        for key in keys:
            heightset = heights[key]
            sum = 0.0
            for height in heightset:
                sum += height
            sum /= len(heightset)
            self.avg_heights.append(sum)

        for key in keys:
            areaset = areas[key]
            sum = 0.0
            for area in areaset:
                sum += area
            sum /= len(areaset)
            self.avg_areas.append(sum)
            
    def get_indices(self):
        return trough1_indices, peak_indices, trough2_indices
    
    def get_areas(self):
        return self.areas

    def get_areas_as_lists(self):
        keys = self.areas.keys()
        keys.sort()
        lists = []
        for i in range(len(self.areas[keys[0]])):
            lists.append([])
        for i in range(len(self.areas[keys[0]])):
            for j in range(len(keys)):
                lists[i].append(self.areas[keys[j]][i])
        return lists
    
    def get_avg_areas(self):
        return self.avg_areas

    def get_weight(self):
        return self.weight
    
    def get_heights(self):
        """ returns a dictionary, where timepoint labels are keys to lists of heights
        """
        return self.heights
    
    def get_avg_heights(self):
        return self.avg_heights
    
    def get_tag_ids(self):
        return self.tag_ids

    def set_tag_ids(self, tags):
        self.tag_ids = tags

class ReactionProfile:
    def __init__(self, fragments, peak, modified_peak=None):
        self.possible_fragment_list = fragments
        self.peak = peak
        self.modified_peak = modified_peak

    def get_timepoint_labels(self):
        """ returns a list of sorted keys """
        keys = self.peak.heights.keys()
        keys.sort()
        return keys
        
    def get_possible_fragments(self):
        return self.possible_fragment_list
    
    def get_peak(self):
        return self.peak

    def get_modified_peak(self):
        return self.modified_peak

    def print_contents(self, filename=None):
        # unmodified heights lines
        u_lines = []
        u_lines.append('unmod wt %s\n'%(self.peak.weight))
        fraginfo = ""
        i = 0
        for frag in self.possible_fragment_list:
            fraginfo += 'pf%s %s %s %s\n'%(i, frag.get_nterm_index(), frag.get_cterm_index(), frag.get_weight())
            i += 1
        u_lines.append(fraginfo)
        
        keys = self.peak.heights.keys()
        keys.sort()

        longest = 0
        for key in keys:
            if len(self.peak.heights[key]) > longest:
                longest = len(self.peak.heights[key])

        for i in range(longest):
            newline = ','
            for key in keys:
                try:
                    newline += '%5.3f, '%(self.peak.heights[key][i])
                except IndexError:
                    newline += '     , '
            u_lines.append(newline + '\n')

        if self.modified_peak:
            m_lines = []
            m_lines.append('mod wt %s\n'%(self.modified_peak.weight))
            for i in range(longest):
                newline2 = ','
                for key in keys:
                    try:
                        newline2 += '%5.3f, '%(self.modified_peak.heights[key][i])
                    except IndexError:
                        newline2 += '     , '
                m_lines.append(newline2 + '\n')

        if filename:
            f = open(filename, 'a')
            f.writelines(u_lines)
            if self.modified_peak:
                f.writelines(m_lines)
            f.close()
        else:
            print u_lines
            if self.modified_peak:
                print m_lines
            print '%s possible fragment associations'%(len(self.possible_fragment_list))

class MiniViewScroller(Frame):
    def __init__(self, frame, parent, width, ht):
        Frame.__init__(self, frame, borderwidth=0, width=width, height=ht, colormap="new", visual='truecolor', bg='white', relief=FLAT)
        self.width = width
        self.height = ht
        self.parent = parent
        self.left_scroll_photo = ImageTk.PhotoImage(file="./Tools/Plotting/left_arrow.ppm")
        # create and pack a canvas, and draw the axes.
        self.scroll_left  = Button(self,
                                   relief=FLAT,
                                   image=self.left_scroll_photo,
                                   height=30,
                                   width=30,
                                   command=self.step_left,
                                   repeatdelay=150,
                                   repeatinterval=100,
                                   borderwidth=0,
                                   bg='white')
        self.right_scroll_photo=ImageTk.PhotoImage(file="./Tools/Plotting/right_arrow.ppm")
        self.scroll_right = Button(self,
                                   relief=FLAT,
                                   borderwidth=0,
                                   image=self.right_scroll_photo,
                                   height=30,
                                   width=30,
                                   command=self.step_right,
                                   repeatdelay=150,
                                   repeatinterval=100,
                                   bg='white')
        
        self.canvas = Canvas(self,
                             bg='white',
                             height=self.height,
                             width=self.width,
                             highlightbackground='white',
                             relief=FLAT,
                             borderwidth=0,
                             selectbackground='white')

        # draw tic marks every thousand, on both top and bottom
        self.tic_count = int(self.parent.max_plot.x[-1]%1000)
        self.canvas.delete('axis')
        for i in range(0, self.tic_count):
            x = int((i*1000)/self.parent.max_plot.x[-1]*self.width)
            y = 0.15*self.height
            self.canvas.create_line(x,0,x,y, fill='black', tags='axis')
            self.canvas.create_line(x,self.height,x,self.height+y, fill='black', tags='axis')

        
        self.scroll_left.bind('<Enter>', self._activate_left)
        self.scroll_left.bind('<Leave>', self._deactivate_left)
        self.scroll_right.bind('<Enter>', self._activate_right)
        self.scroll_right.bind('<Leave>', self._deactivate_right)
        self.scroll_left.pack_propagate(0)
        self.scroll_right.pack_propagate(0)
        self.scroll_left.pack(side=LEFT, expand=YES)
        self.canvas.pack(side=LEFT, expand=NO)
        self.scroll_right.pack(side=LEFT, expand=YES)
        # bind canvas clicks to the scroll function
        self.canvas.bind('<ButtonPress-1>', self.scroll)
        self.canvas.bind('<B1-Motion>', self.scroll)
        self.canvas.bind('<ButtonRelease-1>', self.scroll)
        # get some info for the panel representations
        self.peaks = {}
        self.current_x = 0
        self.last_event = 0
        class fake_event:
            def __init__(self, xval):
                self.x = xval
        event = fake_event(self.last_event)
        self.redraw_scrollbox(event)

    def update_width(self, width):
        self.canvas.config(width=width)
        class fake_event:
            def __init__(self, xval):
                self.x = xval
        event = fake_event((width/float(self.width))*self.last_event)
        self.redraw_scrollbox(event)
        self.width = width
        self.canvas.delete('axis')
        for i in range(0, self.tic_count):
            x = int((i*1000)/self.parent.max_plot.x[-1]*self.width)
            y = 0.15*self.height
            self.canvas.create_line(x,0,x,y, fill='black', tags='axis')
            self.canvas.create_line(x,self.height,x,self.height+y, fill='black', tags='axis')
        self.update()

    def _activate_left(self, event):
        self.left_scroll_photo = PhotoImage(file="./Tools/Plotting/left_arrow_selected.ppm")
        self.scroll_left.config(image=self.left_scroll_photo)

    def _activate_right(self, event):
        self.right_scroll_photo = PhotoImage(file="./Tools/Plotting/right_arrow_selected.ppm")
        self.scroll_right.config(image=self.right_scroll_photo)

    def _deactivate_left(self, event):
        self.left_scroll_photo = PhotoImage(file="./Tools/Plotting/left_arrow.ppm")
        self.scroll_left.config(image=self.left_scroll_photo)

    def _deactivate_right(self, event):
        self.right_scroll_photo = PhotoImage(file="./Tools/Plotting/right_arrow.ppm")
        self.scroll_right.config(image=self.right_scroll_photo)

    def step_left(self):
        class fake_event:
            def __init__(self, xval):
                if xval < 0:
                    xval = 0
                self.x = xval
        event = fake_event(self.current_x - (self.width/100))
        self.scroll(event)

    def step_right(self):
        class fake_event:
            def __init__(self, xval, width):
                if xval > width:
                    xval = width
                self.x = xval
        event = fake_event(self.current_x + (self.width/100), self.width)
        self.scroll(event)

    def scroll(self, event):
        self.current_x = event.x
        percent_move = (event.x/float(self.width)) - ((self.width/2.0)/self.parent.max_plot.x[-1])
        self.parent.scrollXmoveto(percent_move)
        self.redraw_scrollbox(event)

    def scroll_to(self, weight):
        class fake_event:
            def __init__(self, xval):
                self.x = xval
        event = fake_event((weight/self.parent.max_plot.x[-1])*self.width)
        self.scroll(event)

    def redraw_scrollbox(self, event):
        self.canvas.delete('cursorbox')
        self.last_event = event.x
        start = self.width * (event.x/float(self.width)) - (self.width * ((self.width/2.0)/self.parent.max_plot.x[-1]))
        end   = self.width * (event.x/float(self.width)) + (self.width * ((self.width/2.0)/self.parent.max_plot.x[-1]))
        self.canvas.create_line(start-1, 5,           start-1, self.height, fill='black', tags='cursorbox')
        self.canvas.create_line(start-1, 5,           end+1,   5,           fill='black', tags='cursorbox')
        self.canvas.create_line(end+1,   5,           end+1,   self.height, fill='black', tags='cursorbox')
        self.canvas.create_line(start-1, self.height, end+1,   self.height, fill='black', tags='cursorbox')

    def update(self):
        # delete any previous peaks drawn
        keys = ['unmodified', 'int_standards', 'modified', 'uninteresting']
        for key in keys:
            self.canvas.delete(key)

        self.peaks = self.parent.get_peaks()
        
        peaksets = {'uninteresting':'dark grey',
                    'int_standards':'orange',
                    'unmodified':'red',
                    'modified':'green'}

        for name in keys:
            for peak in self.peaks[name]:
                # get the scaled x value
                # find the ratio covered in the large spectra
                x = int(peak.get_weight()/self.parent.max_plot.x[-1]*self.width)
                # draw the line
                lin = self.canvas.create_line(x, 5,x,self.height, fill=peaksets[name], tags=name)
            self.canvas.tag_lower(name)


class PlotWindow(Frame):
    def __init__(self, parent, average_plots, background_plots, controls=None, experiment=None, possible_peaks = [], all_fragments= [], filter_resolution=None, width=720, height=300):
        """ MSWindow acts as a holder for multiple PlotPanel objects and coordinates their activity.
            It also offers tools to integrate controls
        """
        Frame.__init__(self, parent, borderwidth=0, height=height, width=width)
        self.possible_peaks = possible_peaks
        self.all_fragment_weights = all_fragments
        print '%s possible peaks, of %s fragments'%(len(possible_peaks), len(all_fragments))
        if filter_resolution == None:
            self.filter_resolution = 0.005
        else:
            self.filter_resolution = filter_resolution
        self.experiment = experiment
        #hScaleFrame = parent.hScaleFrame
        if experiment:
            self.reaction_type = experiment.get_modifying_reagent()
        else:
            self.reaction_type = None
        self.width=width
        self.height=height
        self.parent=parent
        self.config(height=self.height, width=self.width, bg='white')
        labelfont = ('helvetica', 9, 'bold')
        self.plots = average_plots

        self.max_plot = PlottingTools.max_of_plots(self.plots)

        self.background_plots = background_plots
        self.modified_peak_indices  = []
        self.unmodified_peak_indices = []
        self.fast_plot_type = 'areas'  # decides whether to use areas or heights for quanitfication

        self.spectra_frame = Frame(self)
        
        # heightScaler scales the heights of all ms plots
        self.hscaleframe = Frame(self.spectra_frame)
        Label(self.hscaleframe, text='scale\nall', width=5, font=labelfont).pack(side=TOP)
        self.heightScaler = Scale(self.hscaleframe, from_=10.0, to=0.2, resolution=0.2, orient='vertical', width=15)
        self.heightScaler.set(1.0)
        self.heightScaler.pack(side=TOP,fill=Y, expand=1)
        # fix the callback handler
        self.heightScaler.config(command=self.scaleY)
        self.hscaleframe.pack(side=LEFT, expand=NO, fill=Y)
        Label(self.hscaleframe, text='             ').pack(side=LEFT)
        self.plotCount = len(self.plots)

        # this active message box is below the control panels and shows info on
        # the fragment tags below the x-axis. Look below for the rest of the
        # active messages.
        self.active_message_frame = Frame(self.spectra_frame)
        
        self.position_message = ""
        self.position_text      = StringVar()
        self.position_text.set(self.position_message)
        self.position_label     = Label(self.active_message_frame, textvariable=self.position_text)
        self.position_label.pack(side='left', anchor='w', expand=0, fill='x')

        self.active_peak_message = ""
        self.active_peak_text      = StringVar()
        self.active_peak_text.set(self.active_peak_message)
        self.active_peak_label     = Label(self.active_message_frame, textvariable=self.active_peak_text)
        self.active_peak_label.pack(side='left', anchor='e', expand=0, fill='x')

        self.active_message_frame.pack(side='bottom', expand=0, fill='x')        

        self.plot_panels = []
        self.plots_frame = Frame(self.spectra_frame, height=self.height)
        # add an x-scrollbar that handles all of the MSPanels
        self.xscrollframe = Frame(self.spectra_frame, relief=RAISED, bd=2)
        simple_scroll=0
        if simple_scroll:
            self.xsbar=Scrollbar(self.xscrollframe)
            self.xsbar.config(orient='horizontal', command=self.scrollX)
        else:
            self.xsbar = MiniViewScroller(self.xscrollframe, self, 600, 30)

        self.plots_frame.pack(side=TOP, expand=YES, fill=BOTH)
        self.information_frame = Frame(self, bg='white')
        
        # add a y-axis scroll bar
        self.ysbar=Scrollbar(self.plots_frame)
        self.ysbar.pack(side=RIGHT, fill=Y)
        self.ysbar.config(command=self.scrollY, cursor='arrow')

        # create a scrollable canvas, the same height as the combined set of plots
        self.plots_text             = Text(self.plots_frame, height=4, yscrollcommand=self.ysbar.set)
        self.plots_text.parent      = self
        self.observations_info_text = Text(self.information_frame, height=4, width=17, yscrollcommand=self.ysbar.set, bg='white')

        # and add the mass spec panels
        self.plot_panels = []
        self.info_canvases = []
        labelFont = ('Helvetica', parms.get('label_font_size'))
        # create the plot panels
        for i in range(len(self.plots)):
            avg_plot = self.plots[i]
            bkg_plots = self.background_plots[i]
            panel = PlotPanel(self.plots_text, avg_plot, bkg_plots, possible_peaks, self.filter_resolution)
            self.plot_panels.append(panel)
            canvas = Canvas(self.observations_info_text, height=100, width=100, bg='white')
            self.info_canvases.append(canvas)
            canvas.plot_panel = panel

        # insert them into the canvas
        i = 0
        for panel in self.plot_panels:
            panel.pack(expand=YES, fill=BOTH, side=TOP)
            self.plots_text.window_create('end', window=panel)
            if simple_scroll:
                panel.plot.config(xscrollcommand=self.xsbar.set, state=DISABLED)
            i += 1
            if (i < self.plotCount-1):
                self.plots_text.insert('end', '\n')

        # create information panels for each spectrum
        i = 0
        for canvas in self.info_canvases:
            canvas.pack(expand=NO, fill=NONE, side=TOP)
            self.observations_info_text.window_create('end', window=canvas)
            if (int(self.plots[i].get_label())/60.0)%1.0 == 0:
                label = '%s min.'%(int(int(self.plots[i].get_label())/60.0))
            else:
                label = '%s sec.'%(int(self.plots[i].get_label()))
            widget = Label(canvas, text='%s'%(label), bg='white', fg='black', font=labelFont, width=11+len(label))
            widget.pack(side=TOP, expand=NO)
            canvas.create_window(50, 15, window=widget, tags='names')
            # add an entry box for independent rescaling
            canvas.view_all_frame = Frame(canvas, bg='white')
            canvas.view_all_entry = Pmw.EntryField(canvas.view_all_frame,
                                                   labelpos='w',
                                                   label_text='Ratio:',
                                                   label_width=6,
                                                   label_height=1,
                                                   label_bg='white',
                                                   label_font=labelFont,
                                                   entry_width=5,
                                                   entry_bg='white',
                                                   validate = None,
                                                   hull_bg='white',
                                                   command=None)
            canvas.view_all_entry.setvalue("100.0")
            canvas.view_all_entry.pack(side=LEFT)
            c_lambda = lambda x=self.plot_panels[i], y=canvas.view_all_entry: x.hard_rescale_y(float(y.getvalue()))
            canvas.view_all_button = Button(canvas.view_all_frame,
                                            text='Apply',
                                            width=5,
                                            height=1,
                                            font=labelFont,
                                            command=c_lambda)
            canvas.view_all_button.pack(side=LEFT)
            canvas.view_all_frame.pack(side=TOP)
            canvas.create_window(50, 40, window=canvas.view_all_frame, tags='view_all')

            # a radioselect is here used to visualize the average or all spectra
            d_lambda = lambda x, y, z=self.plot_panels[i]: self.adjust_spectra_view(z, x, y)
            canvas.checkbuttons = Pmw.RadioSelect(canvas,
                                                  buttontype = 'checkbutton',
                                                  orient = 'horizontal',
                                                  hull_bg='white',
                                                  command = d_lambda)
            canvas.checkbuttons.pack(side = 'top', expand=0, fill=NONE)
            for text in ('avg', 'all'):
                canvas.checkbuttons.add(text)
            # this should check without invoking
            canvas.checkbuttons.setvalue(['avg'])
            canvas.checkbuttons.button('avg').config(font=labelFont, bg='white')
            canvas.checkbuttons.button('all').config(font=labelFont, bg='white')
            canvas.create_window(50, 61, window=canvas.checkbuttons, tags='view_all')

            # a button launches a popup for deselecting/reselecting multiple observations
            e_lambda = lambda z=self.plot_panels[i]: self.selection_on_multiple_observations(z)
            canvas.selection_button = Button(canvas, text='Observation Selection', font=labelFont, command=e_lambda)
            canvas.selection_button.pack(side='top', expand=0, fill=None)
            canvas.create_window(50, 85, window=canvas.selection_button, tags='selection_buttons')
            i += 1
            if (i < self.plotCount):
                self.observations_info_text.insert('end', '\n')

        # create the initial view of the panels
        self.draw_PlotPanels()

        # draw the controls panel
        if self.experiment.has_controls():
            self.controls_frame = Frame(self.spectra_frame, height=self.height, relief='flat')
            self.controls_text  = Text(self.controls_frame, height=8)
            self.controls_panel = PlotPanel(self.controls_text, None, self.experiment.get_control_plots())
            p = self.controls_panel
            p.pack(expand=YES, fill=BOTH, side=TOP)
            self.controls_text.window_create('end', window=p)
            if simple_scroll:
                p.plot.config(xscrollcommand=self.xsbar.set)
            self.controls_text.config(state=DISABLED)
            self.controls_panel.draw_coordinates()

        self.fast_plot_canvas = Canvas(self.information_frame, height=70, width=100, bg='white')
        # place some text with the time point
        self.spacer_frame3 = Frame(self.information_frame, height=60, width=100, background='white')
        self.spacer_frame1 = Frame(self.information_frame, height=35, width=100, background='white', relief='raised', bd=4)

        activeLabelFont = ('Helvetica', 7)

        self.active_class_message = "%-s"%("")
        self.active_class_text    = StringVar()
        self.active_class_text.set(self.active_class_message)
        self.active_class_label     = Label(self.spacer_frame3, height=1, textvariable=self.active_class_text, background='white', font=activeLabelFont)
        self.active_class_label.pack(side='top', anchor='w', expand=0, fill='none')

        self.active_weight_message = "%-s"%("")
        self.active_weight_text    = StringVar()
        self.active_weight_text.set(self.active_weight_message)
        self.active_weight_label     = Label(self.spacer_frame3, height=1, textvariable=self.active_weight_text, background='white', font=activeLabelFont)
        self.active_weight_label.pack(side='top', anchor='w', expand=0, fill='none')

        self.active_intensity_message = "%-s"%("")
        self.active_intensity_text    = StringVar()
        self.active_intensity_text.set(self.active_intensity_message)
        self.active_intensity_label     = Label(self.spacer_frame3, height=1, textvariable=self.active_intensity_text, background='white', font=activeLabelFont)
        self.active_intensity_label.pack(side='top', anchor='w', expand=0, fill='none')

        self.plots_text.config(state=DISABLED)

        self.xsbar.pack(side=TOP, fill=BOTH, expand=NO)
        self.xscrollframe.pack(side=TOP, fill=X, expand=NO)

        self.observations_info_text.config(state=DISABLED)
        self.plots_text.pack(side=TOP, expand=YES, fill=BOTH)

        self.information_frame.pack(side=RIGHT, expand=NO, fill=Y, anchor=N)
        self.spectra_frame.pack(side=RIGHT, expand=YES, fill=BOTH, anchor=N)

        self.observations_info_text.pack(side=TOP, expand=YES, fill=Y)
        self.spacer_frame1.pack(side=TOP, anchor='w', expand=NO, fill='x')

        if not self.experiment.has_controls():
            print "what do you mean you don't have any controls?"

        self.spacer_frame3.pack(side=TOP, anchor='w', expand=NO, fill='x')
        if self.experiment.has_controls():
            self.controls_text.pack(expand=YES, fill=BOTH)
            self.controls_frame.pack(side=TOP, expand=NO, fill=X)
        self.fast_plot_canvas.pack(side=TOP, expand=NO, fill=NONE, anchor=SW)

        self.internal_standards = []
        self.reaction_profiles = []
        self.peak_profiles = {}

        self.blocked_internal_standard_sequences = []

        self.centroids_on = 0
        self.unrationalized_peak_indices = []

    def get_mass_range(self):
        return self.plot_panels[0].get_mass_range()

    def set_mass_range(self, parms):
        for panel in self.plot_panels:
            panel.set_mass_range(parms)

    def get_filter_resolution(self):
        return self.filter_resolution

    def set_filter_resolution(self, new_resolution):
        self.filter_resolution = new_resolution
        for panel in self.plot_panels:
            panel.set_filter_resolution(new_resolution)

    def toggle_centroids(self):
        if self.centroids_on:
            for plot in self.plot_panels:
                plot.undraw_centroids()
            self.centroids_on = 0
        else:
            for plot in self.plot_panels:
                plot.draw_centroids()
            self.centroids_on = 1
        

    def get_peaks(self):
        return self.peak_profiles

    def get_reactions(self):
        return self.reaction_profiles

    # selection is initiated in a Plot object, and broadcast to all plots here    
    def _broadcast_select_range_start(self, event):
        for panel in self.plot_panels:
            panel._draw_select_range_start(event)

    def _broadcast_select_range_end(self, event):
        self.active_select_range_end_index = event.x
        for panel in self.plot_panels:
            panel._draw_select_range_end(event)

    def _broadcast_select_range_motion(self, event):
        for panel in self.plot_panels:
            panel._draw_select_range_motion(event.x)

    # extend selection (Ctrl + button-1)
    def _broadcast_extend_range_start(self, event):
        for panel in self.plot_panels:
            panel._draw_extend_range_start(event)

    def _broadcast_extend_range_end(self, event):
        for panel in self.plot_panels:
            panel._draw_extend_range_end(event)

    def _broadcast_extend_range_motion(self, event):
        for panel in self.plot_panels:
            panel._draw_extend_range_motion(event.x)

    def _fastplot_max_in_selected_range(self):
        vals = []
        for panel in self.plot_panels:
            val = panel.get_area_under_selected_regions()
            if val > 0:
                vals.append(panel.get_area_under_selected_regions())
            else:
                vals.append(0.01)
        self._draw_fast_plot(vals)
                  
    def _broadcast_recalculate_baseline(self):
        # first recalculate the baseline for the max_plot
        loq = self.max_plot._recalculate_baseline(self.plot_panels[0].get_last_selection_indices())
        # now the actual
        loqsum = 0.0
        for panel in self.plot_panels:
            loqsum += panel.recalculate_baseline()
        # loq is the average loq from the normal plots
        self.max_plot.recognize_peaks(self.experiment.get_mass_range()[0], loqsum/float(len(self.plot_panels)))
        print 'max plot loq %s gives %s peaks'%(loqsum/float(len(self.plot_panels)), len(self.max_plot.peaks))
        for peak in self.max_plot.peaks:
            print self.max_plot.x[peak],
        self.draw_peaks()

    def _clear_selection(self):
        for panel in self.plot_panels:
            panel._clear_selection()

    def _fastplot_heights(self):
        self.fast_plot_type = 'heights'

    def _fastplot_areas(self):
        self.fast_plot_type = 'areas'

    def _undraw_fast_plot(self):
        self.fast_plot_canvas.delete('plot')

    def _draw_fast_plot(self, heights):
        self.fast_plot_canvas.delete('plot')
        # find the highest point
        highest = 0.0
        for height in heights:
            if height > highest:
                highest = height
            
        self.fast_plot_canvas.create_line(3,3,3,67, tags='plot')
        self.fast_plot_canvas.create_line(3,67,67,67, tags='plot')
        
        increment = 100.0/len(heights)

        i = 0
        for i in range(len(heights)-1):
            x1 = 5+round(i*increment)
            x2 = 5+round((i+1)*increment)
            h1 = 65 - round((heights[i]/float(highest))*60.0)
            h2 = 65 - round((heights[i+1]/float(highest))*60.0)
            self.fast_plot_canvas.create_line(x1, h1, x2, h2, tags='plot')

    def selection_on_multiple_observations(self, plot_panel):
        c_lambda = lambda x, y=plot_panel: self.execute_selection_on_multiple_observations(x, y)
        self.observations_selection_dialog = Pmw.Dialog(self,
                                                        master=self,
                                                        buttons = ('Ok', 'Cancel'),
                                                        defaultbutton = 'OK',
                                                        title = 'Select a subset of spectra',
                                                        command = c_lambda)
        
        self.observations_selection_dialog.radioselect = Pmw.RadioSelect(self.observations_selection_dialog.interior(), 
                                                                         buttontype = 'checkbutton',
                                                                         orient = 'horizontal',
                                                                         hull_bg='white')
        
        self.observations_selection_dialog.radioselect.pack(side = 'top', expand=1, fill=NONE)

        colorset = ['blue', 'red', 'green', 'orange', 'yellow', 'purple']
        bg_specs = plot_panel.background_plots
        bg_specs_selection = plot_panel.background_plots_selected
        i = 0
        to_select = []
        for i in range(len(bg_specs)):
            self.observations_selection_dialog.radioselect.add('Obs. %s'%(i))
            self.observations_selection_dialog.radioselect.button('Obs. %s'%(i)).config(bg='white', fg=colorset[i%len(bg_specs)])
            if bg_specs_selection[i]:
                to_select.append('Obs. %s'%(i))
        self.observations_selection_dialog.radioselect.setvalue(to_select)
        
    def execute_selection_on_multiple_observations(self, result, plot_panel):
        if result == 'Cancel':
            self.observations_selection_dialog.destroy()
            return
        selected = self.observations_selection_dialog.radioselect.getvalue()
        for i in range(len(plot_panel.background_plots)):
            for selection in selected:
                if string.atoi(string.split(selection)[-1]) == i:
                    plot_panel.background_plots_selected[i] = 1
                    break
            else:
                plot_panel.background_plots_selected[i] = 0
        selected_specs = []
        for i in range(len(plot_panel.background_plots)):
            if plot_panel.background_plots_selected[i]:
                selected_specs.append(plot_panel.background_plots[i])
        average_plot = PlottingTools.average_plots(selected_specs)
        plot_panel.plot_object = average_plot
        plot_panel.x = average_plot.x
        plot_panel.y = average_plot.y
        for canvas in self.info_canvases:
            if canvas.plot_panel == plot_panel:
                vals = canvas.checkbuttons.getvalue()
                if 'avg' in vals:
                    plot_panel.draw_coordinates()
                else:
                    plot_panel.undraw_coordinates()
                if 'all' in vals:
                    plot_panel.draw_background_coordinates()
                else:
                    plot_panel.undraw_background_coordinates()
        self.observations_selection_dialog.destroy()

    def adjust_spectra_view(self, plot, key, on_off):
        if key == 'avg':
            if on_off:
                plot.draw_coordinates()
            else:
                plot.undraw_coordinates()
        elif key == 'all':
            if on_off:
                plot.draw_background_coordinates()
            else:
                plot.undraw_background_coordinates()

    def set_x_scale(self, val):
        for panel in self.plot_panels:
            panel.set_x_scale(val)
        
    def calculate_peaks(self):
        protease_specificity = self.experiment.get_protease_specificity()
        reagent_specificity  = self.experiment.get_modifying_reagent_specificity()
        digestion_type       = self.experiment.get_digestion_type()

        # three types of peak recognition:
        # 1. chem mod with reagent specificity in protease specificity
        # 2. chem mod without ""
        # 3. limited proteolysis

        if digestion_type == 'limited':
            self.calculate_limited_proteolysis_peaks()
        else:
            for protease_target_res in protease_specificity:
                if protease_target_res in reagent_specificity:
                    self.calculate_modified_reagent_peaks()
                    break
            else:
                self.calculate_unmodified_reagent_peaks()
                
    def calculate_limited_proteolysis_peaks(self):
        possible_peaks = []
        unmodified_peak_indices = []

        for avg_peak in self.max_plot.peaks:
            for possible_peak in self.possible_peaks:
                if abs(self.max_plot.x[avg_peak]-possible_peak) < (self.filter_resolution*possible_peak):
                    if possible_peak not in possible_peaks:
                        possible_peaks.append(possible_peak)
                        unmodified_peak_indices.append(avg_peak)
                elif abs(self.max_plot.x[avg_peak]-possible_peak) < 2*self.filter_resolution*possible_peak:
                    print 'fragment weight %s nearly matches peak %s but the difference (%s) is greater than %s'%(possible_peak, self.max_plot.x[avg_peak], abs(self.max_plot.x[avg_peak]-possible_peak), self.filter_resolution*possible_peak)
        self._disambiguate_and_collect_indices(possible_peaks, unmodified_peak_indices)
    
    def calculate_modified_reagent_peaks(self):
        # look for unmodified peaks in the zero time point (i.e. NHS with V8)
        # given single-hit possible_peaks from the MassSpecExperiment, identify
        # those present in the zero timepoint
        print 'collecting single-reaction peaks from the modified versions'
        temp_unmodified_peaks = []
        temp_unmodified_peak_indices = []

        reagent_weight = float(self.experiment.get_modifying_reagent_weight())
        reagent_specificity = self.experiment.get_modifying_reagent_specificity()

        # collect modified peaks
        #target_plot = self.plots[-1]
        target_plot = self.max_plot
        for avg_peak in target_plot.peaks:
            for possible_peak in self.possible_peaks:
                # recognize by weight + reagent_weight
                if abs(target_plot.x[avg_peak]-(possible_peak+reagent_weight)) < (self.filter_resolution*(possible_peak+reagent_weight)):
                    if possible_peak not in temp_unmodified_peaks:
                        # collect the weights of unmodified peaks
                        temp_unmodified_peaks.append(possible_peak)

        unmodified_peaks = []
        unmodified_peak_indices = []
        # get their indices
        for i in range(len(temp_unmodified_peaks)):
            unmodified_peak_weight = temp_unmodified_peaks[i]
            for x_ind in range(len(target_plot.x)):
                if target_plot.x[x_ind] > unmodified_peak_weight:
                    if abs(target_plot.x[x_ind-1]-unmodified_peak_weight) < abs(target_plot.x[x_ind]-unmodified_peak_weight):
                        unmodified_peak_indices.append(x_ind-1)
                        unmodified_peaks.append(unmodified_peak_weight)
                    else:
                        unmodified_peak_indices.append(x_ind)
                        unmodified_peaks.append(unmodified_peak_weight)
                    break
            
        self._disambiguate_and_collect_indices(unmodified_peaks, unmodified_peak_indices)

    def calculate_unmodified_reagent_peaks(self):
        # look for unmodified peaks in the zero time point (i.e. NHS with V8)
        # given single-hit possible_peaks from the MassSpecExperiment, identify
        # those present in the zero timepoint
        print 'collecting single-reaction peaks from the unmodified versions'

        unmodified_peaks = []
        modified_peaks   = []
        unmodified_peak_indices = []
        modified_peak_indices = []
        
        loq = self.max_plot._recalculate_baseline(self.plot_panels[0].get_last_selection_indices())
        self.max_plot.recognize_peaks(self.experiment.get_mass_range()[0], loq)
        for avg_peak in self.max_plot.peaks:
            for possible_peak in self.possible_peaks:
                if abs(self.plots[0].x[avg_peak]-possible_peak) < (self.filter_resolution*possible_peak):
                    if possible_peak in unmodified_peaks:
                        for i in range(len(unmodified_peaks)):
                            if unmodified_peaks[i] == possible_peak:
                                if abs(self.plots[0].x[avg_peak]-possible_peak) < abs(self.plots[0].x[unmodified_peak_indices[i]]-possible_peak):
                                    unmodified_peak_indices[i] = avg_peak
                    else:
                        unmodified_peaks.append(possible_peak)
                        unmodified_peak_indices.append(avg_peak)

        self._disambiguate_and_collect_indices(unmodified_peaks, unmodified_peak_indices)

    def _disambiguate_and_collect_indices(self, possible_peaks, unmodified_peak_indices):
        reagent_weight = self.experiment.get_modifying_reagent_weight()
        reagent_specificity = self.experiment.get_modifying_reagent_specificity()
        print 'Collecting peaks for reagent weight %s that hits %s'%(reagent_weight, reagent_specificity)
        print 'disambiguating peaks at',
        for ind in possible_peaks:
            print ind,
        print
        digestion_type = self.experiment.get_digestion_type()
        modified_peaks   = []
        modified_peak_indices = []
        # get all of the possible modified fragments
        i = 0
        new_unmodified_peaks = []
        new_unmodified_peak_indices = []
        for possible_peak in possible_peaks:
            results = self.experiment.get_fragment_for_weight_from_all(possible_peak)
            # get the number of apparent modifications by comparing the first
            # and last time point.
            """
            if len(results) > 1:
                # see how many peaks correspond to reagent multiplets in the last
                # time point, with an arbitrary max of 10
                reaction_peak_count = 0
                for j in range(10):
                    breakout = 0
                    for last_peak_index in self.plots[-1].peaks:
                        last_peak = self.plots[-1].x[last_peak_index]
                        # if it is present in the last plot
                        if abs(last_peak-possible_peak+(j*reagent_weight)) < (self.filter_resolution*(possible_peak+(j*reagent_weight))):
                            # but not present in the first plot
                            for first_peak_index in self.plots[0].peaks:
                                first_peak = self.plots[0].x[first_peak_index]
                                if abs(last_peak-possible_peak+(j*reagent_weight)) < (self.filter_resolution*(possible_peak+(j*reagent_weight))):
                                    # it is present. breakout
                                    breakout=1
                                    break
                            else:
                                reaction_peak_count += 1
                        else:
                            breakout = 1
                            break
                        if breakout:
                            break
                    if breakout:
                        break
                # if there are more than one fragment results, disambiguate
                realistic_results = []
                for result in results:
                    sequence = result['sequence']
                    reaction_site_count = 0
                    for aa in sequence:
                        if aa in reagent_specificity:
                            reaction_site_count += 1
                    if reaction_peak_count <= reaction_site_count:
                        realistic_results.append(result)
                results = realistic_results
            """
            if len(results) > 1 or len(results) == 0:
                if len(results) > 1:
                    print '%s fragments in mass range for peak %s'%(len(results), possible_peak)
                elif len(results) == 0:
                    print 'no fragments in mass range for peak %s'%(possible_peak)
                                                                
                i += 1
                continue
            # count the number of reactive residues in the fragment
            sequence = results[0]['sequence']
            reaction_hits = 0
            r_wt = 0.0
            for aa in sequence[:-1]:
                reagent = self.experiment.get_modifying_reagent()
                if reagent and reagent != 'None':
                    if aa in reagent_specificity:
                        reaction_hits += 1
                        r_wt = string.atof(reagent_weight) * reaction_hits
            # store and create the modified counterpart
            if reaction_hits == 1 or digestion_type == 'limited':
                new_unmodified_peaks.append(possible_peaks[i])
                new_unmodified_peak_indices.append(unmodified_peak_indices[i])
                if digestion_type != 'limited':
                    modified_peaks.append(possible_peak+r_wt)
            i += 1
        possible_peaks = new_unmodified_peaks
        unmodified_peak_indices = new_unmodified_peak_indices

        self.unmodified_peak_indices = []
        self.modified_peak_indices = []
        if digestion_type == 'exhaustive':
            for i in range(len(modified_peaks)):
                # first figure out the appropriate index to associate with the peak
                modified_peak_weight = modified_peaks[i]
                modified_peak_index = 0
                for x_ind in range(len(self.plots[-1].x)):
                    if self.plots[-1].x[x_ind] >= modified_peak_weight:
                        if abs(self.plots[-1].x[x_ind-1]-modified_peak_weight) < abs(self.plots[-1].x[x_ind]-modified_peak_weight):
                            modified_peak_index = x_ind-1
                        else:
                            modified_peak_index = x_ind
                        break
                
                # find the nearest peak to that index in the max_plot
                closest = -1
                dist    = 10000
                for max_peak_index in self.max_plot.peaks:
                    d = abs(self.max_plot.x[max_peak_index]-self.max_plot.x[modified_peak_index]) 
                    if d < dist:
                        dist = d
                        closest = max_peak_index
                # make sure it is close enough
                print 'closest peak to %s at %s'%(modified_peak_weight, self.max_plot.x[closest]),
                if abs(modified_peak_weight - self.max_plot.x[closest]) > (self.filter_resolution)*modified_peak_weight:
                    print ' -  distance %s too big, skipping'%(abs(modified_peak_weight - self.max_plot.x[closest]))
                    continue
                else:
                    print
                """
                    # make sure the peak is not present in the first plot
                    closest = -1
                    dist    = 10000
                    for max_peak_index in self.plots[0].peaks:
                        d = abs(self.plots[0].x[max_peak_index]-self.plots[0].x[modified_peak_index]) 
                        if d < dist:
                            dist = d
                            closest = max_peak_index
                    #if abs(modified_peak_weight - self.plots[0].x[closest]) < (self.filter_resolution)*modified_peak_weight:
                    #    if self.plots[0].y[closest] > 5.0:
                    #        print '%s first plot peak (ht %s) too close (dst %s), skipping'%(modified_peak_weight, self.plots[0].y[closest], abs(modified_peak_weight - self.plots[0].x[closest]))
                    #        continue
                """
                self.modified_peak_indices.append(modified_peak_index)
                self.unmodified_peak_indices.append(unmodified_peak_indices[i])

        else:
            for i in range(len(possible_peaks)):
                # first figure out the appropriate index to associate with the peak
                possible_peak_weight = possible_peaks[i]
                possible_peak_index = 0
                for x_ind in range(len(self.max_plot.x)):
                    if self.max_plot.x[x_ind] >= possible_peak_weight:
                        if abs(self.max_plot.x[x_ind-1]-possible_peak_weight) < abs(self.max_plot.x[x_ind]-possible_peak_weight):
                            possible_peak_index = x_ind-1
                        else:
                            possible_peak_index = x_ind
                        break
                
                # find the nearest peak to that index in the max_plot
                closest = -1
                dist    = 10000
                for max_peak_index in self.max_plot.peaks:
                    d = abs(self.max_plot.x[max_peak_index]-self.max_plot.x[possible_peak_index]) 
                    if d < dist:
                        dist = d
                        closest = max_peak_index
                        
                # make sure it is close enough
                if abs(possible_peak_weight - self.max_plot.x[closest]) > (self.filter_resolution)*possible_peak_weight:
                    print '%s max plot distance %s > %s , skipping'%(possible_peak_weight, abs(possible_peak_weight - self.max_plot.x[closest]), (self.filter_resolution)*possible_peak_weight)
                    continue
            
                # make sure the peak is not present in the first plot
                closest = -1
                dist    = 10000
                for max_peak_index in self.plots[0].peaks:
                    d = abs(self.plots[0].x[max_peak_index]-self.plots[0].x[possible_peak_index]) 
                    if d < dist:
                        dist = d
                        closest = max_peak_index
                if abs(possible_peak_weight - self.plots[0].x[closest]) < (self.filter_resolution)*possible_peak_weight:
                    if self.plots[0].y[closest] > 2.0:
                        print '%s first plot peak (ht %s) too close (dst %s), skipping'%(possible_peak_weight, self.plots[0].y[closest], abs(possible_peak_weight - self.plots[0].x[closest]))
                        continue
                self.unmodified_peak_indices.append(possible_peak_index)

        if len(self.unmodified_peak_indices) > 0:
            print 'unmodified'
            for index in self.unmodified_peak_indices:
                print '%s, '%(self.plots[0].x[index]),
            print '\n'
        if len(self.modified_peak_indices) > 0:
            print 'modified'
            for index in self.modified_peak_indices:
                print '%s, '%(self.plots[0].x[index]),
            print '\n'

    def remove_ambiguous_indices(self, indices):
        # this function locates ambiguous peak assignments from a list of indices
        # and returns the set of unambiguous peaks
        indices_to_remove = []
        for i in range(len(indices)-1):
            ind1 = indices[i]
            # first just remove any internal standards that are
            # close enough to other interanal standards to consider ambiguous
            for j in range(i+1,len(indices)):
                ind2 = indices[j]
                trough1a,peak1,trough1b = self.max_plot.get_trough_indices_for_peak_index(ind1)
                trough2a,peak2,trough2b = self.max_plot.get_trough_indices_for_peak_index(ind2)
                # if they're too close, remove the more distant
                if abs(self.max_plot.x[ind1]-self.max_plot.x[ind2]) < self.filter_resolution*self.max_plot.x[ind2]:
                    if self.max_plot.y[ind1] < self.max_plot.y[ind2]:
                        if ind1 not in indices_to_remove:
                            indices_to_remove.append(ind1)
                    else:
                        if ind2 not in indices_to_remove:
                            indices_to_remove.append(ind2)
                else:
                    # if they occupy the same turf, get the highest peak and see which is
                    # closer to it
                    if trough1a == trough2a and trough1b == trough2b:
                        if self.max_plot.y[peak1] > self.max_plot.y[peak2]:
                            mainpeak = peak1
                        else:
                            mainpeak = peak2
                        if abs(self.max_plot.x[mainpeak]-self.max_plot.x[ind1]) < abs(self.max_plot.x[mainpeak]-self.max_plot.x[ind2]):
                            indices_to_remove.append(ind2)
                        else:
                            indices_to_remove.append(ind1)

        new_indices = []
        for ind in indices:
            if ind not in indices_to_remove:
                new_indices.append(ind)
                
        return new_indices, indices_to_remove
            
    def _create_peak_profiles(self):
        self.peak_profiles = {}
        names = ['uninteresting',
                 'int_standards',
                 'unmodified',
                 'modified']
        for name in names:
            self.peak_profiles[name] = []

        peaklists = {'uninteresting': self.uninteresting_indices,
                     'int_standards': self.internal_standard_indices,
                     'unmodified': self.unmodified_peak_indices,
                     'modified': self.modified_peak_indices}

        umod_indices = self.unmodified_peak_indices
        umod_weights = self.find_weights_for_indices(umod_indices)

        keys = peaklists.keys()
        keys.sort()
        subtract_baseline = self.experiment.get_subtract_baseline_from_areas()
        store_a, store_b, store_c = [], [],[]
        for key in keys:
            data = peaklists[key]
            for i in range(len(data)):
                ind = data[i]
                # first determine the highest intensity example of the peak, through all background
                # plots. This will later be used to assess negative noise with a eg 10% rule. The
                # following two blocks of code both use the assumption that a diagonal baseline removal
                # is optimal.
                highest = 0.00
                for plot_panel in self.plot_panels:
                    a,b,c = plot_panel.get_trough_indices_for_peak_index(ind)
                    label = plot_panel.plot_object.get_label()
                    if label == None:
                        continue
                    for j in range(len(plot_panel.background_plots)):
                        # subtract the angle between the two troughs to remove any baseline
                        plt = plot_panel.background_plots[j]
                        subtract_range = plt.y[c] - plt.y[a]
                        if plot_panel.background_plots_selected[j]:
                            for indx in range(a, c+1):
                                if float(c-a) == 0:
                                    print 'here, %s %s %s (%s %s %s)'%(a, b, c, plt.x[a], plt.x[b], plt.x[c])
                                val = plt.y[indx] - (plt.y[a] + (((indx-a)/float(c-a))*subtract_range))
                                if val > highest:
                                    highest = val
                # first pass through the data. if any values are negative, subtract a constant
                # value (the lower of the two troughs), instead of the diagonal between the two
                # troughs. this is done to eliminate problems observed with peaks that appear
                # concave netagive when compared to a line between the troughs. Note that this shouldn't
                # actually ever happen. If it does, there's a problem with the trough boundary definitions,
                # but this bit of code adds robustness in that event.
                use_straight = 0
                for plot_panel in self.plot_panels:
                    a,b,c = plot_panel.get_trough_indices_for_peak_index(ind)
                    label = plot_panel.plot_object.get_label()
                    if label == None:
                        continue
                    for j in range(len(plot_panel.background_plots)):
                        # subtract the angle between the two troughs to remove any baseline
                        plt = plot_panel.background_plots[j]
                        subtract_range = plt.y[c] - plt.y[a]
                        if plot_panel.background_plots_selected[j]:
                            for indx in range(a, c+1):
                                val = plt.y[indx] - (plt.y[a] + (((indx-a)/float(c-a))*subtract_range))
                                if val < 0 - (0.1 * highest):
                                    use_straight = 1
                                    print 'using non-diagonal baseline for range %s to %s - %s is > 10 percent of %s'%(plt.x[a], plt.x[c], val, highest)
                                    break
                        if use_straight:
                            break
                    if use_straight:
                        break
                use_straight = 0
                weight = self.max_plot.x[ind]
                heights = {}
                weights = {}
                areas   = {}
                for plot_panel in self.plot_panels:
                    a,b,c = plot_panel.get_trough_indices_for_peak_index(ind)
                    store_a.append(a)
                    store_b.append(b)
                    store_c.append(c)
                    label = plot_panel.plot_object.get_label()
                    if label == None:
                        continue
                    heights[label] = []
                    weights[label] = []
                    areas[label]   = []
                    for j in range(len(plot_panel.background_plots)):
                        # subtract the angle between the two troughs to remove any baseline
                        plt = plot_panel.background_plots[j]
                        subtract_range = plt.y[c] - plt.y[a]
                        if plot_panel.background_plots_selected[j]:
                            heights[label].append(plt.y[b])
                            weights[label].append(plt.x[b])
                            area = 0.0
                            for indx in range(a, c+1):
                                if subtract_baseline:
                                    if use_straight:
                                        if plt.y[a] > plt.y[c]:
                                            if plt.y[c] > 0:
                                                val = plt.y[indx] - plt.y[c]
                                                if val > 0:
                                                    area += val
                                            else:
                                                val = plt.y[indx]
                                                if val > 0:
                                                    area += val
                                        else:
                                            if plt.y[a] > 0:
                                                val = plt.y[indx] - plt.y[a]
                                                if val > 0:
                                                    area += plt.y[indx] - plt.y[a]
                                            else:
                                                val = plt.y[indx]
                                                if val > 0:
                                                    area += val
                                    else:
                                        val = plt.y[indx] - (plt.y[a] + (((indx-a)/float(c-a))*subtract_range))
                                        if val > 0:
                                            area += val
                                else:
                                    if plt.y[indx] > 0:
                                        area += plt.y[indx]
                            areas[label].append(area)

                indices_store = [store_a, store_b, store_c]
                u = PeakProfile(weight, heights, areas, indices_store)
                self.peak_profiles[key].append(u)

        for i in range(len(self.unmodified_peak_indices)):
            # collect any associable fragment objects
            frags = []
            u_weight = self.plots[-1].x[self.unmodified_peak_indices[i]]
            for fragment in self.experiment.fragment_objects:
                frag_weight = fragment.get_weight()
                if abs(frag_weight-u_weight) < (self.filter_resolution*u_weight):
                    # create the reaction profile
                    frags.append(fragment)
            u = self.peak_profiles['unmodified'][i]
            if len(self.unmodified_peak_indices) == len(self.modified_peak_indices):
                m = self.peak_profiles['modified'][i]
            else:
                m = None
            r = ReactionProfile(frags, u, m)
            self.reaction_profiles.append(r)
            
    def print_reaction_profiles(self, filename=None):
        for reaction_profile in self.reaction_profiles:
            reaction_profile.print_contents(filename)

    def draw_modified_peaks_for_peak_weight(self, wt):
        i = wt
        if len(self.modified_peak_indices) == 0 and len(self.unmodified_peak_indices) == 0:
            self.calculate_peaks()
        distance = 10000
        closest  = -1
        i = 0
        for unmod in self.unmodified_peak_indices:
            pkwt = self.plots[0].x[unmod]
            if abs(pkwt-wt) < distance:
                distance = abs(pkwt-wt)
                closest = i
            i += 1
        #for panel in self.plot_panels:
        #    panel.draw_peaks([self.modified_peak_indices[closest]], 'modified peaks', 'light green')
            
        
    def block_internal_standard_by_sequence(self, sequence):
        if sequence not in self.blocked_internal_standard_sequences:
            self.blocked_internal_standard_sequences.append(sequence)

    def unblock_internal_standard_by_sequence(self, sequence):
        if sequence in self.blocked_internal_standard_sequences:
            for i in range(len(self.blocked_internal_standard_sequences)):
                if self.blocked_internal_standard_sequences[i] == sequence:
                    del self.blocked_internal_standard_sequences[i]

    def block_internal_standard_by_weight(self, weight):
        frag_results = self.experiment.get_fragment_for_weight_from_all(weight)
        if len(frag_results) > 0:
            self.block_internal_standard_by_sequence(frag_results[0]['sequence'])

    def is_internal_standard_weight_blocked(self, weight):
        int_std_weights = self.find_weights_for_sequences(self.blocked_internal_standard_sequences)
        for wt in int_std_weights:
            if abs(float(weight)-wt) < self.filter_resolution * float(weight):
                return 1
        else:
            return 0

    def is_internal_standard_sequence_blocked(self, sequence):
        if sequence in self.blocked_internal_standard_sequences:
            return 1
        else:
            return 0
        
    def apply_internal_standards(self):
        # offer both internal and external standards
        intstds = self.experiment.get_external_standards()
        print 'got external standards %s'%(intstds)
        keys = intstds.keys()
        external_standards = []
        for key in keys:
            if len(key) > 0:
                if float(intstds[key]) not in external_standards:
                    external_standards.append(float(intstds[key]))

        print 'external_standards %s'%(self.find_indices_for_weights(external_standards))                    
                
        builtin_standards = self.experiment.get_unreactive_fragments()
        print 'unreactive fragments %s'%(self.find_indices_for_weights(builtin_standards))

        internal_standards = self.find_indices_for_weights(external_standards) + self.find_indices_for_weights(builtin_standards)

        int_std_indices,removed_indices = self.remove_ambiguous_indices(internal_standards)

        # make sure none of the indices are in the blocked_internal_standard_sequences list
        #if len(self.blocked_internal_standard_sequences) > 0:
        #    blocked_indices = self.find_indices_for_sequences(self.blocked_internal_standard_sequences)
        #    unblocked_indices = []
        #    for index in int_std_indices:
        #        if index not in blocked_indices:
        #            unblocked_indices.append(index)

        #    int_std_indices = unblocked_indices

        weights = []
        for int_std_index in int_std_indices:
            weights.append(self.max_plot.x[int_std_index])
        sequences = self.find_sequences_for_weights(weights)


        #int_std_indices = store_iss

        max_trough_factor = 0.5
        # first remove any that overlap signficantly with other crap
        # or that are insignificantly high
        good_int_std_indices = []

        print 'int std indices gg %s'%(int_std_indices)
        qweights = self.find_weights_for_indices(int_std_indices)
        print 'int std weights gg %s'%(qweights)

        if len(int_std_indices) > 1:  # be picky only if there are multiple internal standards offered
            for ind in int_std_indices:
                a_avg = 0.0
                b_avg = 0.0
                c_avg = 0.0
                for panel in self.plot_panels:
                    a,b,c = panel.get_trough_indices_for_peak_index(ind)
                    a_avg += panel.y[a]
                    b_avg += panel.y[b]
                    c_avg += panel.y[c]
                a_avg /= len(self.plot_panels)
                b_avg /= len(self.plot_panels)
                c_avg /= len(self.plot_panels)
                if a_avg > max_trough_factor*b_avg or c_avg > max_trough_factor*b_avg:
                    print 'int standard %s < %s over the higher trough'%(panel.x[ind], max_trough_factor*b_avg)
                else:
                    good_int_std_indices.append(ind)
        elif len(int_std_indices) == 1:
            good_int_std_indices = int_std_indices

        if len(good_int_std_indices) == 0:
            print 'no internal standards to use. aborting normalization'
            return []

        good_int_std_indices.sort()
        print 'calling intstd indices good %s'%(good_int_std_indices)
        molecular_weights = []
        for int_std in good_int_std_indices:
            a,b,c = self.max_plot.get_trough_indices_for_peak_index(int_std)
            molecular_weights.append(self.max_plot.x[b])

        ht = []     # heights table
        for index in good_int_std_indices:
            ht.append([])
            height = 0.0
            i = 0
            for panel in self.plot_panels:
                # the commented out lines perform normalization by heights instead of area
                a = panel.plot_object.get_area_for_peak_index(index)
                ht[-1].append(a)
                t1,p,t2 = panel.plot_object.get_trough_indices_for_peak_index(index)
                #ht[-1].append(panel.y[p])
                i += 1

        # calculate noise levels yada yada
        if len(good_int_std_indices) >= 1:
            nht = []
            invert = 1
            normalize_peaks = 1
            if invert:
                width = len(ht[0])
                height = len(ht)
                # create the matrix
                new_ht = []
                for i in range(width):
                    new_ht.append([])
                    for j in range(height):
                        new_ht[-1].append(0.0)
                # invert
                for i in range(height):
                    for j in range(width):
                        new_ht[j][i] = ht[i][j]
                old_ht = copy.deepcopy(ht)
                ht = copy.deepcopy(new_ht)

            # normalize by row
            for peak_obs in ht:
                nht.append([])
                max = 0.0
                for peak in peak_obs:
                    #if peak > max:
                    max += peak
                max = max / len(peak_obs)
                for peak in peak_obs:
                    nht[-1].append(peak/max)

            if invert:
                ht = copy.deepcopy(old_ht)
                # now revert by putting back into old_ht
                for i in range(height):
                    for j in range(width):
                        old_ht[i][j] = nht[j][i]
                nht = copy.deepcopy(old_ht)

        if len(good_int_std_indices) > 1:
            # now calculate averages
            peak_averages = []
            for peak_obs in nht:
                sum = 0.0
                for peak in peak_obs:
                    sum += peak
                peak_averages.append(sum/float(len(peak_obs)))

            if normalize_peaks:
                tmp_nht = []
                new_averages = []
                i = 0
                for peak_obs in nht:
                    tmp_nht.append([])
                    new_averages.append(1.0)
                    for peak in peak_obs:
                        tmp_nht[-1].append(peak/peak_averages[i])
                    i += 1
                nht = copy.deepcopy(tmp_nht)
                peak_averages = copy.deepcopy(new_averages)

            if len(nht[0]) > 1:
                # and SDs
                peak_SDs = []
                i = 0
                for peak_obs in nht:
                    sum = 0.0
                    for peak in peak_obs:
                        sum += (peak-peak_averages[i])**2
                    peak_SDs.append(math.sqrt(sum/(len(peak_obs)-1)))
                    i += 1

                # reprint the table with normalized SD scores
                i = 0
                for peak_obs in nht:
                    print 'wt %5.3f'%(molecular_weights[i]),
                    for peak in peak_obs:
                        print '%5.3f'%(peak),
                    print ': %5.3f, %5.3f, %5.3f'%(peak_averages[i], peak_SDs[i], peak_SDs[i]/peak_averages[i])
                    i += 1

                # now print out information on how many outliers are in each timepoint
                outlier_counts = []
                for tp in range(len(nht[0])):
                    outlier_counts.append(0)
                    for j in range(len(nht)):
                        outlier_counts[tp] += abs(nht[j][tp]-peak_averages[j])/peak_SDs[j]
                    outlier_counts[tp] /= len(nht)

                for tp in range(len(ht[0])):
                    print 'timepoint %s averages %5.3f SDs from the mean'%(tp, outlier_counts[tp])

                # now that the analysis is done, scale all by the lowest RMSD peak
                rmsd_thresh = 0.10
                lowest_indices = []
                while len(lowest_indices) < 3 and len(lowest_indices) < len(nht):
                    rmsd_thresh += 0.05
                    for i in range(len(nht)):
                        if i not in lowest_indices:
                            if peak_SDs[i]/peak_averages[i] < rmsd_thresh:
                                #lowest_rmsd = peak_SDs[i]/peak_averages[i]
                                lowest_indices.append(i)

                avg_peak_heights = []
                for i in range(len(ht)):
                    avg_peak_height = 0.0
                    for h in ht[i]:
                        avg_peak_height += h

                    avg_peak_heights.append(avg_peak_height / len(ht[i]))
            else:
                avg_peak_heights = []
                lowest_indices   = []
                for low_ind in range(len(nht)):
                    avg_peak_height = 0.0
                    lowest_indices.append(low_ind)
                    for h in ht[low_ind]:
                        avg_peak_height += h
                    avg_peak_heights.append(avg_peak_height / len(ht[low_ind]))
            
        elif len(good_int_std_indices) == 1:
            lowest_rmsd = 'NA'
            avg_peak_height = 0.0
            for val in ht[0]:
                avg_peak_height += val
            avg_peak_heights = [avg_peak_height/len(ht[0])]
            lowest_indices = [0]

        """ at this point we have decided how many internal standards might be available
            if none are available, make the user aware of this
            for every external or built-in internal standard located, offer as usable.
            include the following information:
            avg. height, variance after normalization, avg. distance from ideal location
            suggest the set included in lowest_indices.
        """

        # lowest_indices are the recommended subset. this function can be replaced with
        # a simple assignment of lowest_indices to selected_indices for rapid automation.
        automated = 1

        #if len(nht) == 1:
        nht = ht

        if automated:
            print 'passing int_std_indices %s'%(good_int_std_indices)
            indices = self.user_select_internal_standards(good_int_std_indices, lowest_indices, nht, avg_peak_heights)
        else:
            indices = self.finish_apply_internal_standards(good_int_std_indices, lowest_indices, nht, avg_peak_heights)

        return indices
       
    def finish_apply_internal_standards(self, good_int_std_indices, selected_indices, nht, avg_peak_heights):
        print 'in finish_apply'
        print 'heights table'
        for i in nht:
            for j in i:
                print '%s, '%(j),
            print
            
        print 'good indices %s'%(good_int_std_indices)
        print 'selected     %s'%(selected_indices)
        print 'avg hts      %s'%(avg_peak_heights)
        for i in range(len(nht[0])):
            rescales = []
            for ind in range(len(selected_indices)):
                height = nht[selected_indices[ind]][i]
                rescales.append(100.0*(100.0/height))
            print 'rescales %s'%(rescales)

            # now average the rescales
            sum = 0.0
            for rescale in rescales:
                sum += rescale

            real_rescale = sum/float(len(rescales))
            print 'real_rescale %s'%(real_rescale)
            
            self.plot_panels[i].plot_object.hard_rescale_y(real_rescale)
            for bplot in self.plot_panels[i].background_plots:
                bplot.hard_rescale_y(real_rescale)

        # now recalculate the averages and assign to .y
        for panel in self.plot_panels:
            new_ys = []
            for i in range(len(panel.y)):
                new_ys.append(0.0)
            for bkplot in panel.background_plots:
                for i in range(len(panel.y)):
                    new_ys[i] += bkplot.y[i]
            z = len(panel.background_plots)
            for i in range(len(panel.y)):
                new_ys[i] /= z
            panel.y = new_ys
            panel.plot_object.y = new_ys

            panel.draw_coordinates()

        used_standards_indices = []
        for selected in selected_indices:
            used_standards_indices.append(good_int_std_indices[selected])
            
        standards_weights = []
        for ind in selected_indices:
            standards_weights.append(self.max_plot.x[good_int_std_indices[ind]])

        sequences = self.find_sequences_for_weights(standards_weights)
        
        standards = []
        for i in range(len(standards_weights)):
            standards.append([standards_weights[i], sequences[i]])

        print 'setting internal standards'
        print standards
        
        self.experiment.set_normalization_standards(standards)
        
        return used_standards_indices

    def user_select_internal_standards(self, good_indices, suggested_indices, nht, avg_peak_heights):
        good_weights = self.find_weights_for_indices(good_indices)

        peak_height_SDs = []
        for i in range(len(nht)):
            peak_sd = 0.0
            for h in nht[i]:
                peak_sd += (avg_peak_heights[i]-h)**2
            peak_height_SDs.append(math.sqrt(peak_sd/float(len(nht[i])-1)))

        self.select_standards_top = Toplevel(self)
        self.select_standards_top.title('Select Internal Standards')
        self.select_standards_top.wm_transient(self.parent.winfo_toplevel())
        self.info_text  = Pmw.ScrolledText(self.select_standards_top,
                                          labelpos = 'nw',
                                          label_text='                    |   weight   |   avg ht.   |     SD  ')
        values = []
        self.standards_vars = []
        for i in range(len(good_indices)):
            self.standards_vars.append(IntVar())
            c = Checkbutton(self.info_text.component('text'), variable=self.standards_vars[-1])
            c.pack(expand=NO, fill=NONE, side='left')
            self.info_text.window_create('end', window=c)

            c_lambda = lambda x=good_weights[i]: self.xsbar.scroll_to(x)
            b = Button(self.info_text.component('text'), text='View', command=c_lambda)
            self.info_text.window_create('end', window=b)
            self.info_text.insert('end', '| %10.2f | %10.2f | %10.2f\n'%(self.max_plot.x[good_indices[i]], avg_peak_heights[i], peak_height_SDs[i]))

        for i in range(len(good_indices)):
            if i in suggested_indices:
                self.standards_vars[i].set(1)
            else:
                self.standards_vars[i].set(0)

        selected_indices = []
        for i in range(len(values)):
            if values[i]:
                selected_indices.append(i)

        self.buttonBox = Pmw.ButtonBox(self.select_standards_top, orient='horizontal')

        c_lambda = lambda a=good_indices, c=selected_indices, d=nht, e=avg_peak_heights: self.apply_user_select_internal_standards(a,c,d,e)
        
        self.buttonBox.add('Ok', command=c_lambda)
        self.buttonBox.add('Cancel', command=self.select_standards_top.destroy)
        self.info_text.pack(side='left', expand=0, fill='x', anchor='w')
        self.buttonBox.pack(side='top', expand=1, fill='both')

    def apply_user_select_internal_standards(self,good_int_std_indices, selected_indices, nht, avg_peak_heights):
        self.select_standards_top.destroy()

        selected_indices = []
        i = 0
        for var in self.standards_vars:
            if var.get():
                selected_indices.append(i)
            i += 1

        self.finish_apply_internal_standards(good_int_std_indices, selected_indices, nht, avg_peak_heights)

        standards_weights = []
        for selected_index in selected_indices:
            standards_weights.append(self.plot_panels[0].x[good_int_std_indices[selected_index]])

        sequences = self.find_sequences_for_weights(standards_weights)
        
        used_standards_weights = []
        for i in range(len(standards_weights)):
            used_standards_weights.append([standards_weights[i], sequences[i]])

        print 'setting intsts to %s'%(used_standards_weights)
        self.experiment.set_normalization_standards(used_standards_weights)
        print 'normalization standards %s'%(self.experiment.get_normalization_standards())
        
        self.draw_peaks()
        
    def close_user_select_internal_standards(self):
        self.select_standards_top.withdraw()

    def find_indices_for_weights(self, weights):
        # finds the indices of peaks associated with theoretical fragment weights
        indices = {}
        for frag_weight in weights:
            this_frags_peaks = []
            for peak in self.max_plot.peaks:
                if abs(frag_weight-self.max_plot.x[peak]) < frag_weight * self.filter_resolution:
                    this_frags_peaks.append(peak)
            # if there are multiple possibilities, retrieve the closest
            if len(this_frags_peaks) > 1:
                closest_peak_index = 0
                closest_distance   = 10000
                i = 0
                for peak in this_frags_peaks:
                    dist = abs(frag_weight-self.max_plot.x[peak])
                    if dist < closest_distance:
                        closest_distance = dist
                        closest_peak_index = i
                    i += 1
                this_frags_peak = this_frags_peaks[closest_peak_index]
                indices[this_frags_peak] = frag_weight
            elif len(this_frags_peaks) == 1:
                this_frags_peak = this_frags_peaks[0]
                indices[this_frags_peak] = frag_weight

            """
            # this is the one. if the peak is already present,
            if this_frags_peak in indices.keys():
                # figure out which fragment weight is closer and assign
                if abs(frag_weight-self.max_plot.x[this_frags_peak]) < abs(indices[this_frags_peak]-self.max_plot.x[this_frags_peak]):
                    indices[this_frags_peak] = frag_weight
            else:
                indices[this_frags_peak] = frag_weight
            """
        return indices.keys()

    def find_weights_for_sequences(self, sequences):
        # translate to weights and use find_indices_for_weights
        aa_weights = {'A':71.09,  'C':103.15, 'D':115.09, 'E':129.12, 'F':147.18,
                      'G':57.05,  'H':137.14, 'I':113.16, 'K':128.17, 'L':113.16,
                      'M':131.19, 'N':114.11, 'P':97.12,  'Q':128.14, 'R':156.19,
                      'S':87.08,  'T':101.11, 'V':99.14,  'W':186.21, 'Y':163.18}
        frag_weights = []
        for seq in sequences:
            frag_weights.append(18.0)
        i = 0
        for seq in sequences:
            for aa in seq:
                frag_weights[i] += aa_weights[aa]
            i += 1
        return frag_weights

    def find_sequences_for_weights(self, weights):
        sequences = []
        for wt in weights:
            frag_results = self.experiment.get_fragment_for_weight_from_all(wt)
            if len(frag_results) > 0:
                sequences.append(frag_results[0]['sequence'])
            else:
                sequences.append("UNK%s"%(random.randint(0,10000)))
        return sequences

    def find_weights_for_indices(self, indices):
        panel = self.plot_panels[0]
        weights = []
        for ind in indices:
            weights.append(panel.x[ind])
        return weights
        

    def find_indices_for_sequences(self, sequences):
        frag_weights = self.find_weights_for_sequences(sequences)
        return self.find_indices_for_weights(frag_weights)


    def draw_PlotPanels(self, type=0):
        """ draw the coordinates, orange bars for uninteresting frags, red for interesting, and x-axes
            type is 0 to draw the plots themselves
                 is 1 to just draw the axis, fragments, and possible peaks without the plots
                 is 2 to just (re)draw the plots
        """
        for plot_panel in self.plot_panels:
            plot_panel.plot.config(scrollregion=(plot_panel.x_scale*plot_panel.x[0], 0, plot_panel.x_scale*plot_panel.x[-1],self.height))
            if type != 1:
                plot_panel.draw_coordinates()
            if type != 2:
                plot_panel.draw_all_fragment_locations(self.all_fragment_weights)
                plot_panel.draw_peak_rec_possible_peaks()
                plot_panel.draw_x_axes()
                
    def set_all_fragment_weights(self, all_fragment_weights):
        self.all_fragment_weights = all_fragment_weights

    def scaleY(self, value):
        val = self.heightScaler.get()
        for panel in self.plot_panels:
            panel.rescale(panel, val)
            panel.rescale_y(panel, val)
        if self.experiment.has_controls():
            self.controls_panel.rescale(self.controls_panel, val)
            self.controls_panel.rescale_y(self.controls_panel, val)

    def scrollX(self, *args):
        for x in self.plot_panels:
            apply(x.plot.xview, args)
        if self.experiment.has_controls():
            apply(self.controls_panel.plot.xview, args)

    def scrollXmoveto(self, fraction):
        for panel in self.plot_panels:
            panel.plot.xview_moveto(fraction)
        if self.experiment.has_controls():
            self.controls_panel.plot.xview_moveto(fraction)

    def scrollXby(self, number):
        for panel in self.plot_panels:
            panel.plot.xview_scroll(number, "units")
        if self.experiment.has_controls():
            self.controls_panel.plot.xview_scroll(number, "units")

    def scrollY(self, *args):
        apply(self.observations_info_text.yview, args)
        apply(self.plots_text.yview, args)

    def autoscale_intensities(self):
        # not sure why the following two lines were here.
        #if len(self.unmodified_peak_indices) == 0:
        #    self.draw_unmodified_peaks()

        # ranges are collected over the maximum of the interpolated plots
        max_plot = PlottingTools.max_of_plots(self.plots)
        loq = self.max_plot._recalculate_baseline(self.plot_panels[0].get_last_selection_indices())
        max_plot.recognize_peaks(self.experiment.get_mass_range()[0], loq)

        start_at = 0
        for i in range(len(max_plot.x)):
            if max_plot.x[i] > self.experiment.get_mass_range()[0]:
                start_at = i
                break

        ranges = []
        thresh = 0.5
        direction = 0
        for i in range(start_at, len(max_plot.x)):
            if direction == 0:
                if max_plot.y[i] > thresh:
                    new_range = [i]
            if direction == 1:
                if max_plot.y[i] < thresh:
                    new_range.append(i)
                    ranges.append(new_range)
            
        sum_intensities = []
        for plot in self.plot_panels:
            #sum_intensities.append(plot._calculate_peak_height_sum(self.unmodified_peak_indices, self.modified_peak_indices, self.filter_resolution))
            sum_intensities.append(plot._calculate_height_sum(ranges))
        greatest = max(sum_intensities)
        ratios = []
        for intensity in sum_intensities:
            if intensity:
                ratios.append(greatest/intensity)
            else:
                ratios.append(1.0)
                
        for i in range(len(sum_intensities)):
            self.plot_panels[i].hard_rescale_y(ratios[i])
            
        for plot in self.plot_panels:
            plot.draw_coordinates()
            
    def _draw_modified_possible_fragments(self, i, id, wt):
        # i  = fragment_line_id
        # id = the line's actual id
        # wt = the molecular weight
        for plot in self.plot_panels:
            plot.visible_modification_ids[id] = plot._draw_modified_possible_fragments(i)
        #self.draw_modified_peaks_for_peak_weight(id)
            
    def _delete_modified_possible_fragments(self, id):
        for plot in self.plot_panels:
            plot._delete_modified_possible_fragments(id)


    def draw_peaks(self):
        """ draw_peaks generates the peak set including uninteresting_indices, as well as selected
            and unselected internal standards and reactions
        """
        avg_loq = 0.0
        for panel in self.plot_panels:
            panel.plot_object.recognize_peaks(self.experiment.get_mass_range()[0], panel.peak_LOQ)
        #self.max_plot.recognize_peaks(self.experiment.get_mass_range()[0], self.max_plot.peak_LOQ)

        istds = self.experiment.get_normalization_standards()
        keys  = istds.keys()
        int_std_weights = []
        for key in keys:
            int_std_weights.append(istds[key])

        int_std_indices = self.find_indices_for_weights(int_std_weights)
        self.calculate_peaks()

        # first get the indices
        all_indices = self.find_indices_for_weights(self.experiment.get_all_possible_fragments(1))

        real_unmod_indices = self.unmodified_peak_indices
        real_mod_indices =   self.modified_peak_indices
        
        xs = self.max_plot.x

        all_interesting_indices = int_std_indices + real_unmod_indices + real_mod_indices
        
        uninteresting_indices = []

        # collect peak indices that are not near interesting indices
        # and that are not near the fragments associated with interesting indices
        for index in all_indices:
            for comparison_index in all_interesting_indices:
                if abs(xs[index]-xs[comparison_index]) < self.filter_resolution * xs[index]:
                    break
                indices1 = self.find_indices_for_weights([xs[comparison_index]])
                if len(indices1) > 0:
                    uninteresting_peak_index = indices1[0]
                else:
                    continue
                    
                indices2 = self.find_indices_for_weights([xs[index]])
                if len(indices2) > 0:
                    interesting_peak_index   = indices2[0]
                else:
                    continue
                    
                if abs(xs[uninteresting_peak_index]-xs[interesting_peak_index]) < self.filter_resolution * xs[uninteresting_peak_index]:
                    break
            else:
                uninteresting_indices.append(index)
        

        print 'uninteresting indices'
        for ii in uninteresting_indices:
            print '%s, '%(xs[ii]),
        print
        
        self.uninteresting_indices = uninteresting_indices
        self.internal_standard_indices = int_std_indices

        peaklists = {'uninteresting':self.uninteresting_indices,
                     'unmodified':self.unmodified_peak_indices,
                     'modified':self.modified_peak_indices,
                     'int_standards':self.internal_standard_indices}
                     
        # create the peak profiles first so that draw_peaks can associate them with tag ids
        self._create_peak_profiles()
        filename = 'reaction_profiles.csv'
        f = open(filename, 'w')
        line = '%s,'%(len(self.experiment.get_protein_sequences()[0])-1)
        for panel in self.plot_panels:
            line += '%s, '%(panel.plot_object.label)
        
        f.writelines(line + '\n')
        f.close()
        self.print_reaction_profiles(filename)
        
        for key in peaklists.keys():
            tags = {}
            for panel in self.plot_panels:
                # draw_peaks returns the set of tag ids used to reference the peak drawings
                tags[panel.plot_object.label] = panel.draw_peaks(peaklists[key],key)
            invert_tags = {}

            if len(tags) > 0:
                i = 0
                keys = tags.keys()
                keys.sort()
                for i in range(len(tags[keys[0]])):                   # for each peak in the spectrum
                    invert_tags = {}
                    keys2 = tags.keys()
                    keys2.sort()
                    for key2 in keys2:                  # for each set of tag ids (each spectrum)
                        invert_tags[key2] = tags[key2][i]          # collecting the set of tags for the i'th peak
                    self.peak_profiles[key][i].set_tag_ids(invert_tags)   # and submit

        self.xsbar.update()
        
        print 'done drawing peaks\n'

    def popup_data_overview(self):
        self.data_overview_top = Toplevel(self)
        self.data_overview_top.title('Data Overview')
        self.data_overview_top.wm_transient(self.parent.winfo_toplevel())
        n_controls = len(self.experiment.get_control_plots())
        tps = self.experiment.get_timepoints()
        n_timepoints = len(tps)
        keys = tps.keys()
        keys.sort()
        n_files = n_controls
        for key in keys:
            n_files += len(tps[key])

        self.seq_text  = Pmw.ScrolledText(self.data_overview_top)
        #filesd
        ovt = "%s files in %s timepoints and %s controls\n"%(n_files, n_timepoints, n_controls)   # overview text
        # protease, reagent
        ovt += "protease %s cuts %s, reagent %s adds %s Da. to %s\n"%(self.experiment.get_protease(),
                                                                      self.experiment.get_protease_specificity(),
                                                                      self.experiment.get_modifying_reagent(),
                                                                      self.experiment.get_modifying_reagent_weight(),
                                                                      self.experiment.get_modifying_reagent_specificity())
        # internal standards
        ovt += "%s external standards available for calibration\n"%(len(self.experiment.get_external_standards().keys()))
        #itss_weights = self.experiment.get_internal_standards().keys()
        #blockedcnt = 0
        #for itss_weight in itss_weights:
        #    if self.is_internal_standard_weight_blocked(itss_weight):
        #        blockedcnt += 1

        #ovt += "%s of those standards are blocked\n"%(blockedcnt)
        
        # number of fragments in range
        mass_range = self.experiment.get_mass_range()
        total_count = 0
        for f in self.experiment.get_all_possible_fragments():
            if f > mass_range[0] and f < mass_range[1]:
                total_count += 1
        interesting_count = 0
        for f in self.experiment.get_single_reaction_fragment_weights():
            if f > mass_range[0] and f < mass_range[1]:
                interesting_count += 1
        ovt += "%s fragments in range, %s zero-missed cutsites, single reaction\n"%(total_count, interesting_count)
        # how many peaks above LOQ
        labels = {}
        for p in self.plot_panels:
            labels[p.plot_object.label] = [p.plot_object.peak_count, p.peak_LOQ]
        keys = labels.keys()
        keys.sort()
        for key in keys:
            ovt += "\ntp %s: %s peaks, LOQ %s\n"%(key, labels[key][0], labels[key][1])
        ovt += "\n"
        # # of peak/fragment associations per class
        keys = self.peak_profiles.keys()
        keys.sort()
        for key in keys:
            if len(self.peak_profiles[key]) > 0:
                ovt += "%s peak-fragment associations in class %s\n"%(len(self.peak_profiles[key]), key)
        
        self.seq_text.settext(ovt)
        self.seq_text.pack(side='right', fill='both', expand=1, padx=2, pady=2)

    def resize_plots(self, event):
        if event.widget.__class__.__name__ == "Toplevel":
            geom = event.widget.wm_geometry()
            geom = geom.split('+')[0]
            geom = int(geom.split('x')[0])

            for panel in self.plot_panels:
                panel.plot.config(width=geom-230)

            self.xsbar.update_width(geom-230)                

            if self.experiment.has_controls():
                self.controls_panel.plot.config(width=geom-230)

    def set_peak_rec_possible_peaks(self, possible_peaks):
        self.possible_peaks = possible_peaks
        for panel in self.plot_panels:
            panel.set_peak_rec_possible_peaks(possible_peaks)

    def get_unrationalized_peaks(self):
        return_set = []
        peaklists = {'uninteresting': self.uninteresting_indices,
                     'int_standards': self.internal_standard_indices,
                     'unmodified': self.unmodified_peak_indices,
                     'modified': self.modified_peak_indices}
        for peak_index in self.max_plot.peaks:
            val1 = self.max_plot.x[peak_index]
            lobound = self.max_plot.x[peak_index] - (self.experiment.get_filter_resolution() * self.max_plot.x[peak_index])
            hibound = self.max_plot.x[peak_index] + (self.experiment.get_filter_resolution() * self.max_plot.x[peak_index])
            foundit = 0
            for index in self.uninteresting_indices:
                if self.max_plot.x[index] > lobound and self.max_plot.x[index] < hibound:
                    foundit = 1
            for index in self.internal_standard_indices:
                if self.max_plot.x[index] > lobound and self.max_plot.x[index] < hibound:
                    foundit = 1
            for index in self.unmodified_peak_indices:
                if self.max_plot.x[index] > lobound and self.max_plot.x[index] < hibound:
                    foundit = 1
            for index in self.modified_peak_indices:
                if self.max_plot.x[index] > lobound and self.max_plot.x[index] < hibound:
                    foundit = 1
            if not foundit:
                return_set.append(peak_index)
        return return_set

    def get_rationalized_peaks(self):
        return_set = []
        peaklists = {'uninteresting': self.uninteresting_indices,
                     'int_standards': self.internal_standard_indices,
                     'unmodified': self.unmodified_peak_indices,
                     'modified': self.modified_peak_indices}
        for listkey in peaklists.keys():
            for index in peaklists[listkey]:
                return_set.append(index)
        return return_set

class SimplePlotWindow(PlotWindow):
    def __init__(self, parent, average_plots, background_plots, label="", height=100, width=400):
        """ similar to the above without all the bells and whistles
        """
        """ MSWindow acts as a holder for multiple PlotPanel objects and coordinates their activity.
            It also offers tools to integrate controls
        """
        Frame.__init__(self, parent, borderwidth=0, height=height, width=width)
        self.width=width
        self.height=height
        self.parent=parent
        self.config(height=self.height, width=self.width, bg='white')
        labelfont = ('helvetica', 9, 'bold')
        self.plots = average_plots
        self.background_plots = background_plots
        self.spectra_frame = Frame(self)
        
        # heightScaler scales the heights of all ms plots
        self.hscaleframe = Frame(self.spectra_frame)
        Label(self.hscaleframe, text='scale\nall', width=5, font=labelfont).pack(side=TOP)
        self.heightScaler = Scale(self.hscaleframe, from_=10.0, to=0.2, resolution=0.2, orient='vertical', width=15)
        self.heightScaler.set(1.0)
        self.heightScaler.pack(side=TOP,fill=Y, expand=1)
        # fix the callback handler
        self.heightScaler.config(command=self.scaleY)
        self.hscaleframe.pack(side=LEFT, expand=NO, fill=Y)
        Label(self.hscaleframe, text='             ').pack(side=LEFT)
        self.plotCount = len(self.plots)

        # this active message box is below the control panels and shows info on
        # the fragment tags below the x-axis. Look below for the rest of the
        # active messages.
        self.active_message = ""
        self.active_text      = StringVar()
        self.active_text.set(self.active_message)
        self.active_label     = Label(self.spectra_frame, textvariable=self.active_text)
        self.active_label.pack(side='bottom', expand=0, fill='x')

        self.active_message = label
        self.active_text.set(self.active_message)

        self.plot_panels = []
        self.plots_frame = Frame(self.spectra_frame, height=self.height)
        
        self.plots_frame.pack(side=TOP, expand=YES, fill=BOTH)
        self.information_frame = Frame(self, bg='white')

        # add a y-axis scroll bar
        self.ysbar=Scrollbar(self.plots_frame)
        self.ysbar.pack(side=RIGHT, fill=Y)
        self.ysbar.config(command=self.scrollY, cursor='arrow')

        # create a scrollable canvas, the same height as the combined set of plots
        self.plots_text             = Text(self.plots_frame, height=4, yscrollcommand=self.ysbar.set)
        self.plots_text.parent      = self

        # and add the mass spec panels
        self.plot_panels = []
        labelFont = ('Helvetica', parms.get('label_font_size'))
        # create the plot panels
        for i in range(len(self.plots)):
            avg_plot = self.plots[i]
            bkg_plots = self.background_plots[i]
            panel = PlotPanel(self.plots_text, avg_plot, bkg_plots, [], None, (self.width-40)/2.0, 100, 1)
            self.plot_panels.append(panel)

        # insert them into the canvas
        i = 0
        for panel in self.plot_panels:
            panel.pack(expand=YES, fill=BOTH, side=TOP)
            self.plots_text.window_create('end', window=panel)
            # a radioselect is here used to visualize the average or all spectra
            d_lambda = lambda x, y, z=self.plot_panels[i]: self.adjust_spectra_view(z, x, y)
            cbframe = Frame(self.plots_text)
            panel.checkbuttons = Pmw.RadioSelect(cbframe,
                                                  buttontype = 'checkbutton',
                                                  orient = 'vertical',
                                                  hull_bg='white',
                                                  command = d_lambda)
            panel.checkbuttons.pack(side = 'top', expand=0, fill=NONE)
            for text in ('avg', 'all'):
                panel.checkbuttons.add(text)
            panel.checkbuttons.invoke('avg')
            panel.checkbuttons.invoke('all')
            panel.checkbuttons.button('avg').config(font=labelFont, bg='white')
            panel.checkbuttons.button('all').config(font=labelFont, bg='white')
            cbframe.pack(expand=YES, fill=BOTH, side=TOP)
            self.plots_text.window_create('end', window=cbframe)

            i += 1
            if (i < self.plotCount-1):
                self.plots_text.insert('end', '\n')

        # create the initial view of the panels
        self.draw_PlotPanels(1)

        self.plots_text.config(state=DISABLED)
        self.plots_text.pack(expand=YES, fill=BOTH)
        self.spectra_frame.pack(side=LEFT, expand=YES, fill=BOTH, anchor=N)
        self.information_frame.pack(side=LEFT, expand=NO, fill=Y, anchor=N)

        self.internal_standards = []
        self.reaction_profiles = []
        self.peak_profiles = {}

    def scrollX(self, *args):
        for x in self.plot_panels:
            apply(x.plot.xview, args)

    def scrollY(self, *args):
        apply(self.plots_text.yview, args)

    def scaleY(self, value):
        val = self.heightScaler.get()
        for x in self.plot_panels:
            x.rescale(x, val)
            x.rescale_y(x, val)

    def draw_PlotPanels(self, type=0):
        """ draw the coordinates, orange bars for uninteresting frags, red for interesting, and x-axes
            type is 0 to draw the plots themselves
                 is 1 to just draw the axis, fragments, and possible peaks without the plots
                 is 2 to just (re)draw the plots
        """
        for plot_panel in self.plot_panels:
            plot_panel.plot.config(scrollregion=(plot_panel.x_scale*plot_panel.x[0], 0, plot_panel.x_scale*plot_panel.x[-1],self.height))
            if type != 1:
                plot_panel.draw_coordinates()
                plot_panel.draw_background_coordinates()
            if type != 2:
                plot_panel.draw_x_axes()
        

class PlotPanel(Frame):
    def __init__(self, parent, average_plot, bkg_plots, possible_peaks=[], filter_resolution=None, width=600, height=100, autoscale=0):
        # init
        Frame.__init__(self, parent, bg='white')
        self.possible_peaks = possible_peaks
        if filter_resolution == None:
            self.filter_resolution = 0.002
        else:
            self.filter_resolution = filter_resolution
        self.all_my_labels = []
        self.parent=parent

        self.background_plots_visible = 0
        self.average_plot_visible = 0

        self.background_plots = bkg_plots
        self.background_plots_selected = []
        for plot in bkg_plots:
            self.background_plots_selected.append(1)
        has_average = 1
        if average_plot == None:
            has_average = 0
            average_plot = bkg_plots[0]
        self.plot_object = average_plot
        self.y = self.plot_object.y
        self.x = self.plot_object.x
        self.x_scale = 1.0
        self.height = height
        self.width  = width
        # get some global parameters
        self.smoothing_gap     = parms.get('smoothing_gap')
        self.min_weight        = parms.get('min_weight')
        self.height_thresh     = parms.get('height_thresh')
        self.min_peak_height   = parms.get('min_peak_height')
        self.x_axis_mark_every = parms.get('x_axis_mark_every')
        self.y_axis_mark_every = parms.get('y_axis_mark_every')

        if autoscale:
            # get the major heights and widths
            highest = 0.0
            for b in self.background_plots:
                for i in range(len(b.x)):
                    if b.y[i] > highest:
                        highest = b.y[i]
            last = float(b.x[-1])

            self.x_scale = width/last
            
            self.x_axis_mark_every = 50/self.x_scale
            self.y_axis_mark_every = 20
           
        # create a canvas to draw on
        self.plot = Canvas(self, bg='white', relief=SUNKEN, bd=0, cursor='crosshair', height=self.height, width=self.width)
        if has_average:
            self.plot.bind('<Motion>', self.canvas_motion)
        self.plot.bind('<Control-p>', self.plot_to_postscript)
        self.plot.bind('<ButtonPress-1>', self._select_range)
        #self.plot.bind('<ButtonRelease-1>', self._end_select_range)
        self.plot.bind('<Control-ButtonPress-1>', self._extend_range)
        self.plot.bind('<Control-ButtonRelease-1>', self._end_extend_range)
        #self.plot.bind('<ButtonPress-2>', self._selection_middle_callback)
        self.plot.bind('<ButtonPress-3>', self._selection_right_callback)
        self.plot.bind('<Leave>', self.clear_position_text)

        self.plot.focus_set()
        self.plot.config(highlightthickness=1)
        self.yaxis_width = 25
        # create a canvas for the first y-axis
        self.preplot = Canvas(self, relief=SUNKEN, height=self.height, width=self.yaxis_width, bd=0)
        # and one for the second y-axis
        self.postplot = Canvas(self, relief=SUNKEN, height=self.height, width=self.yaxis_width, bd=0)
        self.preplot.bind('<Motion>', self.canvas_motion)
        self.postplot.bind('<Motion>', self.canvas_motion)

        # initialize the plots
        self.plot.delete('plot')
        self.plot.last_mag = 1.0
        # draw the axes
        self.draw_x_axes()
        self.draw_y_axes()
        self.plot.config(scrollregion=(self.x_scale*self.x[0], 0, self.x_scale*self.x[len(self.x)-1], self.height))
        # finally, pack everything up
        self.preplot.pack(expand=NO, fill=NONE, side=LEFT)
        self.postplot.pack(expand=NO, fill=NONE, side=RIGHT)
        self.plot.pack(expand=NO, fill=NONE, side=TOP)
        # used to hold heights for the sum of quantification distributions
        self.dist_y = []
        self.fragment_line_ids = []
        self.fragment_line_weights = []
        self.mod_fragment_line_ids = []
        self.mod_fragment_line_weights = []
        self.peak_line_ids = {}
        self.peak_line_weights = {}
        self.peak_intensities = {}
        self.peak_indices = {}

        self.drawn_peaks = []
        self.visible_modification_ids = {}
        self.all_fragment_line_ids = []
        self.all_fragment_line_weights = []

        self.color_peak_associations = {}
        self.color_peak_associations['uninteresting']  = 'light grey'
        self.color_peak_associations['int_standards']  = 'papaya whip'#'coral'
        self.color_peak_associations['unmodified']     = 'misty rose'#'tomato'
        self.color_peak_associations['modified']       = 'pale green'#'lime green'

        self.active_select_start_indices = []
        self.hilite_tags = []

        self.right_selection_popup = FixedMenu(self.parent.winfo_toplevel(), tearoff=0)
        self.peak_SD  = -1
        self.peak_LOD = -1
        self.peak_LOQ = self.min_peak_height
        self.last_selection_indices = []

    def set_mass_range(self, parms):
        for bkg in self.background_plots:
            bkg.set_mass_range(parms)
        self.plot_object.set_mass_range(parms)

    def get_mass_range(self, parms):
        return self.plot_object.get_mass_range()

    def get_filter_resolution(self):
        return self.filter_resolution

    def set_filter_resolution(self, new_resolution):
        self.filter_resolution = new_resolution
        for plot in self.background_plots:
            plot.set_filter_resolution(new_resolution)
        self.plot_object.set_filter_resolution(new_resolution)

    def draw_centroids(self):
        for i in self.plot_object.peaks:
            self.plot.create_line(self.x_scale*self.x[i], 0, self.x_scale*self.x[i], self.height, width=1, fill='black', tags='centroids')

    def undraw_centroids(self):
        self.plot.delete('centroids')
        

    def get_current_selection_indices(self):
        # returns a list of lists, in the embedded lists, 0 is reserved as start, 1 as end.
        return self.active_select_start_indices
    
    def get_last_selection_indices(self):
        return self.last_selection_indices

    def _select_range(self, event):
        # this is designed to reject selections on control plots
        # its just one way of doing it
        try:
            self.parent.parent._broadcast_select_range_start(event)
        except AttributeError:
            return
        else:
            self.plot.bind('<ButtonRelease-1>', self._end_select_range)

        
    def _end_select_range(self, event):
        # this is designed to reject selections on control plots
        # its just one way of doing it
        try:
            self.parent.parent._broadcast_select_range_end(event)
        except AttributeError:
            return
        else:
            self.plot.bind('<ButtonRelease-1>', None)


    def _extend_range(self, event):
        # this is designed to reject selections on control plots
        # its just one way of doing it
        try:
            self.parent.parent._broadcast_extend_range_start(event)
        except AttributeError:
            return
        else:
            self.plot.bind('<Control-ButtonRelease-1>', self._end_extend_range)


    def _end_extend_range(self, event):
        self.parent.parent._broadcast_extend_range_end(event)
        self.plot.bind('<Control-ButtonRelease-1>', None)
        

    def _draw_select_range_start(self, event):
        # first clear the old select range
        self.plot.delete('select-range-hilite')
        self.active_select_start_indices = [[int(self.plot.canvasx(event.x)),-1]]
        self.plot.bind('<Motion>', self.parent.parent._broadcast_select_range_motion)

    def _draw_extend_range_start(self, event):
        self.active_select_start_indices.append([int(self.plot.canvasx(event.x)),-1])
        self.plot.bind('<Motion>', self.parent.parent._broadcast_extend_range_motion)

    def _draw_select_range_end(self, event):
        self.active_select_start_indices[0][1] = int(self.plot.canvasx(event.x))
        self.plot.unbind_all('<Motion>')
        self.plot.bind('<Motion>', self.canvas_motion)
        self.plot.dtag('select-range-hilite', 'select-range-hilite-current')

    def _draw_extend_range_end(self, event):
        self.active_select_start_indices[-1][1] = int(self.plot.canvasx(event.x))
        self.plot.unbind_all('<Motion>')
        self.plot.bind('<Motion>', self.canvas_motion)
        self.plot.dtag('select-range-hilite', 'select-range-hilite-current')

    def _draw_select_range_motion(self, x):
        x1 = self.active_select_start_indices[0][0]
        x2 = int(self.plot.canvasx(x))
        self.plot.delete('select-range-hilite')
        self.plot.delete('select-range-hilite-current')
        id = self.plot.create_rectangle(self.x_scale * x1, self.height, self.x_scale*x2, 0, width=1, fill='light goldenrod yellow', tags=('select-range-hilite', 'select-range-hilite-current'))
        #self.plot.tag_bind(id, '<Button-3>', self._selection_right_callback)
        self.plot.tag_lower(id)

    def _draw_extend_range_motion(self, x):
        x1 = self.active_select_start_indices[-1][0]
        x2 = self.plot.canvasx(x)
        self.plot.delete('select-range-hilite-current')
        id = self.plot.create_rectangle(self.x_scale * x1, self.height, self.x_scale*x2, 0, width=1, fill='light goldenrod yellow', tags=('select-range-hilite', 'select-range-hilite-current'))
        #self.plot.tag_bind(id, '<Button-3>', self._selection_right_callback)
        self.plot.tag_lower(id)

    def get_area_under_selected_regions(self):
        start_index = 0
        end_index = 0
        started = 0
        for start_stop in self.active_select_start_indices:
            if start_stop[0] < start_stop[1]:
                i = 0
                for xval in self.x:
                    if xval > start_stop[0] and xval < start_stop[1]:
                        if not started:
                            start_index = i
                            started = 1
                        else:
                            end_index = i
                    i += 1
            else:
                i = 0
                for xval in self.x:
                    if xval > start_stop[1] and xval < start_stop[0]:
                        if not started:
                            start_index = i
                            started = 1
                        else:
                            end_index = i
                    i += 1

        subtract_this = 0
        total = 0.0
        if self.y[start_index] > 0 and self.y[end_index] > 0:
            if self.y[start_index] < self.y[end_index]:
                subtract_this = self.y[start_index]
            for j in range(start_index, end_index):
                total += self.y[j]-subtract_this
        else:
            for j in range(start_index, end_index):
                total += self.y[j]
        return total

    def _selection_middle_callback(self, event):
        self._select_right_callback(event)

    def _selection_right_callback(self, event):
        # create a popup at event.x, event.y
        while self.right_selection_popup.index(END) != None:
            self.right_selection_popup.delete_special(self.right_selection_popup.index(END))

        id = event.widget.find_overlapping(self.plot.canvasx(event.x)-1,self.plot.canvasy(event.y)-1,self.plot.canvasx(event.x)+1,self.plot.canvasy(event.y)+1)

        for tag in id:
            print 'got tag %s'%(self.plot.gettags(tag)[0])
            if self.plot.gettags(tag)[0] == "select-range-hilite":
                self.right_selection_popup.add_command(label='use as baseline', command=self.use_selection_as_baseline)
                # if there is a single peak in the hilite range
                selection_list = self.get_current_selection_indices()
                if len(selection_list) != 1:
                    continue
                selection_pair = selection_list[0]
                # use maxplot to see if there is one and only one peak in the hilited area
                peaklist = []
                for peak_index in self.parent.parent.max_plot.peaks:
                    if self.parent.parent.max_plot.x[peak_index] >= selection_pair[0] and self.parent.parent.max_plot.x[peak_index] <= selection_pair[1]:
                        peaklist.append(peak_index)
                # there should be just one peak
                if len(peaklist) != 1:
                    continue
                # see if the peak has already been rationalized
                low_weight = self.parent.parent.max_plot.x[selection_pair[0]]
                high_weight = self.parent.parent.max_plot.x[selection_pair[1]]
                
                candidate_frag_list = []
                for fragment in self.parent.parent.experiment.fragment_objects:
                    if abs(fragment.get_weight() - self.parent.parent.max_plot.x[peaklist[0]]) < fragment.get_weight()*self.filter_resolution:
                        candidate_frag_list.append(fragment)
                if len(candidate_frag_list) == 0:   # if no fragments are within resolution range
                    c_lambda = lambda x=peaklist[0]: self.rationalize_peak(x)
                    self.right_selection_popup.add_command(label='rationalize peak', command=c_lambda)
                    c_lambda = lambda x=peaklist[0]: self.add_peak_to_external_standards(x)
                    self.right_selection_popup.add_command(label='add to internal standards', command=c_lambda)

            elif self.plot.gettags(tag)[0] == 'all_fragments':
                for j in range(len(self.all_fragment_line_ids)):
                    if tag == self.all_fragment_line_ids[j]:
                        fragment = j
                frag_results = self.parent.parent.experiment.get_fragment_for_weight_from_all(self.all_fragment_line_weights[fragment])
                c_lambda = lambda x=frag_results: self.view_fragment(x)
                self.right_selection_popup.add_command(label='view fragment', command=c_lambda) # launch the molecular viewer
                
            elif self.plot.gettags(tag)[0] == 'mod_fragments':
                for j in range(len(self.all_fragment_line_ids)):
                    if tag == self.all_fragment_line_ids[j]:
                        fragment = j
                frag_results = self.parent.parent.experiment.get_fragment_for_weight_from_all(self.fragment_line_weights[fragment])
                c_lambda = lambda x=frag_results: self.view_fragment(x)
                self.right_selection_popup.add_command(label='view fragment', command=c_lambda) # launch the molecular viewer
            elif self.plot.gettags(tag)[0] == 'uninteresting':
                pass
            elif self.plot.gettags(tag)[0] == 'int_standards':
                # updates self.blocked_internal_standard_sequences then redraw the peaks
                
                # add a command to deselect (block) this internal standard
                i = 0
                for j in self.peak_line_ids['int_standards']:
                    if tag == j:
                        fragment = i
                    i += 1
                # first see if the weight corresponds to a blocked internal standard
                wt = self.peak_line_weights['int_standards'][fragment]
                blocked_test = self.parent.parent.is_internal_standard_weight_blocked(wt)
                if blocked_test:
                    frag_results = self.parent.parent.experiment.get_fragment_for_weight_from_all(self.peak_line_weights['int_standards'][fragment])
                    c_lambda = lambda x=frag_results: self.unblock_sequence_from_internal_standards(x)
                    self.right_selection_popup.add_command(label='unblock', command=c_lambda) # launch the molecular viewer
                else:
                    frag_results = self.parent.parent.experiment.get_fragment_for_weight_from_all(self.peak_line_weights['int_standards'][fragment])
                    c_lambda = lambda x=frag_results: self.block_sequence_from_internal_standards(x)
                    self.right_selection_popup.add_command(label='block', command=c_lambda) # launch the molecular viewer

            elif self.plot.gettags(tag)[0] == 'unmodified':
                i = 0
                for j in self.peak_line_ids['unmodified']:
                    if tag == j:
                        fragment = i
                    i += 1
                frag_results = self.parent.parent.experiment.get_fragment_for_weight_from_all(self.peak_line_weights['unmodified'][fragment])
                c_lambda = lambda x=frag_results: self.view_fragment(x)
                self.right_selection_popup.add_command(label='view fragment', command=c_lambda) # launch the molecular viewer

            elif self.plot.gettags(tag)[0] == 'modified':
                i = 0
                for j in self.peak_line_ids['modified']:
                    if tag == j:
                        fragment = i
                    i += 1
                frag_results = self.parent.parent.experiment.get_fragment_for_weight_from_all(self.peak_line_weights['modified'][fragment] - float(self.parent.parent.experiment.get_modifying_reagent_weight()))

                c_lambda = lambda x=frag_results: self.view_fragment(x)
                self.right_selection_popup.add_command(label='view fragment', command=c_lambda) # launch the molecular viewer
            break

        self.right_selection_popup.tk_popup(event.x_root, event.y_root, 0)

    def add_peak_to_external_standards(self, peak_index):
        self.parent.parent.experiment._append_to_external_standards([[peak_index, 'user_defined']])

    def rationalize_peak(self, peak_index):
        reload(PeakRationalizer)
        self.peak_rationalizer = PeakRationalizer.Viewer(self, self.parent.parent.experiment, peak_index)

    def block_sequence_from_internal_standards(self, fragment_results):
        self.parent.parent.block_internal_standard_by_sequence(fragment_results[0]['sequence'])
        self.parent.parent.draw_peaks()

    def unblock_sequence_from_internal_standards(self, fragment_results):
        self.parent.parent.unblock_internal_standard_by_sequence(fragment_results[0]['sequence'])
        self.parent.parent.draw_peaks()

    def view_fragment(self, fragment_results):
        reload(FragmentPeakTablet)
        f = FragmentPeakTablet.Viewer(self, self.parent.parent.reaction_profiles, fragment_results)

    def use_selection_as_baseline(self, event=None):
        self.parent.parent._broadcast_recalculate_baseline()
        self.parent.parent._clear_selection()

    def _clear_selection(self):
        self.plot.delete('select-range-hilite')
        self.active_select_start_indices = [0,0]

    def recalculate_baseline(self):
        self.last_selection_indices = self.get_current_selection_indices()
        self.peak_LOQ = self.plot_object._recalculate_baseline(self.last_selection_indices)
        return self.peak_LOQ
    
    def set_x_scale(self, val=1.0):
        xevery = parms.get('x_axis_mark_every')
        self.x_axis_mark_every = int(float(xevery)/val)
        if self.x_axis_mark_every == 0:
            self.x_axis_mark_every = 1
        self.x_scale = val
        
    def takefocus(self, event):
        self.plot.focus_set()
    
    def draw_coordinates(self):
        self.plot.delete("plot")
        self.average_plot_visible = 1
        # and draw the lines
        for i in range(len(self.y)-1):
            self.plot.create_line(self.x_scale*self.x[i], self.y[i], self.x_scale*self.x[i+1], self.y[i+1], width=1, fill='black', tags="plot")
        # scale the image
        self.plot.scale("plot", 0, 0, 1, -0.8)
        self.plot.move("plot", 0, (self.height * 0.8))
        
    def undraw_coordinates(self):
        self.plot.delete("plot")
        self.average_plot_visible = 0

    def undraw_background_coordinates(self):
        self.plot.delete("background")
        self.background_plots_visible = 0
        
    def draw_background_coordinates(self):    
        self.plot.delete("background")
        self.background_plots_visible = 1
        j = 0
        colorset = ['blue', 'red', 'green', 'orange', 'yellow', 'purple', 'beige', 'cyan', 'magenta']
        for i in range(len(self.background_plots)):
            plot = self.background_plots[i]
            if self.background_plots_selected[i]:
                # and draw the lines
                color = colorset[j%len(self.background_plots)]
                for i in range(len(plot.y)-1):
                    self.plot.create_line(self.x_scale*plot.x[i], plot.y[i], self.x_scale*plot.x[i+1], plot.y[i+1], width=1, fill=color, tags="background")
                j += 1
        # scale the image
        self.plot.scale("background", 0, 0, 1, -0.8)
        self.plot.move("background", 0, (self.height * 0.8))
    
    def get_area_for_peak_index(self, peak):
        last, closest, next = self.get_trough_indices_for_peak_index(peak)
        
    def get_trough_indices_for_peak_index(self, peak):
        #last_trough, closest, next_trough = self.plot_object.get_trough_indices_for_peak_index(peak)
        last_trough, closest, next_trough = self.parent.parent.max_plot.get_trough_indices_for_peak_index(peak)
        return last_trough, closest, next_trough
        
    def draw_peaks(self, peaks=None, type='unmodified'):
        self.plot.delete(type)
        if peaks == None:
            peaks = self.plot_object.peaks
        self.drawn_peaks = peaks

        color = self.color_peak_associations[type]

        self.peak_line_ids[type] = []
        self.peak_line_weights[type] = []
        self.peak_intensities[type] = []
        self.peak_indices[type] = []
        for peak in peaks:
            last_trough, closest, next_trough = self.get_trough_indices_for_peak_index(peak)
            # move the peak drawing 1 index in on each side
            x0 = self.x_scale*(self.x[last_trough])
            y0 = self.height
            x1 = self.x_scale*(self.x[next_trough])
            y1 = 0
            lin = self.plot.create_rectangle(x0,y0,x1,y1, fill=color, outline='white', tags=type)
            self.peak_line_ids[type].append(lin)
            self.peak_line_weights[type].append(self.x[closest])
            self.peak_intensities[type].append(self.y[closest])
            self.peak_indices[type].append(closest)

            # check against the parent's list of blocked internal standards.
            if self.parent.parent.is_internal_standard_weight_blocked(self.x[closest]):
                # draw an x
                lin2 = self.plot.create_line(x0-1,y0,x1+1,y1, fill='red', tags=type)
                lin3 = self.plot.create_line(x1+1,y0,x0-1,y1, fill='red', tags=type)
                self.plot.lower(lin2)
                self.plot.lower(lin3)
                self.plot.scale(lin2, 0, 0, 1, -0.8)
                self.plot.move(lin2, 0, (self.height * 0.8))
                self.plot.scale(lin3, 0, 0, 1, -0.8)
                self.plot.move(lin3, 0, (self.height * 0.8))

            self.plot.lower(lin)
            self.plot.scale(lin, 0, 0, 1, -0.8)
            self.plot.move(lin, 0, (self.height * 0.8))
                    
        return self.peak_line_ids[type]
        
    def draw_y_axes(self):
        self.preplot.delete('y_axis')
        self.postplot.delete('y_axis')
        labelFont = ('Helvetica', parms.get('label_font_size'))
        for zz in range(0, int(self.height), int(round(self.height/5.0))):
            # pre plot first
            widget = Label(self.preplot, text=zz, fg='black', font=labelFont)
            widget.pack()
            self.preplot.create_window(10, zz, window=widget, tags='y_axis')
            self.preplot.create_line(self.yaxis_width-5, zz, self.yaxis_width, zz, fill='red', tags='y_axis')
            # post plot
            widget = Label(self.postplot, text=zz, fg='black', font=labelFont)
            widget.pack()
            self.postplot.create_window(self.yaxis_width-10, zz, window=widget, tags='y_axis')
            self.postplot.create_line(0, zz, 5, zz, fill='red', tags='y_axis')

        self.preplot.scale("y_axis", 0, 0, 1, -0.8)
        self.preplot.move("y_axis", 0, (self.height * 0.8))

        self.postplot.scale("y_axis", 0, 0, 1, -0.8)
        self.postplot.move("y_axis", 0, (self.height * 0.8))

        self.preplot.last_mag = 1.0
        
    def rescale(MSWindow, self, val):
        for type in ("plot", "dist_plot", "peak_plot", "background"):
            self.plot.move(type, 0, -(self.height * 0.8))
            self.plot.scale(type, 0, 0, 1, float(1.0/float(self.plot.last_mag)))
            self.plot.scale(type, 0, 0, 1, val)
            self.plot.move(type, 0, (self.height * 0.8))
        self.plot.last_mag = val

    def rescale_y(MSWindow, self, val):
        """ value should be in (0,1] (i think). Just alters the appearance"""
        self.preplot.delete('y_axis')
        self.postplot.delete('y_axis')
        labelFont = ('Helvetica', parms.get('label_font_size'))
        arg1 = 0
        arg2 = int((1.0/val) * self.height)
        arg3 = int(math.floor(((1.0/val) * self.height)/5.0))
        for zz in range(arg1, arg2, arg3):
            # pre plot first
            widget = Label(self.preplot, text=zz, fg='black', font=labelFont)
            widget.pack()
            self.preplot.create_window(10, val*zz, window=widget, tags='y_axis')
            self.preplot.create_line(self.yaxis_width-5, val*zz, self.yaxis_width, val*zz, fill='red', tags='y_axis')
            # post plot
            widget = Label(self.postplot, text=zz, fg='black', font=labelFont)
            widget.pack()
            self.postplot.create_window(self.yaxis_width-10, val*zz, window=widget, tags='y_axis')
            self.postplot.create_line(0, val*zz, 5, val*zz, fill='red', tags='y_axis')

        self.preplot.scale("y_axis", 0, 0, 1, -0.8)
        self.preplot.move("y_axis", 0, (self.height * 0.8))
        self.postplot.scale("y_axis", 0, 0, 1, -0.8)
        self.postplot.move("y_axis", 0, (self.height * 0.8))

        #for plot in [self.preplot, self.postplot]:
        #    plot.move('y_axis', 0, -(self.height * 0.80))
        #    plot.scale('y_axis', 0, 0, 1, (1.0/self.preplot.last_mag))
        #    plot.scale('y_axis', 0, 0, 1, val)
        #    plot.move('y_axis', 0, (self.height * 0.80))

        self.preplot.last_mag = val

    def hard_rescale_y(self, ratio):
        """ ratio should be out of 100. this function alters the data itself"""
        self.plot_object.hard_rescale_y(ratio)
        for plot in self.background_plots:
            plot.hard_rescale_y(ratio)

        if self.average_plot_visible:
            self.draw_coordinates()

        if self.background_plots_visible:
            self.draw_background_coordinates()
        
        
    def draw_x_axes(self):
        self.plot.delete('x_axis')
        # draw the lines now
        self.plot.create_line(0,0,self.x_scale*self.x[-1],0,width=0.5, fill='black', tags='x_axis')
        labelFont = ('Helvetica', parms.get('label_font_size'))     
        for zz in range(0,int(self.x[-1]), int(self.x_axis_mark_every)):
            a = self.plot.create_line(self.x_scale*zz, -2, self.x_scale*zz, 2, width=0.5, fill='black', tags='x_axis')
            my_label = str(round(zz))
            x = self.plot.create_text(self.x_scale*zz, -10, fill='black', text = my_label, tags = 'x_axis', font=labelFont)
            self.plot.lower(a)
            self.plot.lower(x)
        self.plot.scale("x_axis", 0, 0, 1, -0.8)
        self.plot.move("x_axis", 0, (self.height * 0.8))
        
    def set_peak_rec_possible_peaks(self, possible_peaks):
        self.possible_peaks = possible_peaks


    def set_peak_rec_filter_resolution(self, resolution=0.002):
        self.filter_resolution = resolution

    def draw_peak_rec_possible_peaks(self):
        # show all of the proteolysis fragments
        self.plot.delete("fragments")
        self.fragment_line_ids = []
        self.fragment_line_weights = []
        i = 0
        for frag in self.possible_peaks:
            if frag < self.x[-1]:
                lin = self.plot.create_line(self.x_scale*frag, -4, self.x_scale*frag, 2, width=2, fill='red', tags='fragments')
                self.fragment_line_ids.append(lin)
                self.fragment_line_weights.append(frag)
                self.plot.tkraise(lin)
            i += 1
        self.plot.scale("fragments", 0, 0, 1, -0.8)
        self.plot.move("fragments", 0, (self.height * 0.82))

    def remove_peak_rec_possible_peaks(self):
        self.plot.delete("fragments")
        
    def draw_all_fragment_locations(self, locations):
        # show all of the proteolysis fragments
        self.plot.delete("all_fragments")
        self.all_fragment_line_ids = []
        self.all_fragment_line_weights = []
        i = 0
        for frag in locations:
            if frag < self.x[-1]:
                lin = self.plot.create_line(self.x_scale*frag, -4, self.x_scale*frag, 2, width=2, fill='orange', tags="all_fragments")
                self.all_fragment_line_ids.append(lin)
                self.all_fragment_line_weights.append(frag)
                self.plot.tkraise(lin)
            i += 1
        self.plot.scale("all_fragments", 0, 0, 1, -0.8)
        self.plot.move( "all_fragments", 0, (self.height * 0.82))

    def remove_all_fragment_locations(self):
        self.plot.delete("all_fragments")
        
        
    def canvas_motion(self, event):
        # the following rejects any activity over the control plot
        # there are better ways of doing this but I have to give a demo tomorrow ;)
        try:
            message = 'x:%5s y:%5.2f'%(int(self.plot.canvasx(event.x)),((1.0/self.plot.last_mag)*(1.0/0.8)*((self.height-event.y)-(0.2*self.height))))
            self.parent.parent.position_message = message
            self.parent.parent.position_text.set(self.parent.parent.position_message)
        except AttributeError:
            return

        try:
            self.parent.parent.experiment
        except AttributeError:
            return
        
        id = event.widget.find_overlapping(self.plot.canvasx(event.x)-1,self.plot.canvasy(event.y)-1,self.plot.canvasx(event.x)+1,self.plot.canvasy(event.y)+1)
        if (len(id) > 0):
            #fragset = self.parent.parent.experiment.get_all_possible_fragments_dictionary()
            all_fragments = self.fragment_line_ids+self.all_fragment_line_ids
            all_fragment_line_weights = self.fragment_line_weights+self.all_fragment_line_weights
            # first look for it in the possible fragments list
            for i in range(len(all_fragments)):
                if (all_fragments[i] in id):                        # if this is the line
                    frag_results = self.parent.parent.experiment.get_fragment_for_weight_from_all(all_fragment_line_weights[i])
                    #print 'frag restulst %s'%(all_fragment_line_weights[i])
                    #print frag_results
                    if len(frag_results) > 1:
                        message = 'ambiguous'
                        self.parent.parent.active_peak_message = message
                        self.parent.parent.active_peak_text.set(self.parent.parent.active_peak_message)
                    elif len(frag_results) == 0:
                        message = 'no fragment found?'
                        self.parent.parent.active_peak_message = message
                        self.parent.parent.active_peak_text.set(self.parent.parent.active_peak_message)
                    else:
                        if id[0] not in self.visible_modification_ids.keys():
                            message = 'weight %s fragment %s %s %s'%(all_fragment_line_weights[i], frag_results[0]['site1'], frag_results[0]['sequence'], frag_results[0]['site2'])
                            self.parent.parent.active_peak_message = message
                            self.parent.parent.active_peak_text.set(self.parent.parent.active_peak_message)
                            if self.parent.parent.experiment.get_digestion_type() == 'exhaustive':
                                self.parent.parent._draw_modified_possible_fragments(i, id[0], all_fragment_line_weights[i])
                            break
                        else:
                            self.parent.parent._delete_modified_possible_fragments(id[0])
                            self.parent.parent.active_peak_message = ""
                            self.parent.parent.active_peak_text.set(self.parent.parent.active_peak_message)
            else:
                # now look for it in the peaks list
                keys = self.parent.parent.peak_profiles.keys()
                found_it = 0
                for key in keys:
                    profiles = self.parent.parent.peak_profiles[key]
                    for peak_profile in profiles:
                        tag_ids = peak_profile.get_tag_ids()
                        tag_id = -1
                        if self.plot_object.label in tag_ids.keys():
                            tag_id = tag_ids[self.plot_object.label]
                        if tag_id == id[0]:
                            weight = peak_profile.get_weight()
                            height = peak_profile.get_heights()
                            heights_array = []
                            keys = height.keys()
                            keys.sort()
                            avg_height = 0.0
                            for a in height[self.plot_object.label]:
                                avg_height += a
                            avg_height /= len(height[self.plot_object.label])
                            message = ""
                            self.parent.parent.active_peak_message = message
                            self.parent.parent.active_peak_text.set(self.parent.parent.active_peak_message)
                            tokens = string.split(key)
                            if tokens[0] in ['selected', 'unselected']:
                                keytouse = tokens[1]
                            else:
                                keytouse = tokens[0]
                            self.parent.parent.active_class_message = '%-s  %s'%("class",keytouse)
                            self.parent.parent.active_class_text.set(self.parent.parent.active_class_message)
                            self.parent.parent.active_weight_message = '%-s %5.1f'%("weight", weight)
                            self.parent.parent.active_weight_text.set(self.parent.parent.active_weight_message)
                            self.parent.parent.active_intensity_message = '%-s %5.3f'%("height", avg_height)
                            self.parent.parent.active_intensity_text.set(self.parent.parent.active_intensity_message)
                            if self.parent.parent.fast_plot_type == 'areas':
                                self.parent.parent._draw_fast_plot(peak_profile.get_avg_areas())
                            elif self.parent.parent.fast_plot_type == 'heights':
                                self.parent.parent._draw_fast_plot(peak_profile.get_avg_heights())
                            if key == 'unmodified':
                                pass
                            elif key == 'modified':
                                pass
                            found_it = 1
                            break
                    if found_it:
                        break
                if not found_it:
                    # see if it's a selected region
                    tags = self.plot.find_withtag('select-range-hilite')
                    if id[0] in tags:
                        self.parent.parent._fastplot_max_in_selected_range()

    def clear_position_text(self, event):
        try:
            message = 'x:      y:     '
            self.parent.parent.position_message = message
            self.parent.parent.position_text.set(self.parent.parent.position_message)
        except AttributeError:
            pass

    def _delete_modified_possible_fragments(self, id):
        for tag in self.visible_modification_ids[id]:
            self.plot.delete(tag)
        del self.visible_modification_ids[id]
            
    def _draw_modified_possible_fragments(self, i):
        all_frags = self.fragment_line_weights+self.all_fragment_line_weights
        frag = all_frags[i]
        self.mod_fragment_line_ids = []
        self.mod_fragment_line_weights = []
        collected_tags = []
        if frag < self.x[-1]:
            results = self.parent.parent.experiment.get_fragment_for_weight_from_all(frag)
            if len(results) > 1:
                # disambiguate. get the modifying reagent
                return
            sequence = self.parent.parent.experiment.get_fragment_for_weight_from_all(frag)[0]['sequence']
            protease = self.parent.parent.experiment.get_protease()
            reagent = self.parent.parent.experiment.get_modifying_reagent()
            
            reagent_entry = self.parent.parent.parent.reactions_dict[reagent]
            reagent_specificity = reagent_entry['target_AA']
            reagent_weight      = reagent_entry['added_weight']
            reaction_hits = 0

            if protease == 'trypsin' and reagent == 'NHS':
                sequence = sequence[:-1]

            for aa in sequence[:-1]:
                if aa in reagent_specificity:
                    reaction_hits += 1
                    r_wt = string.atof(reagent_weight) * reaction_hits
                    xval = self.x_scale*(frag+r_wt)
                    lin = self.plot.create_line(xval, -4, xval, 2, width=2, fill='lime green', tags="mod_fragments")
                    collected_tags.append(lin)
                    self.mod_fragment_line_ids.append(lin)
                    self.mod_fragment_line_weights.append(frag+r_wt)
                    self.plot.tkraise(lin)
        for tag in collected_tags:
            self.plot.scale(tag, 0, 0, 1, -0.8)
            self.plot.move(tag, 0, (self.height * 0.82))
        return collected_tags
    
    def print_visible_modifications_report(self):
        # modification_ids are sets of canvas ids for fragments that have been selected
        saved_intensities = []
        for id in self.visible_modification_ids:
            for i in range(len(self.fragment_line_ids)):
                if (id == self.fragment_line_ids[i]):                        # if this is the line
                    weight = self.fragment_line_weights[i]
                    # find the peak that is closest to this weight
                    min_distance = 10000
                    closest_peak = -1
                    for peak_ind in self.plot_object.peaks:
                        if abs(self.plot_object.x[peak_ind] - weight) < min_distance:
                            closest_peak = peak_ind
                            min_distance = abs(self.plot_object.x[peak_ind]-weight)
                    saved_intensities.append([self.plot_object.y[closest_peak]])

                    frag_results = self.parent.parent.experiment.get_fragment_for_weight(self.fragment_line_weights[i])
                    if len(frag_results) == 1:
                        sequence = frag_results[0]['sequence']
                        reagent_entry = self.parent.parent.parent.reactions_dict[self.parent.parent.experiment.get_modifying_reagent()]
                        reagent_specificity = reagent_entry['target_AA']
                        reagent_weight      = reagent_entry['added_weight']
                        reaction_hits = 0
                        for aa in sequence:
                            if aa in reagent_specificity:
                                reaction_hits += 1
                                r_wt = string.atof(reagent_weight) * reaction_hits
                                xval = self.fragment_line_weights[i]+r_wt
                                min_distance = 10000
                                closest_point = -1
                                for plot_ind in range(len(self.plot_object.x)):
                                    if abs(self.plot_object.x[plot_ind] - xval) < min_distance:
                                        closest_point = plot_ind
                                        min_distance = abs(self.plot_object.x[plot_ind]-xval)
                                saved_intensities[-1].append(self.plot_object.y[closest_point])

    def plot_to_postscript(self, event, filename='temp.ps'):
        print 'drawing to file %s'%(filename)
        self.plot.postscript(file=filename)
        
    def _calculate_peak_height_sum(self, unmodified_peak_indices, modified_peak_indices, filter_resolution):
        sum_intensity = 0.0
        for peak in self.plot_object.peaks:
            xval = self.plot_object.x[peak]
            foundit = 0
            for upeak in unmodified_peak_indices:
                xval_u = self.plot_object.x[upeak]
                if abs(xval-xval_u)/xval < filter_resolution:
                    sum_intensity += self.plot_object.y[peak]
                    foundit = 1
                    break
            for mpeak_set in modified_peak_indices:
                for mpeak in mpeak_set:
                    xval_m = self.plot_object.x[mpeak]
                    if abs(xval-xval_m)/xval < filter_resolution:
                        sum_intensity += self.plot_object.y[peak]
        return sum_intensity
    
    def _calculate_height_sum(self, ranges):
        sum_intensity = 0.0
        for r in ranges:
            for i in range(r[0], r[1]):
                #if self.plot_object.y[i] > 1.0:
                sum_intensity += self.plot_object.y[i]
        return sum_intensity
    

        
