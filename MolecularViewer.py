# python imports
import sys
import os
import os.path
sys.path.append(os.getcwd())
import copy
import math
import time
import random
from Tkinter import *
# tool imports
import SPADE
sys.path.append('./Tools/DetectDomains')
import DetectDomains
sys.path.append('./Tools/Selection')
import SystemSelectionDialog
sys.path.append('./Tools/Aligner')
import SequenceAligner
sys.path.append('./Tools/ConservationTools')
import ConservationTools

# handle command line the lame way
homedir = '.'
from sys import argv         
#if len(argv) == 2:
homedir = 'C:\Users\Dude\Desktop\SPADE'

sys.path.append(homedir)

# dependency imports
#from vtk import *
import vtk.tk.vtkTkRenderWidget
import tkFileDialog
#from vtk.tk.vtkTkRenderWindowInteractor import *
from MolecularComponents.classFutamuraHash import FutamuraHash 
import vtk

# internal imports
sys.path.append(os.getcwd())
sys.path.append(os.path.join(homedir, './Dependencies'))
import Pmw
import MolecularSystem
import parms
import SystemMemory
#sys.path.append(os.path.join(homedir, './Tools/SequenceFetcher'))
#import SequenceFetcher

verbose = 0
import string

class MolecularViewer(Frame):             # Molecular Viewer
    def __init__(self, parent, system, ht=400, wd=500, menu=1, skip_update=1, undoredo=1):
        Frame.__init__(self, parent)
        self.system = system
        self.parent = parent
        self.undoredo_toggle = undoredo # if 1, undo/redo activated... slow but useful
        self.axes_visible = 0
        self.axes_actor = 'None'
        self.menuBar = Pmw.MenuBar(self.parent, hull_relief = 'raised',hull_borderwidth = 4)
        self.has_menu = 0
        if menu == 1:
            self.has_menu = 1
            self._build_menu()
        #self.screen = vtk.tk.vtkTkRenderWindowInteractor(parent, width=wd, height=ht)
        self.screen = vtk.tk.vtkTkRenderWidget.vtkTkRenderWidget(parent, width=wd, height=ht)
        win = self.screen.GetRenderWindow()
        # these look interesting, but I didn't see any difference with them
        #win.PointSmoothingOn()
        #win.LineSmoothingOn()
        #win.PolygonSmoothingOn()
        # add a busy signal
        self.busy_text      = StringVar()
        self.busy_text.set('Loading')
        self.busy_label     = Label(parent, textvariable=self.busy_text, width=10, height=1)
        #
        self.renderer = vtk.vtkRenderer()
        self.renderer.GetActiveCamera().GlobalWarningDisplayOff()
        self.screen.GetRenderWindow().AddRenderer(self.renderer)
        self.menuBar.pack(fill='x', side='top', expand=NO)
        self.screen.pack(side=TOP, expand=YES, fill=BOTH)
        self.busy_label.pack(side=TOP, anchor=W, expand=NO)
        start_time = time.clock()
        print 'creating the view  ',
        self.visitor = GraphicsVisitor(None, self)
        self.renderer.ResetCamera()
        self.screen.Render()
        #self._set_lights()
        if self.undoredo_toggle:
            if self.system and self.system != 'None':
                self.graphics_memories = [SystemMemory.GraphicsMemory(self.system)]
                self.selection_memories = [SystemMemory.SelectionMemory(self.system)]
        end_time = time.clock()
        self.graphics_memories = []
        self.selection_memories = []
        self.graphics_memory_index = 0
        self.selection_memory_index = 0
        self.has_icon_file = 0
        print 'done creating view%5.3f seconds'%(end_time - start_time)
        if self.system and self.system != 'None':
            self.loadSystem(self.system, skip_update)
        self.busy_text.set('Ready')
        self.cb = None      # codebox
        self.menu_select_mode = 0
        self.white_background = 0
    
    def _set_lights(self):
        cam = self.renderer.GetActiveCamera()

        light_collection = self.renderer.GetLights()
        light_collection.RemoveAllItems()        

        red_light = vtk.vtkLight()
        red_light.SetColor(1.0,0.5,0.5)
        red_light.SetPosition(1000,25,25)
        red_light.SetFocalPoint(25,25,25)
        red_light.SetIntensity(0.5)

        green_light = vtk.vtkLight()
        green_light.SetColor(0.5,1.0,0.5)
        green_light.SetPosition(25,1000,25)
        green_light.SetFocalPoint(25,25,25)
        green_light.SetIntensity(0.5)

        blue_light = vtk.vtkLight()
        blue_light.SetColor(0.5,0.5,1.0)
        blue_light.SetPosition(25,25,1000)
        blue_light.SetFocalPoint(25,25,25)
        blue_light.SetIntensity(0.5)

        #light1.PositionalOn()
        #light1.SetLightTypeToCameraLight()
        #light1.SetConeAngle(20)

        self.renderer.AddLight(red_light)
        self.renderer.AddLight(green_light)
        self.renderer.AddLight(blue_light)
        #self.renderer.LightFollowCameraOn()

    def close_system(self):
        self.closeSystem()
    
    def closeSystem(self):
        # first remove any actors previously created
        actors = self.renderer.GetActors()
        cnt = actors.GetNumberOfItems()
        for i in range(0,cnt):
            self.renderer.RemoveActor(actors.GetLastActor())
        self.renderer.ResetCamera()
        self.screen.Render()
        
    def load_system(self, system, skip_update=0):
        self.loadSystem(system, skip_update)
        
    def loadSystem(self, system, skip_update=0):
        self.system = system
        # first remove any actors previously created
        actors = self.renderer.GetActors()
        cnt = actors.GetNumberOfItems()
        for i in range(0,cnt):
            self.renderer.RemoveActor(actors.GetLastActor())
        print 'graphics visitor'
        #self.visitor = GraphicsVisitor(system, self)
        self.visitor.load_system(system, 0)
        print 'done visitor'
        self.renderer.ResetCamera()
        # if no icon has been created, create now
        icon_filename = self.system.get_filename_by_extension('jpg')
        if not os.path.exists(icon_filename):
            self.create_icon()
        # the viewer needs a memory for surfaces
        self.graphics_memories = [SystemMemory.GraphicsMemory(self.system)]
        self.selection_memories = [SystemMemory.SelectionMemory(self.system)]
        self.graphics_memory_index = 0
        self.selection_memory_index = 0
        if skip_update == 0:
            self.update_view()
        if self.has_menu:
            self.rebuild_color_menu()
        
    def update_view(self, graphics_state=None):
        if graphics_state:
            graphics_state.restore_system()
        start_time = time.clock()
        self.busy_text.set('Busy')
        self.busy_label.update()
        last_selection_memory = self.selection_memories[self.selection_memory_index]
        last_graphics_memory  = self.graphics_memories[self.graphics_memory_index]
        print 'updating the view  ',
        self.visitor.visit(last_selection_memory, last_graphics_memory)
        if self.undoredo_toggle:
            # undo/redo stuff can be disabled because of its slowness
            self.graphics_memories = self.graphics_memories[:self.graphics_memory_index+1]    # erase any forward memories
            self.selection_memories = self.selection_memories[:self.selection_memory_index+1]
            new_graphics_memory = SystemMemory.GraphicsMemory(self.system)
            new_selection_memory = SystemMemory.SelectionMemory(self.system)
            self.graphics_memories.append(new_graphics_memory)
            self.selection_memories.append(new_selection_memory)
            self.graphics_memory_index += 1
            self.selection_memory_index += 1
        self.screen.Render()
        end_time = time.clock()
        print 'done %5.3f seconds'%(end_time - start_time)
        self.busy_text.set('Ready')
        self.busy_label.update()

    def restore_initial_graphics_state(self):
        self.undo(1)
        
    def undo(self, load_initial=0):
        if self.graphics_memory_index == 0:
            print 'out of do\'s to undo'
            return
        else:
            # don't erase the first memory -- that's the initial state
            if load_initial==1:
                self.graphics_memory_index=1
                self.selection_memory_index=1
            # save the current graphic and selection images
            current_graphics_snap = SystemMemory.GraphicsMemory(self.system)
            current_selection_snap = SystemMemory.SelectionMemory(self.system)
            # select all of the molecular objects so that the graphics image is right
            SystemMemory.ObjectSelectionVisitor(self.system)
            # restore the last graphics image
            reload_graphics_memory = self.graphics_memories[self.graphics_memory_index-1]
            reload_graphics_memory.restore_system()
            reload_selection_memory = self.selection_memories[self.selection_memory_index-1]
            reload_selection_memory.restore_system()
            
            self.visitor.visit(current_selection_snap, current_graphics_snap)
            self.graphics_memory_index -= 1
            self.selection_memory_index -= 1
            self.screen.Render()
            
    def redo(self):
        if len(self.graphics_memories) > self.graphics_memory_index+1:
            current_graphics_snap = SystemMemory.GraphicsMemory(self.system)
            current_selection_snap = SystemMemory.SelectionMemory(self.system)
            reload_graphics_memory = self.graphics_memories[self.graphics_memory_index + 1]
            reload_selection_memory = self.selection_memories[self.selection_memory_index + 1]
            reload_graphics_memory.restore_system()
            reload_selection_memory.restore_system()
            self.visitor.visit(current_selection_snap, current_graphics_snap)
            self.graphics_memory_index += 1
            self.selection_memory_index += 1
        self.screen.Render()
        
    def _build_menu(self):
        self.menuBar.addmenu('Viewer', 'Viewer Controls')
        self.menuBar.addmenuitem('Viewer', 'command', 'Undo', command=self.undo, label='Undo')
        self.menuBar.addmenuitem('Viewer', 'command', 'Redo', command=self.redo, label='Redo')
        c_lambda = lambda: self.undo(1)
        self.menuBar.addmenuitem('Viewer', 'command', 'Restore Initial Memory', command=c_lambda, label='Restore Initial Memory')
        self.menuBar.addmenuitem('Viewer', 'command', 'Save Image', command=self.save_image, label='Save Image')
        self.menuBar.addmenuitem('Viewer', 'command', 'Toggle Axes', command=self.toggle_axes, label='View Axes')
        self.menuBar.addmenuitem('Viewer', 'command', 'Update View', command=self.update_view, label='Update View')
        self.menuBar.addmenuitem('Viewer', 'command', 'Update Icon', command=self.create_icon, label='Update Icon')
        self.menuBar.addmenuitem('Viewer', 'command', 'Toggle Prompt', command=self.toggle_codebox, label='Toggle Prompt')
        self.menuBar.addmenuitem('Viewer', 'command', 'Toggle Undo/Redo Ability', command=self.toggle_undoredo, label='Toggle Undo/Redo Ability')
        self.menuBar.addmenuitem('Viewer', 'command', 'Toggle Background', command=self.toggle_background, label='Toggle Background')
        self.menuBar.addmenuitem('Viewer', 'command', 'Toggle Hydrogens', command=self.toggle_hydrogens, label='Toggle Hydrogens')

        self.menuBar.addmenu('Calculate', 'Calculate Features')
        self.menuBar.addmenuitem('Calculate', 'command', 'Calculate Accessibility', command=self.calculate_accessibility, label='Calculate SA')
        self.menuBar.addmenuitem('Calculate', 'command', 'Calculate Shielding', command=self.calculate_shielding, label='Calculate Shielding')
        self.menuBar.addmenuitem('Calculate', 'command', 'Calculate Atomic Densities', command=self.calculate_atomic_densities, label='Calculate Atomic Densities')
        self.menuBar.addmenuitem('Calculate', 'command', 'Calculate Domains', command=self.calculate_domains, label='Calculate Domains')
        self.menuBar.addmenuitem('Calculate', 'command', 'Calculate Conservation', command=self.calculate_conservation, label='Calculate Conservation')
        self.menuBar.addmenuitem('Calculate', 'command', 'Calculate Interface Stats', command=self.calculate_interface_statistics, label='Calculate Interface Stats')
        self.menuBar.addmenu('Display', 'Change display')
        self.menuBar.addmenuitem('Display', 'command', 'toggle select mode', command=self.toggle_menu_select_mode, label='toggle select mode')
        self.menuBar.addcascademenu('Display', 'Display Atoms')
        c_lambda = lambda: self.display('atoms', 1)
        self.menuBar.addmenuitem('Display Atoms','command','Atoms display on', command=c_lambda, label='on')
        c_lambda = lambda: self.display('atoms', 0)
        self.menuBar.addmenuitem('Display Atoms','command','Atoms display off', command=c_lambda, label='off')
        self.menuBar.addcascademenu('Display', 'Display Trace')
        c_lambda = lambda: self.display('trace', 1)
        self.menuBar.addmenuitem('Display Trace','command','Trace display on', command=c_lambda, label='on')
        c_lambda = lambda: self.display('trace', 0)
        self.menuBar.addmenuitem('Display Trace','command','Trace display off', command=c_lambda, label='off')
        self.menuBar.addcascademenu('Display', 'Display Volumes')
        c_lambda = lambda: self.display('volume', 1)
        self.menuBar.addmenuitem('Display Volumes','command','Volumes display on', command=c_lambda, label='on')
        c_lambda = lambda: self.display('volume', 0)
        self.menuBar.addmenuitem('Display Volumes','command','Volumes display off', command=c_lambda, label='off')
        self.menuBar.addcascademenu('Display', 'Display HBonds')
        c_lambda = lambda: self.display('hbonds', 1)
        self.menuBar.addmenuitem('Display HBonds','command','Hbonds display on', command=c_lambda, label='on')
        c_lambda = lambda: self.display('hbonds', 0)
        self.menuBar.addmenuitem('Display HBonds','command','Hbonds display off', command=c_lambda, label='off')

        self.menuBar.addmenu('Color', 'Apply coloring scheme')
        self._add_color_menu()

        self.menuBar.addmenu('Options', 'Change Options')
        self.menuBar.addmenuitem('Options', 'command', 'toggle select mode', command=self.toggle_menu_select_mode, label='toggle select mode')
        self.menuBar.addcascademenu('Options', 'Atoms Options')
        self.menuBar.addmenuitem('Atoms Options', 'command', 'Show wireframe representation', command=self.represent_atoms_as_wireframe, label='wireframe')
        self.menuBar.addmenuitem('Atoms Options', 'command', 'Show sticks representation', command=self.represent_atoms_as_sticks, label='sticks')
        self.menuBar.addmenuitem('Atoms Options', 'command', 'Show sphere representation', command=self.represent_atoms_as_spheres, label='spheres')
        self.menuBar.addcascademenu('Options', 'Trace Options')
        self.menuBar.addmenuitem('Trace Options', 'command', 'Show line representation', command=self.represent_trace_as_line, label='line')
        self.menuBar.addmenuitem('Trace Options', 'command', 'Show tube representation', command=self.represent_trace_as_tube, label='tube')
        self.menuBar.addmenuitem('Trace Options', 'command', 'Toggle trace transparency', command=self.toggle_trace_transparency, label='transparency')
        self.menuBar.addcascademenu('Options', 'Volume Options')
        self.menuBar.addmenuitem('Volume Options', 'command', 'Show surface representation', command=self.represent_volumes_as_surface, label='surface')
        self.menuBar.addmenuitem('Volume Options', 'command', 'Show wireframe representation', command=self.represent_volumes_as_wireframe, label='wireframe')
        self.menuBar.addmenuitem('Volume Options', 'command', 'Show points representation', command=self.represent_volumes_as_points, label='points')
        self.menuBar.addmenuitem('Volume Options', 'command', 'Toggle volume transparency', command=self.toggle_volume_transparency, label='transparency')
        self.menuBar.addcascademenu('Options', 'H Bond Options')
        self.menuBar.addmenuitem('H Bond Options', 'command', 'Show lines representation', command=self.represent_hbonds_as_lines, label='lines')
        self.menuBar.addmenuitem('H Bond Options', 'command', 'Show tubes representation', command=self.represent_hbonds_as_tubes, label='tubes')
        
    def toggle_hydrogens(self):
        self.visitor.toggle_hydrogens()
        viewer.update_view()
        
    def toggle_background(self):
        if self.white_background:
            self.renderer.SetBackground([0.0,0.0,0.0])
        else:
            self.renderer.SetBackground([1.0,1.0,1.0])
        

    def toggle_menu_select_mode(self):
        if self.menu_select_mode:
            self.menu_select_mode=0
        else:
            self.menu_select_mode=1

    def toggle_trace_transparency(self):
        # define a function to pass to the selection dialog
        def exit_function(viewer, parms, top=None, snap=None):
            print 'executing exit function'
            for pchain in viewer.system.ProteinList:
                for res in pchain.residues:
                    if res.selected:
                        print 'res %s %s selected for opacity'%(res.res_type1, res.res_number)
                        res.vtk_arg_list['trace']['opacity'] = 0.2
                    else:
                        res.vtk_arg_list['trace']['opacity'] = 1.0
            graphics_state = SystemMemory.GraphicsMemory(viewer.system)

            viewer.update_view(graphics_state)
            if snap:
                #snap.restore_system()
                top.destroy()

        if self.menu_select_mode:
            # create a selection dialog to handle selection for transparency
            for pol in self.system.PolymerList:
                for res in pol.residues:
                    res.deselect()
            for pol in self.system.PolymerList:
                for res in pol.residues:
                    if res.vtk_arg_list['trace']['opacity'] < 1.0:
                        print 'opaque residue %s'%(res.res_number)
                        res.select()
                
            selection_top = Toplevel()
            geometry_string = "%dx%d%+d%+d" %(440,300,1,220) # width,height,x-offset,y-offset
            selection_top.geometry(geometry_string)
            selection_top.title('Select to toggle transparency')
            selection_top.wm_transient(self.parent) # transient to the System window
            current_selection_snap = SystemMemory.SelectionMemory(self.system)
            c_lambda = lambda v=self, p=parms, t=selection_top, s=current_selection_snap:exit_function(v, p, t, s)

            selection_dialog = SystemSelectionDialog.SystemSelectionDialog(self,
                                                                           selection_top,
                                                                           self.system,
                                                                           0,
                                                                           [],
                                                                           0,
                                                                           c_lambda)
            selection_dialog.pack(expand='yes', fill='both')
        else:
            for pchain in self.system.ProteinList:
                for res in pchain.residues:
                    if res.vtk_arg_list['trace']['opacity'] == 1.0:
                        res.vtk_arg_list['trace']['opacity'] = 0.2
                    else:
                        res.vtk_arg_list['trace']['opacity'] = 1.0
            self.update_view()

    def toggle_volume_transparency(self):
        # define a function to pass to the selection dialog
        def exit_function(viewer, parms, top=None, snap=None):
            print 'executing exit function'
            for pchain in viewer.system.ProteinList:
                for atom in pchain.atoms:
                    if atom.selected:
                        print 'atom %s %s selected for opacity'%(atom.atom_type, atom.atom_number)
                        atom.vtk_arg_list['volume']['opacity'] = 0.1
                    else:
                        atom.vtk_arg_list['volume']['opacity'] = 1.0
            graphics_state = SystemMemory.GraphicsMemory(viewer.system)

            viewer.update_view(graphics_state)
            if snap:
                #snap.restore_system()
                top.destroy()

        if self.menu_select_mode:
            # create a selection dialog to handle selection for transparency
            for pol in self.system.PolymerList:
                for atom in pol.atoms:
                    atom.deselect()
            for pol in self.system.PolymerList:
                for atom in pol.atoms:
                    if atom.vtk_arg_list['volume']['opacity'] < 1.0:
                        print 'opaque atom %s'%(atom.atom_number)
                        atom.select()
                
            selection_top = Toplevel()
            geometry_string = "%dx%d%+d%+d" %(440,300,1,220) # width,height,x-offset,y-offset
            selection_top.geometry(geometry_string)
            selection_top.title('Select to toggle transparency')
            selection_top.wm_transient(self.parent) # transient to the System window
            current_selection_snap = SystemMemory.SelectionMemory(self.system)
            c_lambda = lambda v=self, p=parms, t=selection_top, s=current_selection_snap:exit_function(v, p, t, s)

            selection_dialog = SystemSelectionDialog.SystemSelectionDialog(self,
                                                                           selection_top,
                                                                           self.system,
                                                                           0,
                                                                           [],
                                                                           0,
                                                                           c_lambda)
            selection_dialog.pack(expand='yes', fill='both')
        else:
            for pchain in self.system.ProteinList:
                for atom in pchain.atoms:
                    if atom.vtk_arg_list['volume']['opacity'] == 1.0:
                        atom.vtk_arg_list['volume']['opacity'] = 0.1
                    else:
                        atom.vtk_arg_list['volume']['opacity'] = 1.0
            self.update_view()
        

    def toggle_codebox(self):
        if self.cb == None:
            self.cb = SPADE.CodeBox(self.parent, self.system, self)
            sys.stdout = self.cb
            sys.stderr = self.cb
        else:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            self.cb.destroy()
            self.cb = None

    def rebuild_color_menu(self):
        self.menuBar.deletemenu('Color Atoms')
        self.menuBar.deletemenu('Color Trace')
        self.menuBar.deletemenu('Color Volumes')
        self._add_color_menu()
        
    def _add_color_menu(self):
        """ Automatically color by value
        All residues or atoms in each protein need to have the value. The value should
        be a digit. If not normalized, a new feature will be added which appends '_normalized'
        to the key """
        print 'adding color menu'
        self.menuBar.addcascademenu('Color', 'Color Atoms'); 
        c_lambda = lambda: self.color_wireframe('cpk');
        self.menuBar.addmenuitem('Color Atoms','command','Color wireframes cpk', command=c_lambda, label='cpk')
        c_lambda = lambda: self.color_wireframe('type');
        self.menuBar.addmenuitem('Color Atoms','command','Color wireframes by type', command=c_lambda, label='type')
        c_lambda = lambda: self.color_wireframe('chain');
        self.menuBar.addmenuitem('Color Atoms','command','color wireframes by chain', command=c_lambda, label='chain')
        c_lambda = lambda: self.color_wireframe('hydrogen_type');
        self.menuBar.addmenuitem('Color Atoms','command','color wireframes by H type', command=c_lambda, label='H Type')
        
        self.menuBar.addcascademenu('Color', 'Color Trace')
        self.menuBar.addmenuitem('Color Trace','command','Color tubes by secondary', command=self.color_trace_by_secondary,label='secondary')
        self.menuBar.addmenuitem('Color Trace','command','Color tubes by type', command=self.color_tubes_type,label='type')
        self.menuBar.addmenuitem('Color Trace','command','Color tubes by chain', command=self.color_tubes_chain,label='chain')

        self.menuBar.addcascademenu('Color', 'Color Volumes')
        self.menuBar.addmenuitem('Color Volumes','command','Color volumes cpk', command=self.color_volumes_cpk,label='cpk')
        self.menuBar.addmenuitem('Color Volumes','command','Color volumes by type', command=self.color_volumes_type,label='type')
        self.menuBar.addmenuitem('Color Volumes','command','Color volumes by chain', command=self.color_volumes_chain,label='chain')

        # create menu items for .features keys for atoms and residues
        if self.system != 'None' and self.system != None:
            key_store = {}
            key_store['atom'] = self.system.ProteinList[0].atoms[0].features.keys()
            key_store['residue'] = self.system.ProteinList[0].residues[0].features.keys()
            for run_type in ['atom', 'residue']:
                broken = 0
                for key in key_store[run_type]:
                    for pol in self.system.ProteinList:
                        if key == 'domain':
                            self.print_domain_info(pol)
                        normalized = 1
                        # if the feature includes non-digits, pass. if it is all digits, see if 
                        # it is normalized
                        if run_type == 'atom':
                            item_list = pol.atoms
                        elif run_type == 'residue':
                            item_list = pol.residues
                        same_val_count = 0
                        try:
                            item_list[0].features[key]
                        except KeyError:
                            continue
                        else:
                            first_val = item_list[0].features[key]
                        for item in item_list:
                            try:
                                feature = item.features[key]
                            except KeyError:
                                print 'key error on %s, breaking'%(key)
                                broken = 1
                                break
                            try:
                                int(feature)
                            except ValueError:
                                print '%s not digit, breaking'%(feature)
                                broken = 1
                                break
                            else:
                                if feature != -1 and (feature < 0.0 or feature > 1.0):
                                    normalized = 0
                            if feature == first_val:
                                same_val_count += 1
                            if same_val_count == len(item_list):
                                print '%s all the same value; breaking'%(key)
                                broken = 1
                                break
                            if key == 'domain':
                                if item.features[key] == 0.0:
                                    item.features[key] = -1
                        else:
                            # if not normalized, make a new key called key+'_normalized', and swap the old
                            # key with the new key to color by it
                            old_key = copy.copy(key)
                            if not normalized and (key+'_normalized' not in item.features.keys()):
                                min_f = 1000000
                                max_f = -1000000
                                for item2 in item_list:
                                    feature = item2.features[key]
                                    if feature != -1:
                                        if feature < min_f:
                                            min_f = feature
                                        if feature > max_f:
                                            max_f = feature
                                key = key + '_normalized'
                                for item2 in item_list:
                                    if item2.features[old_key] != -1.0:
                                        d = (item2.features[old_key]-min_f) / (max_f-min_f+0.0)
                                        item2.features[key] = d
                                    else:
                                        item2.features[key] = -1.0
                            if run_type == 'residue':
                                c_lambda1 = lambda p=pol, k=key: self.color_trace_by_residue_feature(p, k)
                                self.menuBar.addmenuitem('Color Trace','command','Color trace by res '+key, command=c_lambda1, label='%s %s'%(pol.chain_name, key))
                                c_lambda2 = lambda p=pol, k=key: self.color_volume_by_residue_feature(p, k)
                                self.menuBar.addmenuitem('Color Volumes','command','Color volumes by res '+key, command=c_lambda2, label='%s %s'%(pol.chain_name, key))
                                c_lambda3 = lambda p=pol, k=key: self.color_atoms_by_residue_feature(p, k)
                                self.menuBar.addmenuitem('Color Atoms','command','Color atoms by res '+key, command=c_lambda3, label='%s %s'%(pol.chain_name, key))
                            elif run_type == 'atom':
                                c_lambda1 = lambda p=pol, k=key: self.color_trace_by_atom_feature(p, k)
                                self.menuBar.addmenuitem('Color Trace','command','Color trace by atom '+key, command=c_lambda1, label='%s %s'%(pol.chain_name, key))
                                c_lambda2 = lambda p=pol, k=key: self.color_volume_by_atom_feature(p, k)
                                self.menuBar.addmenuitem('Color Volumes','command','Color volumes by atom '+key, command=c_lambda2, label='%s %s'%(pol.chain_name, key))
                                c_lambda3 = lambda p=pol, k=key: self.color_atoms_by_atom_feature(p, k)
                                self.menuBar.addmenuitem('Color Atoms','command','Color atoms by atom '+key, command=c_lambda3, label='%s %s'%(pol.chain_name, key))
                            key = old_key
                            #broken = 1
                            #break
                        if broken:
                            break

    def print_domain_info(self, pchain):
        last_domain = -2
        domains = []
        starts  = []
        ends    = []
        last_res_number = -1
        for res in pchain.residues:
            domain = res.features['domain']
            if domain > 0:
                if domain != last_domain:
                    starts.append(res.res_number)
                    if last_res_number != -1:
                        ends.append(last_res_number)
                    domains.append(domain)
                last_res_number = res.res_number
                last_domain = domain
        ends.append(len(pchain.residues))
        if len(domains) > 1:
            print 'chain %s has %s domains'%(pchain.chain_name, len(domains))
        else:
            print 'chain %s has 1 domain'%(pchain.chain_name)
        for i in range(len(domains)):
            print 'd%s %s - %s'%(domains[i], starts[i], ends[i])
            
    def toggle_undoredo(self):
        if self.undoredo_toggle:
            self.undoredo_toggle = 0
        else:
            self.undoredo_toggle = 1
            
    def draw_label(self, text, x, y, z):
        pass
        
    def color_atoms_by_atom_feature(self, pol, feature):
        for atom in pol.atoms:
            if atom.features[feature] == -1:
                color = [1.0,0.0,0.0]
            else:
                color = [0.0+atom.features[feature],  0.0+atom.features[feature], 1.0]
            atom.vtk_arg_list['atoms']['color'] = color
        self.update_view()

    def color_trace_by_atom_feature(self, pol, feature):
        averages_list = []
        for res in pol.residues:
            sum_feature = 0.0
            num = 0.0
            for atom in res.atoms:
                if atom.features[feature] == -1:
                    continue
                else:
                    sum_feature += atom.features[feature]
                    num += len(res.atoms)
            averages_list.append(sum_feature/num)
        # normalize
        maxval = 0.0
        minval = 100000.0
        for avg in averages_list:
            if avg > maxval:
                maxval = avg
            if avg < minval:
                minval = avg
        for ind in range(len(pol.residues)):
            val = (averages_list[ind] - minval)/(maxval - minval)
            color = [0.0+(val),  0.0+(val), 1.0]
            pol.residues[ind].vtk_arg_list['trace']['color'] = color
        self.update_view()

    def color_volume_by_atom_feature(self, pol, feature):
        for atom in pol.atoms:
            if atom.features[feature] == -1:
                color = [1.0,0.0,0.0]
            else:
                color = [0.0+atom.features[feature],  0.0+atom.features[feature], 1.0]
            atom.vtk_arg_list['volume']['color'] = color
        self.update_view()

    def color_atoms_by_residue_feature(self, pol, feature):
        for res in pol.residues:
            if res.features[feature] == -1:
                color = [1.0,0.0,0.0]
            else:
                color = [0.0+res.features[feature],  0.0+res.features[feature], 1.0]
            for atom in res.atoms:
                atom.vtk_arg_list['atoms']['color'] = color
        self.update_view()
        
    def color_trace_by_residue_feature(self, pol, feature):
        for res in pol.residues:
            if res.features[feature] == -1:
                color = [1.0,0.0,0.0]
            else:
                color = [0.0+res.features[feature],  0.0+res.features[feature], 1.0]
            res.vtk_arg_list['trace']['color'] = color
        self.update_view()

    def color_volume_by_residue_feature(self, pol, feature):
        for res in pol.residues:
            if res.features[feature] == -1:
                color = [1.0,0.0,0.0]
            else:
                color = [0.0+res.features[feature],  0.0+res.features[feature], 1.0]
            for atom in res.atoms:
                atom.vtk_arg_list['volume']['color'] = color
        self.update_view()
        
    def color_wireframe(self, style):
        if style == 'hydrogen_type':
            if len(self.system.HBonds) == 0:
                self.system.create_hbonds(None, None, True)
        chain_colors = parms.get('atom_color_list')
        if style == 'cpk' or style == 'type' or style == 'hydrogen_type':
            for mol in self.system.MoleculeList:
                if mol.selected:
                    for atom in mol.atoms:
                        atom.vtk_arg_list['atoms']['color'] = self.visitor.get_atom_color(atom, style)
        elif style == 'chain':
            for mol in self.system.MoleculeList:
                if mol.selected:
                    color = self.visitor.get_atom_color(mol.atoms[0], 'chain')
                    for atom in mol.atoms:
                        atom.vtk_arg_list['atoms']['color'] = color
        self.update_view()
        
    def color_tubes_type(self):
        for pchain in self.system.ProteinList:
            for res in pchain.residues:
                res.vtk_arg_list['trace']['color'] = self.visitor.get_residue_color(res, 'type')
        self.update_view()

    def color_trace_by_secondary(self):
        for pchain in self.system.ProteinList:
            # first make sure secondary is assigned
            for res in pchain.residues:
                try:
                    res.features['secondary']
                except KeyError:
                    pchain.assign_ss()
                    break
            for res in pchain.residues:
                if res.features['secondary'] == 'A':
                    res.vtk_arg_list['trace']['color'] = [1.0,0.1,0.1]
                elif res.features['secondary'] == 'B':
                    res.vtk_arg_list['trace']['color'] = [0.1,0.1,1.0]
                elif res.features['secondary'] == 'C':
                    res.vtk_arg_list['trace']['color'] = [0.8,0.8,0.8]
            self.update_view()
        
    def color_tubes_chain(self):
        for pchain in self.system.ProteinList:
            color = self.visitor.get_residue_color(pchain.residues[0], 'chain')
            for res in pchain.residues:
                res.vtk_arg_list['trace']['color'] = color
        for nchain in self.system.NucleotideChainList:
            color = self.visitor.get_residue_color(nchain.residues[0], 'chain')
            for res in nchain.residues:
                res.vtk_arg_list['trace']['color'] = color
        self.update_view()

    def color_volumes_chain(self):
        for pol in self.system.PolymerList:
            for atom in pol.atoms:
                atom.vtk_arg_list['volume']['color'] = self.visitor.get_atom_color(atom, 'chain')
        self.update_view()

    def color_volumes_type(self):
        for pol in self.system.PolymerList:
            for atom in pol.atoms:
                atom.vtk_arg_list['volume']['color'] = self.visitor.get_atom_color(atom, 'volume_type')
        self.update_view()

    def color_volumes_cpk(self):
        for pol in self.system.PolymerList:
            for atom in pol.atoms:
                atom.vtk_arg_list['volume']['color'] = self.visitor.get_atom_color(atom, 'cpk')
        self.update_view()

    def display(self, type, onoff):
        # define a function to pass to the selection dialog
        def exit_function(viewer, parms, top=None, snap=None):
            print 'executing exit function'
            listdict = {'system':self.system,
                        'trace':self.system.PolymerList,
                        'atoms':self.system.MoleculeList,
                        'volume':self.system.MoleculeList}
            for x in listdict[type]:
                if x.selected:
                    x.vtk_arg_list[type]['visualize'] = onoff

            graphics_state = SystemMemory.GraphicsMemory(viewer.system)

            viewer.update_view(graphics_state)
            if snap:
                snap.restore_system()
                top.destroy()

        listdict = {'system':self.system,
                    'trace':self.system.PolymerList,
                    'atoms':self.system.MoleculeList,
                    'volume':self.system.MoleculeList}

        if self.menu_select_mode:
            # create a selection dialog to handle selection for transparency
            current_selection_snap = SystemMemory.SelectionMemory(self.system)
            for x in listdict[type]:
                x.deselect()
            for x in listdict[type]:
                if x.vtk_arg_list[type]['visualize'] == 1.0:
                    x.select()
                
            selection_top = Toplevel()
            geometry_string = "%dx%d%+d%+d" %(440,300,1,220) # width,height,x-offset,y-offset
            selection_top.geometry(geometry_string)
            selection_top.title('Select to toggle %s display'%(type))
            selection_top.wm_transient(self.parent) # transient to the System window
            c_lambda = lambda v=self, p=parms, t=selection_top, s=current_selection_snap:exit_function(v, p, t, s)

            selection_dialog = SystemSelectionDialog.SystemSelectionDialog(self,
                                                                           selection_top,
                                                                           self.system,
                                                                           0,
                                                                           [],
                                                                           0,
                                                                           c_lambda)
            selection_dialog.pack(expand='yes', fill='both')
        else:
            listdict = {'hbonds':self.system,
                        'trace':self.system.PolymerList,
                        'atoms':self.system.MoleculeList,
                        'volume':self.system.MoleculeList}
            if type == 'hbonds':
                if len(self.system.HBonds) == 0:
                    self.system.create_hbonds(None, None, True)
                self.system.vtk_arg_list[type]['visualize'] = onoff
            else:
                for x in listdict[type]:
                    if x.selected:
                        x.vtk_arg_list[type]['visualize'] = onoff
            self.update_view()

    def represent_hbonds_as_lines(self):
        # caluclate hbonds if nec.
        if len(self.system.HBonds) == 0:
            self.system.create_hbonds(None, None, True)
        if self.system.selected:
            self.system.vtk_arg_list['hbonds']['representation'] = 'lines'
        self.update_view()

    def represent_hbonds_as_tubes(self):
        # calculate hbonds if nec.
        if len(self.system.HBonds) == 0:
            self.system.create_hbonds(None, None, True)
        if self.system.selected:
            self.system.vtk_arg_list['hbonds']['representation'] = 'tubes'
        self.update_view()
        
    def represent_atoms_as_wireframe(self):
        for mol in self.system.MoleculeList:
            if mol.selected:
                mol.vtk_arg_list['atoms']['representation'] = 'wireframe'
        self.update_view()

    def represent_atoms_as_sticks(self):
        for mol in self.system.MoleculeList:
            if mol.selected:
                mol.vtk_arg_list['atoms']['representation'] = 'sticks'
        self.update_view()
                
    def represent_atoms_as_spheres(self):
        for mol in self.system.MoleculeList:
            if mol.selected:
                mol.vtk_arg_list['atoms']['representation'] = 'spheres'
        self.update_view()

    def represent_trace_as_line(self):
        for mol in self.system.PolymerList:
            if mol.selected:
                mol.vtk_arg_list['trace']['representation'] = 'line'
        self.update_view()
        
    def represent_trace_as_tube(self):
        for mol in self.system.PolymerList:
            if mol.selected:
                mol.vtk_arg_list['trace']['representation'] = 'tube'
        self.update_view()

    def represent_trace_as_ribbon(self):
        for mol in self.system.PolymerList:
            if mol.selected:
                mol.vtk_arg_list['trace']['representation'] = 'ribbon'
        self.update_view()
        
    def represent_volumes_as_surface(self):
        for mol in self.system.MoleculeList:
            if mol.selected:
                mol.vtk_arg_list['volume']['representation'] = 'surface'
        self.update_view()

    def represent_volumes_as_wireframe(self):
        for mol in self.system.MoleculeList:
            if mol.selected:
                mol.vtk_arg_list['volume']['representation'] = 'wireframe'
        self.update_view()

    def represent_volumes_as_points(self):
        for mol in self.system.MoleculeList:
            if mol.selected:
                mol.vtk_arg_list['volume']['representation'] = 'points'
        self.update_view()
        
    
    def save_image(self):
        strg = tkFileDialog.asksaveasfilename(title = 'Save as', defaultextension='.jpg', filetypes=[("JPEG", "*.jpg"),
                                                                                      ("PostScript", "*.ps"),
                                                                                      ("PNG", "*.png"),
                                                                                      ("BMP", "*.bmp"),
                                                                                      ("TIFF", "*.tif"),
                                                                                      ("all files", "*")])
        if len(strg) > 0:
            tokens = string.split(strg, '.')
            if tokens[-1] == "jpg":
                w2i = vtk.vtkWindowToImageFilter()
                w2i.SetInput(self.screen.GetRenderWindow())
                writer = vtk.vtkJPEGWriter()
                writer.SetInput(w2i.GetOutput())
                writer.SetFileName(strg)
                writer.Write()
    def toggle_axes(self):
        if self.axes_visible == 0:
            # first figure out what the bounds are; ignore waters
            bounds = self.system.get_bounds()
            if self.axes_actor == 'None':
                self.axes_actor = vtk.vtkCubeAxesActor2D()
                self.axes_actor.SetCamera(self.renderer.GetActiveCamera())
                self.renderer.AddActor(self.axes_actor)
                self.axes_actor.SetFlyMode(0)
                self.axes_actor.SetInertia(10)
                self.axes_actor.SetNumberOfLabels(5)
                self.axes_actor.ScalingOff()
                self.axes_actor.GetProperty().SetColor([0.2,1.0,0.2])
                self.axes_actor.GetProperty().SetOpacity(0.6)
            self.axes_actor.SetBounds(bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5])
            self.axes_visible = 1
            self.axes_actor.XAxisVisibilityOn()
            self.axes_actor.YAxisVisibilityOn()
            self.axes_actor.ZAxisVisibilityOn()
        else:
            self.axes_visible = 0
            self.axes_actor.XAxisVisibilityOff()
            self.axes_actor.YAxisVisibilityOff()
            self.axes_actor.ZAxisVisibilityOff()
        self.screen.Render()
    def create_icon(self):
        icon_filename = self.system.get_filename_by_extension('jpg')
        w2i = vtk.vtkWindowToImageFilter()
        w2i.SetInput(self.screen.GetRenderWindow())
        writer = vtk.vtkJPEGWriter()
        writer.SetInput(w2i.GetOutput())
        writer.SetFileName(icon_filename)
        writer.Write()
                
    def calculate_accessibility(self):
        # i dont see any use in calculating this for nucleotide chains or small molecules,
        # but similar code will apply
        self.system.calculate_differential_system_asa(1.4, 1000)
        self.rebuild_color_menu()

    def calculate_interface_statistics(self):
        ConservationTools.calculate_conservation(self.system, 'sidechain_asa')
        self.system.calculate_interface_significance('sidechain_asa')

        #a = SequenceAligner.SequenceAligner(0, 'local')
        #a.add_target(self.system.ProteinList[0])
        #a.add_template(self.system.ProteinList[1])
        #pid = a.align_sequences()
        #print '%5.3f percent identity between chains %s and %s'%(pid, self.system.ProteinList[0].chain_name, self.system.ProteinList[1].chain_name)
        
        self.rebuild_color_menu()

    def calculate_shielding(self):
        for pchain in self.system.ProteinList:
            pchain.fill_pseudo_sidechains()
            pchain.fill_neighbors_lists()
        self.rebuild_color_menu()

    def calculate_atomic_densities(self):
        for pchain in self.system.ProteinList:
            pchain.fill_densities()
        self.rebuild_color_menu()

    def calculate_domains(self):
        #fetcher = SequenceFetcher.SequenceFetcher(self.system)
        #fetcher.fetch()
        #fetcher.normalize_domain_labels()
        for pchain in self.system.ProteinList:
            DetectDomains.auto_decompose(pchain, self)
        self.rebuild_color_menu()
        
    def calculate_conservation(self):
        """ currently designed to read Molnir .msq files,
            through the Protein object's calculate_conservation function
        """
        for pchain in self.system.ProteinList:
            ConservationTools.fetch_msq_conservation(pchain)
        ConservationTools.calculate_conservation(self.system)
        self.rebuild_color_menu()

class GraphicsVisitor:
    """ GraphicsVisitor visits the structure, using vtk to create actors """
    def __init__(self, system, parent):
        self.renderer = parent.renderer
        self.screen = parent.screen
        self.parent = parent
        self.polymer_actors        = {}         # append lists of atom_actors, one for each chain
        self.small_molecule_actors = {}
        self.polymer_splines       = {}
        self.system = system
        self.hydrogens_on          = 0
        if system != 'None' and system != None:
            self.load_system(system, 1)
    """ System tools: load, visits """
    def toggle_hydrogens(self):
        if self.hydrogens_on:
            self.hydrogens_on = 0
        else:
            self.hydrogens_on = 1

        for pol in self.system.PolymerList:
            # atoms
            self.renderer.RemoveActor(self.polymer_actors[pol.key]['atoms'])
            pol.vtk_arg_list['atoms']['currently_on'] = 0
            pol.vtk_arg_list['atoms']['visualize'] = 1
            self.polymer_actors[pol.key]['atoms'] = None

            # surface
            if pol.vtk_arg_list['volume']['currently_on']:
                self.renderer.RemoveActor(self.polymer_actors[pol.key]['volume'])
                pol.vtk_arg_list['volume']['currently_on'] = 0
                pol.vtk_arg_list['volume']['visualize'] = 1
                self.polymer_actors[pol.key]['volume'] = None                

    def load_system(self, system, initialize_chain_colors=1):
        self.system=system
        scheme = 'cpk'
        renderer = self.renderer
        # the following lists and dictionaries end up holding all of the actors. The instantiation
        # of the actors is postponed til the first rendering of that type, so as to speed up initial
        # loading
        self.system_actors = {'hbonds':None}
        for pol in system.PolymerList:
            self.polymer_actors[pol.key] = {}
            for style in ['volume','trace','atoms']:
                self.polymer_actors[pol.key][style] = None
            if initialize_chain_colors:
                chain_color = pol.get_chain_color()
                for atom in pol.atoms:
                    atom.vtk_arg_list['volume']['color'] = chain_color
                    atom.vtk_arg_list['atoms']['color'] = chain_color
                for res in pol.residues:
                    res.vtk_arg_list['trace']['color'] = chain_color
        for smol in system.SmallMoleculeList:
            self.small_molecule_actors[smol.key] = {}
            for style in ['atoms', 'volume']:
                self.small_molecule_actors[smol.key][style] = None
        # scale graphics args by the size of the system
        l1 = [12,14,18,24]
        l2 = [12,14,18,24]
        l3 = [0.4,0.6,0.8,1.0]
        l4 = [30,40,50,60]

        if self.system.atom_count > 8000:
            i = 0
        elif self.system.atom_count > 4000:
            i = 1
        elif self.system.atom_count > 2000:
            i = 2
        else:
            i = 3

        self.system.vtk_arg_list['hbonds']['sides']          = l2[i]
        self.system.vtk_arg_list['hbonds']['specular']       = l3[i]
        self.system.vtk_arg_list['hbonds']['specular_power'] = l4[i]
        for pol in self.system.PolymerList:
            pol.vtk_arg_list['trace']['splines']             = l1[i]
            pol.vtk_arg_list['trace']['sides']               = l2[i]
            pol.vtk_arg_list['trace']['specular']            = l3[i]
            pol.vtk_arg_list['trace']['specular_power']      = l4[i]
        self.initial_visit()

    def initial_visit(self):
        current_selection_snap = SystemMemory.SelectionMemory(self.system)
        # select all of the molecular objects so that the graphics image is right
        SystemMemory.ObjectSelectionVisitor(self.system)
        for style in ['hbonds']:
            if self.system.vtk_arg_list[style]['visualize']:
                self.system_actors[style] = self.draw_molecule(style, self.system, self.system_actors[style])
                self.system.vtk_arg_list[style]['currently_on'] = 0
        for pol in self.system.PolymerList:
            for style in ['atoms', 'volume', 'trace']:
                if pol.vtk_arg_list[style]['visualize']:
                    self.polymer_actors[pol.key][style] = self.draw_molecule(style, pol, self.polymer_actors[pol.key][style])
                    pol.vtk_arg_list[style]['currently_on'] = 1
        for smol in self.system.SmallMoleculeList:
            for style in ['atoms', 'volume']:
                if smol.vtk_arg_list[style]['visualize'] == 1:
                    self.small_molecule_actors[smol.key][style] = self.draw_molecule(style, smol, self.small_molecule_actors[smol.key][style])
                    smol.vtk_arg_list[style]['currently_on'] = 1
        # restore the selection state to loaded values
        self.system.visible = 0
        for chain in self.system.PolymerList:
            chain.visible = 1
            for res in chain.residues:
                res.visible = 0
        for smol in self.system.SmallMoleculeList:
            smol.visible = 1
        current_selection_snap.restore_system()

    def visit(self, selection_memory, graphics_memory):
        for style in ['hbonds']:
            if self.system.vtk_arg_list[style]['currently_on'] == 1 and self.system.vtk_arg_list[style]['visualize'] == 0:
                self.system_actors[style].VisibilityOff()
                self.system.vtk_arg_list[style]['currently_on'] = 0
            if self.system.vtk_arg_list[style]['currently_on'] == 0 and self.system.vtk_arg_list[style]['visualize'] == 1:
                if self.system_actors[style] == None:
                    self.system_actors[style] = self.draw_molecule(style, self.system, self.system_actors[style])
                else:
                    self.system_actors[style].VisibilityOn()
                self.system.vtk_arg_list[style]['currently_on'] = 1
            if self.system.vtk_arg_list[style]['currently_on'] == 1 and self.system.vtk_arg_list[style]['visualize'] == 1:
                # already on, just changing some features
                if self.system.vtk_arg_list[style]['representation'] != graphics_memory.get_representation(self.system, style):
                    self.system_actors[style] = self.draw_molecule(style, self.system, self.system_actors[style])
                mapper = self.system_actors[style].GetMapper()
                lu = self.build_lookup_table(self.system, style)
                mapper.SetLookupTable(lu)
                mapper.Update()
            
        for pol in self.system.PolymerList:
            for style in ['atoms','volume','trace']:
                if pol.vtk_arg_list[style]['currently_on'] == 1 and pol.vtk_arg_list[style]['visualize'] == 0:
                    self.polymer_actors[pol.key][style].VisibilityOff()
                    pol.vtk_arg_list[style]['currently_on'] = 0
                if pol.vtk_arg_list[style]['currently_on'] == 0 and pol.vtk_arg_list[style]['visualize'] == 1:
                    if self.polymer_actors[pol.key][style] == None:
                        self.polymer_actors[pol.key][style] = self.draw_molecule(style, pol, self.polymer_actors[pol.key][style])
                    else:
                        self.polymer_actors[pol.key][style].VisibilityOn()
                    pol.vtk_arg_list[style]['currently_on'] = 1
                if pol.vtk_arg_list[style]['currently_on'] == 1 and pol.vtk_arg_list[style]['visualize'] == 1:
                    # already on, just changing some features
                    if pol.vtk_arg_list[style]['representation'] != graphics_memory.get_representation(pol, style):
                        if style in ['atoms','trace']:
                            self.polymer_actors[pol.key][style] = self.draw_molecule(style, pol, self.polymer_actors[pol.key][style])
                        elif style == 'volume':
                            if pol.vtk_arg_list['volume']['representation'] == 'surface':
                                self.polymer_actors[pol.key][style].GetProperty().SetRepresentationToSurface()
                            elif pol.vtk_arg_list['volume']['representation'] == 'wireframe':
                                self.polymer_actors[pol.key][style].GetProperty().SetRepresentationToWireframe()
                            elif pol.vtk_arg_list['volume']['representation'] == 'points':
                                self.polymer_actors[pol.key][style].GetProperty().SetRepresentationToPoints()
                    mapper = self.polymer_actors[pol.key][style].GetMapper()
                    lu = self.build_lookup_table(pol, style)
                    mapper.SetLookupTable(lu)
                    mapper.Update()
                        
        for smol in self.system.SmallMoleculeList:
            if smol.selected:
                for style in ['atoms','volume']:
                    if smol.vtk_arg_list[style]['currently_on'] == 1 and smol.vtk_arg_list[style]['visualize'] == 0:
                        self.small_molecule_actors[smol.key][style].VisibilityOff()
                        smol.vtk_arg_list[style]['currently_on'] = 0
                    if smol.vtk_arg_list[style]['currently_on'] == 0 and smol.vtk_arg_list[style]['visualize'] == 1:
                        if self.small_molecule_actors[smol.key][style] == None:
                            self.small_molecule_actors[smol.key][style] = self.draw_molecule(style, smol, self.small_molecule_actors[smol.key][style])
                        else:
                            self.small_molecule_actors[smol.key][style].VisibilityOn()
                        smol.vtk_arg_list[style]['currently_on'] = 1
                    if smol.vtk_arg_list[style]['currently_on'] == 1 and smol.vtk_arg_list[style]['visualize'] == 1:
                        if smol.vtk_arg_list[style]['representation'] != graphics_memory.get_representation(smol, style):
                            if style == 'atoms':
                                self.small_molecule_actors[smol.key][style] = self.draw_molecule(style, smol, self.small_molecule_actors[smol.key][style])
                            elif style == 'volume':
                                if smol.vtk_arg_list['volume']['representation'] == 'surface':
                                    self.small_molecule_actors[smol.key][style].GetProperty().SetRepresentationToSurface()
                                elif smol.vtk_arg_list['volume']['representation'] == 'wireframe':
                                    self.small_molecule_actors[smol.key][style].GetProperty().SetRepresentationToWireframe()
                                elif smol.vtk_arg_list['volume']['representation'] == 'points':
                                    self.small_molecule_actors[smol.key][style].GetProperty().SetRepresentationToPoints()
                        mapper = self.small_molecule_actors[smol.key][style].GetMapper()
                        lu = self.build_lookup_table(smol, style)
                        mapper.SetLookupTable(lu)
                        mapper.Update()

    """  get color functions  """
    def get_atom_color(self, atom, scheme):
        chain_colors = ['chain1','chain2','chain3','chain4','chain5','chain6',
                       'chain7','chain8','chain9','chain10','chain11','chain12']
        colormap = { 'chain1':[0.25,0.25,0.90], 'chain2':[0.90,0.90,0.25], 'chain3':[0.90,0.25,0.25],
                     'chain4':[0.25,0.90,0.25], 'chain5':[0.25,0.90,0.90], 'chain6':[0.90,0.25,0.90],
                     'chain7':[0.05,0.05,0.50], 'chain8':[0.50,0.50,0.05], 'chain9':[0.50,0.05,0.05],
                     'chain10':[0.05,0.50,0.05],'chain11':[0.05,0.50,0.50],'chain12':[0.50,0.05,0.50],
                     'white':[1.00,1.00,1.00],  'blue':[0.00,0.00,1.00],   'red':[1.00,0.00,0.00],
                     'orange':[1,0.5,0],        'yellow':[1,1,0],          'cyan':[0,1,1],
                     'magenta':[1.0,0.0,0.5],   'grey':[0.8,0.8,0.8],      'green':[0.0,0.8,0.0],
                     'purple':[0.65,0.5,1.0]}
        
        """ return a color for the atom, based on the scheme.
        scheme can be 'cpk', 'type', or 'chain' """
        s = atom.atom_type[0]
        t = atom.atom_type
        r = atom.res_type
        
        
        rn = atom.res_number
        cn = atom.chain_name
        color = 'grey'
        isH = 0
        if scheme == 'cpk':
            if s == 'C' or s == 'H':
                color = 'white'
            elif s == 'N':
                color = 'blue'
            elif s == 'O':
                color = 'red'
            elif s == 'S':
                color = 'yellow'
            elif s == 'P':
                color = 'green'
            elif s == 'Z':
                color = 'magenta'
        elif scheme in ['type','volume_type','hydrogen_type']:               # these 'typ' rules are designed to apply to amino acids
            if t[0] == 'H':
                # find the closest atom from the hydrogen's "parent_molecule" data feature
                parent_molecule = atom.data['parent_molecule']
                closest_nonH = None
                closest_dist = 100000
                for at in parent_molecule.atoms:
                    if at.atom_type[0] != 'H':
                        d = atom.dist(at)
                        if d < closest_dist:
                            closest_dist = d
                            closest_nonH = at
                isH = 1
                color = self.get_atom_color(closest_nonH, scheme)

            if t=='C' or t=='CA':
                if scheme == 'type' or scheme == 'hydrogen_type':
                    color = 'blue'
                else:
                    color = 'cyan'
            elif t=='N' or t=='O':
                if atom.parent.res_type == 'PRO':
                    color='blue'
                if scheme == 'type' or scheme == 'hydrogen_type':
                    if scheme == 'hydrogen_type':
                        if len(self.system.HBonds):
                            for hbond in self.system.HBonds:
                                if atom.atom_number == hbond['accAtom'].atom_number or atom.atom_number == hbond['donorAtom'].atom_number:
                                    color = 'blue'
                                    break
                            else:
                                color = 'orange'
                        else:
                            color = 'blue'
                    else:
                        color = 'blue'
                else:
                    color = 'orange'
            elif t == 'CB':
                color = 'cyan'
            elif t=='CG' or t=='CG1' or t=='CG2':
                if r=='ASP':
                    color = 'red'
                elif r=='ASN':
                    color = 'orange'
                else:
                    color = 'cyan'
            elif t=='OG' or t=='OG1' or t=='OG2':
                color = 'orange'
            elif t=='SG':
                color = 'cyan'
            elif t=='CD' or t=='CD1' or t=='CD2':
                if r=='GLU':
                    color = 'red'
                elif r=='ARG' or r=='HIS':
                    color = 'yellow'
                elif r=='GLN':
                    color = 'orange'
                else:
                    color = 'cyan'
            elif t=='SD':
                color = 'cyan'
            elif t=='OD1' or t=='OD2' or t=='ND1' or t=='ND2':
                if r=='ASP':
                    color = 'red'
                elif r=='ASN' or r=='HIS':
                    color = 'orange'
            elif t=='CE' or t=='CE1' or t=='CE2' or t=='CE3':
                if r=='LYS' or r=='HIS':
                    color = 'yellow'
                else:
                    color = 'cyan'
            elif t=='OE1' or t=='OE2' or t=='NE' or t=='NE1' or t=='NE2':
                if r=='GLU':
                    color = 'red'
                elif r=='ARG' or r=='HIS':
                    color = 'yellow'
                elif r=='GLN' or r=='TRP':
                    color = 'orange'
            elif t=='CZ' or t=='CZ1' or t=='CZ2' or t=='CZ3':
                if r=='ARG':
                    color = 'yellow'
                else:
                    color = 'cyan'
            elif t=='CH' or t=='CH1' or t=='CH2' or t=='CH3':
                color = 'cyan'
            elif t=='NZ' or t=='NH1' or t=='NH2':
                color = 'yellow'
            elif t=='OH':
                color = 'orange'
        elif scheme == 'chain':
            color = 'grey'
            pidx = 0
            for pchain in self.system.ProteinList:
                if atom.chain_name == pchain.chain_name:
                    return colormap[chain_colors[pidx%len(chain_colors)]]
                pidx = pidx + 1
            for nchain in self.system.NucleotideChainList:
                if atom.chain_name == nchain.chain_name:
                    return colormap[chain_colors[pidx%len(chain_colors)]]
                pidx = pidx + 1
            for lig in self.system.LigandList:
                if atom.chain_name == lig.chain_name:
                    return colormap[chain_colors[pidx%len(chain_colors)]]
                pidx = pidx + 1
            for wat in self.system.WaterList:
                if atom.chain_name == wat.chain_name:
                    return colormap[chain_colors[pidx%len(chain_colors)]]
                pidx = pidx + 1
        if color in colormap.keys():
            return colormap[color]
        else:
            return color

    def get_residue_color(self, res, scheme):
        if scheme == 'type':
            buried = ['A','C','F','I','L','M','V','W']
            uncharged = ['G','N','P','Q','S','T','Y']
            chargedminus = ['D','E']
            chargedplus  = ['H','K','R']
            c = res.res_type1
            if c in buried:
                color = [0.5,0.75,1.0]
            elif c in uncharged:
                color = [1.0,0.5,0.0]
            elif c in chargedminus:
                color = [0.9,0.2,0.2]
            elif c in chargedplus:
                color = [1.0,1.0,0.2]
        elif scheme == 'chain':
            color = self.get_atom_color(res.atoms[0], 'chain')
        return color

    def get_polymer_color(self, pol, scheme):
        return self.get_atom_color(pol.atoms[0], 'chain')
    """  draw  """
    def draw_molecule(self, type, obj, actor):
        if actor == None:
            print 'creating new actor'
            actor = vtk.vtkActor()
            if type in ['trace','volume','hbonds']:
                actor.GetProperty().SetSpecular(obj.vtk_arg_list[type]['specular'])
                actor.GetProperty().SetSpecularPower(obj.vtk_arg_list[type]['specular_power'])
            self.renderer.AddActor(actor)
           
        if type == 'trace':
            if obj.vtk_arg_list['trace']['representation'] == 'line':
                return self._draw_lines_for_polymer(obj, actor)
            elif obj.vtk_arg_list['trace']['representation'] == 'tube':
                return self._draw_tube_for_polymer(obj, actor)
            elif obj.vtk_arg_list['trace']['representation'] == 'ribbon':
                return self._draw_ribbon_for_polymer(obj, actor)
        elif type == 'hbonds':
            if obj.vtk_arg_list['hbonds']['representation'] == 'lines':
                return self._draw_lines_for_hbonds(obj, actor)
            elif obj.vtk_arg_list['hbonds']['representation'] == 'tubes':
                print 'here4'
                return self._draw_tubes_for_hbonds(obj, actor)
        elif type == 'atoms':
            if obj.vtk_arg_list['atoms']['representation'] == 'wireframe':
                print 'about to draw wireframe'
                return self._draw_wireframe_for_molecule(obj, actor)
            elif obj.vtk_arg_list['atoms']['representation'] == 'sticks':
                return self._draw_sticks_for_molecule(obj, actor)
            elif obj.vtk_arg_list['atoms']['representation'] == 'spheres':
                return self._draw_spheres_for_molecule(obj, actor)
        elif type == 'volume':
            return self._draw_volume_for_molecule(obj, actor)
        print 'done %s'%(type)
        
    """  lookup tables  """
    def build_lookup_table(self, obj, type):
        if type == 'trace':
            return self._build_polymer_lookup_table(obj)
        elif type == 'atoms':
            return self._build_atom_lookup_table(obj)
        elif type == 'volume':
            return self._build_surface_lookup_table(obj)
        elif type == 'hbonds':
            return self._build_hbonds_lookup_table(obj)

    def _build_hbonds_lookup_table(self, system):
        luTable = vtk.vtkLookupTable()
        luTable.SetNumberOfTableValues(len(system.HBonds))
        for index in range(len(system.HBonds)):
            color = system.HBonds[index]['vtk_arg_list']['color']
            s = system.HBonds[index]['HBondInfo']['strength']
            luTable.SetTableValue(index*2, s, s, 1.0, 1.0)
            luTable.SetTableValue((index*2)+1, s, s, 1.0, 1.0)
        return luTable
    
    def _build_atom_lookup_table(self, molecule, deselected_opacity=None):
        if deselected_opacity==None:
            deselected_opacity = parms.get('deselected_opacity')

        luTable = vtk.vtkLookupTable()
        if not self.hydrogens_on:
            atom_count = 0
            for a in molecule.atoms:
                if a.atom_type[0] != 'H':
                    atom_count += 1
            luTable.SetNumberOfTableValues(atom_count)
        else:
            luTable.SetNumberOfTableValues(len(molecule.atoms))

        i = 0
        for index in range(len(molecule.atoms)):
            white_override = 0
            if not self.hydrogens_on and molecule.atoms[index].atom_type[0] == 'H':
                continue
            if not white_override:
                color = molecule.atoms[index].vtk_arg_list['atoms']['color']
                molecule.atoms[index].visible = 1
                opacity = molecule.atoms[index].vtk_arg_list['atoms']['opacity']
                luTable.SetTableValue(i, color[0], color[1], color[2], opacity)
            else:
                luTable.SetTableValue(i, 1.0,1.0,1.0,1.0)
            i += 1
        return luTable
    
    def _build_lines_lookup_table(self, polymer):
        return self._build_polymer_lookup_table(polymer)
    
    def _build_tubes_lookup_table(self, polymer):
        return self._build_polymer_lookup_table(polymer)

    def _build_ribbons_lookup_table(self, polymer):
        return self._build_polymer_lookup_table(polymer)

    def _build_polymer_lookup_table(self, polymer):
        deselected_opacity = parms.get('deselected_opacity')

        luTable = vtk.vtkLookupTable()
        luTable.SetNumberOfTableValues(len(polymer.residues))

        for index in range(len(polymer.residues)):
            white_override = 0
            if not white_override:
                res = polymer.residues[index]
                color = res.vtk_arg_list['trace']['color']
                res.visible = 1
                opacity = res.vtk_arg_list['trace']['opacity']
                if opacity != 1.0:
                    print 'res %s %s going opaque'%(res.res_type1, res.res_number)
                luTable.SetTableValue(index, color[0], color[1], color[2], opacity)
            else:
                luTable.SetTableValue(index, 1.0,1.0,1.0,1.0)
        return luTable
    
    def _build_surface_lookup_table(self, mol):
        deselected_opacity = parms.get('deselected_opacity')

        luTable = vtk.vtkLookupTable()
        luTable.SetNumberOfTableValues(len(mol.atoms))
        
        for index in range(len(mol.atoms)):
            white_override = 0
            if not white_override:
                color = mol.atoms[index].vtk_arg_list['volume']['color']
                opacity = mol.atoms[index].vtk_arg_list['volume']['opacity']
                luTable.SetTableValue(index, color[0], color[1], color[2], opacity)
            else:
                luTable.SetTableValue(index, 1.0,1.0,1.0,1.0)
        return luTable

    """  atom rendering  """
    def _draw_spheres_for_molecule(self, molecule, actor):
        atom_count = len(molecule.atoms)
        atom_color_list = parms.get('atom_color_list')
        
        # see if bond points need added, by checking the first atom from the molecule for them
        try:
            molecule.atoms[0].bonded_points
        except AttributeError:
            for atom in molecule.atoms:
                atom.bonded_points = atom.parent.get_bonds_list(atom)

        radii  = {'C':1.75,
                  'O':1.4,
                  'N':1.55,
                  'S':1.8,
                  'P':2.0,
                  'H':1.17,
                  'Z':3.0}
        default_distance = 1.8

        spheres = vtk.vtkAppendPolyData()
                
        scalars = vtk.vtkIntArray()
        scalars.SetNumberOfComponents(1)
        
        number_index_map = {}
        i = 0
        for ind in range(len(molecule.atoms)):
            if self.hydrogens_on == 0 and molecule.atoms[ind].atom_type[0] == 'H':
                continue
            else:
                number_index_map[molecule.atoms[ind].atom_number] = i
                i += 1

        i = 0
        points_count = 0
        atoms_used = 0
        for atom in molecule.atoms:
            if self.hydrogens_on == 0 and atom.atom_type[0] == 'H':
                continue
            atoms_used += 1
            sphere = vtk.vtkSphereSource()
            sphere.SetCenter(atom.x,atom.y,atom.z)
            sphere.SetRadius(radii.get(atom.atom_type[0], default_distance))
            sphere.SetThetaResolution(8)
            sphere.SetPhiResolution(8)
            x = sphere.GetOutput()
            x.Update()
            sphere_points = x.GetPoints()
            for point_index in range(sphere_points.GetNumberOfPoints()):
                scalars.InsertTuple1(i, number_index_map[atom.atom_number])            
                i = i + 1
            spheres.AddInput(x)
        x = spheres.GetOutput()
        x.Update()
        x.GetPointData().SetScalars(scalars)
        
        if i > 0:               
            luTable = self._build_atom_lookup_table(molecule)
            
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetLookupTable(luTable)
            mapper.SetInput(x)
            mapper.ScalarVisibilityOn()
            mapper.SetScalarModeToUsePointData()
            mapper.SetColorModeToMapScalars()
            mapper.SetScalarRange(0,atoms_used)

            actor.SetMapper(mapper)
        else:
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInput(x)
            actor.SetMapper(mapper)
        return actor

    def _draw_wireframe_for_molecule(self, molecule, actor):
        atom_count = len(molecule.atoms)
        atom_color_list = parms.get('atom_color_list')
        
        # see if bond points need added, by checking the first atom from the molecule for them
        #try:
        #    molecule.atoms[0].bonded_points
        #except AttributeError:
        # the above try clause doesn't allow dynamic hydrogens on/off
        for atom in molecule.atoms:
            atom.bonded_points = atom.parent.get_bonds_list(atom, self.hydrogens_on)

        bond_line_points = vtk.vtkPoints()
        bond_line_cells = vtk.vtkCellArray()
        
        scalars = vtk.vtkIntArray()
        scalars.SetNumberOfComponents(1)
        
        number_index_map = {}
        i = 0
        for ind in range(len(molecule.atoms)):
            if self.hydrogens_on == 0 and molecule.atoms[ind].atom_type[0] == 'H':
                continue
            else:
                number_index_map[molecule.atoms[ind].atom_number] = i
                i += 1

        i = 0
        points_count = 0
        atoms_used = 0
        for atom in molecule.atoms:
            if self.hydrogens_on == 0 and atom.atom_type[0] == 'H':
                continue
            atoms_used += 1
            points_count = points_count + len(atom.bonded_points)
            point_id = i
            bond_line_points.InsertPoint(i, atom.x, atom.y, atom.z)
            scalars.InsertTuple1(i, number_index_map[atom.atom_number])
            i=i+1
            for point in atom.bonded_points:
                bond_line_points.InsertPoint(i, point[0], point[1], point[2])
                scalars.InsertTuple1(i, number_index_map[atom.atom_number])
                bond_line_cells.InsertNextCell(2)
                bond_line_cells.InsertCellPoint(i)
                bond_line_cells.InsertCellPoint(point_id)
                i = i + 1

        bond_line_profile = vtk.vtkPolyData()
        bond_line_profile.SetPoints(bond_line_points)
        bond_line_profile.SetLines(bond_line_cells)
        bond_line_profile.GetPointData().SetScalars(scalars)
        
        if points_count > 0:               
            luTable = self._build_atom_lookup_table(molecule)
            
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetLookupTable(luTable)
            mapper.SetInput(bond_line_profile)
            mapper.ScalarVisibilityOn()
            mapper.SetScalarModeToUsePointData()
            mapper.SetColorModeToMapScalars()
            mapper.SetScalarRange(0,atoms_used)

            actor.SetMapper(mapper)
        else:
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInput(bond_line_profile)
            actor.SetMapper(mapper)
        return actor

    def _draw_lines_for_hbonds(self, system, actor):
        atom_count = len(system.HBonds)
        
        bond_line_points = vtk.vtkPoints()
        bond_line_cells = vtk.vtkCellArray()
        
        scalars = vtk.vtkIntArray()
        scalars.SetNumberOfComponents(1)

        i = 0
        points_count = 0
        for hbond in system.HBonds:
            points_count = points_count + 2
            point_id = i
            hyd = hbond['hydAtom']
            bond_line_points.InsertPoint(i, hyd.x, hyd.y, hyd.z)
            scalars.InsertTuple1(i, i)
            i=i+1
            don = hbond['accAtom']
            bond_line_points.InsertPoint(i, don.x, don.y, don.z)
            scalars.InsertTuple1(i, i)
            bond_line_cells.InsertNextCell(2)
            bond_line_cells.InsertCellPoint(i)
            bond_line_cells.InsertCellPoint(point_id)
            i = i + 1

        bond_line_profile = vtk.vtkPolyData()
        bond_line_profile.SetPoints(bond_line_points)
        bond_line_profile.SetLines(bond_line_cells)
        bond_line_profile.GetPointData().SetScalars(scalars)
        
        if points_count > 0:               
            luTable = self._build_hbonds_lookup_table(system)
            
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetLookupTable(luTable)
            mapper.SetInput(bond_line_profile)
            mapper.ScalarVisibilityOn()
            mapper.SetScalarModeToUsePointData()
            mapper.SetColorModeToMapScalars()
            mapper.SetScalarRange(0,len(system.HBonds))

            actor.SetMapper(mapper)
        else:
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInput(bond_line_profile)
            actor.SetMapper(mapper)
        return actor
        
    def _draw_tubes_for_hbonds(self, system, actor):
        print 'in draw sticks for hbonds'
        hbond_count = len(system.HBonds)

        hbond_line_points = vtk.vtkPoints()
        hbond_line_cells = vtk.vtkCellArray()

        scalars = vtk.vtkIntArray()
        scalars.SetNumberOfComponents(1)

        i = 0
        points_count = 0
        for hb in system.HBonds:
            points_count = points_count + 2
            point_id = i
            hyd = hb['hydAtom']
            don = hb['accAtom']
            hbond_line_points.InsertPoint(i, hyd.x, hyd.y, hyd.z)
            scalars.InsertTuple1(i, i)
            i=i+1
            hbond_line_points.InsertPoint(i, don.x, don.y, don.z)
            scalars.InsertTuple1(i, i)
            hbond_line_cells.InsertNextCell(2)
            hbond_line_cells.InsertCellPoint(i)
            hbond_line_cells.InsertCellPoint(point_id)
            i = i + 1

        hbond_line_profile = vtk.vtkPolyData()
        hbond_line_profile.SetPoints(hbond_line_points)
        hbond_line_profile.SetLines(hbond_line_cells)
        hbond_line_profile.GetPointData().SetScalars(scalars)
        
        if points_count > 0:               
            profileTubes = vtk.vtkTubeFilter()
            profileTubes.SetNumberOfSides(system.vtk_arg_list['hbonds']['sides'])
            profileTubes.SetInput(hbond_line_profile)
            profileTubes.SetRadius(system.vtk_arg_list['hbonds']['width'])

            luTable = self._build_hbonds_lookup_table(system)
            
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetLookupTable(luTable)
            mapper.SetInput(profileTubes.GetOutput())
            mapper.ScalarVisibilityOn()
            mapper.SetScalarModeToUsePointData()
            mapper.SetColorModeToMapScalars()
            mapper.SetScalarRange(0,len(system.HBonds))

            actor.SetMapper(mapper)
        else:
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInput(hbond_line_profile)
            actor.SetMapper(mapper)
        return actor
        
    def _draw_sticks_for_molecule(self, molecule, actor):
        atom_count = len(molecule.atoms)
        atom_color_list = parms.get('atom_color_list')
        
        # see if bond points need added, by checking the first atom from the molecule for them
        for atom in molecule.atoms:
            atom.bonded_points = atom.parent.get_bonds_list(atom, self.hydrogens_on)

        bond_line_points = vtk.vtkPoints()
        bond_line_cells = vtk.vtkCellArray()
        
        scalars = vtk.vtkIntArray()
        scalars.SetNumberOfComponents(1)
        
        number_index_map = {}
        for ind in range(len(molecule.atoms)):
            number_index_map[molecule.atoms[ind].atom_number] = ind

        i = 0
        points_count = 0
        for atom in molecule.atoms:
            points_count = points_count + len(atom.bonded_points)
            point_id = i
            bond_line_points.InsertPoint(i, atom.x, atom.y, atom.z)
            scalars.InsertTuple1(i, number_index_map[atom.atom_number])
            i=i+1
            for point in atom.bonded_points:
                bond_line_points.InsertPoint(i, point[0], point[1], point[2])
                scalars.InsertTuple1(i, number_index_map[atom.atom_number])
                bond_line_cells.InsertNextCell(2)
                bond_line_cells.InsertCellPoint(i)
                bond_line_cells.InsertCellPoint(point_id)
                i = i + 1

        bond_line_profile = vtk.vtkPolyData()
        bond_line_profile.SetPoints(bond_line_points)
        bond_line_profile.SetLines(bond_line_cells)
        bond_line_profile.GetPointData().SetScalars(scalars)
        
        if points_count > 0:               
            profileTubes = vtk.vtkTubeFilter()
            profileTubes.SetNumberOfSides(molecule.vtk_arg_list['atoms']['sticks_sides'])
            profileTubes.SetInput(bond_line_profile)
            profileTubes.SetRadius(molecule.vtk_arg_list['atoms']['sticks_width'])

            luTable = self._build_atom_lookup_table(molecule)
            
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetLookupTable(luTable)
            mapper.SetInput(profileTubes.GetOutput())
            mapper.ScalarVisibilityOn()
            mapper.SetScalarModeToUsePointData()
            mapper.SetColorModeToMapScalars()
            mapper.SetScalarRange(0,len(molecule.atoms))

            actor.SetMapper(mapper)
        else:
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInput(bond_line_profile)
            actor.SetMapper(mapper)
        return actor
        
    """  trace rendering  """
    def _draw_lines_for_polymer(self, polymer, actor):
        try:
            self.polymer_splines[polymer.chain_name]
        except:
            self._calculate_splines(polymer)

        splines,lines,scalars = self.get_polymer_splines_lines_scalars(polymer)

        data = vtk.vtkPolyData()
        data.SetPoints(splines)
        data.SetLines(lines)
        data.GetPointData().SetScalars(scalars)

        if scalars.GetNumberOfTuples() > 0:
            luTable = self._build_tubes_lookup_table(polymer)
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetLookupTable(luTable)
            mapper.SetInput(data)
            mapper.ScalarVisibilityOn()
            mapper.SetScalarModeToUsePointData()
            mapper.SetColorModeToMapScalars()
            #profileTubes.Update()
            mapper.SetScalarRange(0,len(polymer.residues))
        else:
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInput(data)
        
        actor.SetMapper(mapper)
        return actor


    def _draw_ribbon_for_polymer(self, polymer, actor):
        try:
            self.polymer_splines[polymer.chain_name]
        except:
            self._calculate_splines(polymer)

        splines,lines,scalars = self.get_polymer_splines_lines_scalars(polymer)

        profileData = vtk.vtkPolyData()
        profileData.SetPoints(splines)
        profileData.SetLines(lines)
        profileData.GetPointData().SetScalars(scalars)

        if scalars.GetNumberOfTuples() > 0:
            profileTubes = vtk.vtkRibbonFilter()
            #profileTubes.SetNumberOfSides(polymer.vtk_arg_list['trace']['sides'])
            profileTubes.SetInput(profileData)
            #profileTubes.SetRadius(polymer.vtk_arg_list['trace']['width'])
            
            # not sure why this was here, but it causes crashes now.
            #normals = vtk.vtkPolyDataNormals()
            #normals.SetInput(profileTubes.GetOutput())
            #normals.FlipNormalsOn()

            luTable = self._build_tubes_lookup_table(polymer)
        
            profileMapper = vtk.vtkPolyDataMapper()
            profileMapper.SetLookupTable(luTable)
            profileMapper.SetInput(profileTubes.GetOutput())
            #profileMapper.ScalarVisibilityOn()
            #profileMapper.SetScalarModeToUsePointData()
            #profileMapper.SetColorModeToMapScalars()
            profileTubes.Update()
            profileMapper.SetScalarRange(0,len(polymer.residues))
        else:
            profileMapper = vtk.vtkPolyDataMapper()
            profileMapper.SetInput(profileData)
        
        actor.SetMapper(profileMapper)
        return actor


    def _draw_tube_for_polymer(self, polymer, actor):
        print 'drawing tube for polymer'
        try:
            self.polymer_splines[polymer.chain_name]
        except:
            self._calculate_splines(polymer)

        splines,lines,scalars = self.get_polymer_splines_lines_scalars(polymer)
        profileData = vtk.vtkPolyData()
        profileData.SetPoints(splines)
        profileData.SetLines(lines)
        profileData.GetPointData().SetScalars(scalars)

        if scalars.GetNumberOfTuples() > 0:
            profileTubes = vtk.vtkTubeFilter()
            profileTubes.SetNumberOfSides(polymer.vtk_arg_list['trace']['sides'])
            profileTubes.SetInput(profileData)
            profileTubes.SetRadius(polymer.vtk_arg_list['trace']['width'])
            
            # don't know why I had this here, but it causes crashes now for some reason
            #normals = vtk.vtkPolyDataNormals()
            #normals.SetInput(profileTubes.GetOutput())
            #normals.FlipNormalsOn()
            #normals.Update()

            luTable = self._build_tubes_lookup_table(polymer)
            
            profileMapper = vtk.vtkPolyDataMapper()
            profileMapper.SetLookupTable(luTable)
            profileMapper.SetInput(profileTubes.GetOutput())
            profileMapper.ScalarVisibilityOn()
            profileMapper.SetScalarModeToUsePointData()
            profileMapper.SetColorModeToMapScalars()

            profileTubes.Update()
            profileMapper.SetScalarRange(0,len(polymer.residues))

        else:
            profileMapper = vtk.vtkPolyDataMapper()
            profileMapper.SetInput(profileData)

        actor.SetMapper(profileMapper)
        return actor

    def get_polymer_splines_lines_scalars(self, polymer):
        # if you're looking for the bug that leaves gaps between residues in trace mode,
        # set the number of splines per aminoacid to an even number (currently in GraphicsVisitor.load_system)
        splines_per = polymer.vtk_arg_list['trace']['splines']
        points  = vtk.vtkPoints()
        lines   = vtk.vtkCellArray()
        scalars = vtk.vtkDoubleArray()
        scalars.SetNumberOfComponents(1)
        number_of_central_points = len(polymer.residues)
        number_of_spline_points  = (number_of_central_points-1.0) * splines_per
        frac = (number_of_central_points-1.0)/(number_of_spline_points-1.0)
        index = 0
        ptrcollection = []            
        for res in polymer.residues:
            spline_start  = (splines_per*index)-int(round(splines_per/2))
            spline_middle = (splines_per*index)
            spline_end    = (splines_per*index)+int(round(splines_per/2))
            if res.is_Nterm == 1:
                start = spline_middle
                end   = spline_end
            elif res.is_Cterm == 1:
                start = spline_start
                end   = spline_middle
            else:
                start = spline_start
                end   = spline_end
            for ptcntr in range(start, end):
                t = frac*ptcntr
                points.InsertPoint(ptcntr, self.polymer_splines[polymer.chain_name]['x'].Evaluate(t),
                                           self.polymer_splines[polymer.chain_name]['y'].Evaluate(t),
                                           self.polymer_splines[polymer.chain_name]['z'].Evaluate(t))
                ptrcollection.append(ptcntr)
                scalars.InsertTuple1(ptcntr, index)
            index = index + 1
            
        start = 0
        end = 0
        for i in range(len(ptrcollection)-1):
            if ptrcollection[i] == ptrcollection[i+1]-1:
                end = end + 1
            else:
                lines.InsertNextCell(end-start)
                for j in range(start, end):
                    lines.InsertCellPoint(ptrcollection[j])
                start = end+1
        lines.InsertNextCell(end-start)
        for j in range(start, end):
            lines.InsertCellPoint(ptrcollection[j])

        return points,lines,scalars

    def _calculate_splines(self, polymer):
        central_point_list  = polymer.get_central_point_list()
        self.polymer_splines[polymer.chain_name] = {}
        self.polymer_splines[polymer.chain_name]['x'] = vtk.vtkCardinalSpline()
        self.polymer_splines[polymer.chain_name]['y'] = vtk.vtkCardinalSpline()
        self.polymer_splines[polymer.chain_name]['z'] = vtk.vtkCardinalSpline()
        id_cntr = 0
        for rez in polymer.residues:
            if rez.has_central_pt==1:
                pt = rez.central_atom
                xpt = pt.x#get_coordinate_list()[0]
                ypt = pt.y#get_coordinate_list()[1]
                zpt = pt.z#get_coordinate_list()[2]
                self.polymer_splines[polymer.chain_name]['x'].AddPoint(id_cntr, xpt)
                self.polymer_splines[polymer.chain_name]['y'].AddPoint(id_cntr, ypt)
                self.polymer_splines[polymer.chain_name]['z'].AddPoint(id_cntr, zpt)
                id_cntr = id_cntr + 1
            else:
                print 'missing central point for %s'%(rez.res_number)
        
    """  surface rendering  """
    def _draw_volume_for_molecule(self, mol, actor, forced_rewrite=0):
        """ will accept polymers, ligands, small molecules, but not residues """
        self.solvent_radius = 1.4
        self.grid_spacing = 4*self.solvent_radius
        if self.hydrogens_on:
            filename = mol.parent.get_filename_by_extension('.hms', mol.key)   # spade molecular surface
        else:
            filename = mol.parent.get_filename_by_extension('.sms', mol.key)   # spade molecular surface
        create_new = 0
        if forced_rewrite:
            create_new = 1
        else:
            try:
                bms_file = open(filename)
            except IOError:
                create_new = 1
        if create_new:
            return self._chomp_futamura_surface(mol, actor, 250, 1.4, 50, 1)
        else:
            bms_file.close()
            reader = vtk.vtkPolyDataReader()
            reader.SetFileName(filename)

            rev = vtk.vtkReverseSense()
            rev.SetInput(reader.GetOutput())
            rev.ReverseCellsOff()
            
            normals = vtk.vtkPolyDataNormals()
            normals.SetInput(rev.GetOutput())
            normals.FlipNormalsOn()
            
            luTable = self._build_surface_lookup_table(mol)
            
            map = vtk.vtkPolyDataMapper()
            map.SetLookupTable(luTable)
            map.ScalarVisibilityOn()
            map.SetScalarModeToUsePointData()
            map.SetColorModeToMapScalars()
            map.SetInput(normals.GetOutput())
            map.SetScalarRange(0,len(mol.atoms))
            map.ImmediateModeRenderingOn()

            actor.SetMapper(map)
            self.renderer.AddActor(actor)
            return actor

    def _get_surface_color_scalars(self, mol, solvent_radius, surface_points, smooth_input):
        """ return a list of scalars the same size as surface points, where colors are
        mapped from atoms to surface_points. returns rgb [0-1] triples  """
        grid = FutamuraHash(mol)
        T = grid.T
        radii  = {'C':1.75,
                  'O':1.4,
                  'N':1.55,
                  'S':1.8,
                  'P':2.0,
                  'H':1.17,
                  'Z':3.0}
        default_distance = 1.8
        print 'locating nearest atoms'
        scalars = vtk.vtkIntArray()
        scalars.SetNumberOfComponents(1)
        # now locate the intersections
        number_index_map = {}
        for ind in range(len(mol.atoms)):
            number_index_map[mol.atoms[ind].atom_number] = ind
            
        last_atom = 'None'
        if smooth_input:
            new_points = []
        ptctr = 0
        for point in surface_points:
            x_val = y_val = z_val = 0
            # figure out which bin it goes in
            for x_ind in range(0, grid.volume_count_x):
                if point[0] < grid.volume_indices_x[x_ind]:
                    break
                else:
                    x_val = x_ind
            for y_ind in range(grid.volume_count_y):
                if point[1] < grid.volume_indices_y[y_ind]:
                    break
                else:
                    y_val = y_ind
            for z_ind in range(grid.volume_count_z):
                if point[2] < grid.volume_indices_z[z_ind]:
                    break
                else:
                    z_val = z_ind

            start_array = [0,0,0]
            end_array   = [0,0,0]
            # figure out starts and ends
            counts = [grid.volume_count_x, grid.volume_count_y, grid.volume_count_z]
            keys   = [x_val, y_val, z_val]
            for ind in [0,1,2]:
                if keys[ind] == 0:
                    start_array[ind] = 0
                    end_array[ind]   = 2
                elif keys[ind] == counts[ind] - 1:
                    start_array[ind] = keys[ind]-1
                    end_array[ind]   = keys[ind]+1
                else:
                    start_array[ind] = keys[ind]-1
                    end_array[ind]   = keys[ind]+2
            min_dist = 1000.0
            sec_dist = 1000.0
            id2 = -1
            id = -1
            escape = 0          # turns 1 once the correct atom is found
            if smooth_input == 0:
                identification_distance = 0.1
                # figure out if its in range of the last atom chosen (arbitrary, but tends to speed up the calculations)
                if last_atom != 'None':
                    dist = math.sqrt(pow(point[0]-last_atom.x,2) + pow(point[1]-last_atom.y,2) + pow(point[2]-last_atom.z,2))
                    dif = abs(dist - radii.get(last_atom.atom_type[0], default_distance))
                    if dif < identification_distance:
                        id = last_atom.atom_number             # assume this is it
                        escape = 1
                        
                if not escape:
                    # now look for atoms in the same bin as the last atom
                    ky = '%s %s %s'%(x_val,y_val,z_val)
                    if ky in T.keys():            # first look in this atoms bin
                        for atom in T[ky]:
                            # do not retrieve if type H and protonation is turned off
                            if self.hydrogens_on or ((not self.hydrogens_on) and atom.atom_type[0] != 'H'):
                                dist = math.sqrt(pow(point[0]-atom.x,2) + pow(point[1]-atom.y,2) + pow(point[2]-atom.z,2))
                                if abs(dist - radii.get(atom.atom_type[0], default_distance)) < identification_distance:
                                    id = atom.atom_number             # assume this is it
                                    escape = 1
                                    break
            if not escape:
                for i in range(start_array[0], end_array[0]):
                    for j in range(start_array[1], end_array[1]):
                        for k in range(start_array[2], end_array[2]):
                            key2 = '%s %s %s'%(i,j,k)
                            #if key2 != ky:
                            if key2 in T.keys():
                                for atom in T[key2]:
                                    if self.hydrogens_on or ((not self.hydrogens_on) and atom.atom_type[0] != 'H'):
                                        dist = math.sqrt(pow(point[0]-atom.x,2) + pow(point[1]-atom.y,2) + pow(point[2]-atom.z,2))
                                        if not smooth_input:
                                            if abs(dist - radii.get(atom.atom_type[0], default_distance)) < identification_distance:
                                                id = atom.atom_number
                                                escape = 1
                                                break
                                            elif dist < min_dist:
                                                min_dist = dist
                                                id = atom.atom_number
                                        else:
                                            if dist < min_dist:
                                                sec_dist = min_dist
                                                id2 = id
                                                min_dist = dist
                                                id = atom.atom_number
                            if escape:
                                break
                        if escape:
                            break
                    if escape:
                        break
            # assign the index
            last_atom = mol.atoms[number_index_map[id]]
            scalars.InsertTuple1(ptctr, number_index_map[id])
            # smooth the data
            fitting_back_distance   = 0.2
            if smooth_input:
                x2 = point[0]
                y2 = point[1]
                z2 = point[2]
                if id2 != -1:               # more than one intersection is necessary
                    sec_last_atom = mol.atoms[number_index_map[id2]]
                    if abs(min_dist-radii.get(last_atom.atom_type[0], default_distance)) < fitting_back_distance:   # if this atom is close enough
                        if abs(sec_dist-radii.get(sec_last_atom.atom_type[0], default_distance)) > 0.4:   # if second atom is far enough away
                            r = radii.get(last_atom.atom_type[0], default_distance)
                            d = min_dist
                            x = last_atom.x
                            y = last_atom.y
                            z = last_atom.z
                            x2 = ((r/d)*(point[0]-x)) + x
                            y2 = ((r/d)*(point[1]-y)) + y
                            z2 = ((r/d)*(point[2]-z)) + z
                new_points.append([x2,y2,z2])
                
            ptctr += 1
        if smooth_input:
            return scalars,new_points
        else:
            return scalars
        
    """  surface calculation  """
    def _get_futamura_surface_points(self, mol, sphere_res, solvent_radius):
        if mol.x_table == None:
            mol.build_futamura_intersection_table(solvent_radius)
        x_table = mol.x_table
        radii  = {'C':solvent_radius + 1.75,
                  'O':solvent_radius + 1.4,
                  'N':solvent_radius + 1.55,
                  'S':solvent_radius + 1.8,
                  'P':solvent_radius + 2.0,
                  'H':solvent_radius + 1.17,
                  'Z':solvent_radius + 3.0}
        default_distance = solvent_radius + 1.8
        # create spheres for each atom
        print 'collecting surface points'
        total_points = 0
        surface_points = []
        point_count = sphere_res
        speedup_hits = 0
        for atom in mol.atoms:
            if self.hydrogens_on or ((not self.hydrogens_on) and atom.atom_type != 'H'):
                sphere = []
                radius = radii.get(atom.atom_type[0], default_distance)
                radius_sq = radius**2;
                last_second_atom = 'None'
                for i in range(point_count):
                    angle = random.random() * 2 * 3.141592654
                    z = (random.random() * 2 * radius) - radius
                    z_sq = z**2;
                    x_store = math.sqrt(radius_sq - z_sq) * math.cos(angle) + atom.x
                    y_store = math.sqrt(radius_sq - z_sq) * math.sin(angle) + atom.y
                    z_store = z + atom.z
                    if last_second_atom != 'None':
                        if math.sqrt(pow(x_store-last_second_atom[0],2) + pow(y_store-last_second_atom[1],2) + pow(z_store-last_second_atom[2],2)) <= last_second_atom[3]:
                            speedup_hits += 1
                            continue
                    for second_atom in x_table['%s'%(atom.atom_number)]:
                        if second_atom != last_second_atom:
                            if self.hydrogens_on or ((not self.hydrogens_on) and second_atom[6] != 'H'):
                                if math.sqrt(pow(x_store-second_atom[0],2) + pow(y_store-second_atom[1],2) + pow(z_store-second_atom[2],2)) <= second_atom[3]:
                                    last_second_atom = second_atom
                                    break
                    else:
                        surface_points.append([x_store,y_store,z_store])
                    
        print '%s water centers collected (of %s - %s speedup hits)'%(len(surface_points), len(mol.atoms)*point_count, speedup_hits)
        return surface_points
                
    def _chomp_futamura_surface(self, mol, volume_actor, sphere_res, solvent_radius, surf_res, back_fit_points = 1):
        radii  = {'C':1.75,
                  'O':1.4,
                  'N':1.55,
                  'S':1.8,
                  'P':2.0,
                  'H':1.17,
                  'Z':3.0}
        default_distance = 1.8

        surface_points = self._get_futamura_surface_points(mol, sphere_res, solvent_radius)
        solvent_radius = solvent_radius
        bounds = self.system.get_bounds()
        # expand the bounds a little to avoid boundary conditions
        for i in [0,2,4]:
            bounds[i] = bounds[i] - 2.0
            bounds[i+1] = bounds[i+1] + 2.0
            
        theCream = vtk.vtkImplicitBoolean()
        included_atom_count = 0
        for atom in mol.atoms:
            if self.hydrogens_on or ((not self.hydrogens_on) and atom.atom_type[0] != 'H'):
                included_atom_count += 1
                iceCream = vtk.vtkSphere()
                iceCream.SetCenter(atom.x,atom.y,atom.z)
                iceCream.SetRadius((solvent_radius/2.0)+radii.get(atom.atom_type[0], default_distance))
                #iceCream.SetRadius(2.0)         # the extra should get chomped off
                theCream.AddFunction(iceCream)
            
        # parse the polydata to create the next implicit function
        theChomps = vtk.vtkImplicitBoolean()
        for point in surface_points:
            bite = vtk.vtkSphere()
            bite.SetCenter(point)
            bite.SetRadius(solvent_radius)
            theChomps.AddFunction(bite)

        theSurfaceEqn = vtk.vtkImplicitBoolean()
        theSurfaceEqn.SetOperationTypeToDifference()
        theSurfaceEqn.AddFunction(theCream)
        theSurfaceEqn.AddFunction(theChomps)
        print 'sampling the difference'

        theSurfaceMap = vtk.vtkSampleFunction()
        theSurfaceMap.SetImplicitFunction(theSurfaceEqn)
        theSurfaceMap.SetModelBounds(bounds)
        theSurfaceMap.SetSampleDimensions(surf_res,surf_res,surf_res)
        theSurfaceMap.ComputeNormalsOff()
        theSurfaceMap.Update()

        theSurface = vtk.vtkContourFilter()
        theSurface.SetInput(theSurfaceMap.GetOutput())
        theSurface.SetValue(0,0.0)

        smoother = vtk.vtkWindowedSincPolyDataFilter()
        smoother.SetInput(theSurface.GetOutput())
        smoother.SetNumberOfIterations(10)

        """
        decimate = vtk.vtkDecimatePro()
        decimate.SetInput(theSurface.GetOutput())
        decimate.PreserveTopologyOn()
        decimate.SetTargetReduction(.2)

        """
        x = smoother.GetOutput()
        x.Update()

        final_point_count = x.GetNumberOfPoints()
        final_points = []
        final_copy   = []
        for i in range(final_point_count):
            final_points.append(x.GetPoint(i))
            final_copy.append(x.GetPoint(i))
            
        print '%s final points'%(len(final_points))
        
        luTable = self._build_surface_lookup_table(mol)
        # forces the point back to a precise distance from the closest atom (or sometimes
        # just any atom in range
        if back_fit_points == 0:
            scalars = self._get_surface_color_scalars(mol, 0.8, final_points, 0)
        else:
            scalars,new_points = self._get_surface_color_scalars(mol, 0.8, final_points, 1)
            print '%s final points, %s new_points'%(len(final_points), len(new_points))

            new_vtk_points = vtk.vtkPoints()
            for point_ind in range(len(new_points)):
                new_vtk_points.InsertPoint(point_ind, new_points[point_ind][0], new_points[point_ind][1], new_points[point_ind][2])
            
            x.SetPoints(new_vtk_points)
            
        x.GetPointData().SetScalars(scalars)

        writer = vtk.vtkPolyDataWriter()
        writer.SetInput(x)
        if self.hydrogens_on:
            writer.SetFileName(mol.parent.get_filename_by_extension('.hms', mol.key))   # spade molecular surface
        else:
            writer.SetFileName(mol.parent.get_filename_by_extension('.sms', mol.key))   # spade molecular surface
        writer.SetFileTypeToBinary()
        writer.Write()

        rev = vtk.vtkReverseSense()
        rev.SetInput(x)
        rev.ReverseCellsOff()
        
        normals = vtk.vtkPolyDataNormals()
        normals.SetInput(rev.GetOutput())
        normals.FlipNormalsOn()
        
        map = vtk.vtkPolyDataMapper()
        map.SetLookupTable(luTable)
        map.ScalarVisibilityOn()
        map.SetScalarModeToUsePointData()
        map.SetColorModeToMapScalars()
        map.SetInput(normals.GetOutput())
        map.SetScalarRange(0,len(mol.atoms))
        map.ImmediateModeRenderingOn()

        volume_actor.GetProperty().SetSpecular(mol.vtk_arg_list['volume']['specular'])
        volume_actor.GetProperty().SetSpecularPower(mol.vtk_arg_list['volume']['specular_power'])
        
        volume_actor.SetMapper(map)

        return volume_actor

if __name__ == '__main__':
    from tkFileDialog import *
    def reload_viewer(viewer, type):
        if type == 'pdb':
            new_name = askopenfilename(title = 'Select the PDB', defaultextension='.pdb', filetypes=[("Protein Data Bank", "*.pdb"),("all files", "*")])
            if len(new_name) > 0:
                new_system = MolecularSystem.System(new_name)
                viewer.loadSystem(new_system)
        elif type == 'system':
            new_name = askopenfilename(initialdir=os.path.normpath('./Databases/EC_1.1_done/'),title = 'Select the system', defaultextension='.sps', filetypes=[("SPADE Pickled System", "*.sps"),("all files", "*")])
            if len(new_name) > 0:
                new_system = MolecularSystem.System(new_name)
                viewer.loadSystem(new_system)

    def save_current_system(viewer):
        new_name = asksaveasfilename(title = 'Save as', defaultextension='.sps', filetypes=[("SPADE Pickled System", "*.sps"),("all files", "*")])
        if len(new_name) > 0:
            viewer.system.save_system(new_name)

    def print_system_info(viewer):
        for line in viewer.system.HeaderLines:
            if line[0:5] in ['HEADE', 'TITLE', 'COMPN', 'SOURC', 'KEYWD']:
                print line

    window = Tk()
    viewer = MolecularViewer(window, 'None', 400, 500, 0)
    
    if len(argv) == 2:
        new_name = argv[1]
        if os.path.exists(new_name):
            if argv[1][-3:] == 'pdb':
                new_system = MolecularSystem.System(new_name)
                viewer.loadSystem(new_system)
            elif argv[1][-3:] == 'sps':
                new_system = MolecularSystem.System(new_name)
                viewer.closeSystem()
                viewer.loadSystem(new_system)
        else:
            print 'file %s does not exist'%(new_name)
    
    viewer.menuBar.addmenu('File', 'Load/Unload systems')
    c_lambda = lambda viewer=viewer: reload_viewer(viewer, 'pdb')
    viewer.menuBar.addmenuitem('File', 'command', label='Open PDB', command = c_lambda)
    c_lambda = lambda viewer=viewer: reload_viewer(viewer, 'system')
    viewer.menuBar.addmenuitem('File', 'command', label='Open System', command = c_lambda)
    c_lambda = lambda viewer=viewer: save_current_system(viewer)
    viewer.menuBar.addmenuitem('File', 'command', label='Save System', command = c_lambda)
    c_lambda = lambda viewer=viewer: print_system_info(viewer)
    viewer.menuBar.addmenuitem('File', 'command', label='Print System Info', command = c_lambda)

    viewer._build_menu()
    viewer.has_menu = 1
    window.mainloop()

