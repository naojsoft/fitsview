#
# MoircsAlignPlugin.py -- Ruler plugin for Ginga reference viewer
#
from __future__ import absolute_import
from __future__ import print_function
import os
import threading
import six

import sys
import logging
import time
from time import strftime
import math
import sep
import cv2
import copy
import gc
import numpy as np
import time 
from numpy import ma


import astropy.io.fits as fits
from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy import signal

from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.interpolation import shift 

from ginga.misc.Bunch import Bunch
from ginga import GingaPlugin
from ginga.gw import Widgets
from ginga.gw import Viewers
from ginga.gw import Plot
from ginga.gw import GwHelp
from ginga import AstroImage
from ginga.util import plots


# ginga imports
from ginga.misc.Callback import CallbackError




DIR_MCSRED = '/home/gen2/Procedure/MOIRCS/MCSRED2/'
DIR_PAR_VAR = os.environ['HOME']
PAR_FILENAME = 'mesoffset_parameters.txt'
VAR_FILENAME = 'mesoffset_directories.txt'

#DIR_MCSRED = os.path.join(os.environ['HOME'], '/MOIRCS/MCSRED2/')
NO_SUCH_FILE_ERR = ("No such file or directory: {}\nPlease check your frame "+
                        "numbers and image directory, or run Ginga from a "+
                        "different directory.")
WRONG_CHIP_ERR =   ("{} should be data from chip {}, but is from chip {}. Try "+
                        "a different frame number.")
LOW_ELEV_WARN =   (u"{}MCSA{:08d}.fits has low elevation of {:.1f}\u00B0; the "+
                        "mosaicing database may not be applicable here.")
USER_INTERRUPT_ERR = ("This process was terminated. Please press 'Return to "+
                        "Menu' to start it over.")

SAVE_INTERMEDIATE_FILES = True


class MoircsAlignWindow(GingaPlugin.LocalPlugin):
    """
    Any custom LocalPlugin for ginga that is intended for use as part of the
    MOS Acquisition software for aligning MOIRCS.
    """
    
    HEADER_FONT = GwHelp.get_font('sansFont', 18)
    NORMAL_FONT = GwHelp.get_font('sansFont', 12)
    BODY_FONT   = GwHelp.get_font('sansFont', 10)
    MONO_FONT   = GwHelp.get_font('Monospace', 12)

    def __init__(self, fv, fitsimage):
        """
        Class constructor
        @param fv:
            A reference to the ginga.main.GingaShell object (reference viewer)
        @param fitsimage:
            A reference to the specific ginga.qtw.ImageViewCanvas object
            associated with the channel on which the plugin is being invoked
        """
        # superclass constructor defines self.fv, self.fitsimage, and self.logger:
        super(MoircsAlignWindow, self).__init__(fv, fitsimage)
        
        # now sets up the ginga.canvas.types.layer.DrawingCanvas self.canvas,
        # which is necessary to draw on the image:
        self.dc = fv.get_draw_classes()
        self.canvas = self.dc.DrawingCanvas()
        self.canvas.enable_draw(False)
        self.canvas.set_surface(self.fitsimage)
        self.canvas.register_for_cursor_drawing(self.fitsimage)
        self.canvas.name = type(self).__name__+'-canvas'
        
        
        
    def clear_canvas(self, keep_objects=False,
                           keep_callbacks=False,
                           keep_zoom=False):
        """
        Reset the ImageViewCanvas by deleting objects and callbacks
        @param keep_objects:
            If True, canvas objects will not be deleted
        @param keep_callbacks:
            If True, canvas callbacks will not be cleared
        @param keep_zoom:
            If True, fitsimage zoom level and position will not be reset
        """
        if not keep_objects:
            self.canvas.delete_all_objects()
        if not keep_callbacks:
            for button in ('cursor', 'panset', 'draw'):
                for event in ('-up', '-down'):
                    self.canvas.clear_callback(button+event)
        if not keep_zoom:
            self.fitsimage.zoom_fit()
            self.fitsimage.center_image()
    
    
    def build_gui(self, container, future=None):
        """
        Called when the plugin is invoked; sets up all the components of the GUI
        One of the required LocalPlugin methods
        @param container:
            The widget.Box this GUI must be added into
        """
        # create the outer Box that will hold the GUI and the close button
        out = Widgets.VBox()
        out.set_border_width(4)
        container.add_widget(out, stretch=True)
        
        # create the inner box that will contain the stack of GUIs
        box, box_wrapper, orientation = Widgets.get_oriented_box(container,
                                                                 fill=True)
        box.set_border_width(4)
        box.set_spacing(3)
        out.add_widget(box_wrapper, stretch=True)
        
        # the rest is a stack of GUIs for each step, as decided by the subclass
        stk = Widgets.StackWidget()
        self.stack_guis(stk, orientation)
        box.add_widget(stk, stretch=True)
        self.stack = stk
        
        # end is an HBox that comes at the very end, after the rest of the GUIs
        end = Widgets.HBox()
        end.set_spacing(2)
        out.add_widget(end)
        
        # throw in a close button at the very end, just in case
        btn = Widgets.Button("Close")
        btn.add_callback('activated', lambda w: self.close())
        end.add_widget(btn)
        end.add_widget(Widgets.Label(''), stretch=True)
    
    
    def start(self, future=None):
        """
        Called when the plugin is first invoked, right after build_gui()
        One of the required LocalPlugin methods
        """
        # stick our own canvas on top of the fitsimage canvas
        p_canvas = self.fitsimage.get_canvas()
        if not p_canvas.has_object(self.canvas):
            p_canvas.add(self.canvas, tag='MOSA-canvas')
        
        # clear the canvas
        self.canvas.delete_all_objects()
    
        # Save a reference to the "future" object so we MESInterface
        # can use it later on.
        self.callerInfo = future

    def close(self):
        """
        Called when the plugin is closed
        One of the required LocalPlugin methods
        @returns:
            True. I'm not sure why.
        """
        self.fv.stop_local_plugin(self.chname, str(self))
        return True
    
    
    def pause(self):
        """
        Called when the plugin is unfocused
        One of the required LocalPlugin methods
        """
        self.canvas.ui_setActive(False)
    
    
    def resume(self):
        """
        Called when the plugin is refocused
        One of the required LocalPlugin methods
        """
        self.canvas.ui_setActive(True)
        self.fv.showStatus("Calculate the offset values to align MOIRCS")
    
    
    def stop(self):
        """
        Called when the plugin is stopped
        One of the required LocalPlugin methods
        """
        p_canvas = self.fitsimage.get_canvas()
        try:
            p_canvas.delete_object_by_tag('MOSA-canvas')
        except:
            pass
        self.canvas.ui_setActive(False)
        self.clear_canvas()
        
        # call the termination event, if you have one
        if hasattr(self, 'terminate'):
            self.terminate.set()
    
    
    def redo(self):
        """
        Called whenever a new image is loaded
        One of the required LocalPlugin methods
        """
        pass

    def __str__(self):
        return type(self).__name__


class MoircsAlign(MoircsAlignWindow):
    """
    A custom LocalPlugin for ginga that takes parameters from the user in a
    user-friendly menu, locates a set of calibration objects, asks for users to
    help locate anomolies and artifacts on its images of those objects,
    calculates their centers of masses, graphs some data about some objects'
    positions, and asks the user to modify the data if necessary. Intended for
    use as part of the MOS Acquisition software for aligning MOIRCS.
    """
    SELECTION_MODES = ("Automatic", "Crop", "Mask")
    BOX_COLORS = ('green','red','blue','yellow','magenta','cyan','orange')
    VALUE_NAMES = (("dX","pix"), ("dY","pix"), ("dPA",u"\u00B0"))
    OUTLIER_VALUE = 2.0
    # main menu parameters
    PARAMS_1 = [
        {'name':'star_chip1',
         'label':"Star Frame", 'type':int, 'format':"MCSA{}.fits", 'default':0,
         'desc':"The frame number for the chip1 star FITS image"},
        
        {'name':'sky_chip1',
         'label':"Sky Frame", 'type':int, 'format':"MCSA{}.fits", 'default':0,
         'desc':"The frame number for the chip1 sky FITS image"},

        {'name':'mask_chip1',
         'label':"Mask Frame", 'type':int, 'format':"MCSA{}.fits", 'default':0,
         'desc':"The frame number for the chip1 mask FITS image"},    

        {'name':'rootname',
         'label':"Root Name", 'type':str, 'format':"{}.sbr", 'default':'',
         'desc':"The filename of the mask definition SBR file"},
        
        {'name':'c_file',
         'label':"Config File", 'type':str, 'default':"$DATABASE/ana_apr16.cfg",
         'desc':"The location of the MCSRED configuration file"},
        
        {'name':'img_dir',
         'label':"Image Directory", 'type':str, 'default':"$DATA/",
         'desc':"The directory in which the input FITS images can be found"},
        
        {'name':'exec_mode',
         'label':"Execution Mode", 'type':int, 'options':["Normal","Fine"],
         'desc':"Choose 'Fine' to skip MES Offset 1"},

        {'name':'recalc1',
         'label':"Regenerate Star", 'type':bool,
         'desc':"Do you want to generate new composite star images?"},
        
        {'name':'interact1',
         'label':"Interact Star", 'type':bool,
         'desc':"Do you want to interact with star position measurement?"},
   
        {'name':'debug_mode',
         'label':"Debug", 'type':bool, 'default': False,
         'desc':"Check this box for keeping all by-products"}
    ]

    PARAMS_2 = [
        {'name':'starhole_chip1',
         'label':"Star-Hole Frame", 'type':int, 'format':"MCSA{}.fits",
         'desc':"The frame number for the chip1 star-hole FITS image"},
        
        {'name':'mask_chip1',
         'label':"Mask Frame", 'type':int, 'format':"MCSA{}.fits",
         'desc':"The frame number for the chip1 mask FITS image"},
        
        {'name':'rootname',
         'label':"Root Name", 'type':str, 'format':"{}.sbr",
         'desc':"The filename of the mask definition SBR file"},
        
        {'name':'c_file',
         'label':"Config File", 'type':str, 'default':"$DATABASE/ana_apr16.cfg",
         'desc':"The location of the MCSRED configuration file"},
        
        {'name':'img_dir',
         'label':"Image Directory", 'type':str, 'default':"$DATA/",
         'desc':"The directory in which the raw FITS images can be found"},
        
        {'name':'recalc2',
         'label':"Regenerate", 'type':bool,
         'desc':"Do you want to generate new composite star-hole images?"},
        
        {'name':'interact2',
         'label':"Interact", 'type':bool,
         'desc':"Do you want to interact with star position measurement?"}
    ]
    
    PARAMS_3 = [
        {'name':'starhole_chip1',
         'label':"Star-Hole Frame", 'type':int, 'format':"MCSA{}.fits",
         'desc':"The frame number for the chip1 star-hole FITS image"},
        
        {'name':'sky_chip1',
         'label':"Sky Frame ", 'type':int, 'format':"MCSA{}.fits",
         'desc':"The frame number for the chip1 mask FITS image"},

        {'name':'newmask_chip1',
         'label':"New Mask Frame", 'type':int, 'format':"MCSA{}.fits",
         'desc':"The frame number for the chip1 star-hole FITS image"},
        
        {'name':'rootname',
         'label':"Root Name", 'type':str, 'format':"{}.sbr",
         'desc':"The filename of the mask definition SBR file"},
        
        {'name':'c_file',
         'label':"Config File", 'type':str, 'default':"$DATABASE/ana_apr16.cfg",
         'desc':"The location of the MCSRED configuration file"},
        
        {'name':'img_dir',
         'label':"Image Directory", 'type':str, 'default':"$DATA/",
         'desc':"The directory in which the raw FITS images can be found"},
        
        {'name':'recalc3',
         'label':"Regenerate", 'type':bool,
         'desc':"Do you want to generate new composite star-hole images?"},
        
        {'name':'interact3',
         'label':"Interact", 'type':bool,
         'desc':"Do you want to interact with star position measurement?"}
    ]

    # Establish a null location for MoircsImage class.
    moircsAlignImage = 0

    def __init__(self, fv, fitsimage):
        super(MoircsAlign, self).__init__(fv, fitsimage)
        
        print(self.NORMAL_FONT)
        self.DIR_MCSRED = DIR_MCSRED

        self.variables = self.read_variables()

        self.database = {}   # the variables shared between the departments
        self.stack_idx = {} # the indices of the guis in the stack
        #self.parameter_tabs = 0 

        # the function to run when an image is loaded
        self.image_set_next_step = None
        self.fitsimage.add_callback('image-set', self.image_set_cb)
        
        self.training_dir = None
        self.rootname_logfile = None

        self.moircsAlignImage = MoircsAlignImage()


        self.get_value = []         # getter methods for all parameters
        self.set_value = []         # setter methods for all parameters
        self.resume_mesoffset = {}  # intermediate functions to call after waiting
        self.last_wait_gui = 0      # the last 'wait' gui we were at
        self.variables = self.read_variables()   # defined variables

    def set_params(self, star_chip1, rootname, c_file, img_dir, exec_mode,
                   mcsred_dir, training_dir, work_dir, wait_gui):

        self.PARAMS_1[0]['default'] = int(star_chip1)
        self.PARAMS_1[1]['default'] = int(star_chip1+2)
        self.PARAMS_1[2]['default'] = int(star_chip1+4)
        self.PARAMS_1[3]['default'] = rootname
        self.PARAMS_1[4]['default'] = '$' + c_file
        self.PARAMS_1[5]['default'] = img_dir
        try:
            exec_option = self.PARAMS_1[6]['options'].index(exec_mode.capitalize())
        except ValueError:
            exec_option = 0
        self.PARAMS_1[6]['default'] = exec_option
        
        ##fitsUtils.SAVE_INTERMEDIATE_FILES = self.PARAMS_0[5]['default']
        self.set_params_db(mcsred_dir)
        self.training_dir = training_dir
        
        self.DIR_MCSRED = mcsred_dir + '/'
        self.logger.info('MESOffset self.DIR_MCSRED is %s' % self.DIR_MCSRED)

        self.work_dir = os.path.join(work_dir, 'MESOffset')
        self.logger.info('MESOffset work_dir is %s' % self.work_dir)

        if os.path.isdir(self.work_dir):
            self.logger.info('MESOffset work_dir %s exists' % self.work_dir)
        else:
            self.logger.info('MESOffset creating %s' % self.work_dir)
            os.makedirs(self.work_dir)
        self.rootname_logfile = os.path.join(self.work_dir, rootname+"_log")
        self.wait_gui = wait_gui
        self.logger.info('MESOffset wait_gui input value %s self.wait_gui %s' % (wait_gui, self.wait_gui))
        self.logger.info('MESOffset after PARAMS_1 %s' % self.PARAMS_1)


    def set_params_db(self, mcsred_dir):
        self.variables['DATABASE'] = os.path.join(mcsred_dir, 'DATABASE') 


    def gui_list(self, orientation='vertical'):
        """
        Combine the GUIs necessary for the interface part of this plugin
        Must be implemented for each MESPlugin
        @param orientation:
            Either 'vertical' or 'horizontal', the orientation of this new GUI
        @returns:
            A list of tuples with strings (names) and Widgets (guis)
        """
        # the interface is unique in that it has a tabwidget of similar guis
        self.get_value = []
        self.set_value = []
        tab = Widgets.TabWidget()

        #tab.add_widget(self.make_gui_epar(0, orientation), "MESOffset 0")
        tab.add_widget(self.make_gui_epar(0, orientation), "MESOffset 1")
        tab.add_widget(self.make_gui_epar(1, orientation), "MESOffset 2")
        tab.add_widget(self.make_gui_epar(2, orientation), "MESOffset 3")
        #tab.add_widget(self.make_gui_epar(3, orientation), "MESOffset 3")
        self.logger.info('Tab %s' % tab) 
        self.parameter_tabs = tab
        self.logger.info('Tab %s' % tab) 

        return [('epar',   self.parameter_tabs),
                #('wait 1', self.make_gui_wait(1, orientation)),
                #('wait 3', self.make_gui_wait(3, orientation)),
                ('check',  self.make_gui_look(orientation)),
                ('log',    self.make_gui_log(orientation)),
                ('error',  self.make_gui_err(orientation)),
                ('find',     self.make_gui_find(orientation)),
                ('centroid', self.make_gui_cent(orientation)),
                ('plots',  self.make_gui_plot(orientation))]

        
    def make_gui_epar(self, idx, orientation='vertical'):
        """
        Construct a GUI for the parameter menu, which prepares to launch process
        @param idx:
            The index of this MESOffset process, or None for an ambiguous epar
        @param orientation:
            Either 'vertical' or 'horizontal', the orientation of this new GUI
        @returns:
            A Widgets.Box object containing all necessary buttons, labels, etc.
        """
        name = "MES Offset {}".format(idx+1)
        # start by creating the container
        gui = Widgets.Box(orientation=orientation)
        gui.set_spacing(4)
        
        print("DEFAUT FONT = %s"%self.NORMAL_FONT)
        # fill a text box with brief instructions and put in in an expander
        exp = Widgets.Expander(title="Instructions")
        gui.add_widget(exp)
        txt = Widgets.TextArea(wrap=True, editable=False)
        txt.set_font(super(MoircsAlign, self).NORMAL_FONT)
        txt.set_text("Use the widgets below to specify the parameters for "+
                     name+". Hover over each one to get a description of "+
                     "what it means. When you are finished, press the 'Go!' "+
                     "button, which will begin "+name+".")
        exp.set_widget(txt)
        
        # chose the params
        #if idx == 0:
        #    params = self.PARAMS_1
        if idx == 0:
            params = self.PARAMS_1
        elif idx == 1:
            params = self.PARAMS_2
        elif idx == 2:
            params = self.PARAMS_3
        #elif idx == 3:
        #    params = self.PARAMS_1

        # create a grid to group the different controls
        frm = Widgets.Frame(name)
        gui.add_widget(frm)
        grd, getters, setters = self.build_control_layout(params,
                                                     self.start_process_cb)
        self.get_value.append(getters)
        self.set_value.append(setters)
        frm.set_widget(grd)
        
        # create a box for the defined variables
        frm = Widgets.Frame("Defined Variables")
        gui.add_widget(frm)
        box = self.build_dict_labels(self.variables)
        frm.set_widget(box)
        
        # the go button will go in a box
        box = Widgets.HBox()
        box.set_spacing(3)
        gui.add_widget(box)
        

        # the go button is important
        btn = Widgets.Button("Go!")
        btn.add_callback('activated', self.start_process_cb)
        btn.set_tooltip("Start "+name+" with the given parameters")
        box.add_widget(btn)
        if self.wait_gui:
            btn = Widgets.Button("Done")
            btn.add_callback('activated', self.done_cb)
            btn.set_tooltip("Done with MESOffset")
            box.add_widget(btn)
        box.add_widget(Widgets.Label(''), stretch=True)
        self.logger.info('GO! Button established!')

        # space appropriately and return
        gui.add_widget(Widgets.Label(''), stretch=True)
        return gui

    def make_gui_log(self, orientation='vertical'):
        """
        Construct a GUI for the log: a simple textbox
        @param orientation:
            Either 'vertical' or 'horizontal', the orientation of this new GUI
        @returns:
            A Widget object containing all necessary buttons, labels, etc.
        """
        # there is a box with a single button
        gui = Widgets.Box(orientation=orientation)
        gui.set_spacing(4)
        
        # the 'log' is a gigantic text box
        scr = Widgets.ScrollArea()
        gui.add_widget(scr, stretch=True)
        txt = Widgets.TextArea(wrap=False, editable=False)
        txt.set_font(self.BODY_FONT)
        self.log_textarea = txt
        scr.set_widget(txt)
        
        # underneath it is a 'Stop' button
        box = Widgets.HBox()
        gui.add_widget(box)
        btn = Widgets.Button("Stop")
        btn.set_tooltip("Terminate the current process")
        btn.add_callback('activated', self.terminate_process_cb)
        box.add_widget(btn)
        box.add_widget(Widgets.Label(''), stretch=True)
        
        return gui
    
    def make_gui_find(self, orientation='vertical'):
        """
        Construct a GUI for the first step: finding the objects
        @param orientation:
            Either 'vertical' or 'horizontal', the orientation of this new GUI
        @returns:
            A Widgets.Box object containing all necessary buttons, labels, etc.
        """
        # start by creating the container
        gui = Widgets.Box(orientation=orientation)
        gui.set_spacing(4)

        # fill a text box with brief instructions and put in in an expander
        exp = Widgets.Expander(title="Instructions")
        gui.add_widget(exp)
        txt = Widgets.TextArea(wrap=True, editable=False)
        txt.set_font(super(MoircsAlign, self).NORMAL_FONT)
        txt.set_text("Left click on the object closest to the box labeled "+
                     "'1'. The other objects should appear in the boxes "+
                     "below. Click again to select another position. Click "+
                     "'Next' below or right-click when you are satisfied with "+
                     "your location.\nRemember - bright areas are shown in "+
                     "white.")
        exp.set_widget(txt)
        
        # next, we need the zoomed-in images. This is the grid we put them in
        frm = Widgets.Frame()
        gui.add_widget(frm)
        grd = Widgets.GridBox()
        grd.set_spacing(3)
        frm.set_widget(grd)
        self.viewer_grid = grd

        # create a box to group the primary control buttons together
        box = Widgets.HBox()
        box.set_spacing(3)
        gui.add_widget(box)
        
        # the undo button goes back a click
        btn = Widgets.Button("Undo")
        btn.add_callback('activated', self.undo1_cb)
        btn.set_tooltip("Undo a single click (if a click took place)")
        box.add_widget(btn)

        # the redo button goes forward
        btn = Widgets.Button("Redo")
        btn.add_callback('activated', self.redo1_cb)
        btn.set_tooltip("Undo an undo action (if an undo action took place)")
        box.add_widget(btn)

        # the clear button erases the canvas
        btn = Widgets.Button("Clear")
        btn.add_callback('activated', self.reset_cb)
        btn.set_tooltip("Erase all marks on the canvas")
        box.add_widget(btn)
        
        # the next button moves on to step 2
        btn = Widgets.Button("Next")
        btn.add_callback('activated', self.step2_cb)
        btn.set_tooltip("Accept and proceed to step 2 (right-click)")
        box.add_widget(btn)
        
        # Adding a button to skip all finding process, use system default
        #  for rotation fitting.
        btn = Widgets.Button("Accept all")
        btn.add_callback('activated', self.check_mes_star)
        btn.set_tooltip("Accept and proceed to step 2 (right-click)")
        box.add_widget(btn)
        
        # put in a spacer
        box.add_widget(Widgets.Label(''), stretch=True)
        
        # now add a section for more precise control
        frm = Widgets.Frame()
        gui.add_widget(frm)
        box = Widgets.VBox()
        box.set_spacing(3)
        frm.set_widget(box)
        
        # put in spinboxes for easy precision-alignment
        self.spinboxes = {}
        for var in ("X", "Y"):
            lbl = Widgets.Label(var+" position:")
            box.add_widget(lbl)
            row = Widgets.HBox()
            box.add_widget(row)
            
            num = Widgets.SpinBox(dtype=float)
            num.set_limits(0, 9999, 5)
            num.add_callback('value-changed', self.set_position_cb)
            num.set_tooltip("Use this to fine-tune the "+var+" value")
            row.add_widget(num, stretch=True)
            self.spinboxes[var] = num
            row.add_widget(Widgets.Label(''), stretch=True)
        
        # space appropriately and return
        gui.add_widget(Widgets.Label(''), stretch=True)
        return gui
        
        
    def make_gui_cent(self, orientation='vertical'):
        """
        Construct a GUI for the second step: calculating the centroids
        @param orientation:
            Either 'vertical' or 'horizontal', the orientation of this new GUI
        @returns:
            A Widgets.Box object containing all necessary buttons, labels, etc.
        """
        # start by creating the container
        gui = Widgets.Box(orientation=orientation)
        gui.set_spacing(4)

        # fill a text box with brief instructions and put in in an expander
        exp = Widgets.Expander(title="Instructions")
        gui.add_widget(exp)
        txt = Widgets.TextArea(wrap=True, editable=False)
        txt.set_font(super(MoircsAlign, self).NORMAL_FONT)
        txt.set_text("Help the computer find the centroid of this object. "+
                     "Click and drag to include or exclude regions; "+
                     "left-click will crop to selection and middle-click will "+
                     "mask selection, or you can specify a selection option "+
                     "below. Click 'Next' below or right-click when the "+
                     "centroid has been found.")
        exp.set_widget(txt)
        
        # create a label to display the current object index 
        lbl = Widgets.Label()
        lbl.set_font(self.HEADER_FONT)
        gui.add_widget(lbl)
        self.obj_count_label = lbl
        
        # create a CanvasView for step 2
        viewer = Viewers.CanvasView(logger=self.logger)
        viewer.set_desired_size(420, 420)
        viewer.enable_autozoom('on')
        viewer.enable_autocuts('on')
        viewer.add_callback('button-press', self.viewer_redirect_cb, 'down')
        viewer.add_callback('button-release', self.viewer_redirect_cb, 'up')
        self.step2_viewer = viewer
        
        # put it in a ViewerWidget, and put that in the gui
        frm = Widgets.Frame()
        gui.add_widget(frm)
        pic = Viewers.GingaViewerWidget(viewer=viewer)
        frm.set_widget(pic)
        
        # now make an HBox to hold the primary controls
        box = Widgets.HBox()
        box.set_spacing(3)
        gui.add_widget(box)
        
        # the undo button goes back a crop
        btn = Widgets.Button("Undo")
        btn.add_callback('activated', self.undo2_cb)
        btn.set_tooltip("Undo a single selection (if a selection took place)")
        box.add_widget(btn)

        # the redo button goes forward
        btn = Widgets.Button("Redo")
        btn.add_callback('activated', self.redo2_cb)
        btn.set_tooltip("Undo an undo action (if an undo action took place)")
        box.add_widget(btn)

        # the clear button nullifies all crops
        btn = Widgets.Button("Back")
        btn.add_callback('activated', self.prev_obj_cb)
        btn.set_tooltip("Go back to the previous object")
        box.add_widget(btn)
        
        # the next button moves on to the next object
        btn = Widgets.Button("Next")
        btn.add_callback('activated', self.next_obj_cb)
        btn.set_tooltip("Accept and proceed to the next object (right-click)")
        box.add_widget(btn)
        
        # put in a spacer
        box.add_widget(Widgets.Label(""), stretch=True)
        
        # another HBox holds the skip button, because it doesn't fit on the first line
        box = Widgets.HBox()
        box.set_spacing(3)
        gui.add_widget(box)
        
        # the skip button sets this object to NaN and moves on
        btn = Widgets.Button("Skip")
        btn.add_callback('activated', self.skip_obj_cb)
        btn.set_tooltip("Remove this object from the data set")
        box.add_widget(btn)
        # Save a reference to the Skip button because we need to be
        # able to enable/disable it depending on which object is being
        # displayed. First object cannot be skipped because code uses
        # its coordinates as reference position.
        self.skip_btn = btn
        
        # put in a spacer
        box.add_widget(Widgets.Label(""), stretch=True)
        
        # make a new box for a combobox+label combo
        frm = Widgets.Frame()
        gui.add_widget(frm)
        box = Widgets.VBox()
        box.set_spacing(3)
        frm.set_widget(box)
        lbl = Widgets.Label("Selection Mode:")
        box.add_widget(lbl)
        
        # this is the combobox of selection options
        com = Widgets.ComboBox()
        for text in self.SELECTION_MODES:
            com.append_text(text)
        com.set_index(0)
        com.add_callback('activated', self.choose_select_cb)
        com.set_tooltip("Choose what happens when you click-drag")
        box.add_widget(com)
        
        # space appropriately and return
        gui.add_widget(Widgets.Label(''), stretch=True)
        return gui

    def make_gui_plot(self, orientation='vertical'):
        """
        Construct a GUI for the third step: viewing the plots
        @param orientation:
            Either 'vertical' or 'horizontal', the orientation of this new GUI
        @returns:
            A Widgets.Box object containing all necessary buttons, labels, etc.
        """
        # start by creating the container
        gui = Widgets.Box(orientation=orientation)
        gui.set_spacing(4)

        # fill a text box with brief instructions and put in in an expander
        exp = Widgets.Expander(title="Instructions")
        gui.add_widget(exp)
        txt = Widgets.TextArea(wrap=True, editable=False)
        txt.set_font(super(MoircsAlign, self).NORMAL_FONT)
        txt.set_text("Look at the graphs. Remove any data with residuals "+
                     "greater than 1.0 or less than -1.0. Delete points by "+
                     "right clicking, and restore them by left-clicking. "+
                     "Click 'Next' below when the data is satisfactory.")
        exp.set_widget(txt)
        
        # now a framed vbox to put the plots in
        frm = Widgets.Frame()
        gui.add_widget(frm)
        box = Widgets.VBox()
        box.set_spacing(3)
        frm.set_widget(box)
        
        # add both plots to the frame
        self.plots = []
        for i in range(2):
            self.plots.append(plots.Plot(logger=self.logger))
            fig = Plot.PlotWidget(self.plots[i], width=300, height=300)
            box.add_widget(fig)
        self.plots[0].fig.canvas.mpl_connect("button_press_event",
                                             self.toggle_active_x_cb)
        self.plots[1].fig.canvas.mpl_connect("button_press_event",
                                             self.toggle_active_y_cb)
        
        # fill a text box with brief instructions and put in in an expander
        exp = Widgets.Expander(title="Instructions")
        gui.add_widget(exp)
        txt = Widgets.TextArea(wrap=True, editable=False)
        txt.set_font(super(MoircsAlign, self).NORMAL_FONT)
        txt.set_text("Enter the numbers you see below into the ANA window. dx "+
                     "and dy values are in pixels, and rotation value is in "+
                     "degrees. Values of less than 0.5 pixels and 0.01 "+
                     "degrees have been ignored. Click 'Finish' below when "+
                     "you are done.")
        exp.set_widget(txt)
        
        # now make a frame for the results
        frm = Widgets.Frame()
        gui.add_widget(frm)
        grd = Widgets.GridBox(3, 3)
        grd.set_spacing(3)
        frm.set_widget(grd)

        self.final_displays = {}
        for i, (val, unit) in enumerate(self.VALUE_NAMES):
            lbl = Widgets.Label(val+" =", halign='right')
            lbl.set_font(self.HEADER_FONT)
            grd.add_widget(lbl, i, 0)
            
            txt = Widgets.TextArea(editable=False)
            txt.set_font(self.HEADER_FONT)
            grd.add_widget(txt, i, 1, stretch=True)
            self.final_displays[val] = txt
            
            lbl = Widgets.Label(unit+"\t", halign='left')
            lbl.set_font(self.HEADER_FONT)
            grd.add_widget(lbl, i, 2)


        # now make an HBox to hold the main controls
        box = Widgets.HBox()
        box.set_spacing(1)
        gui.add_widget(box)
        
        #
        #  the next button shows the values
        #
        btn = Widgets.Button("Relocate")
        btn.add_callback('activated', self.relocate_cb)
        btn.set_tooltip("Relocation stars!")
        box.add_widget(btn)
        box.add_widget(Widgets.Label(''), stretch=True)

        btn = Widgets.Button("Finish")
        btn.add_callback('activated', self.finish_cb)
        btn.set_tooltip("Finsh the procedure!")
        box.add_widget(btn)
        box.add_widget(Widgets.Label(''), stretch=True)
        
        box.add_widget(Widgets.Label(""), stretch=True)
        # space appropriately and return
        gui.add_widget(Widgets.Label(''), stretch=True)
        return gui
    
    
 

    def stack_guis(self, stack, orientation='vertical'):
        """
        Get the GUIs from all the different departments of mesoffset, and stacks
        them together neatly in one gw.Widget
        """
        stack.remove_all()
        # all_guis = (self.mes_interface.gui_list(orientation) +
        #             self.mes_locate.gui_list(orientation) +
        #             self.mes_analyze.gui_list(orientation))
        all_guis = (self.gui_list(orientation))
        for i, (name, gui) in enumerate(all_guis):
            stack.add_widget(gui)
            self.stack_idx[name] = i


    def build_control_layout(self, controls, callback=None):
        """
        Build a grid full of labels on the left and input widgets on the right.
        @param controls:
            A list of dictionary where each dictionary has the keys 'name' (the
            name of the parameter), 'type' (str, in, etc.), 'default'
            (the starting value), 'desc' (the tooltip), possibly 'format' (puts
            labels on either side of the input), and possibly 'options' (the
            list of possible values, if type is 'combobox')
        @param callback:
            The function, if any, to be called when 'enter' is pressed
        @returns:
            A Widgets.Box containing controls for all of the layouts,
            and a dictionary whose keys are parameter names and whose values
            are functions to return those parameter values,
            and a dictionary whose keys are parameter names and whose values
            are functions to set those parameter values
        """
        grd = Widgets.GridBox(rows=len(controls), columns=4)
        grd.set_column_spacing(0)
        getters = {}
        setters = {}
        old_pars = self.read_parameters()
        
        # put each of the controls in a row on the grid
        for i, param in enumerate(controls):
            name = param['name']
            # start by labelling the parameter
            lbl = Widgets.Label(param['label']+":  ", halign='right')
            lbl.set_tooltip(param['desc'])
            grd.add_widget(lbl, i, 0)
            
            # create a widget based on type
            if 'options' in param:
                wdg = Widgets.ComboBox()
                for option in param['options']:
                    wdg.append_text(option)
                wdg.set_index(0)
                getters[name] = wdg.get_index
                setters[name] = wdg.set_index
            elif param['type'] == int:
                wdg = Widgets.SpinBox(dtype=int)
                wdg.set_limits(0, 99999999)
                wdg.set_value(0)
                getters[name] = wdg.get_value
                setters[name] = wdg.set_value
            elif param['type'] == str:
                wdg = Widgets.TextEntry(editable=True)
                wdg.set_text("")
                if callback != None:
                    wdg.add_callback('activated', callback)
                getters[name] = wdg.get_text
                setters[name] = wdg.set_text
            elif param['type'] == bool:
                wdg = Widgets.CheckBox()
                wdg.set_state(True)
                getters[name] = wdg.get_state
                setters[name] = wdg.set_state
            else:
                raise TypeError("{} is not a valid par-type.".format(param['type']))
            
            # apply the description and default
            wdg.set_tooltip(param['desc'])
            if old_pars != None and param['name'] in old_pars:
                try:
                    old_value = param['type'](old_pars[param['name']])
                    setters[param['name']](old_value)
                except ValueError:
                    pass
            elif 'default' in param:
                setters[name](param['default'])
            
            # surround the widget with text, if necessary
            if 'format' in param:
                format_str = param['format']
                idx = format_str.index('{}')
                prefix = format_str[:idx]
                suffix = format_str[idx+2:]
                if prefix:
                    grd.add_widget(Widgets.Label(prefix, 'right'), i, 1)
                grd.add_widget(wdg, i, 2)
                if suffix:
                    grd.add_widget(Widgets.Label(suffix, 'left'), i, 3)
            else:
                grd.add_widget(wdg, i, 2)
        
        return grd, getters, setters


    def build_dict_labels(self,dictionary):
        """
        Build a gui that displays the dictionary keys and values
        @param dictionary:
            A dict populated with str keys and str values
        @returns:
            A Widget object that displays the contents of the dictionary
        """
        grd = Widgets.GridBox(rows=len(dictionary), columns=3)
        grd.set_spacing(3)
        for i, key in enumerate(dictionary):
            lbl1 = Widgets.Label("$"+key+":  ",     halign='left')
            lbl2 = Widgets.Label(dictionary[key], halign='left')
            grd.add_widget(lbl1, i, 0)
            grd.add_widget(lbl2, i, 1)
            grd.add_widget(Widgets.Label(''), i, 2, stretch=True)
        return grd

    # def make_gui_wait(self, idx, orientation='vertical'):
    #     """
    #     Construct an intermediate epar GUI, as a break in the middle of a
    #     process
    #     @param orientation:
    #         Either 'vertical' or 'horizontal', the orientation of this new GUI
    #     @returns:
    #         A Widgets.Box object containing all necessary buttons, labels, etc.
    #     """
    #     name = "MES Offset {}".format(idx)
    #     # start by creating the container
    #     gui = Widgets.Box(orientation=orientation)
    #     gui.set_spacing(4)
        
    #     # fill a text box with brief instructions and put in in an expander
    #     exp = Widgets.Expander(title="Instructions")
    #     gui.add_widget(exp)
    #     txt = Widgets.TextArea(wrap=True, editable=False)
    #     txt.set_font(self.NORMAL_FONT)
    #     txt.set_text("Verify parameters in order to continue "+name+", using "+
    #                  "the widgets below. When you are finished and the "+
    #                  "specified files are ready for analysis, press the 'Go!' "+
    #                  "button, which will resume "+name+".")
    #     exp.set_widget(txt)
        
    #     # chose the params
    #     if idx == 1:
    #         params = self.PARAMS_1
    #     elif idx == 3:
    #         params = self.PARAMS_1
        
    #     # create a grid to group the different controls
    #     frm = Widgets.Frame()
    #     gui.add_widget(frm)
    #     grd, getters, setters = self.build_control_layout(params,
    #                                                  self.resume_process_cb)
    #     self.get_value.append(getters)  # NOTE that these getters and setters
    #     self.set_value.append(setters)  # will have different indices than idx
    #     frm.set_widget(grd)
        
    #     # create a box for the defined variables
    #     frm = Widgets.Frame()
    #     gui.add_widget(frm)
    #     box = self.build_dict_labels(self.variables)
    #     frm.set_widget(box)
        
    #     # the go button will go in a box
    #     box = Widgets.HBox()
    #     box.set_spacing(3)
    #     gui.add_widget(box)
        
    #     # the go button is important
    #     btn = Widgets.Button("Go!")
    #     btn.add_callback('activated', self.resume_process_cb)
    #     btn.set_tooltip("Continue "+name+" with the given parameters")
    #     box.add_widget(btn)
    #     if self.wait_gui:
    #         btn = Widgets.Button("Done")
    #         btn.add_callback('activated', self.done_cb)
    #         btn.set_tooltip("Done with MESOffset")
    #         box.add_widget(btn)
    #     box.add_widget(Widgets.Label(''), stretch=True)
        
    #     # space appropriately and return
    #     gui.add_widget(Widgets.Label(''), stretch=True)
    #     return gui

    def make_gui_look(self, orientation='vertical'):
        """
        Construct a GUI for checking results
        @param orientaiton:
            Either 'vertical' or 'horizontal', the orientation of this new GUI
        @returns
            A Widgets.Box object containing all necessary buttons, labels, etc.
        """
        # start by creating the container
        gui = Widgets.Box(orientation=orientation)
        gui.set_spacing(4)
        
        # fill a text box with brief instructions and put them in an expander
        exp = Widgets.Expander(title="Instructions")
        gui.add_widget(exp)
        txt = Widgets.TextArea(wrap=True, editable=False)
        txt.set_font(super(MoircsAlign, self).NORMAL_FONT)
        txt.set_text("Look at the results below. If they seem correct, click "+
                     "'Continue' to proceed to the next step. If they seem "+
                     "inadequate, click 'Try Again', and you will be taken "+
                     "back to the previous step to retake measurements. If "+
                     "you wish to edit the parameters for this process, click "+
                     "'Start Over' to abort and return to the main menu.")
        exp.set_widget(txt)
        
        # put in a label to ask the question:
        lbl = Widgets.Label("Are the results satisfactory?")
        gui.add_widget(lbl)
        
        # now add in the textbox for the results
        txt = Widgets.TextArea(wrap=False, editable=False)
        txt.set_font(self.MONO_FONT)
        gui.add_widget(txt, stretch=True)
        self.results_textarea = txt
        
        # now make an HBox for the controls
        box = Widgets.HBox()
        box.set_spacing(3)
        gui.add_widget(box)
        
        # the Try Again button goes to the last step
        btn = Widgets.Button("Try Again")
        btn.add_callback('activated', self.retake_measurements_cb)
        btn.set_tooltip("Go back and take these measurements again")
        box.add_widget(btn)
        
        # the Start Over button goes to the main menu
        btn = Widgets.Button("Start Over")
        btn.add_callback('activated', self.return_to_menu_cb)
        btn.set_tooltip("Return to the menu to edit your parameters")
        box.add_widget(btn)
        
        # the Continue button goes to the next step
        btn = Widgets.Button("Continue")
        btn.add_callback('activated', self.use_measurements_cb)
        btn.set_tooltip("Proceed to the next part of the process")
        box.add_widget(btn)
        
        # space the buttons
        box.add_widget(Widgets.Label(''), stretch=True)
        
        # space appropriately and return
        gui.add_widget(Widgets.Label(''), stretch=True)
        return gui


    def make_gui_err(self, orientation='vertical'):
        """
        Construct a GUI for errors: there must be a textbox and a back button
        @param orientation:
            Either 'vertical' or 'horizontal', the orientation of this new GUI
        @returns:
            A Widgets.Box object containing all the necessary stuff
        """
        # start by creating the container
        gui = Widgets.Box(orientation=orientation)
        gui.set_spacing(4)
        
        # start with a label to tell the user why they are here
        lbl = Widgets.Label("Sumimasen! :(")
        lbl.set_font(self.HEADER_FONT)
        gui.add_widget(lbl)
        
        # now for the error box itself
        txt = Widgets.TextArea(wrap=True, editable=False)
        txt.set_font(self.BODY_FONT)
        gui.add_widget(txt, stretch=True)
        self.err_textarea = txt
        
        # finish off with a box of important controls at the end
        box = Widgets.HBox()
        gui.add_widget(box)
        
        # in it is a back button
        btn = Widgets.Button("Return to Menu")
        btn.add_callback('activated', self.return_to_menu_cb)
        btn.set_tooltip("Correct the error and try again")
        box.add_widget(btn)
        box.add_widget(Widgets.Label(''), stretch=True)
        
        # space appropriately and return
        gui.add_widget(Widgets.Label(''), stretch=True)
        return gui

    def start_process_cb(self, *args):
        """
        Take the parameters from the gui and begin mesoffset{idx}
        """
        proc_num = self.parameter_tabs.get_index()
        print(self.parameter_tabs)
        print("proc = %d"%proc_num) 

        self.logger.info("Starting MES Offset {}...".format(proc_num))
        try:
            self.update_parameters(self.get_value[proc_num], False)
        except NameError as e:
            self.logger.info("NameError: "+str(e), level='e')
            return
 

        if proc_num == 0:
        #     #print(self.database)
             self.execute_mesoffset0()
        elif proc_num == 1:
             self.begin_mesoffset2()
        elif proc_num == 2:
             self.begin_mesoffset3()
        # elif proc_num == 3:
        #     self.begin_mesoffset3()

    
    def update_parameters(self, getters, write_to_file=False):
        """
        Read parameter values from getters and save them in self.manager, as
        well as in MCSRED2/mesoffset.par
        Scan all strings for variables
        @param getters:
            The dictionary of getter methods for parameter values
        @param write_to_file:
            Whether we should write these parameters to the .par file
        @raises NameError:
            If one of the values contains an undefined variable
        """
        # if DIR_PAR_VAR is a directory, write to the parameter file in it
        #self.logger.info("Updating the information for Offset ".format(proc_num))

        if write_to_file and os.path.isdir(DIR_PAR_VAR):
            par_file = open(os.path.join(DIR_PAR_VAR,PAR_FILENAME), 'w')
        else:
            par_file = None
        

        # now cycle through getters and update the file and manager
        new_params = {}
        for key, get_val in getters.items():
            value = get_val()
            if par_file != None:
                par_file.write('{},{}\n'.format(key,value))
            if type(value) in (str, six.text_type):
                value = self.process_filename(value, self.variables)
            new_params[key] = value
            
        self.logger.info('MESInterface updating database with values %s' % new_params)
        self.database.update(new_params)

    ### ----- MESOFFSET0 FUNCTION ----- ###
    def execute_mesoffset0(self):
        """ Set some stuff up to run mesoffset 1, 2, and 3 continuously """
        # start by using these values to guess at some other values
        self.__dict__.update(self.database)
        #self.database['sky_chip1'] = self.star_chip1 + 2
        #self.database['mask_chip1'] = self.star_chip1 + 4
        self.database['starhole_chip1'] = self.database['mask_chip1'] + 2
        if len(self.img_dir) <= 0 or self.img_dir[-1] != '/':
            self.database['img_dir'] += '/'
        

        # Update Settings for Image Operation
        self.moircsAlignImage.imagepath =  self.database['img_dir']   
        self.moircsAlignImage.star_fits_name = [self.star_chip1, self.star_chip1+1]
        self.moircsAlignImage.sky_fits_name =  [self.database['sky_chip1'], self.database['sky_chip1']+1]
        self.moircsAlignImage.mask_fits_name = [self.database['mask_chip1'], self.database['mask_chip1']+1]
        self.moircsAlignImage.starhole_fits_name =  [self.database['starhole_chip1'], self.database['starhole_chip1']+1]
        self.moircsAlignImage.output_path = self.work_dir            
        self.moircsAlignImage.rootname = self.database['rootname']
        self.moircsAlignImage.MAKEMOSAIC = self.database['recalc1']

        # Open configuration file to get the file name of bad pixel mask 
        # TO-DO read the distortion configuration for the parameters
        print(self.database)
        filename=self.PARAMS_1[4]['default']
        filename=filename.replace("$DATABASE",self.variables['DATABASE'])
        
        print(filename)
        cfg = open(filename, 'r')
        config = []
        line = cfg.readline()
        while line != '':
            if line[0] != '#':
                config.append(line.split()[-1].replace('dir_mcsred$',self.DIR_MCSRED))
            line = cfg.readline()
        cfg.close()
        self.moircsAlignImage.badpix_fits_name = [config[12],config[13]]
        self.moircsAlignImage.dbs_name=[config[2], config[4], config[8], config[10]]

        print(config[8],config[10])
        #set the defaults for all other menus
        for i in (0, 1, 2):
                self.set_defaults(i)
        #self.go_to_mesoffset(1)
        self.process_star_fits()
        

        #self.go_to_mesoffset(1)
        #self.process_star_fits()
         # # next step depends on exec_mode
        #if self.exec_mode == 0:
        #    self.go_to_mesoffset(2)
        #else:
        #    self.go_to_mesoffset(3)

    def begin_mesoffset1(self):
        """
        Start the first, rough, star/hole location, based on the raw data and
        the SBR input file. The first step is processing the star frames
        """
        # modify the database if necessary, and then absorb all values from the database
        img_dir = self.database['img_dir']
        if len(img_dir) <= 0 or img_dir[-1] != '/':
            self.database['img_dir'] += '/'
        self.__dict__.update(self.database)
        self.logger.info('MESOffset1 database %s' % self.database)
        self.logger.info('Setting MESOffset1 debug mode %s' % self.database['debug_mode'])
        #fitsUtils.SAVE_INTERMEDIATE_FILES = self.database['debug_mode']
        self.logger.info('debug mode %s' % self.database['debug_mode'])
        self.process_star_fits()

    def process_star_fits(self, *args):
        """ Use fitsUtils to combine raw data into a usable star mosaic """
        self.process_fits('star', self.recalc1,
                          next_step=self.load_processed_star)
    
    def load_processed_star(self, *args):
        """ Load the star frame FITS image that was processed by fitsUtils """
        self.fv.gui_call(self.displayImage, filename=self.rootname+"_star.fits",
                       next_step=self.mes_star)
        

    def mes_star(self, *args):
        """ Call MESLocate in star mode on the current image """
        self.sbr_data = self.read_sbr_file(os.path.join(self.training_dir,self.rootname+".sbr"), self.logger)
        self.startLocate(self.sbr_data, 'star', self.interact1,
                              next_step=self.check_mes_star)

    def check_mes_star(self, *args):
        """ Review the results from mes_star and give a chance to retry """
        
        # Rearrange positional data to follow the mask setting.
        self.rearrengePositions()

        # Check if there is a attribute for centroid from detailed inspection
        if hasattr(self,'output_data'):
            self.star_locations = self.output_data[:,:2]
        else:
            self.star_locations = self.moircsAlignImage.starMat

        self.logger.info('MoircsAlign check_mes_star: star_locations %s' % self.star_locations)
        self.logger.info('MoircsAlign check_mes_star: hole_locations %s' % self.hole_locations)
        self.check(self.star_locations,
                                 last_step=self.mes_star,
                                 next_step=self.res_viewer_1)
    def res_viewer_1(self, *args):
        """ Call MESAnalyze on the data from mes_star and mes_hole """
        
        # Exam where get the method call
        proc_num = self.parameter_tabs.get_index()

        print(proc_num)
        if proc_num == 0:
            self.startAnalyze(self.star_locations, self.hole_locations,
                               next_step=self.end_mesoffset1)
        if proc_num == 1:
            self.startAnalyze(self.star_locations, self.hole_locations,
                               next_step=self.end_mesoffset2)
        if proc_num == 2:
            self.startAnalyze(self.star_locations, self.hole_locations,
                               next_step=self.end_mesoffset3)

    def end_mesoffset1(self, *args):
        """ Finish off the first, rough, star/hole location """
        self.logwindow("Done with MES Offset 1!\n")
        self.write_to_logfile(self.rootname_logfile,
                                            "MES Offset 1",
                                            self.offset)
        self.logger.info('MESOffset1 offsets are dx %s pix dy %s pix rotate %s deg' % self.offset)
        self.database['starhole_chip1'] =  self.database['mask_chip1'] + 2
        self.go_to_mesoffset(1)

    ### ----- MESOFFSET2 FUNCTIONS ----- ###
    def begin_mesoffset2(self):
        """
        Start the second, star-hole location, based on the raw data. The first
        step is processing the starhole frames
        """
        # modify the database if necessary, and then absorb all values from the database
        img_dir = self.database['img_dir']
        if len(img_dir) <= 0 or img_dir[-1] != '/':
            self.database['img_dir'] += '/'
        self.__dict__.update(self.database)
        self.logger.info('MESOffset2 database %s' % self.database)
        

        self.moircsAlignImage.star_fits_name = [self.database['starhole_chip1'], 
                self.database['starhole_chip1']+1]
        self.moircsAlignImage.sky_fits_name = [self.database['mask_chip1'], self.database['mask_chip1']+1]
        self.moircsAlignImage.MAKEMOSAIC = self.database['recalc2']

        if not hasattr(self, 'hole_locations'):
            self.logwindow("No hole position data found; please run "+
                                   "literally any MESOffset besides 2.",
                                   level='error')
            return

        self.process_starhole_fits()
    
    def process_starhole_fits(self, *args):
        """ Use fitsUtils to combine raw data into a compound starhole image """
        self.process_fits('starhole', self.recalc2,
                          next_step=self.load_processed_starhole)
    
    def load_processed_starhole(self, *args):
        """ Load the finished starhole image into ginga """
        self.fv.gui_call(self.displayImage, filename=self.rootname+"_starhole.fits",
                       next_step=self.mes_starhole)
    
    def mes_starhole(self, *args):
        """ Call MESLocate in starhole mode on the current image """
        self.logger.info('MESOffset mes_starhole calling self.mes_locate.start with self.hole_locations as %s' % self.hole_locations)
        self.startLocate(self.hole_locations, 'starhole', self.interact2,
                              next_step=self.check_mes_starhole)
    
    def check_mes_starhole(self, *args):
        """ Review the results from mes_hole and offer a chance to retry """
        self.star_locations = self.output_data[:,:2]
        self.logger.info('MESOffset check_mes_starhole calling MESInterface.check with star_locations %s' % self.star_locations)
        self.check(self.star_locations,
                                 last_step=self.mes_starhole,
                                 next_step=self.res_viewer_2)
    
    def res_viewer_2(self, *args):
        """ Call MESAnalyze on the data from mes_star and mes_starhole """
        self.startAnalyze(self.star_locations, self.hole_locations,
                               next_step=self.end_mesoffset2)
    
    def end_mesoffset2(self, *args):
        """ Finish off the second star-hole location """
        self.logwindow("Done with MES Offset 2!\n")
        self.write_to_logfile(self.rootname_logfile,
                                            "MES Offset 2",
                                            self.offset)
        self.logger.info('MESOffset2 offsets are dx %s pix dy %s pix rotate %s deg' % self.offset)
        self.database['sky_chip1']=self.sky_chip1
        self.database['starhole_chip1'] = self.starhole_chip1+2
        self.go_to_mesoffset(1)


    def begin_mesoffset3(self):
        """
        Start the third, fine, star-hole location with updated mask locations.
        The first step is processing the mask frames
        """
        img_dir = self.database['img_dir']
        if len(img_dir) <= 0 or img_dir[-1] != '/':
            self.database['img_dir'] += '/'
        self.__dict__.update(self.database)
        

        self.moircsAlignImage.star_fits_name = [self.database['starhole_chip1'], 
                self.database['starhole_chip1']+1]
        self.moircsAlignImage.sky_fits_name = [self.database['sky_chip1'], 
                self.database['sky_chip1']+1]
        self.moircsAlignImage.mask_fits_name = [self.database['newmask_chip1'], 
                self.database['newmask_chip1']+1]
        self.moircsAlignImage.MAKEMOSAIC = self.database['recalc3']
        
        # Set this parameter because we re-use the MESOFFSET 1 code
        self.recalc1=self.database['recalc3']


        self.logger.info('MESOffset3 database %s' % self.database)
        self.process_star_fits()

    def end_mesoffset3(self, *args):
        """ Finish off the third, fine star-hole location with updated masks """
        self.logwindow("Done with MES Offset 3!\n")
        self.write_to_logfile(self.rootname_logfile,
                                            "MES Offset 3",
                                            self.offset)
        self.logger.info('MESOffset3 offsets are dx %s pix dy %s pix rotate %s deg' % self.offset)
        self.database['starhole_chip1'] = self.starhole_chip1 + 2
        self.go_to_mesoffset(0)


    def set_defaults(self, page_num):
        """
        Set the default values for the gui on page page_num
        @param page_num:
            The number for the GUI whose defaults we must set
        """
        setters = self.set_value[page_num]
        for key in setters:
            if key in self.database:
                setters[key](self.database[key])

    def write_to_logfile(self, filename, header, values):
        """
        Write args to a log file
        @param filename:
            The name of the log file
        @param header:
            The string that will serve as the first line of this log file entry
        @param args:
            The tuple of float values to be logged in the form (dx, dy, rotate)
        """
        # write the informaton to the file
        f = open(filename, 'a')
        f.write("="*50+"\n")
        f.write(header+"\n")
        f.write(strftime("%a %b %d %H:%M:%S %Z %Y\n"))
        f.write(("dx = {:7,.2f} (pix) dy = {:7,.2f} (pix) "+
                 "rotate = {:7,.4f} (degree) \n").format(*values))
        f.close()
 

    def resume_process_cb(self, *args):
        """
        Take the parameters from the waiting gui and resume the current process
        """
        try:
            if self.last_wait_gui == 1:
                self.update_parameters(self.get_value[4])
            elif self.last_wait_gui == 3:
                self.update_parameters(self.get_value[5])
        except NameError as e:
            self.log("NameError: "+str(e), level='e')
            return
        self.resume_mesoffset[self.last_wait_gui]()

    def check(self, data, last_step=None, next_step=None):
        """
        Go to the 'check' GUI and let the user review their data, and decide
        whether they want to remeasure those locations or move on
        @param data:
            The location data in the form of a 2-3 column numpy array (x,y[,r])
        @param last_step:
            The function to be executed if these data are unsatisfactory
        @param next_step:
            The function to be executed if these data are deemed good enough
        """
        if data.shape[1] == 2:
            res_string = "   x      y\n"
            fmt_string = "{:5.0f}  {:5.0f}\n"
        elif data.shape[1] == 3:
            res_string = "   x      y      r\n"
            fmt_string = "{:5.0f}  {:5.0f}  {:5.1f}\n"
        
        for row in data:
            res_string += fmt_string.format(*row)
        self.results_textarea.set_text(res_string)
        self.last_step = last_step
        self.next_step = next_step
        self.go_to_gui('check')

    def retake_measurements_cb(self, *args):
        """
        Execute self.last_step
        """
        self.clear_canvas()
        self.last_step()

    def return_to_menu_cb(self, *args):
        """
        Go back to the last parameter menu you were at - either 'epar' or 'wait'
        """
        self.clear_canvas()
        if self.last_wait_gui:
            self.go_to_gui('wait '+str(self.last_wait_gui))
        else:
            self.go_to_gui('epar')

    def go_to_mesoffset(self, idx):
        """
        Go to the 'epar' TabWidget and ask the user for parameters for this
        mesoffset process
        @param idx:
            The index of the process we want parameters for
        """
        self.last_wait_gui = 0
        self.set_defaults(idx)
        self.parameter_tabs.set_index(idx)
        self.go_to_gui('epar')

    def wait(self, idx, next_step=None):
        """
        Go to the 'wait' gui at this index to get more info from the user, and
        prepare to execute the next step.
        @param idx:
            The index of the process that this interrupts
        @param next_step:
            The function to be called when the 'Go!' button is pressed
        """
        self.last_wait_gui = idx
        if idx == 1:
            self.set_defaults(4)
        elif idx == 3:
            self.set_defaults(5)
        self.go_to_gui('wait '+str(idx))
        self.resume_mesoffset[idx] = next_step

    def go_to_gui(self, gui_name):
        """
        Go to the appropriate GUI, and set appropriate callbacks
        @param gui_name:
            The string identifier for this GUI
        """
        if gui_name == 'check':
            self.set_callbacks(right_click=self.use_measurements_cb)
        if gui_name == 'error':
            self.set_callbacks(right_click=self.return_to_menu_cb)
        self.stack.set_index(self.stack_idx[gui_name])
        
    def set_position_cb(self, *args):
        """
        Respond to the spinboxes being used by essentially making a new click
        """
        self.canvas.delete_object_by_tag(self.tag(1, self.click_index))
        self.color_index -= 1
        self.click1_cb(None, None,
                       self.spinboxes['X'].get_value(),
                       self.spinboxes['Y'].get_value())

    def reset_cb(self, *args):
        """
        This function reset the whole canvas to the origingal situation.
        """
        # Deleting all objects 
        self.canvas.delete_all_objects()
        # Now, replot it back to the canvas
        self.color_index = 0
        self.click_index = 0
        self.select_point(self.click_history[0])

    def undo1_cb(self, *args):
        """
        Respond to the undo button in step 1
        by going back one click (if possible)
        """
        if self.click_index > 0:
            self.canvas.delete_object_by_tag(self.tag(1, self.click_index))
            self.click_index -= 1
            self.color_index -= 1
            self.canvas.delete_object_by_tag(self.tag(1, self.click_index))
            self.select_point(self.click_history[self.click_index])
    
    
    def redo1_cb(self, *args):
        """
        Respond to the redo button in step 1
        by going forward one click (if possible)
        """
        if self.click_index < len(self.click_history)-1:
            self.click_index += 1
            self.color_index += 1
            self.select_point(self.click_history[self.click_index])

    def step1_cb(self, *args):
        """
        Respond to back button by returning to step 1
        """
        self.go_to_gui('find')
        self.setLocatecallbacks(step=1)
        self.select_point(self.click_history[self.click_index])
    
    def set_skip_btn_state(self):
        # Disable Skip button for first object. First object is used
        # as reference position and code gets confused if that object
        # is skipped.
        if self.current_obj == 0:
            self.skip_btn.set_enabled(False)
        else:
            self.skip_btn.set_enabled(True)
    
    def step2_cb(self, *args):
        """
        Respond to next button or right click by proceeding to the next step
        """
        # set everything up for the first object of step 2
        self.logger.info('MESLocate mode %s coordinates are %s' % (self.mode, str(self.click_history[self.click_index])))
        self.go_to_gui('centroid')
        self.setLocatecallbacks(step=2)
        # Disable Skip button for first object.
        self.set_skip_btn_state()
        self.select_point(self.click_history[self.click_index], True)
        self.zoom_in_on_current_obj()
        self.mark_current_obj()
        
        # in the event that there are masks here (because of the Back button)...
        for i in range(self.obj_num):
            for j in range(self.drag_index[i]+1):
                self.draw_mask(*self.drag_history[i][j], obj_idx=i, drag_idx=j)
        
        # if interaction is turned off, immediately go to the next object
        if not self.interact:
            self.next_obj_cb()

    def zoom_in_on_current_obj(self):
        """
        Set the position and zoom level on fitsimage such that the user can
        only see the object at index self.current_obj. Also locates and marks it
        """
        x1, y1, x2, y2, r = self.get_current_box()
        self.fitsimage.set_pan((x1+x2)/2, (y1+y2)/2)
        self.fitsimage.zoom_to(350.0/self.square_size)
        self.obj_count_label.set_text("Object {} out of {}".format(
                                            self.current_obj+1, self.obj_num))
        self.mark_current_obj()

    def get_current_box(self, idx=None):
        """
        Calculates the bounds of the box surrounding the current object
        @precondition:
            This method only works in step 2
        @param idx:
            The object index at the instant that we need the box (defaults to
            self.current_obj)
        @returns x1, y1, x2, y2, r:
            The float bounds and radius of the current box
        """
        if idx == None:
            idx = self.current_obj
        xf, yf = self.click_history[self.click_index]
        dx, dy, r = self.obj_arr[idx]
        s = self.square_size
        if math.isnan(r):
            r = 1.42*s
        return (xf+dx-s, yf+dy-s, xf+dx+s, yf+dy+s, r)
    
    
    def mark_current_obj(self, obj=None):
        """
        Puts a point and/or circle on the current object
        @param obj:
            The exact coordinates (x, y, r) of this object. If none are
            provided, they will be calculated using locate_obj
        """
        # if there is already a point for this object, delete it
        t = self.tag(2, self.current_obj, 'pt')
        self.canvas.delete_object_by_tag(t)
        
        # then locate and draw the point (if it exists)
        sq_size = self.square_size
        if obj == None:
            co = self.current_obj
            obj = self.locate_obj(self.get_current_box(),
                             self.drag_history[co][:self.drag_index[co]+1],
                             self.fitsimage.get_image(),
                             min_search_radius=self.exp_obj_size,
                             viewer=self.step2_viewer)
        
        # if any of the coordinates are NaN, then a red x will be drawn in the middle
        if True in [math.isnan(x) for x in obj]:
            x1, y1, x2, y2, r = self.get_current_box()
            self.canvas.add(self.dc.Point((x1+x2)/2, (y1+y2)/2, sq_size/3,
                                          color='red', linewidth=1,
                                          linestyle='dash'),
                            tag=t)
        else:
            self.canvas.add(self.dc.CompoundObject(
                                    self.dc.Circle(obj[0], obj[1], obj[2],
                                                   color='green', linewidth=1),
                                    self.dc.Point(obj[0], obj[1], sq_size/3,
                                                  color='green', linewidth=1)),
                            tag=t)
        
        self.obj_centroids[self.current_obj] = obj

    def finish(self):
        """
        Finishes up and goes to the next step
        """
        # fix up the canvas and clear callbacks
        self.clear_canvas(keep_objects=True)
        
        # tell the manager to do whatever comes next
        self.output_data = np.array(self.obj_centroids)
        if self.next_step != None:
            self.next_step()

    def start_drag_cb(self, c, k, x, y, kind):
        """
        Respond to the mouse finishing a left-drag by finalizing crop selection
        @param kind:
            Either 'mask' or 'crop', depending on the mouse button and the
            current selection mode
        @returns:
            True, in order to prevent the other panset callbacks from going off
        """
        # enable drawing and then start drawing
        self.canvas.enable_draw(True)
        fill = (kind == 'mask')
        self.canvas.set_drawtype(drawtype='rectangle', color='white',
                                 fill=fill, fillcolor='black')
        self.canvas.draw_start(c, k, x, y, self.fitsimage)
        self.drag_start = (x,y)
        return True
    
    
    def end_drag_cb(self, c, k, xf, yf, kind):
        """
        Respond to the mouse finishing a left-drag by finalizing crop selection
        @param xf:
            The final x coordinate of this drag
        @param yf:
            The final y coordinate of this drag
        @param kind:
            Either 'mask' or 'crop', depending on the mouse button and the
            current selection mode
        """
        # if rectangle has area zero, ignore it
        if (xf,yf) == self.drag_start:
            return
        
        # finish the drawing, but make sure nothing is drawn; it won't be visible anyway
        self.canvas.draw_stop(c, k, *self.drag_start, viewer=None)
        self.canvas.enable_draw(False)
        # if anything is ahead of this in drag_history, clear it
        co = self.current_obj
        self.drag_index[co] += 1
        if self.drag_index[co] < len(self.drag_history[co]):
            self.drag_history[co] = self.drag_history[co][:self.drag_index[co]]
        
        # now save that crop
        x1, y1, x2, y2, r = self.get_current_box()
        xi, yi = self.drag_start
        self.drag_history[co].append((min(max(min(xi, xf), x1), x2),
                                      min(max(min(yi, yf), y1), y2),
                                      min(max(max(xi, xf), x1), x2),
                                      min(max(max(yi, yf), y1), y2),
                                      kind))
        
        # shade in the outside areas and remark it
        self.draw_mask(*self.drag_history[self.current_obj][-1])
        self.mark_current_obj()
    

    def viewer_redirect_cb(self, _, button, xv, yv, direction):
        """
        Respond to a click on the step2_viewer by adjusting x and y for the
        position of the viewer, and then making the canvas react as though it
        had been clicked on.
        @param button:
            The key code representing the mouse configuration (1 for cursor,
            2 for panset, and 4 for draw)
        @param xv:
            The x coordinate of the click relative to the viewer
        @param yv:
            The y coordinate of the click relative to the viewer
        @param direction:
            The direction of the event ('up' or 'down')
        """
        x0, y0, x1, y1, r = self.get_current_box()
        xf, yf = xv+int(x0), yv+int(y0)
        if button%2 == 1:
            self.canvas.make_callback('cursor-'+direction, button, xf, yf)
        if button%4/2 == 1:
            self.canvas.make_callback('panset-'+direction, button, xf, yf)
        if button%8/4 == 1:
            self.canvas.make_callback('draw-'+direction,   button, xf, yf)
    
    def undo2_cb(self, *args):
        """
        Respond to the undo button in step 2
        by going back one drag (if possible)
        """
        co = self.current_obj
        if self.drag_index[co] >= 0:
            self.canvas.delete_object_by_tag(self.tag(2, co, self.drag_index[co]))
            self.drag_index[co] -= 1
            self.mark_current_obj()
    
    
    def redo2_cb(self, *args):
        """
        Respond to the redo button in step 2
        by going forward one drag (if possible)
        """
        co = self.current_obj
        if self.drag_index[co] < len(self.drag_history[co])-1:
            self.drag_index[co] += 1
            self.draw_mask(*self.drag_history[co][self.drag_index[co]])
            self.mark_current_obj()

    def terminate_process_cb(self, *args):
        """
        Stop the process currently connected to self.manager.terminate,
        wait for the process to receive the signal, and then go to the main menu
        """
        if hasattr(self, 'terminate'):
            self.terminate.set()
            self.logwindow("Terminating process...")

    def prev_obj_cb(self, *args):
        """
        Respond to back button in step 2 by going back to the last object
        """
        # if there is no previous object, return to step 1
        if self.current_obj > 0:
            self.current_obj -= 1
            # Current object has changed, so set Skip button state
            # based on current object number.
            self.set_skip_btn_state()
            self.zoom_in_on_current_obj()
            self.mark_current_obj()
        else:
            self.step1_cb()
        
    
    def next_obj_cb(self, *args):
        """
        Respond to next button or right click by proceeding to the next object
        """
        # if there is no next object, finish up
        self.current_obj += 1
        if self.current_obj >= self.obj_num:
            self.finish()
            return
        # Current object has changed, so set Skip button state based
        # on current object number.
        self.set_skip_btn_state()

        x1, y1, x2, y2, r = self.get_current_box()
        # Skip over current object if its x or y coordinates are NaN.
        if math.isnan(x1) or math.isnan(y1) or math.isnan(x2) or math.isnan(y2):
            self.obj_centroids[self.current_obj] = (np.nan, np.nan, np.nan)
            self.current_obj += 1

        if self.current_obj >= self.obj_num:
            self.finish()
            return
        self.set_skip_btn_state()
            
        # if there is one, focus in on it
        self.zoom_in_on_current_obj()
        self.mark_current_obj()

        # if interaction is turned off, immediately go to the next object
        if not self.interact:
            self.next_obj_cb()

    def choose_select_cb(self, _, mode_idx):
        """
        Keep track of our selection mode as determined by the combobox
        """
        # update the callbacks to account for this new mode
        self.set_callbacks(step=2, selection_mode=mode_idx, clear=False)    
    
    def skip_obj_cb(self, *args):
        """
        Respond to the skip button by ignoring this star completely
        """
        self.mark_current_obj([float('NaN')]*3)
        self.next_obj_cb()

    
    def use_measurements_cb(self, *args):
        """
        Execute self.next_step
        """
        self.clear_canvas()
        self.next_step()

    def image_set_cb(self, *args):
        """
        Respond to an image being loaded by executing whatever function
        """
        if self.fitsimage.get_image() == None:  # I don't know why this callback sometimes
            return                              # executes for no reason, so I just ignore it.
        if self.image_set_next_step != None:
            self.image_set_next_step()
        self.image_set_next_step = None


    def read_parameters(self):
        """
        Get the last parameters that were used for mesoffset
        @returns:
            A dictionary where keys are parameter names, and values are values,
            or None if the file could not be found or was in the wrong format.
        """
        try:
            output = {}
            par_file = open(os.path.join(DIR_PAR_VAR,PAR_FILENAME), 'r')
            line = par_file.readline().strip()
            while line != "":
                idx = line.index(',')
                output[line[:idx]] = line[idx+1:]
                line = par_file.readline().strip()
            return output
        except Exception as e:
            pass
    
    def _log(self, text, level='i'):
        """
        Print text to the logger TextArea
        @param text:
            The string to be logged
        @param level:
            The level of urgency ('d' for debug, 'i' for info, etc.)
        """
        if level[0].lower() == 'd':
            self.logger.debug(text.strip())
        elif level[0].lower() == 'i':
            self.logger.info(text.strip())
            self.log_textarea.append_text(text+"\n", autoscroll=True)
        elif level[0].lower() == 'w':
            self.logger.warning(text.strip())
            self.log_textarea.append_text("WARN: "+text+"\n", autoscroll=True)
        elif level[0].lower() == 'e':
            self.logger.error(text.strip())
            self.log_textarea.append_text("ERROR: "+text+"\n", autoscroll=True)
            self.err_textarea.set_text(text)
            self.go_to_gui('error')
        else:
            self.logger.critical(text.strip())
            self.log_textarea.append_text("CRIT: "+text+"\n", autoscroll=True)
            self.err_textarea.set_text("CRITICAL!\n"+text)
            self.go_to_gui('error')
    
    def set_callbacks(self, left_click=None, right_click=None):
        """
        Set some basic callbacks
        @param left_click:
            The function to run when the user left clicks
        @param right_click:
            Guess.
        @param enter:
            The function to run when the user hits return
        """
        self.clear_canvas(keep_objects=True)
        self.canvas.clear_callback('cursor-up')
        if left_click != None:
            self.canvas.add_callback('cursor-up', left_click)
        self.canvas.clear_callback('draw-up')
        if right_click != None:
            self.canvas.add_callback('draw-up', right_click)



    def logwindow(self, *args, **kwargs):
        """
        Print text to the logger TextArea from the main thread
        """
        self.fv.gui_do(self._log, *args, **kwargs)

    def read_variables(self):
        """
        Get the defined variable dictionary from mesoffset_directories.txt
        @returns:
            A dictionary where keys are variable names, and values are values
        """
        output = {}
        try:
            var_file = open(os.path.join(DIR_PAR_VAR,VAR_FILENAME), 'r')
            line = var_file.readline()
            while line != "":
                words = line.split()
                output[words[0]] = words[1]
                line = var_file.readline()
            var_file.close()
        except IOError:
            #output["DATABASE"] = "../../MCSRED2/DATABASE"
            output["DATABASE"] = os.path.join(self.DIR_MCSRED, "DATABASE")
        return output

    def process_filename(self,filename, variables=None):
        """
        Take a filename and modifies it to account for any variables
        @param filename:
            The input filename to be processed
        @param variables:
            The dictionary of defined variable names to values
        @returns:
            The updated filename
        @raises NameError:
            If there is an undefined variable
        """
        if variables == None:
            variables = self.read_variables()
        
        # scan the filename for dollar signs
        while "$" in filename:
            ds_idx = filename.find("$")
            sl_idx = filename[ds_idx:].find("/")
            if sl_idx == -1:
                sl_idx = len(filename)
            var_name = filename[ds_idx+1:sl_idx]
            
            # if it is a defined variable, replace it
            if var_name.upper() in variables:
                filename = filename.replace("$"+var_name,
                                            variables[var_name.upper()])
            
            # otherwise, raise an error
            else:
                err_msg = ("$"+var_name+" is not a defined variable. Defined "+
                           "variables are:\n")
                for key in variables:
                    err_msg += "    ${}: {}\n".format(key, variables[key])
                err_msg += "Please check your spelling and try again."
                raise NameError(err_msg)
        
        return filename

    def process_fits(self, mode, recalc=True, next_step=None):
        """
        Plug some values into fitsUtils and start a new thread to create a
        processed FITS image to be loaded and used
        @param mode:
            A string - either 'star', 'mask', or 'starhole'
        @param framenum:
            The chip1 frame number for the first set of input images
        @param recalc:
            Whether new images should even be processed
        @param next_step:
            The function to be called when this is done
        @returns:
            A threading.Event object to be used if the user ever decides to
            terminate the task
        @raises IOError:
            If the specified images cannot be found
        """
        out_filename = os.path.join(self.work_dir, self.rootname+"_"+mode+".fits")
        
        # TO-DO : If this method failed. Adding error log.


        # if regenerate is off, don't bother with any of this
        # if not recalc:
        #     if os.path.isfile(out_filename):
        #         if next_step != None:
        #             next_step()
        #     else:
        #         self.logwindow("No previous image found at "+
        #                                out_filename+". Please change your "+
        #                                "working directory or rootname, or "+
        #                                "enable the 'Regenerate' option",
        #                                level='error')
        
        # # otherwise, start the appropriate process in a new thread
        # else:
        self.stack.set_index(self.stack_idx['log'])
        c, i, w = self.c_file, self.img_dir, self.work_dir
        f = out_filename
        e = threading.Event()


        #l = lambda message: self.fv.gui_call(self.logwindow, message)
        #if mode == 'star':
        #    n1, n2 = int(self.star_chip1), int(self.sky_chip1)
        #elif mode == 'starhole':
        #    n1, n2 = int(self.starhole_chip1), int(self.mask_chip1)
        #elif mode == 'mask':
        #    n1, n2 = int(self.mask_chip1), None
        
        if mode == 'star' or mode == 'mask':
            self.logger.info("Starting to processing star images")
            task = lambda: self.processStarFits(e, self.logwindow, next_step=next_step)
        elif mode == 'starhole':
            self.logger.info("Starting to processing star images")
            task = lambda: self.processStarHoleFits(e,self.logwindow, next_step=next_step)
       
        self.fv.nongui_do(task)
        self.terminate = e

    
    def processStarFits(self, t, log, next_step=None):
        """
        This is an interface for thread operation, telling program to start the
        image processing procedure.
        arguments
        @param mode:
            A string which will determine the method we call - either 'star',
            'mask', or 'starhole'
        """
        self.fv.assert_nongui_thread()
        try:
           
            self.moircsAlignImage.processStarImage(t, log, next_step=next_step)
           
            if t.is_set():
                raise RuntimeError(USER_INTERRUPT_ERR)
        except Exception as e:
            log("{}: {}".format(type(e).__name__, e), level='e')
    
    def processStarHoleFits(self, t, log, next_step=None):
        """
        This is an interface for thread operation, telling program to start the
        image processing procedure.
        arguments
        @param mode:
            A string which will determine the method we call - either 'star',
            'mask', or 'starhole'
        """
        self.fv.assert_nongui_thread()
        try:
           
            self.moircsAlignImage.processStarHoleImage(t, log, next_step=next_step)
           
            if t.is_set():
                raise RuntimeError(USER_INTERRUPT_ERR)
        except Exception as e:
            log("{}: {}".format(type(e).__name__, e), level='e')

    def displayImage(self, filename, next_step=None):
        """
        Open a FITS image and display it in ginga, then call a function
        @param filename:
            The name of the fits file
        @param next_step:
            The function to call once the image has been loaded
        """
        
        full_filename = os.path.join(self.work_dir, filename)
        self.image_set_next_step = next_step
        
        # This is a work around to prevent X-window error
        #fits.writeto(full_filename, mosaic_dataself.moircsAlignImage.mosaic.data)
        try:
            # Before display an image, we have to put it into a frame
            image = AstroImage.AstroImage(logger=self.logger)
            image.set_data(self.moircsAlignImage.star.mosaic.data)
            image.set(name='MOIRCS Guiding')
            self.fitsimage.set_image(image)
            self.fitsimage.gui_do
        except:
            self.logger.info("Can't display image becuase of backend error.")
        else:    
            self.fitsimage.make_callback('drag-drop', [full_filename])

    def __str__(self):
        return type(self).__name__

    def select_point(self, point, draw_circle_masks=True):
        """
        Set a point in step 1 as the current location of object #0,
        draws squares where it thinks all the objects are accordingly,
        and updates all thumbnails and spinboxes
        @param point:
            An int tuple containing the location of object #0
        @param draw_circle_masks:
            Whether we should draw the automatic circular masks
        """
        # define some variables before iterating through the objects
        x, y = point
        src_image = self.fitsimage.get_image()
        color = self.BOX_COLORS[self.color_index%len(self.BOX_COLORS)]   # cycles through all the colors
        shapes = []
        for i, viewer in enumerate(self.thumbnails):
            dx, dy, r = self.obj_arr[i]
            sq_size = self.square_size
            # Skip this object if its dx or dy are NaN (this is an
            # object that was skipped in previous step)
            if math.isnan(dx) or math.isnan(dy):
                continue
        
            # first, draw squares and numbers
            shapes.append(self.dc.SquareBox(x+dx, y+dy, sq_size, color=color))
            shapes.append(self.dc.Text(x+dx+sq_size, y+dy,
                                       str(i+1), color=color))
            
            # draw the circular mask if necessary
            #if draw_circle_masks:
            #    if r <= sq_size:
            #        shapes.append(empty_circle(x+dx, y+dy, r, sq_size, self.dc))
            #        shapes.append(self.dc.Circle(x+dx, y+dy, r, color='white'))

            # then, update the little pictures
            xylimits = (x-sq_size+dx, y-sq_size+dy,
                        x+sq_size+dx, y+sq_size+dy)
            # TODO: cutout_adjust requires integer values. Should we
            # use int function or ceil function or something else to
            # convert float values to integer values?
            cropped_data = src_image.cutout_adjust(*[int(z) for z in xylimits])[0]
            viewer.set_data(cropped_data)
            self.fitsimage.copy_attributes(viewer,
                                           ['transforms','cutlevels','rgbmap'])
        
        # First, check if the number of detected holes is consistent with input
        #  data.
        hole, hole0, dx, dy = self.parse_data(self.sbr_data)
        #print('Input hole number = %d'%len(hole[:,0]))
        #print('Detected holes = %d'%self.moircsAlignImage.mask.gHoleMosaic.shape[1])
        
        # If the number is consistent, plot the hole location    
        if len(hole[:,0]) == self.moircsAlignImage.mask.gHoleMosaic.shape[1]:
            # draw the circular mask if necessary
            if draw_circle_masks:
                if self.moircsAlignImage.mask != 0:
                    circle=self.moircsAlignImage.mask.gHoleMosaic
                    for i in circle[0,:]:
                        shapes.append(self.dc.SquareBox(i[0], i[1], sq_size, color='yellow'))
                        shapes.append(self.dc.Circle(i[0], i[1], i[2], color='yellow'))
        else:
            self.logwindow("incorrect number of detected holes.",level='w')                
            hole[:,0]=hole[:,0]+hole0[0]
            hole[:,1]=hole[:,1]+hole0[1]

            base=self.moircsAlignImage.mask.gHoleMosaic[0,0,:]
            dist=np.sqrt((base[0]-hole[:,0])**2+(base[1]-hole[:,1])**2)
            
            offsetx=float((base[0]-hole[np.where(dist == min(dist)),0])[0][0])
            offsety=float((base[1]-hole[np.where(dist == min(dist)),1])[0][0])
            hole[:,0]=hole[:,0]+offsetx
            hole[:,1]=hole[:,1]+offsety

            # Update the matrix of detected hole with input files. 
            for i in hole:
                dist=np.sqrt(
                     (i[0]-self.moircsAlignImage.mask.gHoleMosaic[0,:,0])**2+
                     (i[1]-self.moircsAlignImage.mask.gHoleMosaic[0,:,1])**2)  
                
                # This case is for the hole detection missed.  The number of 
                if np.min(dist) > 50:
                    ind=np.where(dist == min(dist))
                    mean_radius=np.mean(self.moircsAlignImage.mask.gHoleMosaic[0,:,2])
                    temp=np.vstack(
                                   (self.moircsAlignImage.mask.gHoleMosaic[0,:],
                                   [i[0],i[1],mean_radius])
                                  )            
                    self.moircsAlignImage.mask.gHoleMosaic=np.array([temp])
                if draw_circle_masks:
                    shapes.append(self.dc.SquareBox(i[0], i[1], sq_size, color='yellow'))
                    shapes.append(self.dc.Circle(i[0], i[1], 15, color='yellow'))


         # draw all the squares and numbers to the canvas as one object
        self.canvas.add(self.dc.CompoundObject(*shapes),
                        tag=self.tag(1, self.click_index))
        
        # update the spinboxes
        if self.spinboxes['X'].get_value() != x:
            self.spinboxes['X'].set_value(x)
        if self.spinboxes['Y'].get_value() != y:
            self.spinboxes['Y'].set_value(y)


    def locate_obj(self, bounds, masks, image, viewer=None,
                   min_search_radius=None, thresh=3):
        """
        Finds the center of an object using center of mass calculation
        @param bounds:
            A tuple of floats x1, y1, x2, y2, r. The object should be within
            this box
        @param masks:
            A list of tuples of the form (x1, y1, x2, y2, kind) where kind is either
            'mask' or 'crop' and everything else is floats. Each tuple in masks
            is one drag of the mouse that ommitted either its interior or its
            exterior
        @param image:
            The AstroImage containing the data necessary for this calculation
        @param viewer:
            The viewer object that will display the new data, if desired
        @param min_search_radius:
            The smallest radius that this will search
        @param thresh:
            The number of standard deviations above the mean a data point must
            be to be considered valid
        @returns:
            A tuple of two floats representing the actual location of the object
            or a tuple of NaNs if no star could be found
        """
        # start by getting the raw data from the image matrix

        # TODO: cutout_adjust requires integer values. Should we use int
        # function or ceil function or something else to convert float
        # values to integer values?
        raw, x0,y0,x1,y1 = image.cutout_adjust(*[int(z) for z in bounds[:4]])
        search_radius = bounds[4]
        x_cen, y_cen = raw.shape[0]/2.0, raw.shape[1]/2.0
        yx = np.indices(raw.shape)
        x_arr, y_arr = yx[1], yx[0]
        
        # crop data to circle
        mask_tot = np.hypot(x_arr - x_cen, y_arr - y_cen) > search_radius
        
        # mask data based on masks
        for drag in masks:
            x1, y1, x2, y2, kind = (int(drag[0])-int(x0)+1, int(drag[1])-int(y0)+1,
                                    int(drag[2])-int(x0)+1, int(drag[3])-int(y0)+1,
                                    drag[4])
            mask = np.zeros(raw.shape, dtype=bool)
            mask[y1:y2, x1:x2] = True
            if kind == 'crop':
                mask = np.logical_not(mask)
            mask_tot = np.logical_or(mask_tot, mask)
        
        # apply mask, calculate threshold, normalize, and coerce data positive
        data = ma.masked_array(raw, mask=mask_tot)
        threshold = thresh*ma.std(data) + ma.mean(data)
        data = data - threshold
        data = ma.clip(data, 0, float('inf'))
        
        # display the new data on the viewer, if necessary
        if viewer != None:
            viewer.get_settings().set(autocut_method='minmax')
            viewer.set_data(data)
        
        # exit if the entire screen is masked
        if np.all(mask_tot):
            return (float('NaN'), float('NaN'), float('NaN'))
        
        # iterate over progressively smaller search radii
        if min_search_radius == None:
            min_search_radius = search_radius/2
        has_not_executed_yet = True
        while search_radius >= min_search_radius or has_not_executed_yet:
            has_not_executed_yet = False
            old_x_cen, old_y_cen = float('-inf'), float('-inf')
            
            # repeat the following until you hit an assymptote:
            while np.hypot(x_cen-old_x_cen, y_cen-old_y_cen) >= 0.5:
                # define an array for data constrained to its search radius
                circle_mask = np.hypot(x_arr-x_cen, y_arr-y_cen) > search_radius
                local_data = ma.masked_array(data, mask=circle_mask)
                
                # calculate some moments and stuff
                mom1x = float(ma.sum(local_data*(x_arr)))
                mom1y = float(ma.sum(local_data*(y_arr)))
                mom0 = float(ma.sum(local_data))
                area = float(ma.sum(np.sign(local_data)))
                
                # now try a center-of-mass calculation to find the size and centroid
                try:
                    old_x_cen = x_cen
                    old_y_cen = y_cen
                    x_cen = mom1x/mom0
                    y_cen = mom1y/mom0
                    radius = math.sqrt(area/math.pi)
                except ZeroDivisionError:
                    return (float('NaN'), float('NaN'), float('NaN'))
            
            search_radius = search_radius/2
        
        return (x0 + x_cen - 0.5, y0 + y_cen - 0.5, radius)

    def empty_circle(self, x, y, r, a, dc):
        """
        Create a ginga canvas mixin (whatever that is) composed of a black
        filled square with a circle removed from the center
        @param x:
            The x coordinate of the center of the circle
        @param y:
            The y coordinate of the center of the circle
        @param r:
            The radius of the circle
        @param a:
            The apothem of the square around the circle
        @param dc:
            The drawing classes module
        @returns:
            A canvas.types.layer.CompoundObject, as described above
        """
        # the verticies of the polygon that will approximate this shape
        vertices = [(x+a, y+a), (x+a, y-a), (x-a, y-a), (x-a, y+a), (x+a, y+a)]
        # draw the circle
        for theta in range(45, 406, 10):
            vertices.append((x + r*math.cos(math.radians(theta)),
                             y + r*math.sin(math.radians(theta))))
        # and then fill in the outside
        return dc.Polygon(vertices, color='black', fill=True, fillcolor='black')

    
    def tag(self, step, mod_1, mod_2=None):
        """
        Create a new tag given the step and some modifiers,
        to be used by self.canvas
        @param step:
            Which step we are on
        @param mod_1:
            The primary modifer to distinguish it from other objects
        @param mod_2:
            The secondary modifier, if it is needed for additional distinction
        @returns:
            A 'tag', probably a string, to be passed into CanvasMixin.add
        
        >>> tag(1, 3, 'pt')
        '@1:3:pt'
        """
        if mod_2 == None:
            return '@{}:{}'.format(step, mod_1)
        else:
            return '@{}:{}:{}'.format(step, mod_1, mod_2)

    def rearrengePositions(self):
        """
        This method rearrange the hole positions and star positions in a better 
          format.  
        """
        hole, hole0, dx, dy = self.parse_data(self.sbr_data)
        #print(hole)
        #print(hole0)
        #print(self.star_locations)
        #print(self.hole_locations)
        hole[:,0]=hole[:,0]+hole0[0]
        hole[:,1]=hole[:,1]+hole0[1]
        #print(self.moircsAlignImage.mask.gHoleMosaic)
        #print(self.moircsAlignImage.starMat)
 

        # Check the class has this attribute, if not, create this one. 
        if not hasattr(self,'hole_locations'):
            self.hole_locations = np.reshape(self.moircsAlignImage.mask.gHoleMosaic,(
                    self.moircsAlignImage.mask.gHoleMosaic.shape[1],
                    self.moircsAlignImage.mask.gHoleMosaic.shape[2]))

        # The sequence of detected stars is based on the mask numbering.  
        #   So, rearrange it. 

        tmp_array=np.array([])
        for p in hole[:,:2]:
            d=np.sqrt((p[0]-self.moircsAlignImage.starMat[:,0])**2+
                (p[1]-self.moircsAlignImage.starMat[:,1])**2)
            idx=np.where(d == np.min(d))
            if len(tmp_array) == 0:
                tmp_array = self.moircsAlignImage.starMat[idx]
            else:
                tmp_array = np.vstack((tmp_array,self.moircsAlignImage.starMat[idx]))

        self.moircsAlignImage.starMat = tmp_array
        
        ghole=self.moircsAlignImage.mask.gHoleMosaic
        ghole=np.reshape(ghole,(ghole.shape[1],ghole.shape[2]))
        tmp_array=np.array([])
        for p in hole[:,:2]:
            d=np.sqrt((p[0]-ghole[:,0])**2+(p[1]-ghole[:,1])**2)
            idx=np.where(d == np.min(d))
            if len(tmp_array) == 0:
                tmp_array = ghole[idx]
            else:
                tmp_array = np.vstack((tmp_array,ghole[idx]))
        #self.moircsAlignImage.mask.gHoleMosaic = tmp_array       
        self.hole_locations = tmp_array

        #print(self.moircsAlignImage.starMat)
    
    def startAnalyze(self, star_pos, hole_pos, next_step=None):
        """
        Analyze the data from MESLocate
        @param star_pos:
            A 2-3 column (x,y[,r]) array specifying the star locations and sizes
        @param hole_pos:
            A 2-3 column (x,y[,r]) array specifying the hole locations and sizes
        @param rootname:
            The string that will be used for all temporary filenames
        @param next_step:
            The function to call when this process is done
        """
        # set attributes
        self.data, self.active = self.parseStarHoledata(star_pos, hole_pos)
        self.next_step = next_step
        
        # set the mouse controls
        self.setAnalyzecallbacks()
        
        # adjust the cut levels to make the points easier to see
        self.fitsimage.get_settings().set(autocut_method='stddev')
        
        # initialize the plots
        self.delete_outliers()
        
        # show the GUI
        self.go_to_gui('plots')
        
        self.display_values()

    def parseStarHoledata(self, data1, data2):
        """
        Read the data and return it in a more useful format: a four-columned
        numpy array with the nans removed
        @param data1:
            The first input array: the star locations and/or sizes (ref values)
        @param data2:
            The second input array: the hole locations and/or sizes (in values)
        @returns:
            A four-column array representing the star positions and hole
            positions, and a 1-dimensional array of Trues
        """
        data = np.hstack((data2[:,:2], data1[:,:2]))
        real_idx = np.logical_not(np.any(np.isnan(data), axis=1))
        data = data[np.nonzero(real_idx)]
        
        return data, np.ones(data.shape[0], dtype=bool)

  
    def setAnalyzecallbacks(self, step=3):
        """
        Assign all necessary callbacks to the canvas for the current step
        """
        canvas = self.canvas
        
        # clear all existing callbacks first
        self.clear_canvas()
        
        # the only callbacks are for right-click and left-click
        if step == 3:
            canvas.add_callback('cursor-down', self.set_active_cb, True)
            canvas.add_callback('draw-down', self.set_active_cb, False)
    
    def delete_outliers(self):
        """
        Remove any data points with residuals of absolute values greater than 1.
        Also updates the plots
        """
        active = self.active
        
        xres, yres = self.update_plots()
        residual_mag = np.hypot(xres, yres)*active
        
        # as long as some residuals are out of bounds,
        while np.any(residual_mag > self.OUTLIER_VALUE):
            # delete the point with the worst residual
            idx = np.argmax(residual_mag)
            active[idx] = False
            
            xres, yres = self.update_plots()
            residual_mag = np.hypot(xres, yres)*active

    def transform(self, x, y, trans):
        """
        Applies the given transformation to the given points
        @param x:
            A numpy array of x positions
        @param y:
            A numpy array of y positions
        @param trans:
            A tuple of floats: (x_shift, y_shift, rotation in degrees)
        @returns:
            A tuple of the new x value array and the new y value array
            or NaN, NaN if trans was None
        Algorithm is from IRAF geomap documentation at
        http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?geomap#description
        """
        if trans == None:
            return float('NaN'), float('NaN')
        
        xshift, yshift, thetaR = trans
        b = math.cos(thetaR)
        c = math.sin(thetaR)
        e = -c
        f = b
        newX = xshift + b * x + c * y
        newY = yshift + e * x + f * y

        return newX, newY


    def plot_residual(self, plot, z_observe, z_residual, active, var_name=""):
        """
        Plot the residual of this data against the real value.
        Residual is defined as the difference between the calculated value of
        zref and the observed value of zref.
        @param plot:
            A ginga.util.plots.Plot object, onto which to draw the graph
        @param z_observe:
            A numpy array of the observed values of this variable
        @param z_residual:
            A numpy array of the residuals for this variable
        @param active:
            A numpy array representing which data are active, and which are not
        @param var_name:
            The name of this variable, if it has one
        """
        # separate the active and inactive data
        inactive = np.logical_not(active)
        active_x = z_observe[np.nonzero(active)]
        active_y = z_residual[np.nonzero(active)]
        inactive_x = z_observe[np.nonzero(inactive)]
        inactive_y = z_residual[np.nonzero(inactive)]
        
        # then plot reference values by residual values
        try:
            plot.clear()
        except AttributeError:
            plot.add_axis()
        plot.plot(active_x, active_y,
                  linestyle='None', marker='o', color='blue')
        plot.plot(inactive_x, inactive_y,
                  linestyle='None', marker='o', color='grey',
                  xtitle="{0} Position (pixels)".format(var_name),
                  ytitle="{0} Residual (pixels)".format(var_name),
                  title="{0} Residual by {0}-axis".format(var_name),
                  markerfacecolor='none')
        plot.get_axis().set_ylim([-self.OUTLIER_VALUE-1.5, self.OUTLIER_VALUE+1.5])
        
        plot.xdata = z_observe
        plot.ydata = z_residual
        
        # shade in regions y > 1 and y < -1
        xlimits = plot.get_axis().get_xlim()
        ylimits = plot.get_axis().get_ylim()
        plot.get_axis().fill_between(xlimits, self.OUTLIER_VALUE, ylimits[1]+1,
                                     color='red', alpha=0.3)
        plot.get_axis().fill_between(xlimits, -self.OUTLIER_VALUE, ylimits[0]-1,
                                     color='red', alpha=0.3)
        plot.get_axis().set_xlim(left=xlimits[0], right=xlimits[1])
        plot.get_axis().set_ylim(bottom=ylimits[0], top=ylimits[1])
        
        plot.draw()

    
    # def step4_cb(self, *args):
    #     """
    #     Respond to the Next button in step 3 by displaying the offset values
    #     """
    #     self.setAnalyzecallbacks(step=4)
    #     self.go_to_gui('values')
    #     self.display_values()
    
    
    def set_active_cb(self, _, __, x, y, val):
        """
        Respond to a right or left click on the main ImageViewer by altering the
        datum nearest the cursor
        @param x:
            The x coordinate of the click
        @param y:
            The y coordinate of the click
        @param val:
            The new active value for the point - should be boolean
        """
        distance_from_click = np.hypot(self.data[:,0] - x, self.data[:,1] - y)
        idx = np.argmin(distance_from_click)
        self.active[idx] = val
        self.logger.info('MESAnalyze data point with index %s set to status %s' % (idx, val))
        self.update_plots()
        self.display_values()

    def relocate_cb(self, *args):
            
        proc_num = self.parameter_tabs.get_index()
        if proc_num == 0 or proc == 2:
            self.mes_star
        if proc_num == 1:
            self.mes_starhole

    def finish_cb(self, *args):
        """
        Respond to the 'Finish' button in step 4 by finishing up
        """
        self.clear_canvas()
        if self.next_step != None:
            self.next_step()

    def toggle_active_x_cb(self, e):
        """ Redirect to toggle_active_cb """
        self.toggle_active_cb(e, self.plots[0])
    
    
    def toggle_active_y_cb(self, e):
        """ Redirect to toggle_active_cb """
        self.toggle_active_cb(e, self.plots[1])
    
    
    def toggle_active_cb(self, event, plt):
        """
        Respond to right or left click on one of the plots by altering the datum
        nearest the cursor
        @param event:
            The matplotlib.backend_bases.MouseEvent instance containing all
            important information
        @param plt:
            The MOSPlot that got clicked on
        """
        # first check that the click was on a plot
        x_click, y_click = event.xdata, event.ydata
        if x_click == None or y_click == None:
            return
        
        # extract some info from plt and use it to calculate distances
        fig, x_arr, y_arr = plt.get_data()
        (left, bottom), (right, top) = plt.get_axis().viewLim.get_points()
        dx = (x_arr - x_click) * fig.get_figwidth() / (right - left)
        dy = (y_arr - y_click) * fig.get_figheight() / (top - bottom)
        distance_from_click = np.hypot(dx, dy)
        
        # adjust self.active accordinly
        idx = np.argmin(distance_from_click)
        if event.button == 1:
            self.active[idx] = True
        elif event.button == 3:
            self.active[idx] = False
        self.logger.info('MESAnalyze plot data point with index %s set to status %s' % (idx, self.active[idx]))
        self.update_plots()
        self.display_values()

    def update_plots(self):
        """
        Calculates self.transformation, graph data on all plots, and display it
        @returns:
            The x and y residuals in numpy array form
        """
        data = self.data[np.nonzero(self.active)]

        # calculate the optimal transformation from the input data
        # Algorithm is Kabsch algorithm for calculating the optimal
        # rotation matrix that minimizes the RMSD (root mean squared
        # deviation) between two paired sets of points. This code
        # replaces the IRAF geomap task for the case when
        # fitgeometry=rotate.
        centroid = np.mean(data, axis=0)
        data = data - centroid
        p_i = np.asmatrix(data[:, 0:2])
        p_f = np.asmatrix(data[:, 2:4])
        u, s, v = np.linalg.svd(p_i.T * p_f)
        rot_mat = u * v
        shift =  centroid[2:4] - centroid[0:2]*rot_mat
        try:
            theta = np.mean([math.acos(rot_mat[0,0]), math.asin(rot_mat[1,0])])
        except ValueError:
            theta = 0
        self.transformation = (shift[0,0], shift[0,1], theta)
        
        # use its results to calculate some stuff
        xref = self.data[:, 0]
        yref = self.data[:, 1]
        xin  = self.data[:, 2]
        yin  = self.data[:, 3]
        xcalc, ycalc = self.transform(xref, yref, self.transformation)
        xres = xin - xcalc
        yres = yin - ycalc
        
        # graph residual data on the plots
        self.plot_residual(self.plots[0], xref, xres, self.active, var_name="X")
        self.plot_residual(self.plots[1], yref, yres, self.active, var_name="Y")
        
        # update the vectors on the canvas, as well
        for i in range(0, self.data.shape[0]):
            self.draw_vector_on_canvas(xref[i], yref[i], xres[i], yres[i], i)
        
        return xres, yres

    def draw_vector_on_canvas(self, xref, yref, xres, yres, idx):
        """
        Draws the residual vector at the given point on the canvas
        @param xref:
            The float observed x coordinate of this object
        @param yref:
            The float observed y coordinate of this object
        @param xres:
            The float residual x coordinate of this object
        @param yres:
            The float residual y coordinate of this object
        @param idx:
            The index of this object
        """
        # first calculate the endpoints of the vector
        startX = xref
        startY = yref
        if not math.isnan(xres) and not math.isnan(yres):
            endX = startX + 100*xres
            endY = startY + 100*yres
        else:
            endX = startX
            endY = startY
        magnitude = math.hypot(xres, yres)
        
        # determine the color based on activity and magnitude
        if not self.active[idx]:
            color = 'grey'
        elif magnitude <= 0.5:
            color = 'green'
        elif magnitude <= 1.0:
            color = 'yellow'
        else:
            color = 'red'
        
        # delete the old vector
        self.canvas.delete_object_by_tag(str(idx))
        
        # and draw in the new one
        self.canvas.add(self.dc.CompoundObject(
                                self.dc.Line(startX, startY, endX, endY,
                                             color=color, arrow='end',
                                             showcap=True, alpha = 0.7),
                                self.dc.Text((startX+endX)/2, (startY+endY)/2,
                                             "{:.1f}p".format(magnitude),
                                             color=color)),
                        tag=str(idx))

    def display_values(self):
        """
        Shows the final MES Offset values on the screen, based on
        self.transformation
        """
        # collect values from other sources
        xcenter = 1024.0
        ycenter = 1750.0
        xshift, yshift, thetaR = self.transformation
        thetaD = math.degrees(thetaR)
        
        # calculate dx and dy (no idea what all this math is)
        dx = -yshift + xcenter*math.sin(thetaR) + ycenter*(1-math.cos(thetaR))
        dy = xshift + xcenter*(math.cos(thetaR)-1) + ycenter*math.sin(thetaR)
        # normalize thetaD to the range [-180, 180)
        thetaD = (thetaD+180)%360 - 180
        
        # ignore values with small absolute values
        if abs(dx) < 0.5:
            dx = 0
        if abs(dy) < 0.5:
            dy = 0
        if abs(thetaD) < 0.01:
            thetaD = 0
        
        # then display all values
        self.final_displays["dX"].set_text("{:,.2f}".format(dx))
        self.final_displays["dY"].set_text("{:,.2f}".format(dy))
        self.final_displays["dPA"].set_text(u"{:,.4f}".format(thetaD))
        self.offset = (dx, dy, thetaD)


    def read_sbr_file(self, filename, logger):
        """
        Read the file and return the data within, structured as the position
        of the first star as well as the positions of the other stars
        @param filename:
            The name of the file which contains the data
        @returns:
            A numpy array of two columns, containing the first two data in each
            row of the sbr file
        """
        # open the file
        try:
            sbr = open(filename, 'r')
        except IOError:
            logger.warn('Warning: sbr file %s not found' % filename)
            return np.zeros((1,2))
        
        # declare some variables
        array = []
        line = sbr.readline()
        
        # read and convert the first two variables from each line
        while line != "":
            vals = line.split(', ')
            if vals[0] == 'C':
                sbrX, sbrY = float(vals[1]), float(vals[2])
                array.append(self.imgXY_from_sbrXY(sbrX, sbrY))
            line = sbr.readline()
        
        sbr.close()
        return np.array(array)

    def imgXY_from_sbrXY(self, sbrX, sbrY):
        """
        Converts coordinates from the SBR file to coordinates
        compatible with FITS files
        @param sbr_coords:
            A string or float tuple of x and y read from *.sbr
        @returns:
            A float tuple of x amd u that will work with the rest of Ginga
        """
        # I'm sorry; I have no idea what any of this math is.
        fX = 1078.0 - float(sbrX)*17.57789
        fY = 1857.0 + float(sbrY)*17.57789
        fHoleX = 365.0 + (fX-300.0)
        fHoleY = 2580.0 + (fY-2660.0)
        return (fHoleX, fHoleY)

    def parse_data(self,data):
        """
        Reads the data and returns it in a more useful form
        @param data:
            A numpy array of two or three columns, representing x, y[, and r]
        @returns:
            A numpy array of three columns of floats representing relative locations
            and radii of objects,
            and a single float tuple (absolute location of first object)
        """
        obj_list = []
        
        for row in data:
            # for each line, get the important values and save them in obj_list
            x, y = row[:2]
            if len(row) >= 3:
                r = row[2]
            else:
                r = float('NaN')
            obj_list.append([x, y, r])
        
        # convert obj_list to ndarray, and extract obj0
        obj_list = np.array(obj_list)
        obj0 = (obj_list[0,0], obj_list[0,1])
        
        # Calculate the shift to the best match
        print('Default hole location---\n',obj_list)
        print('Star catalog ----\n',self.moircsAlignImage.starMat)
        min_d=0
        for p in obj_list:
            dist=np.sqrt((p[0]-self.moircsAlignImage.starMat[:,0])**2+
                (p[1]-self.moircsAlignImage.starMat[:,1])**2)
            ind=np.where(dist == np.min(dist))
            dx=np.min(p[0]-self.moircsAlignImage.starMat[ind,0])
            dy=np.min(p[1]-self.moircsAlignImage.starMat[ind,1])
            print(np.min(dist),dx,dy)
            
            # Looking the the offset, remember all the minmum distance should be very close
            #  and the difference should be less than 15 pixels.  Another situation is that
            #  the offset larger than 500 pixels.  This indicated the match goes to wrong 
            #  star and because there are no stars near the hole.
            if min_d == 0 or 0 < min_d - np.min(dist) < 15 or min_d - np.min(dist) > 500:
                min_d = np.min(dist)
                xoffset = dx
                yoffset = dy    
            print(min_d,xoffset,yoffset)        
                
        obj_list[:,0:2] -= obj0
        return obj_list, obj0, xoffset, yoffset

        
    def startLocate(self, initial_data, mode, interact2=True,
              next_step=None):
        """
        Get the positions of a series of objects
        @param initial_data:
            The numpy array containing the approximate positions of the relevant
            objects
        @param mode:
            Either 'star' or 'mask' or 'starhole'; alters the sizes of squares
            and the autocut method
        @param interact2:
            Whether we should give the user a chance to interact with step 2
        @param next_step:
            A function to call when MESLocate is finished
        """
        # read the data
        self.obj_arr, obj0, dx ,dy = self.parse_data(initial_data)
        self.obj_num = self.obj_arr.shape[0]
        
        
        if mode == 'star':
            # Calculate the shift to the best match
            #dist=np.sqrt((obj0[0]-self.moircsAlignImage.starMat[:,0])**2+
            #    (obj0[1]-self.moircsAlignImage.starMat[:,1])**2)
            
            #ind=np.where(dist == np.min(dist))
            
            #dx=np.min(obj0[0]-self.moircsAlignImage.starMat[ind,0])
            #dy=np.min(obj0[1]-self.moircsAlignImage.starMat[ind,1])
            #print(self.obj_arr)
            print(obj0,dx,dy)
            obj0 = (obj0[0]-float(dx),obj0[1]-float(dy))
            #obj0[0] = obj0[0]-float(dx)
            #obj0[1] = obj0[1]-float(dy)



        # define some attributes
        self.click_history = []     # places we've clicked
        self.click_index = -1       # index of the last click
        self.color_index = -1       # index of the current color
        self.current_obj = 0        # index of the current object
        self.drag_history = [[]]    # places we've click-dragged*
        self.drag_index = [-1]      # index of the current drag for each object
        self.drag_start = None      # the place where we most recently began to drag
        self.obj_centroids = np.zeros(self.obj_arr.shape)    # the new obj_arr based on user input and calculations
        self.square_size =  {'star':30, 'mask':60, 'starhole':20}[mode]  # the apothem of the search regions
        self.exp_obj_size = {'star':4,  'mask':20, 'starhole':4}[mode]  # the maximum expected radius of the objects
        self.interact = interact2    # whether we should interact in step 2
        self.next_step = next_step  # what to do when we're done
        
        # set some values based on mode
        self.mode = mode
        if mode == 'star':
            autocut_method = 'stddev'   # fitsimage autocut method
            self.exp_obj_size = 4       # the maximum expected radius of objects
            self.square_size = 30       # the apothem of the search regions
        elif mode == 'mask':
            autocut_method = 'minmax'
            self.exp_obj_size = 20
            self.square_size = 60
        elif mode == 'starhole':
            autocut_method = 'stddev'
            self.exp_obj_size = 4
            # Use np.nanmax to ignore any objects with radius set to
            # NaN (those are objects that were skipped over in the
            # mask image step)
            self.square_size = np.nanmax(self.obj_arr[:,2])
        
        # creates the list of thumbnails that will go in the GUI
        self.fitsimage.get_settings().set(autocut_method=autocut_method)
        self.thumbnails = self.create_viewer_list(self.obj_num, self.logger)
        self.viewer_grid.remove_all()
        for row in range(int(math.ceil(self.obj_num/2.0))):
            for col in range(2):
                i = 2*row + col
                if i < len(self.thumbnails):
                    pic = Viewers.GingaViewerWidget(viewer=self.thumbnails[i])
                    self.viewer_grid.add_widget(pic, row, col)
        
        # set the mouse controls and automatically start if this is starhole mode
        self.setLocatecallbacks()
        self.click1_cb(self.canvas, 1, *obj0)
        self.go_to_gui('find')
        if mode == 'starhole':
            self.step2_cb()

    def create_viewer_list(self, n, logger=None, width=120, height=120):    
        # 147x147 is approximately the size it will set it to, but I have to 
        #  set it manually because otherwise it will either be too small or scale 
        #  to the wrong size at certain points in the program. Why does it do this? 
        #  Why can't it seem to figure out how big the window actually is when it zooms? #
        #  I don't have a clue! It just randomly decides sometime after my plugin's last 
        #  init method and before its first callback method, hey, guess what, the window 
        #  is 194x111 now - should I zoom_fit again to match the new size? 
        #  Nah, that would be TOO EASY. And of course I don't even know where or when or 
        #  why the widget size is changing because it DOESN'T EVEN HAPPEN IN GINGA! 
        #  It happens in PyQt4 or PyQt 5 or, who knows, maybe even Pyside. Obviously. OBVIOUSLY. GAGFAGLAHOIFHAOWHOUHOUH~~!!!!!
        """
        Create a list of n viewers with certain properties
        @param n:
            An integer - the length of the desired list
        @param width:
            The desired width of each viewer
        @param height:
            The desired height of each veiwer
        @param logger:
            A Logger object to pass into the new Viewers
        @returns:
            A list of Viewers.CanvasView objects
        """
        output = []
        for i in range(n):
            viewer = Viewers.CanvasView(logger=logger)
            viewer.set_desired_size(width, height)
            viewer.enable_autozoom('on')
            viewer.enable_autocuts('on')
            output.append(viewer)
        return output

    def setLocatecallbacks(self, step=1, selection_mode=0, clear=True):
        """
        Assign all necessary callbacks to the canvas for the current step
        @param step:
            The number of this step - 1 for finding and 2 for centroid-getting
        @param selection_mode:
            0 for 'Automatic', 1 for 'Crop', or 2 for 'Mask'
        @param clear:
            Whether the canvas should be completely cleared
        """
        canvas = self.canvas
        
        # clear all existing callbacks first
        self.clear_canvas(keep_objects=not clear, keep_zoom=not clear)
        
        # for step one, the only callbacks are for right-click and left-click
        if step == 1:
            canvas.add_callback('cursor-up', self.click1_cb)
            canvas.add_callback('draw-up', self.step2_cb)
        
        # for step two, you need callbacks for left-drag and middle-drag, too
        elif step == 2:
            if selection_mode == 2:
                canvas.add_callback('cursor-down', self.start_drag_cb, 'mask')
                canvas.add_callback('cursor-up', self.end_drag_cb, 'mask')
            else:
                canvas.add_callback('cursor-down', self.start_drag_cb, 'crop')
                canvas.add_callback('cursor-up', self.end_drag_cb, 'crop')
            if selection_mode == 1:
                canvas.add_callback('panset-down', self.start_drag_cb, 'crop')
                canvas.add_callback('panset-up', self.end_drag_cb, 'crop')
            else:
                canvas.add_callback('panset-down', self.start_drag_cb, 'mask')
                canvas.add_callback('panset-up', self.end_drag_cb, 'mask')
            canvas.add_callback('draw-up', self.next_obj_cb)

    def click1_cb(self, _, __, x, y):
        """
        Respond to a left click on the screen in step 1
        @param x:
            The x coordinate of the click
        @param y:
            The y coordiante of the click
        """
        # increment the index
        self.click_index += 1
        self.color_index += 1
        
        # if there are things saved ahead of this index (because of undo), clear them
        if self.click_index < len(self.click_history):
            self.click_history = self.click_history[:self.click_index]
        self.click_history.append((x,y))
        self.select_point(self.click_history[self.click_index], draw_circle_masks=True)
        
        # also clear step2 variables
        self.current_obj = 0
        self.drag_history = [[]]*self.obj_num
        self.drag_index = [-1]*self.obj_num
        return False



class MoircsAlignImage(object):
    
    
    # Keyword for setting re-calculate mosaic map
    MAKEMOSAIC = True

    # Establishing list for needed files
    star_fits_name=[]
    sky_fits_name=[]
    mask_fits_name=[]
    starhole_fits_name=[]
    badpix_fits_name=[]

    dbs_name=[]

    output_path=''
    rootname=''

    star = 0
    badpix = 0
    mask = 0

    starCat = 0
    starMat = 0
    
    def nothing(*args, **kwargs):
        """A placeholder function for log"""
        pass

    def processStarHoleImage(self, terminate, log=nothing, next_step=None):
        self.loadImage()
        
        log("Processing star-hole frames...")
        # Check the altitude of the image.
        if self.star.raw1.header['ALTITUDE'] < 45:
            log(LOW_ELEV_WARN.format("Image file: ",self.star_fits_name[0],
                            44.1), level='warning')
        
        if self.star.raw2.header['ALTITUDE'] < 45:
            log(LOW_ELEV_WARN.format("Image file: ",self.star_fits_name[1],
                            44.1), level='warning')
        if self.MAKEMOSAIC is True:
            if terminate.is_set():  return
            self.removeBackground()
        
        
            log("Making distortion correction...")
            if terminate.is_set():  return
            self.transformImage('remap1','detrend1',self.dcc1)
        
            if terminate.is_set():  return
            self.transformImage('remap2','detrend2',self.dcc2)
        
            log("Making mosaic image...")
            if terminate.is_set():  return
            self.mosaicField()

            # apply gaussian blur
            log("Blurring image...")
            if terminate.is_set():  return
            self.star.mosaic.data = gaussian_filter(self.star.mosaic.data, 1.0)


            out_filename = os.path.join(self.output_path, self.rootname+"_starhole.fits")

            log("Writing image file...")
            #fits.writeto(out)
            self.star.mosaic.writeto(out_filename,overwrite=True)

            out_filename = os.path.join(self.output_path,
                    'star_MCSA%08d_mosaic.fits'% self.star_fits_name[0])
            self.star.mosaic.writeto(out_filename,overwrite=True)

        else:
            log("Reading file from disk")    
            filename = os.path.join(self.output_path,
                    'star_MCSA%08d_mosaic.fits'% self.star_fits_name[0])
            self.star.mosaic=fits.open(filename,memmap=False)[0]
            
            log("Swap the byte order...") 
            self.star.mosaic.data=self.star.mosaic.data.byteswap().newbyteorder()    

        log("Imgage operation finished...")
        if terminate.is_set():  return

        if next_step != None:
            next_step()



    def processStarImage(self, terminate, log=nothing, next_step=None):
        self.loadImage()

        log("Processing star frames...")
        # Check the altitude of the image.
        if self.star.raw1.header['ALTITUDE'] < 45:
            log(LOW_ELEV_WARN.format("Image file: ",self.star_fits_name[0],
                            44.1), level='warning')
        
        if self.star.raw2.header['ALTITUDE'] < 45:
            log(LOW_ELEV_WARN.format("Image file: ",self.star_fits_name[1],
                            44.1), level='warning')

        log("Looking for mask holes...")
        if terminate.is_set():  return
        self.findGuidingHole('ch1',self.dcc1)
        self.findGuidingHole('ch2',self.dcc2)

        log("Combine mask locations...")
        if terminate.is_set():  return
        self.combineHoleLacation()
        
        if self.MAKEMOSAIC is True:

            log("Removing sky backgound...")   
            if terminate.is_set():  return
            self.removeBackground()
            
            log("Making distortion correction...")
            if terminate.is_set():  return
            self.transformImage('remap1','detrend1',self.dcc1)
            
            if terminate.is_set():  return
            self.transformImage('remap2','detrend2',self.dcc2)
            
            log("Making mosaic image...")
            if terminate.is_set():  return
            self.mosaicField()
            
            
            # apply gaussian blur
            log("Blurring image...")
            if terminate.is_set():  return
            self.star.mosaic.data = gaussian_filter(self.star.mosaic.data, 1.0)

            out_filename = os.path.join(self.output_path, self.rootname+"_star.fits")

            print(out_filename)
            log("Writing image file...")
            #fits.writeto(out)
            #self.star.mosaic.writeto(out_filename,overwrite=True)

            out_filename = os.path.join(self.output_path,
                    'star_MCSA%08d_mosaic.fits'% self.star_fits_name[0])
            #self.star.mosaic.writeto(out_filename,overwrite=True)

        else:
            log("Reading file from disk")    
            filename = os.path.join(self.output_path,
                    'star_MCSA%08d_mosaic.fits'% self.star_fits_name[0])
            self.star.mosaic=fits.open(filename,memmap=False)[0]
            
            log("Swap the byte order...") 
            self.star.mosaic.data=self.star.mosaic.data.byteswap().newbyteorder()


        log("Extracting stars on image...")
        if terminate.is_set():  return
        self.extractStar()

        log("Looking for best stars for hole locations...")
        if terminate.is_set():  return
        self.matchingStar()



        log("Imgage operation finished...")
        if terminate.is_set():  return

        if next_step != None:
            next_step()

    def findGuidingHole(self, maskKey, dc):
        img=np.copy(self.mask[maskKey].data)
        inx=np.where(img > 150)
        ind=np.where(img < 0)

        img[inx]=150
        img[ind]=0

        img=img.astype('uint8')
        img = cv2.medianBlur(img, 5)

        #np.clip(img,0,256)
        
        circles = cv2.HoughCircles(img,\
              cv2.HOUGH_GRADIENT,1,100,param1=30,param2=30,minRadius=11,maxRadius=30)
        
        
        if maskKey == 'ch1':
            self.mask.ch1_gHoles=circles
        if maskKey == 'ch2':
            self.mask.ch2_gHoles=circles
        
        cor_circles=circles
        
        cor_circles[0,:,0],cor_circles[0,:,1]=\
            self.reverseTransformLocation([circles[0,:,0],circles[0,:,1]],dc)
        
        if maskKey == 'ch1':
            self.mask.ch1_gHolesCorr=cor_circles
        if maskKey == 'ch2':
            self.mask.ch2_gHolesCorr=cor_circles

        del(img)

    def combineHoleLacation(self):

        '''Remap the ch2 XY location to mosaic image plan'''
        x,y=(self.mask.ch2_gHolesCorr[0,:,0],self.mask.ch2_gHolesCorr[0,:,1])
        
        '''Loading mosaic parameter'''
        mc2=self.mc2
        (xfit,yfit)=np.linalg.solve(np.array([[self.mc2.b,self.mc2.c],
                [self.mc2.e,self.mc2.f]]),np.array([x-self.mc2.a,y-self.mc2.b]))

        
        (self.mask.ch2_gHolesCorr[0,:,0],self.mask.ch2_gHolesCorr[0,:,1])=(xfit,yfit)
        
        
        self.mask.gHoleMosaic=np.zeros((1,self.mask.ch1_gHolesCorr.shape[1]+self.mask.ch2_gHolesCorr.shape[1],3))
        
        
        
        self.mask.gHoleMosaic[0,:,0]=np.append(np.subtract(2046,self.mask.ch1_gHolesCorr[0,:,1]),\
                                                    np.subtract(2086.5,self.mask.ch2_gHolesCorr[0,:,1]))
        self.mask.gHoleMosaic[0,:,1]=np.subtract(np.append(self.mask.ch1_gHolesCorr[0,:,0],\
                                                    self.mask.ch2_gHolesCorr[0,:,0]),33.5)
        self.mask.gHoleMosaic[0,:,2]=np.append(self.mask.ch1_gHolesCorr[0,:,2],
                                                    self.mask.ch2_gHolesCorr[0,:,2])
        
        self.mask.gHoleMosaic=self.mask.gHoleMosaic


    def transformImage(self, remapKey, detrendKey, dc):
        
        #img=self.star[detrendKey].data[:,:]
        
        y,x=np.mgrid[0:self.star[detrendKey].data.shape[0],0:self.star[detrendKey].data.shape[1]]
        indx,indy=self.transformLocation([x,y],dc)
        
        #remap=self.star[detrendKey].data[indy,indx]
        
        self.star[remapKey]=fits.PrimaryHDU(data=self.star[detrendKey].data[indy,indx])

    def reverseTransformLocation(self,uncorrectXY,dc):
        
        y,x=np.mgrid[0:2048,0:2048]
        xx,yy=(self.transformLocation((x,y),dc))
        xd,yd=uncorrectXY

        xd=np.uint16(np.around(xd))
        yd=np.uint16(np.around(yd))
        

        indx=np.copy(xd)
        indy=np.copy(yd)
        
        
        #print(len(xd))
        for i in range(0,len(xd)):
            xn=np.unique(x[np.where(xx==xd[i])])
            yn=np.unique(y[np.where(yy==yd[i])])
            for j in range(0,len(xn)):
                for k in range(0,len(yn)):
                    xf,yf=self.transformLocation(np.array([xn[j],yn[k]]),dc)
                    if xf == xd[i] and yf == yd[i]:
                        indx[i]=xn[j]
                        indy[i]=yn[k]
                        
        return indx,indy


    def transformLocation(self,uncorrectXY,dc):
        
        x,y=uncorrectXY
        
        xfit1=dc.a+dc.b*x+dc.c*y
        yfit1=dc.d+dc.e*x+dc.f*y

        xfit2=dc.xc00+(dc.xc10*x)+(dc.xc01*y)+(dc.xc20*x**2)+(dc.xc30*x**3)+(dc.xc11*x*y)+(dc.xc21*x**2*y)\
                +(dc.xc02*y**2)+(dc.xc12*x*y**2)+(dc.xc03*y**3)
            
        yfit2=dc.yc00+(dc.yc10*x)+(dc.yc20*x**2)+(dc.yc30*x**3)+(dc.yc01*y)+(dc.yc11*x*y)+(dc.yc21*x**2*y)\
                +(dc.yc02*y**2)+(dc.yc12*x*y**2)+(dc.yc03*y**3)

        xfit=xfit1+xfit2
        yfit=yfit1+yfit2

        if type(xfit) is np.ndarray:
            xfit[np.where(xfit < 0)]=0
            xfit[np.where(xfit > 2047)]=2047
        else:
           np.clip(xfit,0,2047)
        
        if type(yfit) is np.ndarray:
            yfit[np.where(yfit < 0)]=0
            yfit[np.where(yfit > 2047)]=2047
        else:
            np.clip(yfit,0,2047)

                
        indx=np.uint(np.round(xfit))
        indy=np.uint(np.round(yfit))

        return indx,indy

    def mosaicField(self):
        ''''''
                
        ''' Prepare a 2d array '''
        temp1=np.zeros((2048,3636))
        
        ''' Shift the channel 1 image '''
        temp1[:,0:2048]=self.star.remap1.data[:,:]

        '''Operate on channel2'''
        temp2=np.zeros((2048,3636))
        temp2[:,0:2048]=self.star.remap2.data
        
        y,x=np.mgrid[0:temp2.shape[0],0:temp2.shape[1]]
        
        '''Loading mosaic parameter'''
        #mc2=MosaicParameterCh2()
        mc2=self.mc2
        xfit=mc2.a+mc2.b*x+mc2.c*y
        yfit=mc2.d+mc2.e*x+mc2.f*y
        
        ''' Arranging indexes in correct range'''
        xfit[np.where(xfit < 0)]=0
        yfit[np.where(yfit < 0)]=0
        xfit[np.where(xfit > 2047)]=2047
        yfit[np.where(yfit > 2047)]=2047
        
        ''' Round the float numbers to integer'''
        indx=np.uint(np.round(xfit))
        indy=np.uint(np.round(yfit))
        
        remap=temp2[indy,indx]
        
        img=temp1*np.subtract(1,self.badpix.ch1.data)+remap*np.subtract(1,self.badpix.ch2.data)
        
        avg=np.median(img)
        img=np.add(0.5*avg,0.5*img)
        
        self.star.mosaic=fits.PrimaryHDU(data=shift(np.rot90(img, k=3),[33.5,0],order=0)[67:,:],
            header=self.star.raw1.header)


    def removeBackground(self):
        
        self.star.detrend1=fits.PrimaryHDU(data=self.star.raw1.data-self.star.bg1.data)
        self.star.detrend2=fits.PrimaryHDU(data=self.star.raw2.data-self.star.bg2.data)


    def loadImage(self):

        self.star.raw1 = fits.open(self.imagepath+'MCSA%08d.fits'% self.star_fits_name[0],memmap=False)[0]
        self.star.raw2 = fits.open(self.imagepath+'MCSA%08d.fits'% self.star_fits_name[1],memmap=False)[0]

        self.star.bg1 = fits.open(self.imagepath+'MCSA%08d.fits'% self.sky_fits_name[0],memmap=False)[0]
        self.star.bg2 = fits.open(self.imagepath+'MCSA%08d.fits'% self.sky_fits_name[1],memmap=False)[0]
        
        self.mask.ch1=fits.open(self.imagepath+'MCSA%08d.fits'% self.mask_fits_name[0],memmap=False)[0]
        self.mask.ch2=fits.open(self.imagepath+'MCSA%08d.fits'% self.mask_fits_name[1],memmap=False)[0]
        
        self.badpix.ch1=fits.open(self.badpix_fits_name[0]+'.fits',memmap=False)[0]
        self.badpix.ch2=fits.open(self.badpix_fits_name[1]+'.fits',memmap=False)[0]

        self.loadDistorCoeff()

    def loadDistorCoeff(self):

        try:        
            #self.dcc1 = MoircsAlignConfig().dcc1 #DistortionMap(self.dbs_name[0])
            self.dcc1 = DistortionMap(self.dbs_name[0])
            self.dcc2 = DistortionMap(self.dbs_name[1])
            self.mc1 = MosaicPrarameter(self.dbs_name[2])
            self.mc2 = MosaicPrarameter(self.dbs_name[3])

        except IOError as err:
            raise IOError
        
        attrs = vars(self.mc2)
        print(', '.join("%s: %s" % item for item in attrs.items()))
        #print(self.dcc2.a,self.dcc2.xc00)

    def extractStar(self):
        
        data=self.star.mosaic.data
        m, s = np.mean(data), np.std(data)
        bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
        objs = sep.extract(data-bkg, 2.5, err=bkg.globalrms,minarea=20)


        aper_radius=8.0

        # Calculate the Kron Radius
        kronrad, krflag = sep.kron_radius(data, objs['x'], objs['y'], \
            objs['a'], objs['b'], objs['theta'], aper_radius)
        
        r_min = 4
        use_circle = kronrad * np.sqrt(objs['a'] * objs['b'])
        cinx=np.where(use_circle <= r_min)
        einx=np.where(use_circle > r_min)

        # Calculate the equivalent of FLUX_AUTO
        flux, fluxerr, flag = sep.sum_ellipse(data, objs['x'][einx], objs['y'][einx], \
            objs['a'][einx], objs['b'][einx], objs['theta'][einx], 2.5*kronrad[einx],subpix=1)      
        
        cflux, cfluxerr, cflag = sep.sum_circle(data, objs['x'][cinx], objs['y'][cinx],
                                        objs['a'][cinx], subpix=1)

        objs['flux'][einx]=flux
        objs['flux'][cinx]=cflux


        r, flag = sep.flux_radius(data, objs['x'], objs['y'], \
            aper_radius*objs['a'], 0.5,normflux=objs['flux'], subpix=5)
        
        flag |= krflag
        
        #objs.append({'r':r})
        objs['flag'][einx]=flag
        objs['flag'][cinx]=cflag

        #eliminate unwanted sources
        index=np.logical_and(objs['b']/objs['a'] > 0.1,
                        np.logical_and(flag >= 0,np.logical_or(\
            np.logical_and(objs['y'] > 120, objs['y'] < 1800),\
            np.logical_and(objs['y'] > 1860, objs['y'] < 3530))))

    
        objects=objs[:][np.where(index == True)]
        
        self.starCat=objects[:]
        

        # fig, ax = plt.subplots()
        # im = ax.imshow(np.flipud(data), interpolation='nearest', cmap='gray',
        #         vmin=0, vmax=200, origin='lower')

        # # plot an ellipse for each object
        # for i in range(len(objects)):
        #     ##print (objects['x'][i], objects['y'][i],objects['npix'][i])
        #     ##if objects['b'][i]/objects['a'][i] > 0.5:
        #     e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
        #                 width=6*objects['a'][i],
        #                 height=6*objects['b'][i],
        #                 angle=objects['theta'][i] * 180. / np.pi)
        #     e.set_facecolor('none')
        #     e.set_edgecolor('red')
        #     ax.add_artist(e)
        
        # for x in self.mask.gHoleMosaic[0,:]:
        #     ##draw the outer circle
        #     r=plt.Circle((x[0],x[1]),x[2], \
        #             edgecolor='yellow', facecolor='none')
        #     ax.add_artist(r)

        
        # plt.show()
        # plt.close()
        
    def matchingStar(self):        
        
        
        # The limit of N nearest stars around the holes.
        nlimit=5
       
        darray=np.array([])
        for i in self.mask.gHoleMosaic[0,:]:
            d = np.sqrt( (i[0] - np.array([self.starCat['x']]))**2 + 
                (i[1] - np.array([self.starCat['y']]))**2 )
            if darray.shape[0] == 0:
                darray=d
            else:
                darray=np.vstack((darray,d))
                
        
        # Extracting the nearest N stars for searching
        nearest_dis=np.array([])
        nearest_ind=np.array([])
        for i in range(len(self.mask.gHoleMosaic[0,:,0])):
            if nearest_dis.shape[0] == 0:
                nearest_dis = darray[i,darray[i,0:].argsort()[0:nlimit]]
            else:
                nearest_dis=np.vstack((nearest_dis,darray[i,darray[i,0:].argsort()[0:nlimit]]))
            
            if nearest_ind.shape[0] == 0:
                nearest_ind = darray[i,0:].argsort()[0:nlimit]
            else:
                nearest_ind=np.vstack((nearest_ind,darray[i,0:].argsort()[0:nlimit]))
        
        print(nearest_dis)
        print(nearest_ind)
            
        # It will be better if we select this value from a data set that  
        #  stars are close to the hole
        si,sj=np.where(nearest_dis == np.min(nearest_dis))
    
        
        # Now looking for matched stars by checking distances, go through 
        #  all N-th nearest stars and find their friends.
        match_set=Bunch()
        weight_set=Bunch()
        for i in range(nlimit):
            
            # Using this distance for friend-looking
            dmin_template=nearest_dis[si,i]
            
            #print('templat from ',nearest_dis[si,i])
            
            # subtract this distance 
            matrix_dis=np.abs(np.subtract(nearest_dis,dmin_template))
            
            # Initiate two arrays for the storage of the pari indexes and 
            #  distances
            sub_match_set=np.array([])
            sub_distance=np.array([])
            
            # Pick the star with minimum distance to a hole from each row
            for x in range(matrix_dis.shape[0]):
                temp_dis = np.array(matrix_dis[x:x+1,:])
                temp_ind = np.array(nearest_ind[x:x+1,:])
                
                # Set the selection condition for stars
                index=np.logical_and(temp_dis < 20, \
                    temp_dis == np.min(temp_dis))
                
                #print(temp_dis)
                #print(temp_ind[np.where(index == True)])
                if len(temp_ind[np.where(index == True)]) == 0:
                    sub_match_set=np.append(sub_match_set,np.nan)
                    sub_distance=np.append(sub_distance,np.nan)
                else:    
                    sub_match_set=np.append(sub_match_set,\
                        temp_ind[np.where(index == True)])
                    sub_distance=np.append(sub_distance,\
                        temp_dis[np.where(index == True)])
                
                #print(sub_match_set)
                #print('-------')
            
            # When storing the match set, remove all NaN elements.  The NaN are 
            #  used to calculate the weighting, so that we do not want to change
            #  it in original array.  The weighting is based on 1) the standard 
            #  deviation, 2) the number of associated stars and 3) the distances.
            match_set['set'+str(i)]=np.uint16(
                sub_match_set[np.where(~np.isnan(sub_match_set) == True)]
                )
            weight_set['set'+str(i)]=\
                np.std(sub_distance[np.where(~np.isnan(sub_distance) == True)])*\
                float(dmin_template[0])*\
                (np.count_nonzero(np.isnan(sub_distance))+1)**3
            print(np.std(sub_distance[np.where(~np.isnan(sub_distance) == True)]))
            print(np.mean(sub_distance[np.where(~np.isnan(sub_distance) == True)]))
            print((np.count_nonzero(np.isnan(sub_distance))+1))
            print(dmin_template)
            
        # When there is only one star detected in the hole, the STD=0 and 
        #  the rest weighting will not working.  This is the fix
        weight_min=0
        for i in weight_set.items():
            if weight_min == 0:
                weight_min = i[1]
                stars=match_set[i[0]]
            elif (i[1] < weight_min and i[1] != 0):
                weight_min = i[1]
                stars=match_set[i[0]]
        
                
        print(match_set)
        print(weight_set)
        print(stars)
        #stars=match_set[min(weight_set,key=weight_set.get)]
        
        #print(self.starCat['x'][stars], self.starCat['y'][stars])
        for i in range(0,len(self.starCat['x'][stars])):
            print([self.starCat['x'][stars[i]],self.starCat['y'][stars[i]]])
            if i == 0:
                self.starMat= np.array([self.starCat['x'][stars[i]],self.starCat['y'][stars[i]]])
            else:
                self.starMat=np.vstack((self.starMat,[self.starCat['x'][stars[i]],self.starCat['y'][stars[i]]]))
       
        # plt.close('all')
        # fig, ax = plt.subplots()
        # plt.imshow(np.flipud(self.star.mosaic.data), \
        #     cmap=plt.get_cmap('gray'),vmin=0, vmax=2000, origin='lower')
        # ## plot an ellipse for each object
        # for i in range(len(stars)):
        #     #if self.scat['flag'][i] == 99:
        #     if np.isnan(stars[i]) != True:
        #         e = Ellipse(xy=(self.starCat['x'][stars[i]], self.starCat['y'][stars[i]]),
        #                     width=6*self.starCat['a'][stars[i]],
        #                     height=6*self.starCat['b'][stars[i]],
        #                     angle=self.starCat['theta'][stars[i]] * 180. / np.pi)
        #         e.set_facecolor('none')
        #         e.set_edgecolor('red')
        #         ax.add_artist(e)

        # for x in self.mask.gHoleMosaic[0,:]:
        #     ###draw the outer circle
        #     r=plt.Circle((x[0],x[1]),x[2], \
        #             edgecolor='yellow', facecolor='none')
        #     ax.add_artist(r)

        # plt.show()
        # plt.close()


    def closeImage(self):
       
 
        if self.star != 0:
            del self.star
        
        if self.mask != 0:
            del self.mask
        
        if self.badpix !=0:
            del self.badpix
        

    def __init__(self):
        #super(MoircsAlignConfig, self).__init__()
        #self.imagepath=MoircsAlignConfig().imagepath
        #self.dcc1=MoircsAlignConfig().dcc1
        #self.dcc2=MoircsAlignConfig().dcc2
        #self.mc2=MoircsAlignConfig().mc2
        
        self.starCat=0
        self.starMat
        self.badpix=Bunch(dict(ch1=0,ch2=0))

        self.mask=Bunch(dict(ch1=0,ch2=0,mosaic=0,ch1_gHoles=0,
                    ch2_gHoles=0,ch1_gHolesCorr=0,ch2_gHolesCorr=0,gHoleMosaic=0))

        self.star=Bunch(dict(raw1=0,raw2=0,bg1=0,bg2=0,
            detrend1=0,detrend2=0,rempa1=0,remap2=0,mosaic=0))



    def __del__(self):
        self.closeImage()
        gc.collect()

class DistortionMap(object):
    
    filename=""
    
    def __init__ (self,filename):   
        self.filename=filename
        self.readDataFile()

    def readDataFile(self):    
        data=ascii.read(self.filename)
        
        field0 = []
        for x in data['begin'].data :
            try:
                field0.append(float(x))
            except ValueError:
                field0.append(x)
        field1 = []

        for x in data.field(1).data :
            try:
                field1.append(float(x))
            except ValueError:
                field1.append(x)

        np.sum(data['begin'] == 'begin')
        if np.sum(data['begin'] == 'begin') != 0:
            ind=np.array(np.where((data['begin'] == 'begin') == True))
        else:
            ind = np.array([-1])

        self.a=field0[int(ind[-1])+24]
        self.b=field0[int(ind[-1])+25]
        self.c=field0[int(ind[-1])+26]
        self.d=field1[int(ind[-1])+24]
        self.e=field1[int(ind[-1])+25]
        self.f=field1[int(ind[-1])+26]

        self.xc00=field0[int(ind[-1])+36]
        self.xc10=field0[int(ind[-1])+37]
        self.xc20=field0[int(ind[-1])+38]
        self.xc30=field0[int(ind[-1])+39]
        self.xc01=field0[int(ind[-1])+40]
        self.xc11=field0[int(ind[-1])+41]
        self.xc21=field0[int(ind[-1])+42]
        self.xc02=field0[int(ind[-1])+43]
        self.xc12=field0[int(ind[-1])+44]
        self.xc03=field0[int(ind[-1])+45]

        self.yc00=field1[int(ind[-1])+36]
        self.yc10=field1[int(ind[-1])+37]
        self.yc20=field1[int(ind[-1])+38]
        self.yc30=field1[int(ind[-1])+39]
        self.yc01=field1[int(ind[-1])+40]
        self.yc11=field1[int(ind[-1])+41]
        self.yc21=field1[int(ind[-1])+42]
        self.yc02=field1[int(ind[-1])+43]
        self.yc12=field1[int(ind[-1])+44]
        self.yc03=field1[int(ind[-1])+45]

class MosaicPrarameter(object):
    
    filename=""
    
    def __init__ (self,filename):   
        self.filename=filename
        self.readDataFile()

    def readDataFile(self):    
        data=ascii.read(self.filename)
        
        field0 = []
        for x in data['begin'].data :
            try:
                field0.append(float(x))
            except ValueError:
                field0.append(x)
        field1 = []

        for x in data.field(1).data :
            try:
                field1.append(float(x))
            except ValueError:
                field1.append(x)

        np.sum(data['begin'] == 'begin')
        if np.sum(data['begin'] == 'begin') != 0:
            ind=np.array(np.where((data['begin'] == 'begin') == True))
        else:
            ind = np.array([-1])

        self.a=field0[int(ind[-1])+24]
        self.b=field0[int(ind[-1])+25]
        self.c=field0[int(ind[-1])+26]
        self.d=field1[int(ind[-1])+24]
        self.e=field1[int(ind[-1])+25]
        self.f=field1[int(ind[-1])+26]

        self.a=-1602.189
        self.b=0.9993759
        self.c=0.009615554
        self.d=40.30523
        self.e=-0.009611026
        self.f=0.9998468


# class DistortionCoeffCh1:
#     a=-8.635143
#     b=1.009153
#     c=2.394803E-4
#     d=-10.43176
#     e=2.192026E-4
#     f=1.00965

#     xc00=-16.37997
#     xc10=0.04216358
#     xc20=-3.961311E-5
#     xc30=1.453919E-8
#     xc01=0.02687643
#     xc11=-2.919222E-5
#     xc21=-6.896180E-11
#     xc02=-1.334916E-5
#     xc12=1.432681E-8
#     xc03=6.087694E-11

#     yc00=-16.98225
#     yc10=0.02702253
#     yc20=-1.542926E-5
#     yc30=2.064173E-10
#     yc01=0.04707864
#     yc11=-2.632263E-5
#     yc21=1.465591E-8
#     yc02=-4.374318E-5
#     yc12=-2.114975E-12
#     yc03=1.430626E-8


# class DistortionCoeffCh2:
#     a=-5.023962
#     b=1.004285
#     c=1.277659E-4
#     d=-5.48291
#     e=7.625483E-5
#     f=1.004831

#     xc00=-27.69665
#     xc10=0.06082649
#     xc20=-4.797236E-5
#     xc30=1.420412E-8
#     xc01=0.03267656
#     xc11=-2.884639E-5
#     xc21=-9.848598E-11
#     xc02=-1.614414E-5
#     xc12=1.410305E-8
#     xc03=9.066495E-11

#     yc00=-24.74531
#     yc10=0.03398627
#     yc20=-1.568090E-5
#     yc30=2.707538E-10
#     yc01=0.05459096
#     yc11=-3.218881E-5
#     yc21=1.433751E-8
#     yc02=-4.413517E-5
#     yc12=-4.351142E-11
#     yc03=1.435618E-8


# class MosaicParameterCh2:
#     a=-1602.189
#     b=0.9993759
#     c=0.009615554
#     d=40.30523
#     e=-0.009611026
#     f=0.9998468

