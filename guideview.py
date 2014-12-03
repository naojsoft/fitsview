#!/usr/bin/env python
#
# guideview.py -- Simple FITS viewer/display server.
#
# Eric Jeschke (eric@naoj.org)
#
"""
guideview.py implements a simple FITS viewer/display server to display FITS
images in GTK widgets.

Usage:
    guideview.py 
"""

# stdlib imports
import sys, os
import threading
import logging
import ssdlog
import traceback

moduleHome = os.path.split(sys.modules[__name__].__file__)[0]
sys.path.insert(0, moduleHome)
pluginHome = os.path.join(moduleHome, 'plugins')
sys.path.insert(0, pluginHome)

# Subaru python stdlib imports
import remoteObjects as ro
import remoteObjects.Monitor as Monitor
import Gen2.soundsink as SoundSink

# Ginga imports
from ginga.misc import ModuleManager, Datasrc, Settings
from ginga.misc.Bunch import Bunch
from ginga.Control import GingaControl, GuiLogHandler
from ginga import toolkit
toolkit.use('gtk2')
from ginga.gtkw.GingaGtk import GingaView

from ginga.util import wcsmod
#wcsmod.use('astropy')
wcsmod.use('kapteyn')

# Local application imports
from util import Receive

serviceName = 'guideview'
version = "20131213.0"

default_layout = ['seq', {},
                   ['vbox', dict(name='top', width=1800, height=1050),
                    dict(row=['hbox', dict(name='menu')],
                         stretch=0),
                    dict(row=['hpanel', {},
                     ['ws', dict(name='left', width=300),
                      # (tabname, layout), ...
                      [("Info", ['vpanel', {},
                                 ['ws', dict(name='uleft', height=300,
                                             show_tabs=False, group=3)],
                                 ['ws', dict(name='lleft', height=430,
                                             show_tabs=False, group=3)],
                                 ]
                        )]
                      ], 
                     ['vpanel', {},
                      ['hpanel', dict(height=400),
                       ['vbox', dict(name='main', width=700),
                        dict(row=['ws', dict(name='channels', group=1)], stretch=1)],
                       ['ws', dict(name='right', width=350, group=2),
                        # (tabname, layout), ...
                        [("Dialogs", ['ws', dict(name='dialogs', group=2)
                                      ]
                          )]
                        ],
                       ],
                      ['hpanel', {},
                       ['ws', dict(name='sub1', width=720, height=520,
                                   group=1)],
                       ['ws', dict(name='sub2', width=720, group=1)],
                       ],
                      ],
                     ], stretch=1),
                    dict(row=['ws', dict(name='toolbar', height=40,
                                             show_tabs=False, group=2)],
                         stretch=0),
                    dict(row=['hbox', dict(name='status')], stretch=0),
                    ]]

global_plugins = [
    Bunch(module='Toolbar', tab='Toolbar', ws='toolbar'),
    Bunch(module='Pan', tab='_pan', ws='uleft', raisekey=None),
    Bunch(module='Info', tab='_info', ws='lleft', raisekey=None),
    Bunch(module='Header', tab='Header', ws='left', raisekey='H'),
    Bunch(module='Zoom', tab='Zoom', ws='left', raisekey='Z'),
    Bunch(module='Thumbs', tab='Thumbs', ws='right', raisekey='T'),
    Bunch(module='Contents', tab='Contents', ws='right', raisekey='C'),
    Bunch(module='WBrowser', tab='Help', ws='right', raisekey='?'),
    Bunch(module='Errors', tab='Errors', ws='right'),
    Bunch(module='Log', tab='Log', ws='right', start=False),
    Bunch(module='Debug', tab='Debug', ws='right', start=False),
    ]

local_plugins = [
    Bunch(module='Pick', ws='dialogs', shortkey='f1'),
    Bunch(module='Ruler', ws='dialogs', shortkey='f2'),
    Bunch(module='MultiDim', ws='dialogs', shortkey='f4'), 
    Bunch(module='Cuts', ws='dialogs', shortkey='f5'),
    Bunch(module='Histogram', ws='dialogs', shortkey='f6'),
    Bunch(module='PixTable', ws='dialogs', shortkey='f7'),
    Bunch(module='Preferences', ws='dialogs', shortkey='f9'),
    Bunch(module='Catalogs', ws='dialogs', shortkey='f10'),
    Bunch(module='Drawing', ws='dialogs', shortkey='f11'),
    Bunch(module='FBrowser', ws='dialogs', shortkey='f12'), 
    ]

default_channels = [('AG', 'sub1'), ('SV', 'sub1'), ('HSCSCAG', 'channels'),
                    ('QDAS_VGW', 'sub2'), ('DSS', 'sub2'),
                    ('SH', 'channels'),
                    ('HSCSHAG', 'channels'), ('HSCSH', 'channels'),
                    ('FMOS', 'channels'), ]

extra_modules = [
    Bunch(module='VGW'),
    ]

extra_plugins = [
    Bunch(module='Region_Selection', ws='dialogs', hidden=True),
    Bunch(module='Sv_Drive', ws='dialogs', hidden=True),
    Bunch(module='AgAutoSelect', ws='dialogs', hidden=True),
    Bunch(module='AgAreaSelection', ws='dialogs', hidden=True),
    ]


class DisplayFITS(GingaControl, GingaView):
    """This class manages the creation and handling of a FITS viewer GUI.
    The class is constructed with a data source and it reads images from the
    source and displays them.
    """
     
    def __init__(self, logger, threadPool, module_manager, preferences,
                 soundsink, ev_quit=None):

        self.controller = None
        self.soundsink = soundsink
        
        GingaView.__init__(self, logger, ev_quit=ev_quit)
        GingaControl.__init__(self, logger, threadPool, module_manager,
                              preferences, ev_quit=ev_quit)


    def load_file(self, filepath, chname=None, wait=True,
                  image_loader=None):
        """Loads a command file from _filepath_ into the commands window.
        """
        try:
            filepath = self.get_filepath(filepath)
            # <-- filepath should now be a real file in the filesystem
            self.logger.debug("filepath=%s" % (filepath))

            image = self.controller.open_fits(filepath, channel=chname,
                                              wait=wait,
                                              image_loader=image_loader)
            return image

        except Exception as e:
            errmsg = "Unable to open '%s': %s" % (
                filepath, str(e))
            self.show_error(errmsg)
            return ro.ERROR


    def play_soundfile(self, filepath, format=None, priority=20):
        self.soundsink.playFile(filepath, format=format,
                                priority=priority)
        
    def gui_load_file(self):
        """Runs dialog to read in a command file into the command window.
        """
        initialdir = os.environ['DATAHOME']
        
        super(DisplayFITS, self).gui_load_file(initialdir=initialdir)
        
        
def main(options, args):
    """Implements the display server.  Creates a DisplayFITS object
    (the GUI), a ReceiveFITS object (ro server) and a datasrc that links
    them together.  It runs until a ^C is used to terminate the server.
    """

    # default of 1000 is a little too small
    sys.setrecursionlimit(2000)
    
    # Create top level logger.
    svcname = options.svcname
    logger = ssdlog.make_logger(svcname, options)

    # Initialize remote objects subsystem.
    try:
        ro.init()

    except ro.remoteObjectError as e:
        logger.error("Error initializing remote objects subsystem: %s" % \
                     str(e))
        sys.exit(1)

        # Otherwise, assume we want to be a GUI/display server
        server(options, args)

    ev_quit = threading.Event()

    # make a name for our monitor
    myMonName = '%s-%s-%d.mon' % (svcname, ro.get_myhost(short=True),
                                  options.monport)

    # monitor channels we are interested in
    channels = []
    if len(options.monchannels) > 0:
        channels = options.monchannels.split(',')

    # Create a local pub sub instance
    mymon = Monitor.Monitor(myMonName, logger,
                            numthreads=options.numthreads)

    threadPool = mymon.get_threadPool()

    sndsink = SoundSink.SoundSource(monitor=mymon, logger=logger,
                                    channels=['sound'])
    
    # Get preferences folder
    if os.environ.has_key('CONFHOME'):
        basedir = os.path.join(os.environ['CONFHOME'], serviceName)
    else:
        basedir = os.path.join(os.environ['HOME'], '.' + serviceName)
    if not os.path.exists(basedir):
        try:
            os.mkdir(basedir)
        except OSError as e:
            logger.warn("Couldn't create ginga settings area (%s): %s" % (
                basedir, str(e)))
            logger.warn("Preferences will not be able to be saved")

    sys.path.insert(0, basedir)
    prefs = Settings.Preferences(basefolder=basedir, logger=logger)

    t_ = prefs.createCategory('general')
    t_.setDefaults(shareReadout=False, useMatplotlibColormaps=False,
                   widgetSet='choose',
                   WCSpkg='astropy', FITSpkg='astropy')
    t_.load(onError='silent')

    # TEMP: ginga needs to find its plugins
    gingaHome = os.path.split(sys.modules['ginga'].__file__)[0]
    childDir = os.path.join(gingaHome, 'gtkw', 'plugins')
    sys.path.insert(0, childDir)
    childDir = os.path.join(gingaHome, 'misc', 'plugins')
    sys.path.insert(0, childDir)

    childDir = os.path.join(basedir, 'plugins')
    sys.path.insert(0, childDir)

    mm = ModuleManager.ModuleManager(logger)
    
    # Start up the display engine
    ginga = DisplayFITS(logger, threadPool, mm, prefs,
                        sndsink, ev_quit=ev_quit)
    ginga.set_layout(default_layout)
    ginga.followFocus(False)

    # User configuration (custom star catalogs, etc.)
    try:
        import ginga_config

        ginga_config.pre_gui_config(ginga)
    except Exception as e:
        try:
            (type, value, tb) = sys.exc_info()
            tb_str = "\n".join(traceback.format_tb(tb))

        except Exception:
            tb_str = "Traceback information unavailable."

        logger.error("Error importing Ginga config file: %s" % (
            str(e)))
        logger.error("Traceback:\n%s" % (tb_str))

    # Build desired layout
    ginga.build_toplevel()
    # TEMP: FIX ME!
    ginga.gpmon.ds = ginga.ds

    # Did user specify geometry
    if options.geometry:
        ginga.setGeometry(options.geometry)

    # Add desired global plugins
    for spec in global_plugins:
        ginga.add_global_plugin(spec)

    # Add GUI log handler (for "Log" global plugin)
    guiHdlr = GuiLogHandler(ginga)
    guiHdlr.setLevel(logging.WARN)
    guiHdlr.setFormatter(ssdlog.get_formatter())
    logger.addHandler(guiHdlr)

    # Add any custom modules
    for spec in extra_modules:
        ginga.add_global_plugin(spec)

    ginga.update_pending()

    # TEMP?
    ginga.ds.raise_tab('Info')
    ginga.ds.raise_tab('Thumbs')

    # Load modules for "local" (per-channel) plug ins
    for spec in local_plugins:
        ginga.add_local_plugin(spec)

    # Add any custom plugins
    for spec in extra_plugins:
        ginga.add_local_plugin(spec)

    # Add custom fitsviewer channels
    for chname, wsname in default_channels:
        datasrc = Datasrc.Datasrc(length=options.bufsize)
        ginga.add_channel(chname, datasrc, workspace=wsname)
    ginga.change_channel('QDAS_VGW')

    receiver = Receive.ReceiveFITS(ginga, mymon, logger)
    ginga.controller = receiver

    # User configuration (custom star catalogs, etc.)
    try:
        ginga_config.post_gui_config(ginga)
    except Exception as e:
        try:
            (type, value, tb) = sys.exc_info()
            tb_str = "\n".join(traceback.format_tb(tb))

        except Exception:
            tb_str = "Traceback information unavailable."

        logger.error("Error processing Ginga config file: %s" % (
            str(e)))
        logger.error("Traceback:\n%s" % (tb_str))

    server_started = False

    # Create receiver and start it
    try:
        # Startup monitor threadpool
        mymon.start(wait=True)

        # Configure logger for logging via our monitor
        ## if options.logmon:
        ##     mymon.logmon(logger, options.logmon, ['logs'])

        # start_server is necessary if we are subscribing, but not if only
        # publishing
        mymon.start_server(wait=True, port=options.monport)
        server_started = True

        if options.monitor:
            # subscribe our monitor to the central monitor hub
            if len(channels) > 0:
                mymon.subscribe_remote(options.monitor, channels, {})
            # publishing for remote command executions
            mymon.publish_to(options.monitor, ['sound', 'g2task'], {})

        # Register local fits info subscription callback
        mymon.subscribe_cb(receiver.arr_fitsinfo, ['fits'])
        mymon.subscribe_cb(receiver.arr_taskinfo, ['taskmgr'])

        # Create our remote service object
        viewsvc = ro.remoteObjectServer(svcname=options.svcname,
                                        obj=receiver,
                                        logger=logger, ev_quit=ev_quit,
                                        port=options.port,
                                        usethread=True,
                                        threadPool=threadPool)


        # Assume remaining arguments are fits files and load them into
        # the datasrc.
        for arg in args:
            receiver.display_fitsfile(arg)
        
        logger.info("Starting fits viewing service.")
        viewsvc.ro_start()

        try:
            # Main loop to handle GTK events
            ginga.mainloop(timeout=0.001)

        except KeyboardInterrupt:
            logger.error("Received keyboard interrupt!")

        except Exception as e:
            logger.error("Received exception: %s" % (str(e)))

    finally:
        logger.info("Shutting down...")

        ginga.stop()
        viewsvc.ro_stop(wait=True)
        mymon.stop_server(wait=True)
        mymon.stop(wait=True)

    sys.exit(0)
        

if __name__ == "__main__":
   
    # Parse command line options with nifty new optparse module
    from optparse import OptionParser

    usage = "usage: %prog [options] cmd [args]"
    optprs = OptionParser(usage=usage, version=('%%prog %s' % version))
    
    optprs.add_option("--bufsize", dest="bufsize", metavar="NUM",
                      type="int", default=25,
                      help="Buffer length to NUM")
    optprs.add_option("--channels", dest="channels", default="Image",
                      help="Specify list of channels to create")
    optprs.add_option("--debug", dest="debug", default=False, action="store_true",
                      help="Enter the pdb debugger on main()")
    optprs.add_option("--display", dest="display", metavar="HOST:N",
                      help="Use X display on HOST:N")
    optprs.add_option("--fits", dest="fitsfile", metavar="FILE",
                      help="Send FITS file FILE to the server")
    optprs.add_option("--fitspath", dest="fitspath",
                      help="Send a FITS pathname to the server")
    optprs.add_option("-g", "--geometry", dest="geometry",
                      metavar="GEOM", default="+20+100",
                      help="X geometry for initial size and placement")
    optprs.add_option("--modules", dest="modules", metavar="NAMES",
                      help="Specify additional modules to load")
    optprs.add_option("--monitor", dest="monitor", metavar="NAME",
                      default='monitor',
                      help="Synchronize from monitor named NAME")
    optprs.add_option("--monchannels", dest="monchannels", 
                      default='', metavar="NAMES",
                      help="Specify monitor channels to subscribe to")
    optprs.add_option("--monport", dest="monport", type="int",
                      help="Register monitor using PORT", metavar="PORT")
    optprs.add_option("--numthreads", dest="numthreads", type="int",
                      default=200,
                      help="Start NUM threads in thread pool", metavar="NUM")
    optprs.add_option("--plugins", dest="plugins", metavar="NAMES",
                      help="Specify additional plugins to load")
    optprs.add_option("--port", dest="port", type="int", default=None,
                      help="Register using PORT", metavar="PORT")
    optprs.add_option("--profile", dest="profile", action="store_true",
                      default=False,
                      help="Run the profiler on main()")
    optprs.add_option("--svcname", dest="svcname", metavar="NAME",
                      default=serviceName,
                      help="Register using NAME as service name")
    ssdlog.addlogopts(optprs)

    (options, args) = optprs.parse_args(sys.argv[1:])

    if options.display:
        os.environ['DISPLAY'] = options.display

    # Are we debugging this?
    if options.debug:
        import pdb

        pdb.run('main(options, args)')

    # Are we profiling this?
    elif options.profile:
        import profile

        print ("%s profile:" % sys.argv[0])
        profile.run('main(options, args)')


    else:
        main(options, args)

# END
