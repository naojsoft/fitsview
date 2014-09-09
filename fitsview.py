#!/usr/bin/env python
#
# fitsview.py -- Simple FITS viewer/display server.
#
# Eric Jeschke (eric@naoj.org)
#
"""
fitsview.py implements a simple FITS viewer/display server to display FITS
images in GTK widgets.

Usage:
    fitsview.py 
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

from ginga.misc import ModuleManager, Datasrc, Settings
from ginga.misc.Bunch import Bunch
from ginga.Control import GingaControl, GuiLogHandler
import ginga.toolkit as ginga_toolkit

# Local application imports
from util import Receive

defaultServiceName = 'fitsview'
version = "20140905.0"

default_layout = ['seq', {},
                   ['vbox', dict(name='top', width=1600, height=1000),
                    dict(row=['hbox', dict(name='menu')],
                         stretch=0),
                    dict(row=['hpanel', dict(name='hpnl'),
                     ['ws', dict(name='left', width=300, group=2),
                      # (tabname, layout), ...
                      [("Info", ['vpanel', {},
                                 ['ws', dict(name='uleft', height=300,
                                             show_tabs=False, group=3)],
                                 ['ws', dict(name='lleft', height=430,
                                             show_tabs=False, group=3)],
                                 ]
                        )]],
                     ['vbox', dict(name='main', width=700),
                      dict(row=['ws', dict(name='channels', group=1)], stretch=1)],
                     ['ws', dict(name='right', width=350, group=2),
                      # (tabname, layout), ...
                      [("Dialogs", ['ws', dict(name='dialogs', group=2)
                                    ]
                        )]
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
    Bunch(module='Contents', tab='Contents', ws='right', raisekey='c'),
    Bunch(module='WBrowser', tab='Help', ws='channels', raisekey='?', start=False),
    Bunch(module='Errors', tab='Errors', ws='right', start=True),
    Bunch(module='RC', tab='RC', ws='right', start=False),
    Bunch(module='SAMP', tab='SAMP', ws='right', start=False),
    Bunch(module='IRAF', tab='IRAF', ws='right', start=False),
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
    Bunch(module='Mosaic', ws='dialogs'),
    Bunch(module='Pipeline', ws='dialogs'),
    Bunch(module='Catalogs', ws='dialogs', shortkey='f10'),
    Bunch(module='Drawing', ws='dialogs', shortkey='f11'),
    Bunch(module='FBrowser', ws='dialogs', shortkey='f12'),
    Bunch(module='SPCAM', ws='dialogs'),
    Bunch(module='HSC', ws='dialogs'),
    ]

def get_displayfits(viewKlass):
    class DisplayFITS(GingaControl, viewKlass):
        """This class manages the creation and handling of a FITS viewer GUI.
        The class is constructed with a data source and it reads images from the
        source and displays them.
        """

        def __init__(self, logger, threadPool, module_manager, preferences,
                     soundsink, ev_quit=None):

            self.controller = None
            self.soundsink = soundsink

            viewKlass.__init__(self, logger, ev_quit=ev_quit)
            GingaControl.__init__(self, logger, threadPool, module_manager,
                                  preferences, ev_quit=ev_quit)


        def load_file(self, filepath, chname=None, wait=True,
                      image_loader=None):
            """Loads a command file from _filepath_ into the commands window.
            """
            try:
                # TODO: what to do about image_loader parameter?
                
                filepath = self.get_filepath(filepath)
                # <-- filepath should now be a real file in the filesystem
                self.logger.debug("filepath=%s" % (filepath))
        
                image = self.controller.open_fits(filepath, channel=chname,
                                                  wait=wait)
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

    return DisplayFITS

        
def main(options, args):
    """Implements the display server.  Creates a DisplayFITS object
    (the GUI), a ReceiveFITS object (ro server) and a datasrc that links
    them together.  It runs until a ^C is used to terminate the server.
    """

    # default of 1000 is a little too small
    sys.setrecursionlimit(2000)
    
    # Create top level logger.
    logger = ssdlog.make_logger(options.svcname, options)

    # Initialize remote objects subsystem.
    rohosts = options.rohosts.split(',')
    try:
        ro.init(rohosts)

    except ro.remoteObjectError as e:
        logger.error("Error initializing remote objects subsystem: %s" % \
                     str(e))
        sys.exit(1)

    ev_quit = threading.Event()

    # make a name for our monitor
    myMonName = 'fitsview-%s-%d.mon' % (ro.get_myhost(short=True),
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
        basedir = os.path.join(os.environ['CONFHOME'], 'fitsview')
    else:
        basedir = os.path.join(os.environ['HOME'], '.fitsview')
    if not os.path.exists(basedir):
        try:
            os.mkdir(basedir)
        except OSError as e:
            logger.warn("Couldn't create ginga settings area (%s): %s" % (
                basedir, str(e)))
            logger.warn("Preferences will not be able to be saved")

    sys.path.insert(0, basedir)
    prefs = Settings.Preferences(basefolder=basedir, logger=logger)
    settings = prefs.createCategory('general')
    settings.load(onError='silent')
    settings.setDefaults(useMatplotlibColormaps=False,
                         widgetSet='choose',
                         WCSpkg='kapteyn', FITSpkg='astropy')

    # Choose a toolkit
    if options.toolkit:
        toolkit = options.toolkit
    else:
        toolkit = settings.get('widgetSet', 'choose')

    ginga_toolkit.use(toolkit)
    tkname = ginga_toolkit.get_family()
    
    if tkname == 'gtk':
        from ginga.gtkw.GingaGtk import GingaView
    elif tkname == 'qt':
        from ginga.qtw.GingaQt import GingaView
    else:
        try:
            from ginga.qtw.GingaQt import GingaView
        except ImportError:
            try:
                from ginga.gtkw.GingaGtk import GingaView
            except ImportError:
                print("You need python-gtk or python-qt4 to run Ginga!")
                sys.exit(1)

    # TEMP: ginga needs to find its plugins
    gingaHome = os.path.split(sys.modules['ginga'].__file__)[0]
    ## widgetDir = tkname + 'w'
    ## childDir = os.path.join(gingaHome, widgetDir, 'plugins')
    ## sys.path.insert(0, childDir)
    childDir = os.path.join(gingaHome, 'misc', 'plugins')
    sys.path.insert(0, childDir)

    childDir = os.path.join(basedir, 'plugins')
    sys.path.insert(0, childDir)

    mm = ModuleManager.ModuleManager(logger)
    
    # Start up the display engine
    disp_klass = get_displayfits(GingaView)
    ginga = disp_klass(logger, threadPool, mm, prefs,
                       sndsink, ev_quit=ev_quit)
    ginga.set_layout(default_layout)
    ginga.followFocus(True)

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

    # Did user specify a particular geometry?
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

    # Load any custom modules
    if options.modules:
        modules = options.modules.split(',')
        for pluginName in modules:
            spec = Bunch(name=pluginName, module=pluginName,
                         tab=pluginName, ws='right')
            ginga.add_global_plugin(spec)

    ginga.update_pending()

    # TEMP?
    ginga.ds.raise_tab('Info')
    ginga.ds.raise_tab('Thumbs')

    # Load modules for "local" (per-channel) plug ins
    for spec in local_plugins:
        ginga.add_local_plugin(spec)

    # Load any custom plugins
    if options.plugins:
        plugins = options.plugins.split(',')
        for pluginName in plugins:
            spec = Bunch(module=pluginName, ws='dialogs',
                         hidden=True)
            ginga.add_local_plugin(spec)

    # Add custom fitsviewer channels
    instruments = options.channels.split(',')
    for chname in instruments:
        datasrc = Datasrc.Datasrc(length=options.bufsize)
        ginga.add_channel(chname, datasrc)
    ginga.change_channel(instruments[0])

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

    # Parse command line options with nifty optparse module
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
                      default=30,
                      help="Start NUM threads in thread pool", metavar="NUM")
    optprs.add_option("--plugins", dest="plugins", metavar="NAMES",
                      help="Specify additional plugins to load")
    optprs.add_option("--port", dest="port", type="int", default=None,
                      help="Register using PORT", metavar="PORT")
    optprs.add_option("--profile", dest="profile", action="store_true",
                      default=False,
                      help="Run the profiler on main()")
    optprs.add_option("--rohosts", dest="rohosts", metavar="NAME",
                      default='localhost',
                      help="List of Gen2 remote object hosts")
    optprs.add_option("--svcname", dest="svcname", metavar="NAME",
                      default=defaultServiceName,
                      help="Register using NAME as service name")
    optprs.add_option("-t", "--toolkit", dest="toolkit", metavar="NAME",
                      default=None,
                      help="Prefer GUI toolkit (gtk|qt)")
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
