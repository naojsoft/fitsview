#!/usr/bin/env python
#
# guideview.py -- Gen2 guide image viewer/controller.
#
# Eric Jeschke (eric@naoj.org)
#
"""
guideview.py implements a guide image display and control interface.

Usage:
    guideview.py --monport=NNNNN --loglevel=20 --stderr
"""
from __future__ import print_function

# stdlib imports
import sys, os
import threading
import logging
import traceback
from six.moves import map
from six.moves import zip

# g2base imports
from g2base import ssdlog
from g2base.remoteObjects import remoteObjects as ro
from g2base.remoteObjects import Monitor

# Gen2 imports
import Gen2.soundsink as SoundSink

# Ginga imports
from ginga.misc import ModuleManager, Datasrc, Settings
from ginga.misc.Bunch import Bunch
import ginga.toolkit as ginga_toolkit

# Local application imports
from Gen2.fitsview.util import Receive

moduleHome = os.path.split(sys.modules[__name__].__file__)[0]
sys.path.insert(0, moduleHome)
pluginHome = os.path.join(moduleHome, 'plugins')
sys.path.insert(0, pluginHome)

serviceName = 'guideview'
version = "20161130.0"

default_layout = ['seq', {},
                   ['vbox', dict(name='top', width=1600, height=1100),
                    dict(row=['hbox', dict(name='menu')],
                         stretch=0),
                    dict(row=['hpanel', {},
                     ['ws', dict(name='left', width=300),
                      # (tabname, layout), ...
                      [("Info", ['vpanel', {},
                                 ['ws', dict(name='uleft', height=300,
                                             show_tabs=False, group=3)],
                                 ['ws', dict(name='lleft', height=430,
                                             show_tabs=True, group=3)],
                                 ]
                        )]
                      ],
                     ['vpanel', dict(width=1300),
                      ['hpanel', dict(height=450),
                       ['vbox', dict(name='main', width=600),
                        dict(row=['ws', dict(name='channels', group=1)], stretch=1),
                        dict(row=['ws', dict(wstype='stack', name='cbar',
                                             group=99)], stretch=0),
                        dict(row=['ws', dict(wstype='stack', name='readout',
                                           group=99)], stretch=0),
                        dict(row=['ws', dict(wstype='stack', name='operations',
                                             group=99)], stretch=0),
                        ],
                       ['ws', dict(name='right', width=600, group=2),
                        # (tabname, layout), ...
                        [("Dialogs", ['ws', dict(name='dialogs', group=2)
                                      ]
                          )]
                        ],
                       ],
                      ['hpanel', dict(height=550),
                       ['ws', dict(name='sub1', width=600, height=420,
                                   group=1)],
                       ['ws', dict(name='sub2', width=600, group=1)],
                       ],
                      ],
                     ], stretch=1),
                    dict(row=['ws', dict(name='toolbar', height=40,
                                             show_tabs=False, group=2)],
                         stretch=0),
                    dict(row=['hbox', dict(name='status')], stretch=0),
                    ]]

global_plugins = [
    Bunch(module='Operations', workspace='operations', start=True,
          hidden=True, category='system'),
    Bunch(module='Toolbar', workspace='toolbar', start=True,
          hidden=True, category='system'),
    Bunch(module='Pan', workspace='uleft', start=True,
          hidden=True, category='system'),
    Bunch(module='Info', tab='Synopsis', workspace='lleft', start=True,
          hidden=True, category='system'),
    Bunch(module='Header', tab='Header', workspace='left', start=True,
          hidden=True, category='system'),
    Bunch(module='Zoom', tab='Zoom', workspace='left', start=True,
          hidden=True, category='system'),
    ## Bunch(module='Thumbs', tab='Thumbs', workspace='right', start=True,
    ##       hidden=True, category='system'),
    ## Bunch(module='Contents', tab='Contents', workspace='right', start=True,
    ##       hidden=True, category='system'),
    Bunch(module='Colorbar', workspace='cbar', start=True,
          hidden=True, category='system'),
    Bunch(module='Cursor', workspace='readout', start=True,
          hidden=True, category='system'),
    Bunch(module='Errors', tab='Errors', workspace='right', start=True,
          hidden=True, category='system'),
    Bunch(module='Command', tab='Command', workspace='right', start=False,
          category='Global'),
    Bunch(module='Log', tab='Log', workspace='right', start=False,
          category='Global'),
    Bunch(module='WBrowser', tab='Help', workspace='right', start=False,
          category='Global'),
    ## Bunch(module='FBrowser', tab='Open File', workspace='right', start=False,
    ##       category='Global'),
    ## Bunch(module='Blink', tab='Blink Channels', workspace='right', start=False,
    ##       category='Global'),
    Bunch(module='ColorMapPicker', tab='Color Map Picker', workspace='right',
          start=False, category='Global'),
    ]

local_plugins = [
    Bunch(module='Pick', workspace='dialogs', category=None),
    Bunch(module='Ruler', workspace='dialogs', category=None),
    ## Bunch(module='MultiDim', workspace='lleft', category=None),
    Bunch(module='Cuts', workspace='dialogs', category=None),
    Bunch(module='Histogram', workspace='dialogs', category=None),
    Bunch(module='Crosshair', workspace='dialogs', category=None),
    Bunch(module='Overlays', workspace='dialogs', category=None),
    ## Bunch(module='Blink', workspace='dialogs', category=None),
    Bunch(module='PixTable', workspace='dialogs', category=None),
    Bunch(module='Preferences', workspace='dialogs', category=None),
    Bunch(module='Catalogs', workspace='dialogs', category=None),
    ## Bunch(module='Mosaic', workspace='dialogs', category=None),
    Bunch(module='Drawing', workspace='dialogs', category=None),
    Bunch(module='FBrowser', workspace='dialogs', category=None),
    ## Bunch(module='Compose', workspace='dialogs', category=None),
    Bunch(module='ScreenShot', workspace='dialogs', category=None),
    ]

default_channels = [('AG', 'channels'), ('SV', 'sub1'), ('HSCSCAG', 'sub1'),
                    ('QDAS_VGW', 'sub2'), ('DSS', 'sub2'),
                    ('SH', 'channels'),
                    ('HSCSHAG', 'channels'), ('HSCSH', 'channels'),
                    ('FMOS', 'channels'), ]

extra_modules = [
    Bunch(module='VGW', hidden=True),
    ]

extra_plugins = [
    Bunch(module='Region_Selection', ws='dialogs', hidden=True),
    Bunch(module='Sv_Drive', ws='dialogs', hidden=True),
    Bunch(module='AgAutoSelect', ws='dialogs', hidden=True),
    Bunch(module='AgAreaSelection', ws='dialogs', hidden=True),
    ]


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
    if 'CONFHOME' in os.environ:
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

    settings = prefs.create_category('general')
    settings.set_defaults(share_readout=False, useMatplotlibColormaps=False,
                          widgetSet='choose',
                          pixel_coords_offset=1.0,
                          WCSpkg='kapteyn', FITSpkg='astropy')
    settings.load(onError='silent')

    # Choose a toolkit
    if options.toolkit:
        toolkit = options.toolkit
    else:
        toolkit = settings.get('widgetSet', 'choose')

    ginga_toolkit.use(toolkit)
    tkname = ginga_toolkit.get_family()

    # TEMP: ginga needs to find its plugins
    gingaHome = os.path.split(sys.modules['ginga'].__file__)[0]
    widgetDir = tkname + 'w'
    childDir = os.path.join(gingaHome, widgetDir, 'plugins')
    sys.path.insert(0, childDir)
    childDir = os.path.join(gingaHome, 'rv', 'plugins')
    sys.path.insert(0, childDir)

    childDir = os.path.join(basedir, 'plugins')
    sys.path.insert(0, childDir)

    mm = ModuleManager.ModuleManager(logger)

    # Start up the display engine
    from Gen2.fitsview.g2viewer import Gen2FITSViewer

    ginga = Gen2FITSViewer(logger, threadPool, mm, prefs,
                           sndsink, ev_quit=ev_quit)
    ginga.set_layout(default_layout)
    ginga.follow_focus(False)

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
        ginga.set_geometry(options.geometry)

    # Add desired global plugins
    for spec in global_plugins:
        ginga.add_global_plugin(spec)

    # Add GUI log handler (for "Log" global plugin)
    from ginga.rv.Control import GuiLogHandler
    guiHdlr = GuiLogHandler(ginga)
    guiHdlr.setLevel(logging.WARN)
    guiHdlr.setFormatter(ssdlog.get_formatter())
    logger.addHandler(guiHdlr)

    # Add any custom modules
    for spec in extra_modules:
        ginga.add_global_plugin(spec)

    ginga.update_pending()

    # TEMP?
    tab_names = list(map(lambda name: name.lower(),
                         ginga.ds.get_tabnames(group=None)))
    if 'info' in tab_names:
        ginga.ds.raise_tab('Info')
    if 'thumbs' in tab_names:
        ginga.ds.raise_tab('Thumbs')

    # Load modules for "local" (per-channel) plug ins
    for spec in local_plugins:
        ginga.add_local_plugin(spec)

    # Add any custom plugins
    for spec in extra_plugins:
        ginga.add_local_plugin(spec)

    # Add custom fitsviewer channels
    for chname, wsname in default_channels:
        ginga.add_channel(chname, workspace=wsname)
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
            # Main loop to handle window events
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
    optprs.add_option("-t", "--toolkit", dest="toolkit", metavar="NAME",
                      default='gtk',
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

        print(("%s profile:" % sys.argv[0]))
        profile.run('main(options, args)')


    else:
        main(options, args)

# END
