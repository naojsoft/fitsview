#!/usr/bin/env python
#
# fitsview.py -- Gen2 Quick Look viewer/display server.
#
# Eric Jeschke (eric@naoj.org)
#
"""
fitsview.py implements a FITS viewer/display server for Gen2 quick look
  activities.

Usage:
    fitsview.py --monport=NNNNN --loglevel=20
"""
from __future__ import print_function

# stdlib imports
import sys, os
import threading
import logging
import traceback

# ginga imports
from ginga.misc import ModuleManager, Datasrc, Settings
from ginga.misc.Bunch import Bunch
import ginga.toolkit as ginga_toolkit

# g2base imports
from g2base import ssdlog
from g2base.remoteObjects import remoteObjects as ro
from g2base.remoteObjects import Monitor

# Gen2 imports
import g2client.soundsink as SoundSink

# Local application imports
from Gen2.fitsview.util import Receive

moduleHome = os.path.split(sys.modules[__name__].__file__)[0]
sys.path.insert(0, moduleHome)
pluginHome = os.path.join(moduleHome, 'plugins')
sys.path.insert(0, pluginHome)

defaultServiceName = 'fitsview'
version = "20171220.0"

default_layout = ['seq', {},
                   ['vbox', dict(name='top', width=1600, height=900),
                    dict(row=['hbox', dict(name='menu')],
                         stretch=0),
                    dict(row=['hpanel', dict(name='hpnl'),
                     ['ws', dict(name='left', wstype='tabs',
                                 width=340, group=2),
                      # (tabname, layout), ...
                      [("Info", ['vpanel', {},
                                 ['ws', dict(name='uleft', wstype='stack',
                                             height=300, group=3)],
                                 ['ws', dict(name='lleft', wstype='tabs',
                                             height=530, group=3)],
                                 ]
                        )]],
                     ['vbox', dict(name='main', width=760),
                      dict(row=['ws', dict(wstype='tabs', name='channels',
                                           group=1, use_toolbar=True)],
                           stretch=1),
                      dict(row=['ws', dict(wstype='stack', name='cbar',
                                           group=99)], stretch=0),
                      dict(row=['ws', dict(wstype='stack', name='readout',
                                           group=99)], stretch=0),
                      dict(row=['ws', dict(wstype='stack', name='operations',
                                           group=99)], stretch=0),
                      ],
                     ['ws', dict(name='right', width=500, group=2),
                      # (tabname, layout), ...
                      [("Dialogs", ['ws', dict(name='dialogs',
                                                wstype='tabs', group=2)
                                    ]
                        )]
                      ],
                     ], stretch=1),
                    dict(row=['ws', dict(name='toolbar', wstype='stack',
                                         height=40,
                                         show_tabs=False, group=2)],
                         stretch=0),
                    dict(row=['hbox', dict(name='status')], stretch=0),
                    ]]

plugins = [
    # hidden plugins, started at program initialization
    Bunch(module='Operations', workspace='operations', start=True,
          hidden=True, category='System', ptype='global'),
    Bunch(module='Toolbar', workspace='toolbar', start=True,
          hidden=True, category='System', ptype='global'),
    Bunch(module='Pan', workspace='uleft', start=True,
          hidden=True, category='System', ptype='global'),
    Bunch(module='Info', tab='Synopsis', workspace='lleft', start=True,
          hidden=True, category='System', ptype='global'),
    Bunch(module='Thumbs', tab='Thumbs', workspace='right', start=True,
          hidden=True, category='System', ptype='global'),
    Bunch(module='Contents', tab='Contents', workspace='right', start=True,
          hidden=True, category='System', ptype='global'),
    Bunch(module='Colorbar', workspace='cbar', start=True,
          hidden=True, category='System', ptype='global'),
    Bunch(module='Cursor', workspace='readout', start=True,
          hidden=True, category='System', ptype='global'),
    Bunch(module='Errors', tab='Errors', workspace='right', start=True,
          hidden=True, category='System', ptype='global'),

    # optional, user-started plugins
    ## Bunch(module='Blink', tab='Blink Channels', workspace='right', start=False,
    ##       menu="Blink Channels [G]", category='Analysis', ptype='global'),
    ## Bunch(module='Blink', workspace='dialogs', menu='Blink Images',
    ##       category='Analysis', ptype='local'),
    Bunch(module='Cuts', workspace='dialogs', category='Analysis',
          ptype='local'),
    Bunch(module='LineProfile', workspace='dialogs',
          category='Analysis.Datacube', ptype='local'),
    Bunch(module='Histogram', workspace='dialogs', category='Analysis',
          ptype='local'),
    Bunch(module='Overlays', workspace='dialogs', category='Analysis',
          ptype='local'),
    Bunch(module='Pick', workspace='dialogs', category='Analysis',
          ptype='local'),
    Bunch(module='PixTable', workspace='dialogs', category='Analysis',
          ptype='local'),
    ## Bunch(module='TVMark', workspace='dialogs', category='Analysis',
    ##       ptype='local'),
    ## Bunch(module='TVMask', workspace='dialogs', category='Analysis',
    ##       ptype='local'),
    ## Bunch(module='WCSMatch', tab='WCSMatch', workspace='right', start=False,
    ##       menu="WCS Match [G]", category='Analysis', ptype='global'),
    Bunch(module='Command', tab='Command', workspace='lleft', start=False,
          menu="Command Line [G]", category='Debug', ptype='global'),
    Bunch(module='Log', tab='Log', workspace='right', start=False,
          menu="Logger Info [G]", category='Debug', ptype='global'),
    Bunch(module='MultiDim', workspace='lleft', category='Navigation',
          ptype='local'),
    ## Bunch(module='IRAF', tab='IRAF', workspace='right', start=False,
    ##       menu="IRAF Interface [G]", category='Remote', ptype='global'),
    ## Bunch(module='RC', tab='RC', workspace='right', start=False,
    ##       menu="Remote Control [G]", category='Remote', ptype='global'),
    ## Bunch(module='SAMP', tab='SAMP', workspace='right', start=False,
    ##       menu="SAMP Client [G]", category='Remote', ptype='global'),
    ## Bunch(module='Compose', workspace='dialogs', category='RGB', ptype='local'),
    Bunch(module='ScreenShot', workspace='dialogs', category='RGB',
          ptype='local'),
    Bunch(module='ColorMapPicker', tab='ColorMapPicker',
          menu="Set Color Map [G]", workspace='right', start=False,
          category='RGB', ptype='global'),
    ## Bunch(module='PlotTable', workspace='dialogs', category='Table',
    ##       ptype='local'),
    Bunch(module='Catalogs', workspace='dialogs', category='Utils',
          ptype='local'),
    Bunch(module='Crosshair', workspace='dialogs', category='Utils',
          ptype='local'),
    Bunch(module='Drawing', workspace='dialogs', category='Utils',
          ptype='local'),
    Bunch(module='FBrowser', workspace='dialogs', category='Utils',
          ptype='local'),
    ## Bunch(module='ChangeHistory', tab='History', workspace='right',
    ##       menu="History [G]", start=False, category='Utils', ptype='global'),
    ## Bunch(module='Mosaic', workspace='dialogs', category='Utils', ptype='local'),
    Bunch(module='FBrowser', tab='Open File', workspace='right',
          menu="Open File [G]", start=False, category='Utils', ptype='global'),
    Bunch(module='Preferences', workspace='dialogs', category='Utils',
          ptype='local'),
    Bunch(module='Ruler', workspace='dialogs', category='Analysis', ptype='local'),
    # TODO: Add SaveImage to File menu.
    ## Bunch(module='SaveImage', tab='SaveImage', workspace='right',
    ##       menu="Save File [G]", start=False, category='Utils', ptype='global'),
    Bunch(module='WCSAxes', workspace='dialogs', category='Utils',
          ptype='local'),
    Bunch(module='WBrowser', tab='Help', workspace='channels', start=False,
          menu="Help [G]", category='Help', ptype='global'),
    Bunch(module='Header', tab='Header', workspace='left', start=True,
          menu="Header [G]", hidden=False, category='Utils', ptype='global',
          optray=False),
    Bunch(module='Zoom', tab='Zoom', workspace='left', start=False,
          menu="Zoom [G]", category='Utils', ptype='global',
          optray=False),

    # Subaru-specific plugins
    Bunch(module='GView', tab='GView', ws='right', ptype='global', start=False,
          category='Subaru'),
    Bunch(module='CHARIS', ws='right', ptype='global', start=False,
          category='Subaru'),

]

def main(options, args):
    """Implements the display server.  Creates a DisplayFITS object
    (the GUI), a ReceiveFITS object (ro server) and a datasrc that links
    them together.  It runs until a ^C is used to terminate the server.
    """
    global plugins

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

    status_srv = ro.remoteObjectProxy('status')

    # Get preferences folder
    if 'CONFHOME' in os.environ:
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
    settings = prefs.create_category('general')
    settings.load(onError='silent')
    settings.set_defaults(useMatplotlibColormaps=False,
                          recursion_limit=2000,
                          widgetSet='choose',
                          WCSpkg='astropy', FITSpkg='astropy',
                          font_scaling_factor=None)

    # default of 1000 is a little too small
    sys.setrecursionlimit(settings.get('recursion_limit'))

    # Choose a toolkit
    if options.toolkit:
        toolkit = options.toolkit
    else:
        toolkit = settings.get('widgetSet', 'choose')

    ginga_toolkit.use(toolkit)
    tkname = ginga_toolkit.get_family()

    logger.info("Chosen toolkit (%s) family is '%s'" % (
        ginga_toolkit.toolkit, tkname))

    if settings.get('useMatplotlibColormaps', False):
        # Add matplotlib color maps if matplotlib is installed
        try:
            from ginga import cmap
            cmap.add_matplotlib_cmaps()
        except Exception as e:
            logger.warning("failed to load matplotlib colormaps: %s" % (str(e)))

    # User wants to customize the WCS package?
    wcspkg = settings.get('WCSpkg', 'choose')
    try:
        from ginga.util import wcsmod
        if wcspkg != 'choose':
            assert wcsmod.use(wcspkg) is True
    except Exception as e:
        logger.warning("failed to set WCS package preference: %s" % (str(e)))

    # User wants to customize the FITS package?
    fitspkg = settings.get('FITSpkg', 'choose')
    try:
        from ginga.util import io_fits
        if wcspkg != 'choose':
            assert io_fits.use(fitspkg) is True
    except Exception as e:
        logger.warning("failed to set FITS package preference: %s" % (str(e)))

    # Check whether user wants to use OpenCv
    use_opencv = settings.get('use_opencv', False)
    if use_opencv:
        from ginga import trcalc
        try:
            trcalc.use('opencv')
        except Exception as e:
            logger.warning("failed to set OpenCv preference: %s" % (str(e)))

    # Check whether user wants to use OpenCL
    use_opencl = settings.get('use_opencl', False)
    if use_opencl:
        from ginga import trcalc
        try:
            trcalc.use('opencl')
        except Exception as e:
            logger.warning("failed to set OpenCL preference: %s" % (str(e)))

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

    ginga_shell = Gen2FITSViewer(logger, threadPool, mm, prefs,
                                 sndsink, status_srv, ev_quit=ev_quit)
    ginga_shell.set_layout(default_layout)
    ginga_shell.follow_focus(True)

    # user wants to set font scaling.
    # NOTE: this happens *after* creation of shell object, since
    # Application object constructor will also set this
    font_scaling = settings.get('font_scaling_factor', None)
    if font_scaling is not None:
        logger.debug("font_scaling_factor={}".format(font_scaling))
        from ginga.fonts import font_asst
        font_asst.default_scaling_factor = font_scaling

    # User configuration (custom star catalogs, etc.)
    try:
        import ginga_config

        ginga_config.pre_gui_config(ginga_shell)
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
    ginga_shell.build_toplevel()

    # Did user specify a particular geometry?
    if options.geometry:
        ginga_shell.set_geometry(options.geometry)

    # Add GUI log handler (for "Log" global plugin)
    from ginga.rv.Control import GuiLogHandler
    guiHdlr = GuiLogHandler(ginga_shell)
    guiHdlr.setLevel(logging.WARN)
    guiHdlr.setFormatter(ssdlog.get_formatter())
    logger.addHandler(guiHdlr)

    # make the list of disabled plugins
    disabled_plugins = []
    if not (options.disable_plugins is None):
        disabled_plugins = options.disable_plugins.lower().split(',')

    # Load built in plugins
    for spec in plugins:
        if spec.module.lower() not in disabled_plugins:
            ginga_shell.add_plugin(spec)

    # Load any custom modules
    if options.modules:
        modules = options.modules.split(',')
        for pluginName in modules:
            spec = Bunch(name=pluginName, module=pluginName,
                         tab=pluginName, ws='right', hidden=True,
                         category='Custom', ptype='global')
            ginga_shell.add_plugin(spec)

    # Load any custom plugins
    if options.plugins:
        plugins = options.plugins.split(',')
        for pluginName in plugins:
            spec = Bunch(module=pluginName, ws='dialogs', ptype='local',
                         hidden=True, category='Custom')
            ginga_shell.add_plugin(spec)

    # start any plugins that have start=True
    ginga_shell.boot_plugins()
    ginga_shell.update_pending()

    # TEMP?
    tab_names = [name.lower()
                 for name in ginga_shell.ds.get_tabnames(group=None)]
    if 'info' in tab_names:
        ginga_shell.ds.raise_tab('Info')
    if 'thumbs' in tab_names:
        ginga_shell.ds.raise_tab('Thumbs')

    # Add custom fitsviewer channels
    instruments = options.channels.split(',')
    for chname in instruments:
        ginga_shell.add_channel(chname)
    ginga_shell.change_channel(instruments[0])

    receiver = Receive.ReceiveFITS(ginga_shell, mymon, logger)
    ginga_shell.controller = receiver

    # User configuration (custom star catalogs, etc.)
    try:
        ginga_config.post_gui_config(ginga_shell)
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
            ginga_shell.mainloop(timeout=0.001)

        except KeyboardInterrupt:
            logger.error("Received keyboard interrupt!")

        except Exception as e:
            logger.error("Received exception: %s" % (str(e)))

    finally:
        logger.info("Shutting down...")

        ginga_shell.stop()
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
    optprs.add_option("--disable-plugins", dest="disable_plugins",
                      metavar="NAMES",
                      help="Specify plugins that should be disabled")
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
