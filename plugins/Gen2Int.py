#
# Gen2.py -- Gen2 plugin for Ginga FITS viewer
#
# Eric Jeschke (eric@naoj.org)
#
import sys, os
import re
import time
import traceback
from io import BytesIO

from astropy.io import fits as pyfits
import numpy as np

from ginga import GingaPlugin
from ginga import AstroImage, BaseImage
from ginga.util import loader, iohelper
from ginga.misc import Future, Bunch

# g2cam imports
from g2base.remoteObjects import remoteObjects as ro
from g2base.remoteObjects import Monitor
from g2base.astro.frame import Frame
from g2cam.INS import INSdata

# Gen2 imports
import g2client.soundsink as SoundSink


# Regex used to discover/parse frame info
regex_frame = re.compile(r'^mon\.frame\.(\w+)\.Archiver$')


class Gen2Error(Exception):
    pass

class Gen2Int(GingaPlugin.GlobalPlugin):
    """
    NOTE: *** All these methods are running as the GUI thread, unless
    otherwise noted. Do not block!! ***
    """

    def __init__(self, fv):
        # superclass defines some variables for us, like logger
        super(Gen2Int, self).__init__(fv)

        # get plugin settings
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_Gen2Int')
        self.settings.add_defaults(svcname='fitsview', port=13088,
                                   monport=13089, monitor_name='monitor',
                                   subscribe_fits=False, disp_cache=None)
        self.settings.load()

        self.svcname = self.settings.get('svcname')
        self.port = self.settings.get('port')
        self.monport = self.settings.get('monport')
        self.monitor_name = self.settings.get('monitor_name')
        # points to a possible display cache of images
        self.disp_cache = self.settings.get('disp_cache')

        self.host = ro.get_myhost(short=True)
        #self.host = 'localhost'

        # monitor channels we are interested in
        self.channels = ['fits']

        ro.init([self.host])

        # for looking up instrument names
        self.insconfig = INSdata()

        # make a name for our monitor
        mymonname = '{}-{}-{}.mon'.format(self.svcname,
                                          self.host,
                                          self.monport)

        threadPool = self.fv.get_threadPool()

        # Create a local pub sub instance
        self.monitor = Monitor.Monitor(mymonname, self.logger,
                                       threadPool=threadPool,
                                       ev_quit=self.fv.ev_quit)

        self.sndsink = SoundSink.SoundSource(monitor=self.monitor,
                                             logger=self.logger,
                                             channels=['sound'])

        self.status_srv = ro.remoteObjectProxy('status')

        # some of the other plugins expect this handle to be available
        # via fv
        self.fv.controller = self

    def start(self):
        self.logger.info("Gen2 plugin starting up...")
        # Startup monitor threadpool
        self.monitor.start(wait=True)

        # Configure logger for logging via our monitor
        ## mymon.logmon(logger, 'monlog', ['logs'])

        # start_server is necessary if we are subscribing, but not if only
        # publishing
        self.logger.info("starting monitor on port {}".format(self.monport))
        self.monitor.start_server(wait=True, port=self.monport)

        # subscribe our monitor to the central monitor hub
        self.monitor.subscribe_remote(self.monitor_name, self.channels, {})
        # publishing for remote command executions
        self.monitor.publish_to('monitor', ['sound', 'g2task'], {})

        if self.settings.get('subscribe_fits', False):
            # Register local fits info subscription callback
            self.monitor.subscribe_cb(self.arr_fitsinfo, ['fits'])
        #self.monitor.subscribe_cb(self.arr_taskinfo, ['taskmgr'])

        # Create our remote service object
        threadPool = self.fv.get_threadPool()
        # methods that can be called from outside via our service
        method_list = ['callGlobalPlugin', 'callGlobalPlugin2',
                       'display_fitsbuf2', 'display_fitsbuf3']
        self.viewsvc = ro.remoteObjectServer(svcname=self.svcname,
                                             obj=self,
                                             logger=self.logger,
                                             ev_quit=self.fv.ev_quit,
                                             port=self.port,
                                             usethread=True,
                                             threadPool=threadPool,
                                             method_list=method_list)

        self.logger.info("starting fits viewing service on port {}".format(self.port))
        self.viewsvc.ro_start()
        self.logger.info("Gen2 plugin startup completed.")

    def stop(self):
        self.logger.info("Gen2 plugin shutting down...")
        self.viewsvc.ro_stop(wait=True)
        self.monitor.stop_server(wait=True)
        self.monitor.stop(wait=True)
        self.logger.info("Gen2 plugin stopped.")


    #############################################################
    #    Here come the Gen2 commands
    #############################################################

    def load_image(self, filepath, idx=None):
        image = self.fv.load_image(filepath, idx=idx)

        # override the name, suppressing "[0]" at the end
        # because this messes up Gen2 operations, which index
        # the image by frameid
        name = image.get('name', None)
        if name is None:
            (dir, filename) = os.path.split(filepath)
            (name, ext) = os.path.splitext(filename)

        match = re.match(r'^(.+)\[0\]$', name)
        if match:
            name = match.group(1)

        image.set(name=name)
        return image

    def load_file(self, filepath, chname=None, wait=True,
                  image_loader=None):
        """Loads a command file from _filepath_ into the commands window.
        """
        try:
            info = iohelper.get_fileinfo(filepath)
            # <-- filepath should now be a real file in the filesystem
            self.logger.debug("fileinfo=%s" % (str(info)))

            image = self.open_fits(info.filepath, channel=chname,
                                   wait=wait,
                                   image_loader=image_loader)
            return image

        except Exception as e:
            errmsg = "Unable to open '%s': %s" % (
                filepath, str(e))
            self.fv.show_error(errmsg)
            return ro.ERROR

    def open_fits(self, filepath, frameid=None, channel=None, wait=False,
                  image_loader=None, display_image=True):

        dirname, filename = os.path.split(filepath)

        if image_loader is None:
            image_loader = self.load_image

        try:
            info = iohelper.get_fileinfo(filepath)
            filepath = info.filepath

            kwdargs = {}
            idx = None
            if info.numhdu is not None:
                idx = max(0, info.numhdu)
                kwdargs['idx'] = idx

            image = image_loader(filepath, **kwdargs)

            assert isinstance(image, BaseImage.BaseImage), \
                   ValueError("Loader did not produce a loadable image: %s" % (
                str(type(image))))

        except Exception as e:
            errmsg = "Failed to load file '%s': %s" % (
                filepath, str(e))
            self.logger.error(errmsg)
            try:
                (etyp, value, tb) = sys.exc_info()
                tb_str = "\n".join(traceback.format_tb(tb))
            except Exception as e:
                tb_str = "Traceback information unavailable."

            self.fv.gui_do(self.fv.show_error, errmsg + '\n' + tb_str)
            return None

        future = Future.Future()
        future.freeze(image_loader, filepath, **kwdargs)

        # Save a future for this image to reload it later if we
        # have to remove it from memory
        image.set(loader=image_loader, image_future=future)

        if image.get('path', None) is None:
            image.set(path=filepath)

        header = image.get_header()

        try:
            if frameid:
                name = frameid
            else:
                (name, ext) = os.path.splitext(filename)
                frameid = header['FRAMEID'].strip()

            chname = self.insconfig.getNameByFrameId(frameid)

        except Exception as e:
            self.logger.warn("Error getting FRAMEID from fits header of %s: %s" % (
                filename, str(e)))
            chname = 'Image'

        if not channel:
            channel = chname

        image.set(name=name, chname=channel)

        # Display image.  If the wait parameter is False then don't wait
        # for the image to load into the viewer
        if wait:
            self.fv.gui_call(self.fv.add_image, name, image, chname=channel)
        else:
            self.fv.gui_do(self.fv.add_image, name, image, chname=channel)

        return image

    ## def gui_load_file(self):
    ##     """Runs dialog to read in a command file into the command window.
    ##     """
    ##     initialdir = os.environ['DATAHOME']

    ##     self.fv.gui_load_file(initialdir=initialdir)


    #############################################################
    # methods that are called from Gen2 observation control
    #############################################################

    def display_fitsbuf2(self, imname, chname, img_buf, shape, dtype,
                         header, metadata):
        """Display a FITS image buffer.  Parameters:
        _imname_: name of the image
        _chname_: name of the channel to display the image
        _img_buf_: ascii encoded numpy containing image data
        _shape_: numpy image shape as a tuple
        _dtype_: numpy data type as a str
        _header_: fits file header as a dictionary
        _metadata_: metadata about image to attach to image
        """

        self.logger.info("received image data: name=%s len=%d shape=%s dtype=%s" % (
            imname, len(img_buf), str(shape), dtype))

        # Decode binary data
        img_buf = ro.binary_decode(img_buf)

        if dtype == '':
            dtype = np.float32
        else:
            # string to actual type
            #dtype = getattr(numpy, dtype)
            pass

        try:
            if metadata.get('compressed', False):
                img_buf = ro.uncompress(img_buf)

            byteswap = metadata.get('byteswap', False)

            # Create image container
            image = AstroImage.AstroImage(logger=self.logger)
            image.load_buffer(img_buf, shape, dtype, byteswap=byteswap,
                              metadata=metadata)

            image.set(name=imname)
            image.update_keywords(header)

        except Exception as e:
            # Some kind of error decoding the value
            self.logger.error("Error creating image data for '%s': %s" % (
                imname, str(e)))
            return ro.ERROR

        # Enqueue image to display datasrc
        self.fv.gui_do(self.fv.add_image, imname, image, chname=chname)

        return ro.OK

    def display_fitsbuf3(self, imname, chname, img_buf, num_hdu,
                         metadata):
        # Decode binary data
        img_buf = ro.binary_decode(img_buf)
        self.logger.info("received image data '%s': len=%d" % (
            imname, len(img_buf)))

        if metadata.get('compressed', False):
            img_buf = ro.uncompress(img_buf)

        in_f = BytesIO(img_buf)

        filepath = None
        if not metadata.get('nocache', False) and self.disp_cache is not None:
            # write out buffer to a file for reloadability
            filepath = os.path.join(self.disp_cache, imname)
            if not filepath.lower().endswith('.fits'):
                filepath = filepath + '.fits'
            self.fv.nongui_do(self._write_cache, filepath, img_buf)

        # Create image container
        self.logger.info("opening buffer with FITS reader")
        with pyfits.open(in_f, 'readonly') as fits_f:

            # Seems to be needed otherwise we sometimes get "unparsable card"
            # errors for broken FITS files
            fits_f.verify('fix')

            image = AstroImage.AstroImage(metadata=metadata,
                                          logger=self.logger)
            image.load_hdu(fits_f[num_hdu], fobj=fits_f)
            image.set(name=imname, path=filepath)
            self.logger.info("image name is '%s'" % image.get('name'))

        # Enqueue image to display datasrc
        self.fv.gui_do(self.fv.add_image, imname, image, chname=chname)

        return ro.OK

    def callGlobalPlugin(self, tag, pluginName, methodName, args, kwdargs):

        self.logger.debug("Command received: plugin=%s method=%s args=%s kwdargs=%s tag=%s" % (
                pluginName, methodName, str(args), str(kwdargs), tag))

        # Get object associated with plugin
        obj = self.fv.gpmon.getPlugin(pluginName)

        # Get method we should call
        if not hasattr(obj, methodName):
            raise Gen2Error("No such method '%s' in plugin object %s" % (
                methodName, pluginName))
        method = getattr(obj, methodName)

        ## # Make a future that will be resolved by the GUI thread
        ## future = Future.Future(data=Bunch.Bunch(tag=tag))
        ## future.freeze(method, *args, **kwdargs)
        ## self.fv.gui_do_future(future)

        ## # Wait for result
        ## return future.wait()
        #return self.fv.gui_call(method, *args, **kwdargs)
        return method(*args, **kwdargs)

    def callGlobalPlugin2(self, tag, pluginName, methodName, args, kwdargs):

        self.logger.debug("Command received: plugin=%s method=%s args=%s kwdargs=%s tag=%s" % (
                pluginName, methodName, str(args), str(kwdargs), tag))

        # Get object associated with plugin
        obj = self.fv.gpmon.getPlugin(pluginName)

        # Get method we should call
        if not hasattr(obj, methodName):
            raise Gen2Error("No such method '%s' in plugin object %s" % (
                methodName, pluginName))
        method = getattr(obj, methodName)

        # Make a future that will be resolved by the GUI thread
        p = Bunch.Bunch()
        future = Future.Future(data=p)
        #future.freeze(method, *args, **kwdargs)
        newargs = [tag, future]
        newargs.extend(list(args))
        future.add_callback('resolved',
                            lambda f: self.fv.nongui_do(self.result_cb2,
                                                        future, tag, p))

        future = Future.Future(data=p)
        future.freeze(method, *newargs, **kwdargs)
        future.add_callback('resolved',
                            lambda f: self.fv.nongui_do(self.result_cb1,
                                                        future, tag, p))

        self.fv.gui_do_future(future)
        return ro.OK

    def result_cb1(self, future, tag, data):
        res = future.get_value(suppress_exception=True)
        if isinstance(res, Exception):
            errmsg = str(res)
            resdata = { 'gui_done': time.time(), 'result': 'error',
                        'errmsg': errmsg }

            self._rpc_cleanse(resdata)

            self.logger.error("Command (%s) terminated by exception: %s" % (
                tag, errmsg))
            self.monitor.setvals(['g2task'], tag, **resdata)
        else:
            self.logger.debug("Command made it to the GUI interaction: %s" % (
                tag))

    def result_cb2(self, future, tag, data):
        self.logger.debug("Command termination: %s" % (tag))
        res = future.get_value(suppress_exception=True)

        if isinstance(res, Exception):
            errmsg = str(res)
            data.update({ 'errmsg': errmsg, 'result': 'error' })
            self.logger.error("Command (%s) terminated by exception: %s" % (
                tag, errmsg))
        else:
            self.logger.debug("Command completed GUI interaction: %s" % (
                tag))

        resdata = { 'gui_done': time.time() }
        resdata.update(data)

        self._rpc_cleanse(resdata)
        self.logger.debug("Result is: %s" % (str(resdata)))

        self.monitor.setvals(['g2task'], tag, **resdata)


    def _rpc_cleanse(self, d):
        # RPC cleansing of return data dictionaries
        for key, val in list(d.items()):
            if isinstance(val, np.integer):
                d[key] = int(val)
            elif isinstance(val, np.floating):
                d[key] = float(val)

    def _write_cache(self, filepath, img_buf):
        self.logger.info("writing cache file: {}".format(filepath))
        try:
            with open(filepath, 'wb') as out_f:
                out_f.write(img_buf)
        except Exception as e:
            self.logger.error("error writing cache file: {}".format(e))

    #############################################################
    # Methods called when we are notified via monitor of a new
    # FITS file.
    #############################################################

    def display_fitsfile(self, fitspath, frameid=None):
        """Display a FITS image file.  Parameters:
        _fitspath_: path of a FITS file.
        """

        try:
            self.logger.debug("Attempting to display '%s'" % (fitspath))
            image = self.open_fits(fitspath, frameid=frameid)

        except IOError as e:
            self.logger.error("Error opening FITS file '%s': %s" % (
                    fitspath, str(e)))
            return ro.ERROR

        return ro.OK

    def file_notify(self, filepath):
        """
        This gets called when we are being notified of a new file.
        """
        self.logger.info("Notified of new file: %s" % (filepath))
        frame = Frame(path=filepath)

        if frame.inscode in ('HSC', 'SUP'):
            # Don't display raw HSC frames; mosaic plugin will display them
            return

        frameid = str(frame)

        try:
            #with self._rlock:
            #self.fv.nongui_do(self.display_fitsfile, filepath, frameid=frameid)
            self.display_fitsfile(filepath, frameid=frameid)

        except Exception as e:
            self.logger.error("Error displaying '%s': %s" % (
                filepath, str(e)))

    def arr_fitsinfo(self, payload, name, channels):
        """Called via the monitor if new information becomes available
        about fits images."""

        self.logger.debug("received values '%s'" % str(payload))

        try:
            bnch = Monitor.unpack_payload(payload)

        except Monitor.MonitorError as e:
            self.logger.error("malformed packet '%s': %s" % (
                str(payload), str(e)))
            return

        # Find out the source of this information by examining the path
        match = regex_frame.match(bnch.path)
        if match:
            frameid = match.group(1)
            vals = bnch.value

            try:
                fitspath = vals['fitspath']
            except KeyError:
                return

            self.fv.nongui_do(self.file_notify, fitspath)
            #self.fv.gui_do(self.file_notify, fitspath)


    def arr_taskinfo(self, payload, name, channels):
        self.logger.debug("received values '%s'" % str(payload))
        try:
            bnch = Monitor.unpack_payload(payload)

        except Monitor.MonitorError as e:
            self.logger.error("malformed packet '%s': %s" % (
                str(payload), str(e)))
            return

        if 'value' not in bnch:
            # delete (vaccuum) packet
            return
        vals = bnch.value

        if 'task_code' in vals:
            res = vals['task_code']
            # Interpret task results:
            #   task_code == 0 --> OK   task_code != 0 --> ERROR
            if res == 3:
                self.logger.info("Task cancelled (%s)" % bnch.path)
                self.cancel_commands(bnch.path)

    #############################################################
    # Methods called by other Gen2 plugins
    #############################################################

    def play_soundfile(self, filepath, format=None, priority=20):
        """Method called by some Gen2 plugins to play sound through
        the Gen2 sound system.
        """
        self.sndsink.playFile(filepath, format=format, priority=priority)

    def get_gen2_status(self, stat_dict):
        return self.status_srv.fetch(stat_dict)

    def __str__(self):
        return 'gen2int'

#END
