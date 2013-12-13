#
# Receive.py 
#
# Eric Jeschke (eric@naoj.org)
#
import sys, traceback
import os, re
import time
import pyfits
import numpy
import threading

import remoteObjects as ro
import remoteObjects.Monitor as Monitor
import cfg.INS as INSconfig

from ginga.misc import Bunch, Future
from ginga import AstroImage, RGBImage
from ginga.util import wcs

# Regex used to discover/parse frame info
regex_frame = re.compile(r'^mon\.frame\.(\w+)\.Archiver$')


class ReceiverError(Exception):
    pass

class ReceiveFITS(object):
    """This class implements a simple remoteObject-based fits image receiver.
    Two methods are served: display_fitsfile() and display_fitsbuf()
    """

    def __init__(self, fv, monitor, logger):
        self.fv = fv
        self.monitor = monitor
        self.logger = logger
        # for looking up things about instruments
        self.insconfig = INSconfig.INSdata()

        self._rlock = threading.RLock()
        

    def open_fits(self, filepath, frameid=None, channel=None, wait=False):

        dirname, filename = os.path.split(filepath)
        
        # Create an image.  Assume type to be an AstroImage unless
        # the MIME association says it is something different.
        try:
            typ, subtyp = self.fv.guess_filetype(filepath)
        except Exception:
            # Can't determine file type: assume and attempt FITS
            typ, subtyp = 'image', 'fits'
        
        if (typ == 'image') and (subtyp != 'fits'):
            image = RGBImage.RGBImage(logger=self.logger)
        else:
            image = AstroImage.AstroImage(logger=self.logger)

        try:
            #self.fv.showStatus("Loading %s" % (filename))
            self.logger.debug("Loading file '%s'" % (filename))
            image.load_file(filepath)

        except Exception, e:
            errmsg = "Failed to load file '%s': %s" % (
                filepath, str(e))
            self.logger.error(errmsg)
            try:
                (type, value, tb) = sys.exc_info()
                tb_str = "\n".join(traceback.format_tb(tb))
            except Exception, e:
                tb_str = "Traceback information unavailable."

            self.gui_do(self.fv.show_error, errmsg + '\n' + tb_str)
            return
            
        header = image.get_header()

        try:
            if frameid:
                name = frameid
            else:
                (name, ext) = os.path.splitext(filename)
                frameid = header['FRAMEID'].strip()

            chname = self.insconfig.getNameByFrameId(frameid)
                
        except KeyError:
            self.logger.error("Error getting FRAMEID from fits header: %s" % (
                filename))
            chname = 'Image'

        if not channel:
            channel = chname
            
        image.set(name=name, path=filepath, chname=channel)
        path = image.get('path', 'NO PATH')
        #print "receive: path=%s" % (path)

        # Display image.  If the wait parameter is False then don't wait
        # for the image to load into the viewer
        if wait:
            self.fv.gui_call(self.fv.add_image, name, image, chname=channel)
        else:
            self.fv.gui_do(self.fv.add_image, name, image, chname=channel)

        return image


    def display_fitsbuf(self, fitsname, chname, data, width, height, na_type,
                        header, metadata):
        """Display a FITS image buffer.  Parameters:
        _fitsfile_: name of the file or title for the title bar
        _data_: ascii encoded numpy containing image data
        _width_, _height_: image dimensions in pixels
        _na_type_: numpy data type (currently ignored)
        _header_: fits file header as a dictionary
        _metadata_: metadata about image to attach to image
        """

        # Decode binary data
        data = ro.binary_decode(data)
        print "Received data: len=%d width=%d height=%d" % \
              (len(data), width, height)

        # TODO: recreate pyfits object
        
        # Temporarily coerce numpy type  TODO: fix
        try:
            na_type = numpy.float32
            data = numpy.fromstring(data, dtype=na_type)
            data.byteswap(True)
            data = data.reshape((height, width))
            #print data

        except Exception, e:
            # Some kind of error decoding the value
            self.logger.error("Error creating image data for '%s': %s" % (
                fitsname, str(e)))
            return ro.ERROR

        # Create image container
        image = AstroImage.AstroImage(data, metadata=metadata,
                                      wcsclass=wcs.WCS,
                                      logger=self.logger)
                                      #wcsclass=wcs.BareBonesWCS)
        image.set(name=fitsname)
        image.update_keywords(header)
        
        # Enqueue image to display datasrc
        self.fv.gui_do(self.fv.add_image, fitsname, image,
                            chname=chname)
        #self.fv.gui_do(self.update_image, fitsname, chname, image,
        #               metadata)

        return ro.OK


    def display_fitsfile(self, fitspath, frameid=None):
        """Display a FITS image file.  Parameters:
        _fitspath_: path of a FITS file.
        """

        try:
            image = self.open_fits(fitspath, frameid=frameid)

        except IOError, e:
            self.logger.error("Error opening FITS file '%s': %s" % (
                    fitspath, str(e)))
            return ro.ERROR

        return ro.OK


    ## def executeCmd(self, subsys, tag, cmdName, args, kwdargs):

    ##     self.logger.debug("Command received: subsys=%s command=%s args=%s kwdargs=%s tag=%s" % (
    ##             subsys, cmdName, str(args), str(kwdargs), tag))

    ##     try:
    ##         # Try to look up the start method
    ##         method = getattr(self.fv, cmdName)

    ##     except AttributeError, e:
    ##         ack_msg = "ERROR: No such method: %s" % (cmdName)
    ##         self.logger.error(ack_msg)
    ##         raise Exception(ack_msg)

    ##     future = Future.Future(data=Bunch.Bunch(cmd=cmdName, tag=tag))
    ##     future.freeze(None, *args, **kwdargs)
    ##     future.add_callback('resolved', self.result_cb)
        
    ##     future = self.fv.gui_do(method, *args)
    ##     return ro.OK


    def pluginCmd(self, tag, chname, cmdName, args, kwdargs):

        self.logger.debug("Command received: chname=%s command=%s args=%s kwdargs=%s tag=%s" % (
                chname, cmdName, str(args), str(kwdargs), tag))

        future = Future.Future(data=Bunch.Bunch(tag=tag))
        future.freeze(None, *args, **kwdargs)
        future.add_callback('resolved', self.result_cb)
        
        self.fv.gui_do(self.fv.start_operation_channel,
                            chname, cmdName, future)
        # Wait for result
        return future.wait()

    def play_soundfile(self, filepath, format=None, priority=20):
        self.fv.play_soundfile(filepath, format=format,
                                priority=priority)

    def reloadLocalPlugin(self, plname):
        self.fv.mm.loadModule(plname)
        for chname in self.fv.get_channelNames():
            chinfo = self.fv.get_channelInfo(chname)
            chinfo.opmon.reloadPlugin(plname, chinfo=chinfo)
            
    def reloadGlobalPlugin(self, plname):
        self.fv.mm.loadModule(plname)
        self.fv.gpmon.reloadPlugin(plname)

    def callGlobalPlugin(self, tag, pluginName, methodName, args, kwdargs):

        self.logger.debug("Command received: plugin=%s method=%s args=%s kwdargs=%s tag=%s" % (
                pluginName, methodName, str(args), str(kwdargs), tag))

        # Get object associated with plugin
        obj = self.fv.gpmon.getPlugin(pluginName)

        # Get method we should call
        if not hasattr(obj, methodName):
            raise ReceiverError("No such method '%s' in plugin object %s" % (
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
            raise ReceiverError("No such method '%s' in plugin object %s" % (
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

    def result_cb1(self, future, tag, data):
        res = future.get_value(suppress_exception=True)
        if isinstance(res, Exception):
            errmsg = str(res)
            resdata = { 'gui_done': time.time(), 'result': 'error',
                        'errmsg': errmsg }
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
        self.logger.debug("Result is: %s" % (str(resdata)))

        self.monitor.setvals(['g2task'], tag, **resdata)


    def quit(self):
        self.ev_quit.set()

        return ro.OK
        
        
    def arr_fitsinfo(self, payload, name, channels):
        """Called via the monitor if new information becomes available
        about fits images."""
        
        self.logger.debug("received values '%s'" % str(payload))

        try:
            bnch = Monitor.unpack_payload(payload)

        except Monitor.MonitorError:
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

            try:
                with self._rlock:
                    self.logger.debug("Attempting to display '%s'" % (
                        fitspath))
                    self.display_fitsfile(fitspath, frameid=frameid)

            except Exception, e:
                self.logger.error("Error displaying '%s': %s" % (
                    fitspath, str(e)))

    def arr_taskinfo(self, payload, name, channels):
        self.logger.debug("received values '%s'" % str(payload))
        try:
            bnch = Monitor.unpack_payload(payload)

        except Monitor.MonitorError:
            self.logger.error("malformed packet '%s': %s" % (
                str(payload), str(e)))
            return

        if not bnch.has_key('value'):
            # delete (vaccuum) packet
            return
        vals = bnch.value

        if vals.has_key('task_code'):
            res = vals['task_code']
            # Interpret task results:
            #   task_code == 0 --> OK   task_code != 0 --> ERROR
            if res == 3:
                self.logger.info("Task cancelled (%s)" % bnch.path)
                self.cancel_commands(bnch.path)


#END
