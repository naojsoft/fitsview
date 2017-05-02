import os

from ginga.rv.Control import GingaShell

class Gen2FITSViewer(GingaShell):
    """This class manages the creation and handling of a FITS viewer GUI.
    The class is constructed with a data source and it reads images from the
    source and displays them.
    """

    def __init__(self, logger, threadPool, module_manager, preferences,
                 soundsink, ev_quit=None):

        self.controller = None
        self.soundsink = soundsink

        GingaShell.__init__(self, logger, threadPool, module_manager,
                            preferences, ev_quit=ev_quit)


    def load_file(self, filepath, chname=None, wait=True,
                  image_loader=None):
        """Loads a command file from _filepath_ into the commands window.
        """
        try:
            info = self.get_fileinfo(filepath)
            # <-- filepath should now be a real file in the filesystem
            self.logger.debug("fileinfo=%s" % (str(info)))

            image = self.controller.open_fits(info.filepath, channel=chname,
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

        super(Gen2FITSViewer, self).gui_load_file(initialdir=initialdir)
