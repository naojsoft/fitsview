FitsView
========

ABOUT
-----
This is the FitsView data monitoring/interaction GUI for Subaru Telescope.

USAGE
-----
  $ fitsview --loglevel=20 --log=path/to/fitsview.log

  or

  $ guideview --loglevel=20 --log=path/to/guideview.log

add --stderr if you want to see logging output to the terminal as well.

REQUIRED PACKAGES
-----------------

- pyqt5
- eclipse (qualsize)
- esolib (iqe)
- inotify
- g2cam (pip install git+https://github.com/naojsoft/g2cam)
- g2client (pip install git+https://github.com/naojsoft/g2client)

Packages that will be installed by the installer if not already
installed:

- pyyaml
- matplotlib
- ginga

INSTALLATION
------------

Install all prerequisites, then:

  $ pip install .

AUTHORS
-------
T. Inagaki
R. Kackley
E. Jeschke


