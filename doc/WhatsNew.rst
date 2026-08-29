++++++++++
What's New
++++++++++

Since v0.9.0 (unreleased)
=========================
- ``LeastSquareFits.iqe()`` prefers Ginga's ``IQCalc.iqe()``, which was
  added in Ginga 7.5, and falls back to the ESO ``iqe()`` wrapped by
  ``esolib`` when the installed Ginga predates it.  If neither is
  available it raises, naming both missing pieces.  Both fit a rotated 2D
  elliptical Gaussian and return the major and minor axis FWHM, so the
  measurement is the same either way; the Ginga one is more accurate,
  where the ESO routine reads about 1% low at a FWHM of 4 pixels and 6%
  low at 2.
- ``esolib`` is now imported defensively.  It is not a declared
  dependency, so importing it unconditionally made
  ``fitsview.util.curvefit`` fail to import wherever it had not been
  built -- including for the ``qualsize`` algorithms, which never needed
  it.
- The fallback passes the region dimensions to the ESO ``iqe()`` in the
  order it expects, ``(mx, my)``.  They had been passed the other way
  round, which silently reinterprets the pixel buffer: a 32x48 region
  holding an 8.0 by 3.0 pixel object measured 13.6 by 6.6, and some
  shapes failed outright.  Square regions were unaffected, which is why
  this went unnoticed.
- The ``qualsize`` algorithms report the per-axis FWHM rather than the
  combined value twice.  ``buildDataPoints()`` pairs the two numbers with
  ``CDELT1`` and ``CDELT2`` in ``starsize()``, so passing the combined
  value -- the quadratic mean of the two axes -- for both over-reported
  an elongated object, by about 2.7% at 8:5 and 5% at 8:4.  Round objects
  are unaffected.  The older eclipse-based ``qualsize`` measures only the
  one combined value, so it continues to use it for both.
- ``g2calc`` imports ``eclipse`` leniently, and ``IQCalc.qualsize_old()``
  measures with ``IQCalc.qualsize()`` where that compiled extension is not
  built.  ``eclipse`` is not a declared dependency, so importing it
  unconditionally stopped ``g2calc`` -- and with it the ``VGW``, ``QDAS``,
  ``AgAreaSelection``, ``Region_Selection`` and ``Sv_Drive`` plugins, none
  of which need it for anything else -- from importing at all wherever it
  was missing.  The result carries the same ``x``, ``y``, ``fwhm``,
  ``brightness``, ``objx`` and ``objy`` either way, which is everything
  those plugins read.
- ``qualsize_old()`` takes ``minfwhm``, ``maxfwhm``, ``minelipse`` and
  ``edgew``, defaulting to the limits the SOSS algorithm has compiled in
  (now exported as ``g2calc.MINFWHM`` and friends).  Without this the
  fallback silently applied ``IQCalc``'s much stricter defaults and
  rejected candidates the SOSS routine measures happily -- a 1.5 pixel
  star, and anything more elongated than 2:1, failed outright.  The
  limits reach only the fallback; the SOSS routine cannot be told
  anything but its compiled-in values.
- The ``Sv_Drive`` plugin passes its FWHM, ellipticity and edge controls
  to ``qualsize_old()``, and those controls now default to the SOSS
  limits rather than to ``IQCalc``'s.  The old algorithm's branch had
  ignored them entirely.  Note this makes the new algorithm's branch,
  which shares the same controls, start out permissive as well.
- ``qualsize_old()`` reports ``fwhm_x`` and ``fwhm_y`` on the SOSS path as
  well, so that the result has the same shape whichever algorithm ran.
  ``AgAreaSelection`` and ``Region_Selection`` read those two directly and
  raised ``AttributeError`` -- caught, so the FWHM X, FWHM Y and Star Size
  readouts simply stayed empty -- whenever the old algorithm ran with
  ``eclipse`` installed.  The SOSS routine measures a single combined
  value and does not resolve the two axes, so both fields carry it; for
  non-square pixels that yields the mean of the two plate scales, and the
  IQCalc algorithm is the one that measures the axes properly.
