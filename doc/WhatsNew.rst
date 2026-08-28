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
