#
# T. Inagaki
#
import os
import sys
import numpy as np

from g2base import ssdlog

class QuadraticError(Exception):
    pass


class QuadraticFunction(object):

    def __init__(self, logger=None):
        self.a = 0
        self.b = 0
        self.c = 0
        self.coeffs = [0, 0, 0]
        self.logger = logger

    def coefficient(self, x_points, y_points, degree=2):

        try:
            self.coeffs = np.polyfit(x_points, y_points, degree)

        except Exception as e:
            self.logger.error('error: failed to calculate coefficient. %s" %e')
            raise QuadraticError("failed to calculate coefficient. %s" %e )

        else:
            self.a = self.coeffs[0]
            self.b = self.coeffs[1]
            self.c = self.coeffs[2]
            self.logger.debug('coefficient a=%f b=%f c=%f' % (
                self.a, self.b, self.c))

    def _calc_vertex(self):

        try:
            y = self.c - (self.b**2 / (4.0*self.a))
            x = -self.b / (2.0*self.a)
            self.logger.debug("vertex x=%f y=%f" %(x,y))

        except Exception as e:
            self.logger.error('error: failed to calculate vertex. %s" %e')
            raise QuadraticError("failed to calculate vertex. %s" %e )

        else:
            return (x, y)


    def _is_zero(self):

        is_zero = False

        zero = 0.0000000001
        try:
            assert not (-zero <= self.a <= zero)

        except AssertionError as e:
            self.logger.error('error: coefficient A is 0')
            is_zero = True

        return is_zero

    def max_vertex(self):

        self.logger.debug('finding max vertex...')

        if self.a >= 0 or self._is_zero():
            self.logger.error('error: coefficient A is equal or greater than 0.')
            raise QuadraticError("coefficient A is equal or greater than 0.")

        else:
            res = self._calc_vertex()

        return res

    def min_vertex(self):

        self.logger.debug('finding min vertex...')

        #try:
        #    assert (self.a > 0)
        #except Exception as e:
        if self.a <= 0 or self._is_zero():
            self.logger.error('error: coefficient A is equal or less than 0.')
            raise QuadraticError("coefficient A is equal or greater than 0.")

        else:
            res = self._calc_vertex()

        return res

    def quadratic(self):
        try:
            f = np.poly1d(self.coeffs)

        except Exception as e:
            self.logger.error('error: calculate quadratic function. %s' %e)
            raise QuadraticError("failed to calculate quadratic function. %s" %e)

        else:
            self.logger.debug('quadratic function=%s' %(str(f)))
            return f


def main(options,args):

    logname = 'polynomial'
    # Create top level logger.
    logger = ssdlog.make_logger(logname, options)

    # test data points
    points = np.array([(1, -1), (1, 0), (1, 0.1),
                       (2, 2), (2, 1.85), (2, 2.01),
                       (3, 4.5), (3, 5.0), (3, 5.52),
                       (4, 8), (4, 7.85), (4, 8.11),
                       (5, 7), (5, 10.05), (5, 9.11),
                       (6, 8.2), (6, 7.95), (6, 6.01),
                       (7, 4.9), (7, 5.0), (7, 6.02),
                       (8, 2), (8, 1.95), (8, 1.31),
                       (9, 1.9), (9, 2.75), (9, 2.001),])
    # get x and y vectors
    x = points[:,0]
    y = points[:,1]

    qf = QuadraticFunction(logger=logger)
    qf.coefficient(x_points=x, y_points=y)
    func = qf.quadratic()

    vertex = qf.max_vertex()
    logger.debug('max vertex=%s' %(str(vertex)))

    x_new = np.linspace(x[0], x[-1], 100)
    y_new = func(x_new)

    import matplotlib.pyplot as plt
    plt.plot(x, y, 'o', x_new, y_new, '-')
    #plt.plot(x, y, 'o', x_new, y_new, '-', x_new, y_est, '--')
    plt.xlim([x[0]-1, x[-1] + 1 ])
    #plt.plot(x, y, 'o-', x_new, y_est)
    plt.show()


if __name__ == "__main__":

    # Parse command line options with nifty new optparse module
    from optparse import OptionParser

    usage = "usage: %prog [options] [file] ..."
    optprs = OptionParser(usage=usage, version=('%prog'))

    optprs.add_option("--debug", dest="debug", default=False,
                      action="store_true",
                      help="Enter the pdb debugger on main()")

    optprs.add_option("--profile", dest="profile", action="store_true",
                      default=False,
                      help="Run the profiler on main()")

    ssdlog.addlogopts(optprs)
    (options, args) = optprs.parse_args(sys.argv[1:])

    #if len(args) > 0:
    #    optprs.error("incorrect number of arguments")

    # Are we debugging this?
    if options.debug:
        import pdb

        pdb.run('main(options, args)')

    # Are we profiling this?
    elif options.profile:
        import profile

        print("%s profile:" % sys.argv[0])
        profile.run('main(options, args)')

    else:
        main(options, args)
