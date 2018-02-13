import warnings
import time
import functools


def deprecated(func):
    """
    This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emmitted
    when the function is used.
    """
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.simplefilter('always', DeprecationWarning)  # turn off filter
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning,
                      stacklevel=2)
        warnings.simplefilter('default', DeprecationWarning)  # reset filter
        return func(*args, **kwargs)
    return new_func


def timed(func):
    """
    This is a decorator which can be used to time functions.
    When the function has completed it will show the time it
    took in seconds.
    """
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        start = time.clock()
        result = func(*args, **kwargs)
        print "Calculating " + func.__name__ + "() took", time.clock() - start, "seconds"
        return result
    return new_func