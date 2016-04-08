import nose
import sys 
from numpy.testing import assert_

def run(*args, **kwargs):
    """Test-running function that loads plugins, sets up arguments, and calls `nose.run_exit()`"""
    try:
        kwargs['argv'] = sys.argv + kwargs['argv'] #sys.argv takes precedence
    except KeyError:
        kwargs['argv'] = sys.argv
    return nose.run_exit(*args, **kwargs)
