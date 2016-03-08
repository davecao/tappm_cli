# -*- coding:utf-8 -*-
from __future__ import division, print_function, absolute_import

import sys
import os
import sysconfig

from jinja2 import Environment, FileSystemLoader

__version__ = '1.0.0'

if sys.version_info[:2] < (2, 7):
    raise Exception('Python version less than 2.7')

try:
    import numpy as np
except ImportError:
    raise ImportError('Numpy is a required package')
else:
    if tuple(map(int, np.__version__.split('.')[:2])) < (1, 8):
        raise ImportError('Numpy v1.8 or later is required.')


_PY3K = PY3K = sys.version_info[0] > 2
PY2K = not PY3K
VERSION = __version__


def is_float(n):
    return isinstance(n, float)

TPLPATH = "{}/templates/".format(os.path.dirname(__file__))
jinja2_ENV = Environment(loader=FileSystemLoader(TPLPATH))
jinja2_ENV.add_extension('jinja2.ext.loopcontrols')
jinja2_ENV.tests['Float'] = is_float
MODELPATH = "{}/models/".format(os.path.dirname(__file__))

data = "{}{}{}{}{}".format(
            sysconfig.get_path('data'),
            os.sep, 'share', os.sep, 'tappm', os.sep, 'data')

__all__ = []

from .dataset_maker import *
from .method_hmm import *

from .utils import *

from .hmm import *
from .io import *
