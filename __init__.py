__all__ = ["pytools"]
print("__init__ invoked")
# from analysator.pytools import *
import warnings
warnings.filterwarnings("once", category=DeprecationWarning)
warnings.filterwarnings("once", category=PendingDeprecationWarning)
warnings.filterwarnings("once", category=FutureWarning)


# from pytools import vslvfile

from os import path as __path
root = __path.dirname(__file__)
with open(__path.join(root,'pytools.py'),'r') as f:
    source = f.read()
    exec(source)
    f.close()
    del f
