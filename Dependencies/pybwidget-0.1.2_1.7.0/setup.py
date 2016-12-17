from glob import glob
from distutils.core import setup
import os

version = "0.1.2_1.7.0"

_datadir = os.path.join("share", "pybwidget")
prefix = "BWidget-1.7.0"
def files(dir, pat):
    return os.path.join(_datadir, dir), glob(os.path.join(prefix, dir, pat))

setup(
    name="pybwidget", version=version,
    description="BWidget for Tkinter",
    author="Jeff Epler", author_email="jepler@unpy.net",
    packages=['bwidget'],
    data_files = [
          files(".", "*.tcl"),
          files(".", "LICENSE.txt"),
          files("images", "*.gif"),
          files("images", "*.xbm"),
          files("lang", "*.rc"),
    ],
    url="http://tkinter.unpythonic.net/bwidget",
    license="Free to use and distribute (see LICENSE.txt)"
)
