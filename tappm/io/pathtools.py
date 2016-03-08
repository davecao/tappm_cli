# -*- coding: utf-8 -*-

import os
import sys
import zipfile
import os.path
import platform
from os import sep as pathsep
from os.path import isfile, isdir
from os.path import getsize, isabs, exists, abspath
from shutil import copy

PLATFORM = platform.system()

__all__ = []

major, minor = sys.version_info[:2]
if major > 2 and minor <= 3:
    import gzip
    from gzip import GzipFile
    import io
    import bz2
    from bz2 import BZ2File

    class TextIOWrapper(io.TextIOWrapper):

        def _getlines(self):

            try:
                lines = self._lines
            except AttributeError:
                self._lines = None

            if self._lines is None:
                self._lines = self.read().split('\n')
            return self._lines

        def readline(self, *args):

            lines = self._getlines()
            if lines:
                return lines.pop(0)
            else:
                return ''

        def readlines(self, size=None):

            lines = self._getlines()
            if size is None:
                self._lines = []
                return lines
            else:
                self._lines = lines[size:]
                return lines[:size]

        def __del__(self):

            self.close()

    def gzip_open(filename, mode="rb", compresslevel=9,
                  encoding=None, errors=None, newline=None):
        """Open a gzip-compressed file in binary or text mode.

        The filename argument can be an actual filename (a str or bytes
        object), or an existing file object to read from or write to.

        The mode argument can be "r", "rb", "w", "wb", "a" or "ab" for binary
        mode, or "rt", "wt" or "at" for text mode. The default mode is "rb",
        and the default compresslevel is 9.

        For binary mode, this function is equivalent to the GzipFile
        constructor:

        GzipFile(filename, mode, compresslevel). In this case, the encoding,
        errors and newline arguments must not be provided.

        For text mode, a GzipFile object is created, and wrapped in an
        io.TextIOWrapper instance with the specified encoding, error handling
        behavior, and line ending(s).

        """
        if "t" in mode:
            if "b" in mode:
                raise ValueError("Invalid mode: %r" % (mode,))
        else:
            if encoding is not None:
                raise ValueError(
                    "Argument 'encoding' not supported in binary mode")
            if errors is not None:
                raise ValueError(
                    "Argument 'errors' not supported in binary mode")
            if newline is not None:
                raise ValueError(
                    "Argument 'newline' not supported in binary mode")

        gz_mode = mode.replace("t", "")
        if isinstance(filename, (str, bytes)):
            binary_file = GzipFile(filename, gz_mode, compresslevel)
        elif hasattr(filename, "read") or hasattr(filename, "write"):
            binary_file = GzipFile(None, gz_mode, compresslevel, filename)
        else:
            raise TypeError(
                "filename must be a str or bytes object, or a file")

        if "t" in mode:
            return TextIOWrapper(binary_file, encoding, errors, newline)
        else:
            return binary_file

    def bz2_open(filename, mode="rb", compresslevel=9,
                 encoding=None, errors=None, newline=None):
        """ Open bz2 compressed file """
        if "t" in mode:
            if "b" in mode:
                raise ValueError("Invalid mode: %r" % (mode,))
        else:
            if encoding is not None:
                raise ValueError(
                    "Argument 'encoding' not supported in binary mode")
            if errors is not None:
                raise ValueError(
                    "Argument 'errors' not supported in binary mode")
            if newline is not None:
                raise ValueError(
                    "Argument 'newline' not supported in binary mode")
        bz_mode = mode.replace("t", "")
        if isinstance(filename, (str, bytes)):
            binary_file = BZ2File(filename, bz_mode, compresslevel)
        elif hasattr(filename, "read") or hasattr(filename, "write"):
            binary_file = BZ2File(None, bz_mode, compresslevel, filename)
        else:
            raise TypeError(
                "filename must be a str or bytes object, or a file")

        if "t" in mode:
            return TextIOWrapper(binary_file, encoding, errors, newline)
        else:
            return binary_file
else:
    import gzip
    import bz2

    def gzip_open(filename, *args, **kwargs):
        if args and "t" in args[0]:
            args = (args[0].replace("t", ""), ) + args[1:]
        return gzip.GzipFile(filename, *args, **kwargs)

    def bz2_open(filename, *args, **kwargs):
        if args and "t" in args[0]:
            args = (args[0].replace("t", ""), ) + args[1:]
        return bz2.BZ2File(filename, *args, **kwargs)

if (major, minor) >= (3, 2):
    from gzip import compress as gzip_compress
    from gzip import decompress as gzip_decompress

OPEN = {
    'bz2': bz2_open,
    'gz': gzip_open,
    'zip': zipfile.ZipFile,
}


def isExecutable(path):
    """Return true if *path* is an executable."""

    return (isinstance(path, str) and exists(path) and
            os.access(path, os.X_OK))


def isReadable(path):
    """Return true if *path* is readable by the user."""

    return (isinstance(path, str) and exists(path) and
            os.access(path, os.R_OK))


def isWritable(path):
    """Return true if *path* is writable by the user."""

    return (isinstance(path, str) and exists(path) and
            os.access(path, os.W_OK))


def relpath(path):
    """Return *path* on Windows, and relative path elsewhere."""

    if PLATFORM == 'Windows':
        return path
    else:
        return os.path.relpath(path)


def sympath(path, beg=2, end=1, ellipsis='...'):
    """Return a symbolic path for a long *path*, by replacing folder names
    in the middle with *ellipsis*.  *beg* and *end* specified how many folder
    (or file) names to include from the beginning and end of the path."""

    abs_items = abspath(path).split(pathsep)
    rel_items = relpath(path).split(pathsep)
    if len(abs_items) <= len(rel_items):
        items = abs_items
    else:
        items = rel_items
    if len(items) <= beg + end:
        return pathsep.join(items)
    else:
        return pathsep.join(items[:beg+1] + [ellipsis] + items[-end:])


def makePath(path):
    """Make all directories that does not exist in a given *path*."""

    if not isdir(path):
        dirs = path.split(pathsep)
        for i, dirname in enumerate(dirs):
            if not dirname:
                continue
            dirname = pathsep.join(dirs[:i+1])
            try:
                if not isdir(dirname):
                    os.mkdir(dirname)
            except OSError:
                raise OSError('{0} could not be created, please '
                              'specify another path'.format(path))
    return path


def which(program):
    """This function is based on the example in:
    http://stackoverflow.com/questions/377017/"""

    fpath, fname = os.path.split(program)
    if fpath and isExecutable(program):
        return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = os.path.join(path, program)
            if isExecutable(path):
                return path
    return None
