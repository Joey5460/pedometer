#!/usr/bin/env python
"""
setup.py file for SWIG em_math
"""

from distutils.core import setup, Extension


em_math_module = Extension('_em_math',
                           include_dirs = ['../em_lib/'],
                           sources = ['em_math_wrap.c', '../em_lib/em_math.c'],
                           )

setup (name = 'em_math',
       version = '0.1',
       author      = "SWIG Docs",
       description = """em_math swig from docs""",
       ext_modules = [em_math_module],
       py_modules = ["em_math"],
       )
