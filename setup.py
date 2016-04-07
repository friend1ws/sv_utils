#!/usr/bin/env python

from distutils.core import setup

setup(name='sv_utils',
      version='0.3.0beta',
      description='Python programs for analyzing GenomonSV results',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/friend1ws/sv_utils',
      package_dir = {'': 'lib'},
      packages=['sv_utils'],
      scripts=['sv_utils'],
      license='GPL-3'
     )

