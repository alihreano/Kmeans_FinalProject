from setuptools import Extension, setup

setup(name='spkmeansmodule',
      version='1.0',
      description='final project ali and jad',
      ext_modules=[Extension('spkmeansmodule', sources=['kmeanspp.c','spkmeansmodule.c', 'utils.c', 'spkmeans.c'])])