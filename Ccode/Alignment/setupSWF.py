from distutils.core import setup, Extension

module1 = Extension('PycSWF', sources = ['PyWatermanFun.c','WatermanFun.c']);

setup(name='PycSWF', version='1.0',  \
      ext_modules=[module1])

