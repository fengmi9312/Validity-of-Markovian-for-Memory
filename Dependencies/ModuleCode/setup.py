from distutils.core import setup, Extension

module = Extension('metapopulation_simulation',
                    sources = ['metapopulation_simulation.cpp', 'metapopulation.cpp', 'model_elem.cpp', 'spreading_func.cpp'],
                    extra_compile_args = ['-std=c++17'])

setup (name = 'PackageName',
       version = '1.0',
       description = 'This is a package for the metapopulation simulation',
       ext_modules = [module,])