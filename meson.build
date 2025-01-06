project('quesadilla',
        ['c', 'fortran'],
        meson_version: '>= 1.1.0',
        default_options : [
                            'warning_level=1',
                            'buildtype=release'
                          ])

# Get Fortran compiler
fc = meson.get_compiler('fortran')

# Setup python
py = import('python').find_installation(pure: false)
py_dep = py.dependency()

# Get numpy and f2py include dirs and setup deps
incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()
inc_np = include_directories(incdir_numpy)
np_dep = declare_dependency(include_directories: inc_np)
inc_np = include_directories(incdir_numpy, incdir_f2py)
incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

# Prepare the f2py files (module and wrappers)
run_command(py, 'tools/run_f2py.py', meson.current_build_dir() ,check : true)
f2py_root = meson.current_build_dir() / 'f2py'

# Install python source files
py.install_sources(['quesadilla/__init__.py',
                    'quesadilla/dynmat.py', 
                    'quesadilla/symmetries.py', 
                    'quesadilla/supercells.py'],
                    subdir : 'quesadilla')

# Install fortran module
py.extension_module('espresso_symm',
                     [
                     'quesadilla/FModules_QE/transformations.f90',
                     'quesadilla/FModules_QE/symmetrize_fcq.f90',
                     'quesadilla/FModules_QE/star_q.f90',
                     'quesadilla/FModules_QE/little_group.f90',
                     'quesadilla/FModules_QE/symm_base.f90',
                     'quesadilla/FModules_QE/eqvect.f90',
                     'quesadilla/FModules_QE/invmat.f90',
                     'quesadilla/FModules_QE/io_global.f90',
                     'quesadilla/FModules_QE/error_handler.f90',
                     f2py_root / 'espresso_symmmodule.c',
                     f2py_root / 'espresso_symm-f2pywrappers2.f90',
                     f2py_root / 'espresso_symm-f2pywrappers.f',
                     fortranobject_c
                     ],
                     include_directories: [inc_np],
                     subdir : 'quesadilla',
                     dependencies : [py_dep],
                     install : true)
