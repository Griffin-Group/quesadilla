project('quesadilla',
        ['c', 'fortran'],
        meson_version: '>= 1.1.0',
        default_options : [
                            'warning_level=1',
                            'buildtype=release'
                          ])

fc = meson.get_compiler('fortran')

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

# Prepare the f2py files
run_command(py, 'tools/run_f2py.py', meson.current_build_dir() ,check : true)
f2py_root = meson.current_build_dir() / 'f2py'

inc_np = include_directories(incdir_numpy)
np_dep = declare_dependency(include_directories: inc_np)

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

inc_np = include_directories(incdir_numpy, incdir_f2py)

py.install_sources(['quesadilla/__init__.py',
                    'quesadilla/dynmat.py', 
                    'quesadilla/symmetries.py', 
                    'quesadilla/supercells.py'],
                    subdir : 'quesadilla')

py.extension_module('espresso_symm',
                     [
                     'quesadilla/FModules_QE/q2qstar_out.f90',
                     'quesadilla/FModules_QE/trntnsc.f90',
                     'quesadilla/FModules_QE/rotate_and_add_dyn.f90',
                     'quesadilla/FModules_QE/star_q.f90',
                     'quesadilla/FModules_QE/symm_base.f90',
                     'quesadilla/FModules_QE/eqvect.f90',
                     'quesadilla/FModules_QE/invmat.f90',
                     'quesadilla/FModules_QE/sgam_ph.f90',
                     'quesadilla/FModules_QE/smallgq.f90',
                     'quesadilla/FModules_QE/cryst_to_cart.f90',
                     'quesadilla/FModules_QE/io_global.f90',
                     'quesadilla/FModules_QE/error_handler.f90',
                     f2py_root / 'espresso_symmmodule.c',
                     f2py_root / 'espresso_symm-f2pywrappers2.f90',
                     f2py_root / 'espresso_symm-f2pywrappers.f',
                     fortranobject_c
                     ],
                     include_directories: [inc_np],
                     dependencies : [py_dep],
                     subdir: 'quesadilla',
                     install : true)
