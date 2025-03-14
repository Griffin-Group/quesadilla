import os
import sys

import numpy.f2py


def main():
    if len(sys.argv) != 2:
        print("Usage: python build_f2py.py <build_dir>")
        sys.exit(1)

    # Get the build directory from the command-line argument
    build_dir = os.path.join(sys.argv[1], "f2py")
    if not os.path.exists(build_dir):
        os.makedirs(build_dir)

    # Define the Fortran source files
    fortran_sources = [
        "quesadilla/FModules_QE/transformations.f90",
        "quesadilla/FModules_QE/symmetrize_fcq.f90",
        "quesadilla/FModules_QE/star_q.f90",
        "quesadilla/FModules_QE/symm_base.f90",
        "quesadilla/FModules_QE/eqvect.f90",
        "quesadilla/FModules_QE/invmat.f90",
        "quesadilla/FModules_QE/little_group.f90",
        "quesadilla/FModules_QE/io.f90",
    ]

    # Prepare the command line arguments for f2py
    comline_list = [
        "--verbose",
        "--lower",
        "--build-dir",
        build_dir,
        "--backend",
        "meson",
        "-m",
        "espresso_symm",
        "--dep",
        "blas",
    ] + fortran_sources

    # Run f2py
    try:
        result = numpy.f2py.run_main(comline_list)
        print(f"f2py completed successfully. Result: {result}")
    except Exception as e:
        print(f"Error during f2py execution: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
