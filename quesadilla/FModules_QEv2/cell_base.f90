MODULE cell_base
    !------------------------------------------------------------------------------!
    !! Cell parameters and initialization.
    !
    !  direct and reciprocal lattice primitive vectors
    double precision :: at(3,3)
    !! The lattice vectors of the simulation cell, a_i, in alat units:
    !! a_i(:) = at(:,i)/alat
    double precision :: bg(3,3)
    !! The reciprocal lattice vectors, b_i, in tpiba=2pi/alat units:
    !! b_i(:) = bg(:,i)/tpiba
END MODULE cell_base