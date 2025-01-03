!
! Copyright (C) 2010-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------------
!
module symm_base
  !
  !! This module contains the variables needed to describe the symmetry properties
  !! and the routines to find crystal symmetries.
  !
  !
  ! ... these are acceptance criteria
  !
  double precision, parameter :: eps1 = 1.0d-6, eps2 = 1.0d-5
  !
  save
  !
  private
  !
  ! ... Exported variables
  !
  public :: s, sr, sname, ft, nrot, nsym, nsym_ns, nsym_na, t_rev, &
            no_t_rev, time_reversal, irt, invs, invsym, d1, d2, d3, &
            allfrac, nofrac, nosym, nosym_evc, fft_fact, spacegroup, &
            chem_symb
  !
  integer :: s(3, 3, 48)
  !! symmetry matrices, in crystal axis
  integer :: invs(48)
  !! index of inverse operation: S^{-1}_i=S(invs(i))
  integer :: fft_fact(3) = 1
  !! FFT dimensions must be multiple of fft_fact
  integer :: nrot
  !! number of bravais lattice symmetries
  integer :: spacegroup = 0
  !! space group index, as read from input
  integer :: nsym = 1
  !! total number of crystal symmetries
  integer :: nsym_ns = 0
  !! nonsymmorphic (fractional translation) symms
  integer :: nsym_na = 0
  !! excluded nonsymmorphic symmetries because
  !! fract. transl. is noncommensurate with FFT grid
  double precision :: ft(3, 48)
  !! fractional translations, in crystal axis
  double precision :: sr(3, 3, 48)
  !! symmetry matrices, in cartesian axis
  double precision :: accep = 1.0d-5
  !! initial value of the acceptance threshold
  !! for position comparison by eqvect in checksym
  !
  character(LEN=45) :: sname(48)
  !! name of the symmetries
  integer :: t_rev(48) = 0
  !! time reversal flag, for noncolinear magnetism
  integer, allocatable :: irt(:, :)
  !! symmetric atom for each atom and sym.op.
  logical :: time_reversal = .true.
  !! if .TRUE. the system has time reversal symmetry
  logical :: invsym
  !! if .TRUE. the system has inversion symmetry
  logical :: nofrac = .false.
  !! if .TRUE. fract. translations are not allowed
  logical :: allfrac = .false.
  !! if .TRUE. all fractionary translations allowed,
  !! even those not commensurate with FFT grid
  logical :: nosym = .false.
  !! if .TRUE. no symmetry is used
  logical :: nosym_evc = .false.
  !! if .TRUE. symmetry is used only to symmetrize
  !! k points
  logical :: no_t_rev = .false.
  !! if .TRUE. remove the symmetries that
  !! require time reversal
  double precision, target :: d1(3, 3, 48)
  !! matrix d1 for rotation of spherical harmonics (d1 for l=1, ...)
  double precision, target :: d2(5, 5, 48)
  !! matrix d2 for rotation of spherical harmonics
  double precision, target :: d3(7, 7, 48)
  !! matrix d3 for rotation of spherical harmonics

  !!!!!!!!!!!!!!!!
  integer, parameter :: ntypx = 10
  !! max number of different types of atom
  character(LEN=6)      :: atm(ntypx)
  !     atm( j )  = name of the type of the j-th atomic specie
  integer :: colin_mag = -1
  !! equal to 0 if the system does not have a collinear magnetism
  !! equal to -1 if the collinearity is not checked.
  !! larger than 0 if the system has a collinear magnetism (nspin_mag = 2)
  !! equal to 1 if the symmetries with time-reversal is not detected
  !! equal to 2 if the symmetries with time-reversal is detected
  !!!!!!!!!!!!!!!!!

  !
  ! ... Exported routines
  !
  public ::  find_sym, inverse_s, copy_sym, checkallsym, &
            s_axis_to_cart, set_sym, set_sym_bl, check_grid_sym
  public ::  find_sym_ifc ! FIXME: should be merged with find_sym
  public ::  remove_sym   ! FIXME: is this still useful?

  !
contains
  !
  !-----------------------------------------------------------------------
  subroutine inverse_s()
    !-----------------------------------------------------------------------
     !! Locate index of \(S^{-1}\).
    !
    implicit none
    !
    integer :: isym, jsym, ss(3, 3)
    logical :: found
    !
    do isym = 1, nsym
      found = .false.
      do jsym = 1, nsym
        !
        ss = matmul(s(:, :, jsym), s(:, :, isym))
        ! s(:,:,1) is the identity
        if (all(s(:, :, 1) == ss(:, :))) then
          invs(isym) = jsym
          found = .true.
        end if
        !
      end do
      if (.not. found) call errore('inverse_s', ' Not a group', 1)
    end do
    !
  end subroutine inverse_s
  !
  !
  !-----------------------------------------------------------------------
  subroutine set_sym_bl(at)
    !---------------------------------------------------------------------
     !! Provides symmetry operations for all bravais lattices.
     !! Tests the 24 proper rotations for the cubic lattice first, then
     !! the 8 rotations specific for the hexagonal axis (special axis c),
     !! then inversion is added.
    !
    !use matrix_inversion
    !
    implicit none
    !
    double precision, intent(IN) :: at(3, 3)
    character(LEN=6), external :: int_to_char
    !
    ! sin3 = sin(pi/3), cos3 = cos(pi/3), msin3 = -sin(pi/3), mcos3 = -cos(pi/3)
    !
    double precision, parameter :: sin3 = 0.866025403784438597d0, cos3 = 0.5d0, &
      msin3 = -0.866025403784438597d0, mcos3 = -0.5d0
    !
    double precision :: s0(3, 3, 32), overlap(3, 3), rat(3), rot(3, 3), value
    ! s0: the s matrices in cartesian axis
    ! overlap: inverse overlap matrix between direct lattice
    ! rat: the rotated of a direct vector ( cartesian )
    ! rot: the rotated of a direct vector ( crystal axis )
    ! value: component of the s matrix in axis basis
    integer :: jpol, kpol, mpol, irot, imat(32)
    ! counters over the polarizations and the rotations
    !
    character(LEN=45) :: s0name(64)
    ! full name of the rotational part of each symmetry operation
    !
    data s0/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, &
      -1.d0, 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, 0.d0, 0.d0, 1.d0, &
      -1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, -1.d0, &
      1.d0, 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, 0.d0, 0.d0, -1.d0, &
      0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, -1.d0, &
      0.d0, -1.d0, 0.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, -1.d0, &
      0.d0, -1.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, &
      0.d0, 1.d0, 0.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, &
      0.d0, 0.d0, 1.d0, 0.d0, -1.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
      0.d0, 0.d0, -1.d0, 0.d0, -1.d0, 0.d0, -1.d0, 0.d0, 0.d0, &
      0.d0, 0.d0, -1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
      0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, -1.d0, 0.d0, 0.d0, &
      -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, &
      -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, -1.d0, 0.d0, &
      1.d0, 0.d0, 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, 1.d0, 0.d0, &
      1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, -1.d0, 0.d0, &
      0.d0, 0.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, &
      0.d0, 0.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, &
      0.d0, 0.d0, -1.d0, 1.d0, 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, &
      0.d0, 0.d0, 1.d0, -1.d0, 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, &
      0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 1.d0, 0.d0, 0.d0, &
      0.d0, -1.d0, 0.d0, 0.d0, 0.d0, -1.d0, 1.d0, 0.d0, 0.d0, &
      0.d0, -1.d0, 0.d0, 0.d0, 0.d0, 1.d0, -1.d0, 0.d0, 0.d0, &
      0.d0, 1.d0, 0.d0, 0.d0, 0.d0, -1.d0, -1.d0, 0.d0, 0.d0, &
      cos3, sin3, 0.d0, msin3, cos3, 0.d0, 0.d0, 0.d0, 1.d0, &
      cos3, msin3, 0.d0, sin3, cos3, 0.d0, 0.d0, 0.d0, 1.d0, &
      mcos3, sin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0, 1.d0, &
      mcos3, msin3, 0.d0, sin3, mcos3, 0.d0, 0.d0, 0.d0, 1.d0, &
      cos3, msin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
      cos3, sin3, 0.d0, sin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
      mcos3, msin3, 0.d0, msin3, cos3, 0.d0, 0.d0, 0.d0, -1.d0, &
      mcos3, sin3, 0.d0, sin3, cos3, 0.d0, 0.d0, 0.d0, -1.d0/
    !
    data s0name/'identity                                     ', &
      '180 deg rotation - cart. axis [0,0,1]        ', &
      '180 deg rotation - cart. axis [0,1,0]        ', &
      '180 deg rotation - cart. axis [1,0,0]        ', &
      '180 deg rotation - cart. axis [1,1,0]        ', &
      '180 deg rotation - cart. axis [1,-1,0]       ', &
      ' 90 deg rotation - cart. axis [0,0,-1]       ', &
      ' 90 deg rotation - cart. axis [0,0,1]        ', &
      '180 deg rotation - cart. axis [1,0,1]        ', &
      '180 deg rotation - cart. axis [-1,0,1]       ', &
      ' 90 deg rotation - cart. axis [0,1,0]        ', &
      ' 90 deg rotation - cart. axis [0,-1,0]       ', &
      '180 deg rotation - cart. axis [0,1,1]        ', &
      '180 deg rotation - cart. axis [0,1,-1]       ', &
      ' 90 deg rotation - cart. axis [-1,0,0]       ', &
      ' 90 deg rotation - cart. axis [1,0,0]        ', &
      '120 deg rotation - cart. axis [-1,-1,-1]     ', &
      '120 deg rotation - cart. axis [-1,1,1]       ', &
      '120 deg rotation - cart. axis [1,1,-1]       ', &
      '120 deg rotation - cart. axis [1,-1,1]       ', &
      '120 deg rotation - cart. axis [1,1,1]        ', &
      '120 deg rotation - cart. axis [-1,1,-1]      ', &
      '120 deg rotation - cart. axis [1,-1,-1]      ', &
      '120 deg rotation - cart. axis [-1,-1,1]      ', &
      ' 60 deg rotation - cryst. axis [0,0,1]       ', &
      ' 60 deg rotation - cryst. axis [0,0,-1]      ', &
      '120 deg rotation - cryst. axis [0,0,1]       ', &
      '120 deg rotation - cryst. axis [0,0,-1]      ', &
      '180 deg rotation - cryst. axis [1,-1,0]      ', &
      '180 deg rotation - cryst. axis [2,1,0]       ', &
      '180 deg rotation - cryst. axis [0,1,0]       ', &
      '180 deg rotation - cryst. axis [1,1,0]       ', &
      'inversion                                    ', &
      'inv. 180 deg rotation - cart. axis [0,0,1]   ', &
      'inv. 180 deg rotation - cart. axis [0,1,0]   ', &
      'inv. 180 deg rotation - cart. axis [1,0,0]   ', &
      'inv. 180 deg rotation - cart. axis [1,1,0]   ', &
      'inv. 180 deg rotation - cart. axis [1,-1,0]  ', &
      'inv.  90 deg rotation - cart. axis [0,0,-1]  ', &
      'inv.  90 deg rotation - cart. axis [0,0,1]   ', &
      'inv. 180 deg rotation - cart. axis [1,0,1]   ', &
      'inv. 180 deg rotation - cart. axis [-1,0,1]  ', &
      'inv.  90 deg rotation - cart. axis [0,1,0]   ', &
      'inv.  90 deg rotation - cart. axis [0,-1,0]  ', &
      'inv. 180 deg rotation - cart. axis [0,1,1]   ', &
      'inv. 180 deg rotation - cart. axis [0,1,-1]  ', &
      'inv.  90 deg rotation - cart. axis [-1,0,0]  ', &
      'inv.  90 deg rotation - cart. axis [1,0,0]   ', &
      'inv. 120 deg rotation - cart. axis [-1,-1,-1]', &
      'inv. 120 deg rotation - cart. axis [-1,1,1]  ', &
      'inv. 120 deg rotation - cart. axis [1,1,-1]  ', &
      'inv. 120 deg rotation - cart. axis [1,-1,1]  ', &
      'inv. 120 deg rotation - cart. axis [1,1,1]   ', &
      'inv. 120 deg rotation - cart. axis [-1,1,-1] ', &
      'inv. 120 deg rotation - cart. axis [1,-1,-1] ', &
      'inv. 120 deg rotation - cart. axis [-1,-1,1] ', &
      'inv.  60 deg rotation - cryst. axis [0,0,1]  ', &
      'inv.  60 deg rotation - cryst. axis [0,0,-1] ', &
      'inv. 120 deg rotation - cryst. axis [0,0,1]  ', &
      'inv. 120 deg rotation - cryst. axis [0,0,-1] ', &
      'inv. 180 deg rotation - cryst. axis [1,-1,0] ', &
      'inv. 180 deg rotation - cryst. axis [2,1,0]  ', &
      'inv. 180 deg rotation - cryst. axis [0,1,0]  ', &
      'inv. 180 deg rotation - cryst. axis [1,1,0]  '/
    !
    ! ... compute the overlap matrix for crystal axis
    do jpol = 1, 3
      do kpol = 1, 3
        rot(kpol, jpol) = at(1, kpol)*at(1, jpol) + &
                          at(2, kpol)*at(2, jpol) + &
                          at(3, kpol)*at(3, jpol)
      end do
    end do
    !
    ! ... then its inverse (rot is used as work space)
    call invmat(3, rot, overlap)
    !
    nrot = 1
    !
    do irot = 1, 32
      !
      ! ... for each possible symmetry
      do jpol = 1, 3
        do mpol = 1, 3
          !
          ! ... compute, in cartesian coordinates the rotated vector
          rat(mpol) = s0(mpol, 1, irot)*at(1, jpol) + &
                      s0(mpol, 2, irot)*at(2, jpol) + &
                      s0(mpol, 3, irot)*at(3, jpol)
        end do

        do kpol = 1, 3
          !
          ! ... the rotated vector is projected on the direct lattice
          rot(kpol, jpol) = at(1, kpol)*rat(1) + &
                            at(2, kpol)*rat(2) + &
                            at(3, kpol)*rat(3)
        end do
      end do
      !
      ! ... and the inverse of the overlap matrix is applied
      do jpol = 1, 3
        do kpol = 1, 3
          value = overlap(jpol, 1)*rot(1, kpol) + &
                  overlap(jpol, 2)*rot(2, kpol) + &
                  overlap(jpol, 3)*rot(3, kpol)
          !
          if (abs(dble(nint(value)) - value) > eps1) then
            !
            ! ... if a noninteger is obtained, this implies that this operation
            ! is not a symmetry operation for the given lattice
            !
            goto 10
          end if
          !
          s(kpol, jpol, nrot) = nint(value)
        end do
      end do
      !
      sname(nrot) = s0name(irot)
      imat(nrot) = irot
      nrot = nrot + 1
      !
10    continue
      !
    end do
    !
    nrot = nrot - 1
    !
    if (nrot /= 1 .and. nrot /= 2 .and. nrot /= 4 .and. nrot /= 6 .and. &
        nrot /= 8 .and. nrot /= 12 .and. nrot /= 24) then
      write (*, '(80("-"),/,"NOTICE: Bravais lattice has wrong number (",&
     & i2,") of symmetries - symmetries are disabled",/,80("-"))') nrot
      nrot = 1
    end if
    !
    ! ... set the inversion symmetry (Bravais lattices have always inversion symmetry)
    do irot = 1, nrot
      sname(irot + nrot) = s0name(imat(irot) + 32)
      do kpol = 1, 3
        do jpol = 1, 3
          s(kpol, jpol, irot + nrot) = -s(kpol, jpol, irot)
        end do
      end do
    end do
    !
    nrot = 2*nrot
    !
    ! ... reset fractional translations to zero before checking the group
    ft(:, :) = 0.0d0
    if (.not. is_group(nrot)) then
      ! ... This happens for instance for an hexagonal lattice with one axis
      ! oriented at 15 degrees from the x axis, the other along (-1,1,0)
      ! --> TEMPORARILY COMMENTED OUT BY OMAR
      CALL infomsg( 'set_sym_bl', 'NOTICE: Symmetry group for Bravais lattice &
                   &is not a group (' // TRIM(int_to_char(nrot)) // &
                    ') - symmetries are disabled' )
      nrot = 1
    end if
    !
    return
    !
  end subroutine set_sym_bl
  !
  !
  !-----------------------------------------------------------------------
  subroutine find_sym(at, bg, nat, tau, ityp, magnetic_sym, m_loc, no_z_inv)
    !-----------------------------------------------------------------------
     !! This routine finds the point group of the crystal, by eliminating
     !! the symmetries of the Bravais lattice which are not allowed
     !! by the atomic positions (or by the magnetization if present).
    !
    implicit none
    !
    double precision, intent(IN) :: at(3, 3), bg(3, 3)
    integer, intent(IN) :: nat
     !! total number of atoms of all species
    integer, intent(IN) :: ityp(nat)
     !! the type of each i-th atom in stdin
    double precision, intent(IN) :: tau(3, nat)
     !! atomic positions
    double precision, intent(IN) :: m_loc(3, nat)
     !! local integrated magnetization
    logical, intent(IN) :: magnetic_sym
     !! magnetic_sym = noncolin .AND. domag
    logical, intent(IN), optional :: no_z_inv
     !! if .TRUE., disable symmetries sending z into -z.
     !! Some calculations (e.g. gate fields) require this.
    !
    ! ... local variables
    !
    integer :: i
    logical :: sym(48)
    ! if true the corresponding operation is a symmetry operation
    !
    if (.not. allocated(irt)) allocate (irt(48, nat))
    irt(:, :) = 0
    !
    !    Here we find the true symmetries of the crystal
    !
    symm: do i = 1, 3 !emine: if it is not resolved in 3 steps it is sth else?
      if (present(no_z_inv)) then
        call sgam_at(at, bg, nat, tau, ityp, sym, no_z_inv)
      else
        call sgam_at(at, bg, nat, tau, ityp, sym)
      end if
      !
      ! ... Here we check for magnetic symmetries
      if (magnetic_sym) then
        call sgam_at_mag(at, bg, nat, m_loc, sym)
        ! ... Here we check for time reversal symmetries for collinear systems
        ! NOTE: This check should be performed in the consistent way as in setup.f90
        ! However, we temporarily use this way not to change the interface
        ! until the structure of the code is fixed.
      else if (colin_mag >= 1) then
        call sgam_at_collin(nat, m_loc, sym)
        ! ... If nosym_evc is true from now on we do not use the symmetry any more
      end if
      !
      if (nosym_evc) then
        sym = .false.
        sym(1) = .true.
      end if
      !
      ! ... Here we re-order all rotations in such a way that true sym.ops
      ! are the first nsym; rotations that are not sym.ops. follow
      nsym = copy_sym(nrot, sym)
      !
      if (.not. is_group(nsym)) then
        if (i == 1) call infomsg('find_sym', &
                                 'Not a group! Trying with lower acceptance parameter...')
        accep = accep*0.5d0
        if (i == 3) then
          call infomsg('find_sym', 'Still not a group! symmetry disabled')
          nsym = 1
        end if
        cycle symm
      else
        if (i > 1) call infomsg('find_sym', 'Symmetry operations form a group')
        exit symm
      end if
    end do symm
    !
    ! ... check if inversion (I) is a symmetry.
    ! If so, it should be the (nsym/2+1)-th operation of the group
    !
    invsym = all(s(:, :, nsym/2 + 1) == -s(:, :, 1))
    !
    call inverse_s()
    !
    call s_axis_to_cart(at, bg)
    !
    return
    !
  end subroutine find_sym
  !
  !
  !-----------------------------------------------------------------------
  subroutine sgam_at(at, bg, nat, tau, ityp, sym, no_z_inv)
    !-----------------------------------------------------------------------
     !! Given the point group of the Bravais lattice, this routine finds
     !! the subgroup which is the point group of the considered crystal.
     !! Non symmorphic groups are allowed, provided that fractional
     !! translations are allowed (nofrac=.false) and that the unit cell
     !! is not a supercell.
    !
     !! On output, the array sym is set to .TRUE.. for each operation
     !! of the original point group that is also a symmetry operation
     !! of the crystal symmetry point group.
    !
    implicit none
    !
    double precision, intent(IN) :: at(3, 3), bg(3, 3)
    integer, intent(IN) :: nat
     !! number of atoms in the unit cell
    integer, intent(IN) :: ityp(nat)
     !! species of each atom in the unit cell
    double precision, intent(IN) :: tau(3, nat)
     !! cartesian coordinates of the atoms
    logical, intent(IN), optional :: no_z_inv
     !! if .TRUE., disable symmetry operations sending z into -z.
     !! Some calculations (e.g. gate fields) require this
    logical, intent(OUT) :: sym(48)
     !! flag indicating if sym.op. isym in the parent group
     !! is a true symmetry operation of the crystal.
    !
    ! ... local variables
    !
    integer :: na, kpol, nb, irot, i, j
    ! counters
    double precision, allocatable :: xau(:, :), rau(:, :)
    ! atomic coordinates in crystal axis
    logical :: fractional_translations
    integer :: nfrac
    double precision :: ft_(3), ftaux(3)
    !
    allocate (xau(3, nat))
    allocate (rau(3, nat))
    !
    ! ... Compute the coordinates of each atom in the basis of
    ! the direct lattice vectors
    do na = 1, nat
      xau(:, na) = bg(1, :)*tau(1, na) + bg(2, :)*tau(2, na) + bg(3, :)*tau(3, na)
    end do
    !
    ! ... check if the identity has fractional translations (this means
    ! that the cell is actually a supercell). When this happens, fractional
    ! translations are disabled, because there is no guarantee that the
    ! generated sym.ops. form a group.
    !
    nb = 1
    irot = 1
    !
    fractional_translations = .not. nofrac
    !
    if (fractional_translations) then
      do na = 2, nat
        if ((colin_mag >= 0 .and. chem_symb(atm(ityp(nb))) == chem_symb(atm(ityp(na)))) &
            .or. (colin_mag < 0 .and. ityp(nb) == ityp(na))) then
          !IF ( ityp(nb) == ityp(na) ) THEN
          !
          ft_(:) = xau(:, na) - xau(:, nb) - nint(xau(:, na) - xau(:, nb))
          sym(irot) = checksym(irot, nat, ityp, xau, xau, ft_)
          if (sym(irot)) then
            fractional_translations = .false.
            write (*, '(5x,"Found identity + (",&
           &   3f8.4, ") symmetry",/,5x,"This is a supercell,", &
           &   " fractional translations are disabled")') ft_
            goto 10
          end if
          !
        end if
      end do
    end if
    !
10  continue
    !
    nsym_ns = 0
    fft_fact(:) = 1
    !
    do irot = 1, nrot
      !
      do na = 1, nat
        ! rau = rotated atom coordinates
        rau(:, na) = s(1, :, irot)*xau(1, na) + &
                     s(2, :, irot)*xau(2, na) + &
                     s(3, :, irot)*xau(3, na)
      end do
      !
      ! ... first attempt: no fractional translation
      ft(:, irot) = 0
      ft_(:) = 0.d0
      !
      sym(irot) = checksym(irot, nat, ityp, xau, rau, ft_)
      !
      if (.not. sym(irot) .and. fractional_translations) then
        nb = 1
        do na = 1, nat
          if ((colin_mag >= 0 .and. chem_symb(atm(ityp(nb))) == chem_symb(atm(ityp(na)))) &
              .or. (colin_mag < 0 .and. ityp(nb) == ityp(na))) then
            !IF ( ityp(nb) == ityp(na) ) THEN
            !
            ! ... second attempt: check all possible fractional translations
            ft_(:) = rau(:, na) - xau(:, nb) - nint(rau(:, na) - xau(:, nb))
            !
            ! ... ft_ is in crystal axis and is a valid fractional translation
            ! only if ft_(i)=0 or ft_(i)=1/n, with n=2,3,4,6
            !
            do i = 1, 3
              if (abs(ft_(i)) > eps2) then
                ftaux(i) = abs(1.0d0/ft_(i) - nint(1.0d0/ft_(i)))
                nfrac = nint(1.0d0/abs(ft_(i)))
                if (ftaux(i) < eps2 .and. nfrac /= 2 .and. &
                    nfrac /= 3 .and. nfrac /= 4 .and. nfrac /= 6) &
                  ftaux(i) = 2*eps2
              else
                ftaux(i) = 0.0d0
              end if
            end do
            !
            if (any(ftaux(:) > eps2)) cycle
            !
            sym(irot) = checksym(irot, nat, ityp, xau, rau, ft_)
            !
            if (sym(irot)) then
              nsym_ns = nsym_ns + 1
              ft(:, irot) = ft_(:)
              !
              ! ... Find factors that must be present in FFT grid dimensions
              ! in order to ensure that fractional translations are
              ! commensurate with FFT grids.
              do i = 1, 3
                if (abs(ft_(i)) > eps2) then
                  nfrac = nint(1.0d0/abs(ft_(i)))
                else
                  nfrac = 0
                end if
                fft_fact(i) = mcm(fft_fact(i), nfrac)
              end do
              !
              goto 20
            end if
          end if
        end do
        !
      end if
      !
20    continue
      !
    end do
    !
    ! ... disable all symmetries z -> -z
    if (present(no_z_inv)) then
      if (no_z_inv) then
        do irot = 1, nrot
          if (s(3, 3, irot) == -1) sym(irot) = .false.
        end do
      end if
    end if
    !
    ! ... deallocate work space
    deallocate (rau)
    deallocate (xau)
    !
    return
    !
  end subroutine sgam_at
  !
  !
  !-----------------------------------------------------------------------
  subroutine sgam_at_mag(at, bg, nat, m_loc, sym)
    !-----------------------------------------------------------------------
     !! Find magnetic symmetries, i.e. point-group symmetries that are
     !! also symmetries of the local magnetization - including
     !! rotation + time reversal operations.
    !
    implicit none
    !
    double precision, intent(IN) :: at(3, 3), bg(3, 3)
    integer, intent(IN) :: nat
     !! numbero of atoms of all species
    double precision, intent(IN) :: m_loc(3, nat)
     !! local magnetization, must be invariant under the sym.op.
    logical, intent(INOUT) :: sym(48)
     !! .TRUE. if rotation isym is a sym.op. of the crystal
     !! (i.e. not of the bravais lattice only)
    !
    ! ... local variables
    !
    integer :: na, nb, irot
    logical :: t1, t2
    double precision, allocatable ::  mxau(:, :), mrau(:, :)
    ! magnetization and rotated magnetization in crystal axis
    !
    allocate (mxau(3, nat), mrau(3, nat))
    !
    ! ... Compute the local magnetization of each atom in the basis of
    ! the direct lattice vectors
    do na = 1, nat
      mxau(:, na) = bg(1, :)*m_loc(1, na) + &
                    bg(2, :)*m_loc(2, na) + &
                    bg(3, :)*m_loc(3, na)
    end do
    !
    do irot = 1, nrot
      !
      t_rev(irot) = 0
      !
      if (sym(irot)) then
        !
        ! ... mrau = rotated local magnetization
        do na = 1, nat
          mrau(:, na) = s(1, :, irot)*mxau(1, na) + &
                        s(2, :, irot)*mxau(2, na) + &
                        s(3, :, irot)*mxau(3, na)
        end do
        !
        if (sname(irot) (1:3) == 'inv') mrau = -mrau
        !
        ! ... check if this a magnetic symmetry
        t1 = .true.
        t2 = .true.
        !
        do na = 1, nat
          !
          nb = irt(irot, na)
          if (nb < 1 .or. nb > nat) call errore('check_mag_sym', &
                                                'internal error: out-of-bound atomic index', na)
          !
          t1 = (abs(mrau(1, na) - mxau(1, nb)) + &
                abs(mrau(2, na) - mxau(2, nb)) + &
                abs(mrau(3, na) - mxau(3, nb)) < eps2) .and. t1
          t2 = (abs(mrau(1, na) + mxau(1, nb)) + &
                abs(mrau(2, na) + mxau(2, nb)) + &
                abs(mrau(3, na) + mxau(3, nb)) < eps2) .and. t2
          !
        end do
        !
        if (.not. t1 .and. .not. t2) then
          ! not a magnetic symmetry
          sym(irot) = .false.
        elseif (t2 .and. .not. t1) then
          ! magnetic symmetry with time reversal, if allowed
          if (no_t_rev) then
            sym(irot) = .false.
          else
            t_rev(irot) = 1
          end if
        end if
        if ((.not. sym(irot)) .and. &
            (abs(ft(1, irot)) > eps2 .or. &
             abs(ft(2, irot)) > eps2 .or. &
             abs(ft(3, irot)) > eps2)) nsym_ns = nsym_ns - 1
        !
      end if
      !
    end do
    !
    ! ... deallocate work space
    deallocate (mrau, mxau)
    !
    return
    !
  end subroutine sgam_at_mag
  !
  !
  !-----------------------------------------------------------------------
  subroutine sgam_at_collin(nat, m_loc, sym)
    !-----------------------------------------------------------------------
      !! Find spin-space-group symmetries of the collinear system, i.e.
      !! the pair of the point-group symmetries and the spin operations
      !! that are symmetries of the atomic configurations and the local magnetization.
      !! The spin operations include the identity and the time reversal.
    !
    implicit none
    !
    integer, intent(IN) :: nat
      !! numbero of atoms of all species
    double precision, intent(IN) :: m_loc(3, nat)
      !! local magnetization, must be invariant under the sym.op.
    logical, intent(INOUT) :: sym(48)
      !! .TRUE. if rotation isym is a sym.op. of the crystal
      !! (i.e. not of the bravais lattice only)
    !
    ! ... local variables
    !
    integer :: na, nb, irot
    logical :: t1, t2
    double precision, allocatable ::  m_org(:), m_op(:)
    ! magnetization and rotated magnetization in crystal axis
    !
    allocate (m_org(nat), m_op(nat))
    !
    ! Set original magnetization
    do na = 1, nat
      m_org(na) = m_loc(3, na)
    end do

    ! Check for time reversal
    do irot = 1, nrot
      !
      t_rev(irot) = 0
      !
      if (sym(irot)) then
        do na = 1, nat
          nb = irt(irot, na)
          if (nb < 1 .or. nb > nat) call errore('check_mag_sym', &
                                                'internal error: out-of-bound atomic index', na)

          m_op(nb) = m_org(na)

        end do

        if (all(abs(m_op - m_org) < 1.0d-6)) then
          ! the operation is a symmetry without time-reversal
          t_rev(irot) = 0
        else if (all(abs(m_op + m_org) < 1.0d-6)) then
          if (colin_mag == 1) then
            ! discard symmteries with time-reversal
            sym(irot) = .false.
          else ! IF ( colin_mag == 2) THEN
            ! the operation is a symmetry with time-reversal
            t_rev(irot) = 1
          end if
        else
          ! the operation is not a symmetry
          sym(irot) = .false.
        end if

        if ((.not. sym(irot)) .and. &
            (abs(ft(1, irot)) > eps2 .or. &
             abs(ft(2, irot)) > eps2 .or. &
             abs(ft(3, irot)) > eps2)) nsym_ns = nsym_ns - 1

      end if
    end do
    ! ... deallocate work space
    deallocate (m_op, m_org)
    !
    return
    !
  end subroutine sgam_at_collin
  !
  !
  !-------------------------------------------------------------------------
  subroutine set_sym(at, bg, nat, tau, ityp, nspin_mag, m_loc)
    !-----------------------------------------------------------------------
     !! This routine receives as input atomic types and positions, if there
     !! is noncollinear magnetism and the initial magnetic moments
     !! it sets the symmetry elements of this module.
     !! Note that \(at\) and \(bg\) are those in \(\textrm{cell_base}\). It sets nrot, nsym, s,
     !! sname, sr, invs, ft, irt, t_rev,  time_reversal, and invsym.
    !
    implicit none
    !
    double precision, intent(IN) :: at(3, 3), bg(3, 3)
    integer, intent(IN) :: nat
     !! number of atoms in the unit cell
    integer, intent(IN) :: ityp(nat)
     !! species of each atom in the unit cell
    integer, intent(IN) :: nspin_mag
     !! =1 when nspin=1,4 (domag=.false.), =2 when
     !! nspin=2, =4 nspin=4 (domag=.true.)
    double precision, intent(IN) :: tau(3, nat)
     !! cartesian coordinates of the atoms
    double precision, intent(IN) :: m_loc(3, nat)
     !! local magnetization, must be invariant under the sym.op.
    !
    time_reversal = (nspin_mag /= 4)
    t_rev(:) = 0
    !
    call set_sym_bl(at)
    call find_sym(at, bg, nat, tau, ityp,.not. time_reversal, m_loc)
    !
    return
    !
  end subroutine set_sym
  !
  !
  !-----------------------------------------------------------------------
  integer function copy_sym(nrot_, sym)
    !-----------------------------------------------------------------------
     !! Copy symmetry operations in sequential order so that:
    !
     !! * \(s(i,j,\text{irot})\), with \(\text{irot} \leq \text{nsym}\) are the symmetry
     !!   operations of the crystal;
     !! * \(s(i,j,\text{irot})\), with \(\text{nsym}+1<\text{irot}\leq \text{nrot}\) are
     !!   the symmetry operations of the lattice.
    !
     !! On exit \(\textrm{copy_sym}\) returns nsym.
    !
    implicit none
    !
    integer, intent(IN) :: nrot_
     !! number of rotations
    logical, intent(INOUT) :: sym(48)
     !! .TRUE. if rotation isym is a sym.op. of the crystal
     !! (i.e. not of the bravais lattice only)
    !
    ! ... local variables
    !
    integer :: stemp(3, 3), ftemp(3), ttemp, irot, jrot
    double precision :: ft_(3)
    integer, allocatable :: irtemp(:)
    character(LEN=45) :: nametemp
    !
    !
    allocate (irtemp(size(irt, 2)))
    !
    jrot = 0
    !
    do irot = 1, nrot_
      if (sym(irot)) then
        jrot = jrot + 1
        if (irot > jrot) then
          stemp = s(:, :, jrot)
          s(:, :, jrot) = s(:, :, irot)
          s(:, :, irot) = stemp
          ft_(:) = ft(:, jrot)
          ft(:, jrot) = ft(:, irot)
          ft(:, irot) = ft_(:)
          irtemp(:) = irt(jrot, :)
          irt(jrot, :) = irt(irot, :)
          irt(irot, :) = irtemp(:)
          nametemp = sname(jrot)
          sname(jrot) = sname(irot)
          sname(irot) = nametemp
          ttemp = t_rev(jrot)
          t_rev(jrot) = t_rev(irot)
          t_rev(irot) = ttemp
        end if
      end if
    end do
    !
    sym(1:jrot) = .true.
    sym(jrot + 1:nrot_) = .false.
    !
    deallocate (irtemp)
    !
    copy_sym = jrot
    !
    return
    !
  end function copy_sym
  !
  !
  !-----------------------------------------------------------------------
  logical function is_group(nsym_)
    !-----------------------------------------------------------------------
     !! Checks that {S} is a group.
    !
    implicit none
    !
    integer, intent(IN) :: nsym_
    integer :: isym, jsym, ksym, ss(3, 3)
    double precision :: st(3), dt(3)
    logical :: found
    !
    do isym = 1, nsym_
      do jsym = 1, nsym_
        !
        ss = matmul(s(:, :, isym), s(:, :, jsym))
        st(:) = ft(:, jsym) + s(1, :, jsym)*ft(1, isym) + &
                s(2, :, jsym)*ft(2, isym) + &
                s(3, :, jsym)*ft(3, isym)
        !
        ! ... here we check that the input matrices really form a group:
        ! S(k) = S(i)*S(j)
        ! ftau_k = S(j)*ftau_i+ftau_j (modulo a lattice vector)
        !
        found = .false.
        !
        do ksym = 1, nsym_
          dt(:) = ft(:, ksym) - st(:) - nint(ft(:, ksym) - st(:))
          if (all(s(:, :, ksym) == ss(:, :)) .and. &
              (abs(dt(1)) < eps2) .and. &
              (abs(dt(2)) < eps2) .and. &
              (abs(dt(3)) < eps2)) then
            if (found) then
              is_group = .false.
              return
            end if
            found = .true.
          end if
        end do
        !
        if (.not. found) then
          is_group = .false.
          return
        end if
        !
      end do
    end do
    !
    is_group = .true.
    !
    return
    !
  end function is_group
  !
  !
  !-----------------------------------------------------------------------
  logical function checksym(irot, nat, ityp, xau, rau, ft_)
    !-----------------------------------------------------------------------
     !! This function receives as input all the atomic positions xau,
     !! and the rotated rau by the symmetry operation ir. It returns
     !! .TRUE. if, for each atom na, it is possible to find an atom nb
     !! which is of the same type of na, and coincides with it after the
     !! symmetry operation. Fractional translations are allowed.
    !
    implicit none
    !
    integer, intent(IN) :: nat
     !! number of atoms
    integer, intent(IN) :: ityp(nat)
     !! the type of each atom
    integer, intent(IN) :: irot
     !! rotation index
    double precision, intent(IN) :: xau(3, nat)
     !! the initial vectors (in crystal coordinates)
    double precision, intent(IN) :: rau(3, nat)
     !! the rotated vectors (as above)
    double precision, intent(IN) :: ft_(3)
     !! fractionary translation (as above)
    !
    ! ... local variables
    !
    integer :: na, nb
    logical, external :: eqvect
    ! the testing function
    !
    do na = 1, nat
      do nb = 1, nat
        !
        if ((colin_mag >= 0 .and. chem_symb(atm(ityp(nb))) == chem_symb(atm(ityp(na)))) &
            .or. (colin_mag < 0 .and. ityp(nb) == ityp(na))) then
          !IF ( ityp(nb) == ityp(na) ) THEN
          checksym = eqvect(rau(1, na), xau(1, nb), ft_, accep)
          if (checksym) then
            !
            ! ... the rotated atom does coincide with one of the like atoms
            ! keep track of which atom the rotated atom coincides with
            irt(irot, na) = nb
            goto 10
            !
          end if
        end if
        !
      end do
      !
      ! ... the rotated atom does not coincide with any of the like atoms
      ! s(ir) + ft is not a symmetry operation
      checksym = .false.
      return
      !
10    continue
    end do
    !
    ! ... s(ir) + ft is a symmetry operation
    !
    return
    !
  end function checksym
  !
  !
  !-----------------------------------------------------------------------
  subroutine checkallsym(at, bg, nat, tau, ityp)
    !-----------------------------------------------------------------------
     !! Given a crystal group this routine checks that the actual atomic
     !! positions and bravais lattice vectors are compatible with it.
     !! Used in relaxation/MD runs to check that atomic motion is
     !! consistent with assumed symmetry.
    !
    implicit none
    !
    double precision, intent(IN) :: at(3, 3), bg(3, 3)
    integer, intent(IN) :: nat
     !! number of atoms
    integer, intent(IN) :: ityp(nat)
     !! the type of each atom
    double precision, intent(IN) :: tau(3, nat)
     !! postion of each atom
    !
    ! ... local variables
    !
    integer :: na, kpol, isym, i, j, k, l
    logical :: loksym(48)
    double precision :: sx(3, 3), sy(3, 3)
    double precision, allocatable :: xau(:, :), rau(:, :)
    !
    allocate (xau(3, nat))
    allocate (rau(3, nat))
    !
    ! ... check that s(i,j, isym) is an orthogonal operation
    do isym = 1, nsym
      sx = dble(s(:, :, isym))
      sy = matmul(bg, sx)
      sx = matmul(sy, transpose(at))
      ! sx is s in cartesian axis
      sy = matmul(transpose(sx), sx)
      ! sy = s*TRANSPOSE(s) = I
      do i = 1, 3
        sy(i, i) = sy(i, i) - 1.0d0
      end do
      if (any(abs(sy) > eps1)) &
        call errore('checkallsym', 'not orthogonal operation', isym)
    end do
    !
    ! ... Compute the coordinates of each atom in the basis of the lattice
    do na = 1, nat
      do kpol = 1, 3
        xau(kpol, na) = bg(1, kpol)*tau(1, na) + &
                        bg(2, kpol)*tau(2, na) + &
                        bg(3, kpol)*tau(3, na)
      end do
    end do
    !
    ! ... Generate the coordinates of the rotated atoms
    do isym = 1, nsym
      do na = 1, nat
        do kpol = 1, 3
          rau(kpol, na) = s(1, kpol, isym)*xau(1, na) + &
                          s(2, kpol, isym)*xau(2, na) + &
                          s(3, kpol, isym)*xau(3, na)
        end do
      end do
      !
      loksym(isym) = checksym(isym, nat, ityp, xau, rau, ft(1, isym))
    end do
    !
    ! ... deallocate work space
    !
    deallocate (rau)
    deallocate (xau)
    !
    do isym = 1, nsym
      if (.not. loksym(isym)) call errore('checkallsym', &
                                          'the following symmetry operation is not satisfied  ', -isym)
    end do
    !
    if (any(.not. loksym(1:nsym))) then
      !call symmetrize_at (nsym, s, invs, ft, irt, nat, tau, at, bg, &
      !                    alat, omega)
      call errore('checkallsym', &
                  'some of the original symmetry operations not satisfied ', 1)
    end if
    !
    return
    !
  end subroutine checkallsym
  !
  !
  !----------------------------------------------------------------------
  subroutine s_axis_to_cart(at, bg)
    !----------------------------------------------------------------------
     !! This routine transforms symmetry matrices expressed in the
     !! basis of the crystal axis into rotations in cartesian axis.
    !
    implicit none
    double precision, intent(IN) :: at(3, 3), bg(3, 3)
    !
    integer :: isym
    double precision :: sa(3, 3), sb(3, 3)
    !
    do isym = 1, nsym
      sa(:, :) = dble(s(:, :, isym))
      sb = matmul(bg, sa)
      sr(:, :, isym) = matmul(at, transpose(sb))
    end do
    !
  end subroutine s_axis_to_cart
  !
  !
  !-----------------------------------------------------------------------
  subroutine find_sym_ifc(at, bg, nat, tau, ityp)
    !-----------------------------------------------------------------------
      !! This routine finds the point group of the crystal, by eliminating
      !! the symmetries of the Bravais lattice which are not allowed
      !! by the atomic positions (for use in the FD package).
    !
    implicit none
    !
    double precision, intent(IN) :: at(3, 3), bg(3, 3)
    integer, intent(IN) :: nat
      !! number of atoms
    integer, intent(IN) :: ityp(nat)
      !! the type of each atom
    double precision, intent(IN) :: tau(3, nat)
      !! postion of each atom
    !
    ! ... local variables
    !
    integer :: i
    logical :: sym(48)
    ! if true the corresponding operation is a symmetry operation
    !
    if (.not. allocated(irt)) allocate (irt(48, nat))
    irt(:, :) = 0
    !
    ! ... Here we find the true symmetries of the crystal
    call sgam_at_ifc(at, bg, nat, tau, ityp, sym)
    !
    ! ... Here we re-order all rotations in such a way that true sym.ops
    ! are the first nsym; rotations that are not sym.ops. follow
    nsym = copy_sym(nrot, sym)
    !
    ! ... check if inversion (I) is a symmetry.
    ! If so, it should be the (nsym/2+1)-th operation of the group
    invsym = all(s(:, :, nsym/2 + 1) == -s(:, :, 1))
    !
    call inverse_s()
    !
    call s_axis_to_cart(at, bg)
    !
    return
    !
  end subroutine find_sym_ifc
  !
  !
  !-----------------------------------------------------------------------
  subroutine sgam_at_ifc(at, bg, nat, tau, ityp, sym)
    !-----------------------------------------------------------------------
      !! Given the point group of the Bravais lattice, this routine finds
      !! the subgroup which is the point group of the considered crystal.
      !! Non symmorphic groups are allowed, provided that fractional
      !! translations are allowed (nofrac=.false), that the unit cell is
      !! not a supercell.
    !
      !! On output, the array sym is set to .TRUE.. for each operation
      !! of the original point group that is also a symmetry operation
      !! of the crystal symmetry point group.
    !
    implicit none
    double precision, intent(IN) :: at(3, 3), bg(3, 3)
    !
    integer, intent(IN) :: nat
      !! number of atoms in the unit cell
    integer, intent(IN) :: ityp(nat)
      !! species of each atom in the unit cell
    double precision, intent(IN) :: tau(3, nat)
      !! cartesian coordinates of the atoms
    logical, intent(OUT) :: sym(48)
      !! flag indicating if sym.op. isym in the parent group
      !! is a true symmetry operation of the crystal
    !
    ! ... local variables
    !
    integer :: na, kpol, nb, irot, i, j
    ! counters
    double precision, allocatable :: xau(:, :), rau(:, :)
    ! atomic coordinates in crystal axis
    logical :: fractional_translations
    double precision :: ft_(3)
    !
    allocate (xau(3, nat))
    allocate (rau(3, nat))
    !
    ! ... Compute the coordinates of each atom in the basis of
    ! the direct lattice vectors.
    !
    do na = 1, nat
      xau(:, na) = bg(1, :)*tau(1, na) + bg(2, :)*tau(2, na) + bg(3, :)*tau(3, na)
    end do
    !
    ! ... check if the identity has fractional translations
    ! (this means that the cell is actually a supercell).
    ! When this happens, fractional translations are disabled,
    ! because there is no guarantee that the generated sym.ops.
    ! form a group.
    !
    nb = 1
    irot = 1
    !
    fractional_translations = .not. nofrac
    !
    do na = 2, nat
      if (fractional_translations) then
        if ((colin_mag >= 0 .and. chem_symb(atm(ityp(nb))) == chem_symb(atm(ityp(na)))) &
            .or. (colin_mag < 0 .and. ityp(nb) == ityp(na))) then

          !IF ( ityp(nb) == ityp(na) ) THEN
          ft_(:) = xau(:, na) - xau(:, nb) - nint(xau(:, na) - xau(:, nb))
          !
          sym(irot) = checksym(irot, nat, ityp, xau, xau, ft_)
          !
          if (sym(irot) .and. &
              (abs(ft_(1)**2 + ft_(2)**2 + ft_(3)**2) < 1.d-8)) &
            call errore('sgam_at_ifc', 'overlapping atoms', na)
        end if
      end if
    end do
    !
    nsym_ns = 0
    !
    do irot = 1, nrot
      do na = 1, nat
        ! rau = rotated atom coordinates
        rau(:, na) = s(1, :, irot)*xau(1, na) + &
                     s(2, :, irot)*xau(2, na) + &
                     s(3, :, irot)*xau(3, na)
      end do
      !
      ! ... first attempt: no fractional translation
      ft(:, irot) = 0
      ft_(:) = 0.d0
      !
      sym(irot) = checksym(irot, nat, ityp, xau, rau, ft_)
      !
      if (.not. sym(irot) .and. fractional_translations) then
        nb = 1
        !
        do na = 1, nat
          if ((colin_mag >= 0 .and. chem_symb(atm(ityp(nb))) == chem_symb(atm(ityp(na)))) &
              .or. (colin_mag < 0 .and. ityp(nb) == ityp(na))) then
            !IF ( ityp(nb) == ityp(na) ) THEN
            !
            !      second attempt: check all possible fractional translations
            !
            ft_(:) = rau(:, na) - xau(:, nb) - nint(rau(:, na) - xau(:, nb))
            !
            sym(irot) = checksym(irot, nat, ityp, xau, rau, ft_)
            !
            if (sym(irot)) then
              nsym_ns = nsym_ns + 1
              ft(:, irot) = ft_(:)
              goto 100
            end if
          end if
        end do
        !
      end if
      !
100   continue
      !
    end do
    !
    deallocate (rau)
    deallocate (xau)
    !
    return
    !
  end subroutine sgam_at_ifc
  !
  !-----------------------------------------------------------------------
  function check_grid_sym(nr1, nr2, nr3) result(compatible)
    !---------------------------------------------------------------------
      !! Check that symmetry operations and FFT grid are compatible
      !! Needed to prevent trouble with real-space symmetrization
    !
    implicit none
    !
    integer, intent(IN) :: nr1, nr2, nr3
    logical :: compatible, bad
    integer :: isym, i, j
    !
    compatible = .true.
    do isym = 1, nsym
      !
      bad = (mod(s(2, 1, isym)*nr1, nr2) /= 0 .or. &
             mod(s(3, 1, isym)*nr1, nr3) /= 0 .or. &
             mod(s(1, 2, isym)*nr2, nr1) /= 0 .or. &
             mod(s(3, 2, isym)*nr2, nr3) /= 0 .or. &
             mod(s(1, 3, isym)*nr3, nr1) /= 0 .or. &
             mod(s(2, 3, isym)*nr3, nr2) /= 0)
      if (bad) then
        write (*, '(5x,"warning: symmetry operation # ",i2, &
             &         " not compatible with FFT grid. ")') isym
        write (*, '(3i4)') ((s(i, j, isym), j=1, 3), i=1, 3)
        compatible = .false.
      end if
      !
    end do
    !
  end function check_grid_sym
  !
  !-----------------------------------------------------------------------
  subroutine remove_sym(at, bg, nr1, nr2, nr3)
    !---------------------------------------------------------------------
      !! Compute ftau used for symmetrization in real space (phonon, exx)
      !! ensure that they are commensurated with the FFT grid.
    !
    implicit none
    !
    double precision, intent(IN) :: at(3, 3), bg(3, 3)
    integer, intent(IN) :: nr1, nr2, nr3
    !
    ! ... local variables
    !
    logical :: sym(48)
    integer :: isym, nsym_, i, j
    double precision :: ftaux(3)
    !
    nsym_ = nsym
    sym(1:nsym_) = .true.
    nsym_na = 0
    !
    do isym = 1, nsym_
      !
      ! check that the grid is compatible with the S rotation
      !
      if (mod(s(2, 1, isym)*nr1, nr2) /= 0 .or. &
          mod(s(3, 1, isym)*nr1, nr3) /= 0 .or. &
          mod(s(1, 2, isym)*nr2, nr1) /= 0 .or. &
          mod(s(3, 2, isym)*nr2, nr3) /= 0 .or. &
          mod(s(1, 3, isym)*nr3, nr1) /= 0 .or. &
          mod(s(2, 3, isym)*nr3, nr2) /= 0) then
        sym(isym) = .false.
        write (*, '(5x,"warning: symmetry operation # ",i2, &
             &         " not compatible with FFT grid. ")') isym
        write (*, '(3i4)') ((s(i, j, isym), j=1, 3), i=1, 3)
        sym(isym) = .false.
        if (abs(ft(1, isym)) > eps2 .or. &
            abs(ft(2, isym)) > eps2 .or. &
            abs(ft(3, isym)) > eps2) nsym_ns = nsym_ns - 1
      end if
      !
      ! convert ft to FFT coordinates, check if compatible with FFT grid
      ! for real-space symmetrization
      !
      ftaux(1) = ft(1, isym)*nr1
      ftaux(2) = ft(2, isym)*nr2
      ftaux(3) = ft(3, isym)*nr3
      ! check if the fractional translations are commensurate
      ! with the FFT grid, discard sym.op. if not
      ! (needed because ph.x symmetrizes in real space)
      if (abs(ftaux(1) - nint(ftaux(1)))/nr1 > eps2 .or. &
          abs(ftaux(2) - nint(ftaux(2)))/nr2 > eps2 .or. &
          abs(ftaux(3) - nint(ftaux(3)))/nr3 > eps2) then
        !     WRITE( *, '(5x,"warning: symmetry operation", &
        !          &     " # ",i2," not allowed.   fractional ", &
        !          &     "translation:"/5x,3f11.7,"  in crystal", &
        !          &     " coordinates")') isym, ft_
        sym(isym) = .false.
        nsym_na = nsym_na + 1
        nsym_ns = nsym_ns - 1
      end if
      !
    end do
    !
    ! ... count symmetries, reorder them exactly as in "find_sym"
    !
    nsym = copy_sym(nsym_, sym)
    invsym = all(s(:, :, nsym/2 + 1) == -s(:, :, 1))
    !
    call inverse_s()
    call s_axis_to_cart(at, bg)
    !
  end subroutine remove_sym
  !
  !
  !--------------------------------------------------------------------------
  integer function mcm(i, j)
    !------------------------------------------------------------------------
      !! Returns minimum common multiple of two integers
      !! if i=0, returns j, and vice versa; if i<0 or j<0, returns -1.
    !
    integer, intent(IN) :: i, j
    integer :: n1, n2, k
    !
    if (i < 0 .or. j < 0) then
      mcm = -1
    elseif (i == 0 .and. j == 0) then
      mcm = 0
    else
      n1 = min(i, j)
      n2 = max(i, j)
      do k = 1, n1
        mcm = k*n2
        if (mod(mcm, n1) == 0) return
      end do
      mcm = n2
    end if
    !
  end function mcm
  !
  !
  !--------------------------------------------------------------------------
  character function chem_symb(symbol)
    !------------------------------------------------------------------------
      !! Returns the chemical symbol used to identify the symmetry
    !
    implicit none
    !
    character(LEN=*), intent(IN) :: symbol
    !
    if (scan(symbol, "0123456789") == 0) then
      chem_symb = symbol
    else
      chem_symb = symbol(1:scan(symbol, "0123456789_-") - 1)
    end if
    !
  end function chem_symb
  !
  !

end module symm_base
