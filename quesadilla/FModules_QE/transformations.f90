!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine trntnsc(phi, at, bg, iflg)
  !-----------------------------------------------------------------------
  !
  ! trasforms a COMPLEX tensor (like the dynamical matrix)
  ! from crystal to cartesian axis (iflg >=  1) or viceversa (iflg <= -1)
  !
  implicit none

  integer :: iflg
  ! input: gives the versus of the trans.

  double complex :: phi(3, 3)
  ! inp/out: the matrix to transform

  double precision :: at(3, 3), bg(3, 3)
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice

  integer :: i, j, k, l
  !
  !  counters on polarizations
  ! /
  !/

  double complex :: wrk(3, 3)
  ! a working array
  if (iflg .gt. 0) then
    !
    ! forward transformation (crystal to cartesian axis)
    !

    call zcopy(9, phi, 1, wrk, 1)
    do i = 1, 3
      do j = 1, 3
        phi(i, j) = (0.d0, 0.d0)
        do k = 1, 3
          do l = 1, 3
            phi(i, j) = phi(i, j) + wrk(k, l)*bg(i, k)*bg(j, l)
          end do
        end do
      end do
    end do
  else
    !
    ! backward transformation (cartesian to crystal axis)
    !
    do i = 1, 3
      do j = 1, 3
        wrk(i, j) = (0.d0, 0.d0)
        do k = 1, 3
          do l = 1, 3
            wrk(i, j) = wrk(i, j) + phi(k, l)*at(k, i)*at(l, j)
          end do
        end do
      end do
    end do
    call zcopy(9, wrk, 1, phi, 1)
  end if
  return
end subroutine trntnsc

!-----------------------------------------------------------------------
subroutine compact_dyn(nat, fcq, phi)
  !-----------------------------------------------------------------------
   !! This routine writes the dynamical matrix from a 3,3,nat,nat array
   !! to a 3*nat,3*nat array.
   !! The versions included with this package are different from the
   !! original ones as they are set up to deal with the way
   !! the force constants are stored in the python interface.
  implicit none
  integer, intent(IN) :: nat
  double complex, intent(IN) :: phi(3, 3, nat, nat)
  double complex, intent(OUT) :: fcq(3*nat, 3*nat)

  integer :: na, nb, i, j, icart, jcart

  do na = 1, nat
    do nb = 1, nat
      do icart = 1, 3
        i = 3*(na - 1) + icart
        do jcart = 1, 3
          j = 3*(nb - 1) + jcart
          fcq(i, j) = phi(icart, jcart, na, nb)
        end do
      end do
    end do
  end do
  return
end subroutine compact_dyn

!-----------------------------------------------------------------------
subroutine scompact_dyn(nat, fcq, phi)
  !-----------------------------------------------------------------------
   !! This routine writes the dynamical matrix from a 3*nat,3*nat array
   !! to a 3,3,nat,nat array.
   !! The versions included with this package are different from the
   !! original ones as they are set up to deal with the way
   !! the force constants are stored in the python interface.

  implicit none
  integer, intent(IN) :: nat
  double complex, intent(OUT) :: phi(3, 3, nat, nat)
  double complex, intent(IN) :: fcq(3*nat, 3*nat)

  integer :: na, nb, i, j, icart, jcart

  do i = 1, 3*nat
    na = (i - 1)/3 + 1
    icart = i - 3*(na - 1)
    do j = 1, 3*nat
      nb = (j - 1)/3 + 1
      jcart = j - 3*(nb - 1)
      phi(icart, jcart, na, nb) = fcq(i, j)
    end do
  end do
  return
end subroutine scompact_dyn
!-----------------------------------------------------------------------
subroutine cryst_to_cart(nvec, vec, trmat, iflag)
  !-----------------------------------------------------------------------
  !
  !     This routine transforms the atomic positions or the k-point
  !     components from crystallographic to cartesian coordinates
  !     ( iflag=1 ) and viceversa ( iflag=-1 ).
  !     Output cartesian coordinates are stored in the input ('vec') array
  !
  !
  implicit none
  !
  integer, intent(in) :: nvec, iflag
  ! nvec:  number of vectors (atomic positions or k-points)
  !        to be transformed from crystal to cartesian and vice versa
  ! iflag: gives the direction of the transformation
  double precision, intent(in) :: trmat(3, 3)
  ! trmat: transformation matrix
  ! if iflag=1:
  !    trmat = at ,  basis of the real-space lattice,       for atoms   or
  !          = bg ,  basis of the reciprocal-space lattice, for k-points
  ! if iflag=-1: the opposite
  double precision, intent(inout) :: vec(3, nvec)
  ! coordinates of the vector (atomic positions or k-points) to be
  ! transformed - overwritten on output
  !
  !    local variables
  !
  integer :: nv, kpol
  ! counter on vectors
  ! counter on polarizations
  double precision :: vau(3)
  ! workspace
  !
  !     Compute the cartesian coordinates of each vectors
  !     (atomic positions or k-points components)
  !
  do nv = 1, nvec
    if (iflag .eq. 1) then
      do kpol = 1, 3
        vau(kpol) = trmat(kpol, 1)*vec(1, nv) + trmat(kpol, 2) &
                    *vec(2, nv) + trmat(kpol, 3)*vec(3, nv)
      end do
    else
      do kpol = 1, 3
        vau(kpol) = trmat(1, kpol)*vec(1, nv) + trmat(2, kpol) &
                    *vec(2, nv) + trmat(3, kpol)*vec(3, nv)
      end do
    end if
    do kpol = 1, 3
      vec(kpol, nv) = vau(kpol)
    end do
  end do
  !
  return
end subroutine cryst_to_cart
