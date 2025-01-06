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

! !-----------------------------------------------------------------------
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
