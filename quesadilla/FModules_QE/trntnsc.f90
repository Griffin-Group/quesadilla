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
!SUBROUTINE compact_dyn(nat, dyn, phi)
!   !-----------------------------------------------------------------------
!   !! This routine writes the dynamical matrix from a 3,3,nat,nat array
!   !! to a 3*nat,3*nat array.
!   !
!   IMPLICIT NONE
!   INTEGER, INTENT(IN) :: nat
!   double complex, INTENT(IN) :: phi(3,3,nat,nat)
!   double complex, INTENT(OUT) :: dyn(3*nat, 3*nat)
!
!   INTEGER :: na, nb, icart, jcart, imode, jmode
!
!   DO na = 1, nat
!      DO icart = 1, 3
!         imode = 3 * ( na - 1 ) + icart
!         DO nb = 1, nat
!            DO jcart = 1, 3
!               jmode = 3 * ( nb - 1 ) + jcart
!               dyn (imode, jmode) = phi (icart, jcart, na, nb)
!            END DO
!         END DO
!      END DO
!   END DO
!   RETURN
!   END SUBROUTINE compact_dyn
! !
! !-----------------------------------------------------------------------
!SUBROUTINE scompact_dyn(nat, dyn, phi)
!   !-----------------------------------------------------------------------
!   !! This routine writes the dynamical matrix from a 3*nat,3*nat array
!   !! to a 3,3,nat,nat array.
!   !
!
!   IMPLICIT NONE
!   INTEGER, INTENT(IN) :: nat
!   double complex, INTENT(OUT) :: phi(3,3,nat,nat)
!   double complex, INTENT(IN) :: dyn(3*nat, 3*nat)
!
!   INTEGER :: na, nb, icart, jcart, imode, jmode
!
!   DO na = 1, nat
!      DO icart = 1, 3
!         imode = 3 * ( na - 1 ) + icart
!         DO nb = 1, nat
!            DO jcart = 1, 3
!               jmode = 3 * ( nb - 1 ) + jcart
!               phi (icart, jcart, na, nb) = dyn (imode, jmode)
!            END DO
!         END DO
!      END DO
!   END DO
!   RETURN
!END SUBROUTINE scompact_dyn
