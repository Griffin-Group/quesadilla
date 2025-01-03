!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine q2qstar_ph(fcq, at, bg, nat, nsym, s, invs, irt, rtau, &
                      nq, sxq, isq, imq, nq_tot, fcqstar)
  !-----------------------------------------------------------------------
  !! Generates the dynamical matrices for the star of q and writes them on
  !! disk for later use.
  !! If there is a symmetry operation such that \(q \rightarrow -q+G \) then
  !! imposes on dynamical matrix those conditions related to time reversal
  !! symmetry.
  !
  implicit none
  !
  integer, intent(in) :: nat
  !! number of atoms in the unit cell
  integer, intent(in) :: nsym
  !! number of symmetry operations
  integer, intent(in) :: s(3, 3, 48)
  !! the symmetry operations
  integer, intent(in) :: invs(48)
  !! index of the inverse operations
  integer, intent(in) :: irt(48, nat)
  !! index of the rotated atom
  integer, intent(in) :: nq
  !! degeneracy of the star of q
  integer, intent(in) :: isq(48)
  !! symmetry op. giving the rotated q
  integer, intent(in) :: imq
  !! index of -q in the star (0 if non present)
  double complex, intent(in) :: fcq(3*nat, 3*nat)
  !! the input dynamical matrix. If \(\text{imq}\) different
  !! from 0 the output matrix is symmetrized w.r.t. time-reversal
  double precision, intent(in) :: at(3, 3)
  !! direct lattice vectors
  double precision, intent(in) :: bg(3, 3)
  !! reciprocal lattice vectors
  double precision, intent(in) :: rtau(3, 48, nat)
  !! for each atom and rotation gives the R vector involved
  double precision, intent(in) :: sxq(3, 48)
  !! list of q in the star
  integer, intent(in) :: nq_tot
  !! Total number of q accounting for the -q if it's not already present
  double complex, intent(out) :: fcqstar(nq_tot, 3, 3, nat, nat)
  !! the output dynamical matrix at each q in the star

  !
  ! ... local variables
  !
  integer :: na, nb, iq, nsq, isym, i, j, counter, icar, jcar
  ! counters
  ! nsq: number of sym.op. giving each q in the list

  double complex :: phi(3, 3, nat, nat), phi2(3, 3, nat, nat), phiqstar(nq_tot, 3, 3, nat, nat)
  ! work space
  counter = 0
  !
  ! Sets number of symmetry operations giving each q in the list
  !
  nsq = nsym/nq
  if (nsq*nq /= nsym) call errore('q2star_ph', 'wrong degeneracy', 1)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writes dyn.mat. dyn(3*nat,3*nat) on the 4-index array phi(3,3,nat,nat)
  !
  !CALL scompact_dyn(nat, fcq, phi)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !phi = fcq
  do i = 1, 3*nat
    na = (i - 1)/3 + 1
    icar = i - 3*(na - 1)
    do j = 1, 3*nat
      nb = (j - 1)/3 + 1
      jcar = j - 3*(nb - 1)
      phi(icar, jcar, na, nb) = fcq(i, j)
      !print *, "PHI1:", icar, jcar, na,nb, phi(icar, jcar, na, nb)
    end do
  end do
  !
  ! Go from Cartesian to Crystal coordinates
  !
  do na = 1, nat
    do nb = 1, nat
      call trntnsc(phi(:, :, na, nb), at, bg, -1)
    end do
  end do
  !
  ! If -q is in the list impose first of all the conditions coming from
  ! time reversal symmetry
  !
  if (imq /= 0) then
    phi2(:, :, :, :) = (0.d0, 0.d0)
    isym = 1
    do while (isq(isym) /= imq)
      isym = isym + 1
    end do
    call rotate_and_add_dyn(phi, phi2, nat, isym, s, invs, irt, &
                            rtau, sxq(1, imq))
    do na = 1, nat
      do nb = 1, nat
        do i = 1, 3
          do j = 1, 3
            phi(i, j, na, nb) = 0.5d0*(phi(i, j, na, nb) + &
                                       conjg(phi2(i, j, na, nb)))
          end do
        end do
      end do
    end do
    phi2(:, :, :, :) = phi(:, :, :, :)
    !
    ! Go from Crystal to Cartesian
    !
    do na = 1, nat
      do nb = 1, nat
        call trntnsc(phi2(1, 1, na, nb), at, bg, +1)
      end do
    end do
     !!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Saves 4-index array phi2(3,3,nat,nat) on the dyn.mat. dyn(3*nat,3*nat)
    !
    !CALL compact_dyn(nat, dyn, phi2)
     !!!!!!!!!!!!!!!
  end if
  !
  ! For each q of the star rotates phi with the appropriate sym.op. -> phi
  !
  do iq = 1, nq
    phi2(:, :, :, :) = (0.d0, 0.d0)
    do isym = 1, nsym
      if (isq(isym) == iq) then
        call rotate_and_add_dyn(phi, phi2, nat, isym, s, invs, irt, &
                                rtau, sxq(1, iq))
      end if
    end do
    phi2(:, :, :, :) = phi2(:, :, :, :)/dble(nsq)
    !
    ! Back to cartesian coordinates
    !
    do na = 1, nat
      do nb = 1, nat
        call trntnsc(phi2(1, 1, na, nb), at, bg, +1)
      end do
    end do
    !
    ! Writes the dynamical matrix in cartesian coordinates on file
    !
    counter = counter + 1
    phiqstar(counter, :, :, :, :) = phi2
    if (imq == 0) then
      !
      ! if -q is not in the star recovers its matrix by time reversal
      !
      do na = 1, nat
        do nb = 1, nat
          do i = 1, 3
            do j = 1, 3
              phi2(i, j, na, nb) = conjg(phi2(i, j, na, nb))
            end do
          end do
        end do
      end do
      !
      ! and writes it (changing temporarily sign to q)
      !
      phiqstar(counter + nq, :, :, :, :) = phi2
    end if
  end do
  fcqstar = phiqstar
  !DO iq = 1, nq_tot
  !  !CALL compact_dyn(nat, fcqstar(iq, :, :), phiqstar(iq, :, :, :, :))
  !  DO i = 1, 3*nat
  !    na = (i - 1)/3 + 1
  !    icar = i - 3*(na - 1)
  !    DO j = 1, 3*nat
  !      nb = (j - 1)/3 + 1
  !      jcar = j - 3*(nb - 1)
  !      fcqstar(iq, i, j) = phiqstar(iq, icar, jcar, na, nb)
  !    END DO
  !  END DO
  !END DO
  !
  return
end subroutine q2qstar_ph
