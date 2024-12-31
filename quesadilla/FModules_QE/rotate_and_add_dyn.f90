!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine rotate_and_add_dyn(phi, phi2, nat, isym, s, invs, irt, &
                              rtau, sxq)
  !-----------------------------------------------------------------------

  !  Rotates a dynamical matrix (phi) in crystal coordinates according
  !  to the specified symmetry operation and add the rotated matrix
  !  to phi2.   phi is left unmodified.
  !
  implicit none
  ! input variables

  integer, intent(in) :: nat, isym, s(3, 3, 48), invs(48), irt(48, nat)
  ! number of atoms in the unit cell
  ! index of the symm.op.
  ! the symmetry operations
  ! index of the inverse operations
  ! index of the rotated atom

  double complex, intent(in) :: phi(3, 3, nat, nat)
  double complex, intent(out) :: phi2(3, 3, nat, nat)
  ! the input dyn.mat. in crystal coordinates
  ! the rotated dyn.mat. in crystal coordinates

  double precision, intent(in) :: rtau(3, 48, nat), sxq(3)
  ! for eaxh atom and rotation gives the R vector
  !involved
  ! the rotated q involved in this sym.op.
  !  local variables
  integer :: na, nb, sna, snb, ism1, i, j, k, l
  ! counters on atoms
  ! indices of rotated atoms
  ! index of the inverse symm.op.
  ! generic counters
  double precision :: arg
  ! argument of the phase
  double complex :: phase, work
  double precision, parameter :: tpi = 6.283185307179586

  ism1 = invs(isym)
  do na = 1, nat
    do nb = 1, nat
      sna = irt(isym, na)
      snb = irt(isym, nb)
      arg = (sxq(1)*(rtau(1, isym, na) - rtau(1, isym, nb)) &
             + sxq(2)*(rtau(2, isym, na) - rtau(2, isym, nb)) + sxq(3) &
             *(rtau(3, isym, na) - rtau(3, isym, nb)))*tpi
      phase = DCMPLX(cos(arg), -sin(arg))
      do i = 1, 3
        do j = 1, 3
          work = DCMPLX(0.d0, 0.d0)
          do k = 1, 3
            do l = 1, 3
              work = work + s(i, k, ism1)*s(j, l, ism1)*phi(k, l, na, nb) &
                     *phase
            end do
          end do
          phi2(i, j, sna, snb) = phi2(i, j, sna, snb) + work
        end do
      end do
    end do
  end do
  !
  return
end subroutine rotate_and_add_dyn

!-----------------------------------------------------------------------
subroutine symdynph_gq_new(xq, at, bg, fcq, s, invs, rtau, irt, nsymq, &
                           nat, irotmq, minus_q, fcqsymm)
!-----------------------------------------------------------------------
  !! This routine receives as input an unsymmetrized dynamical
  !! matrix expressed on the crystal axes and imposes the symmetry
  !! of the small group of q. Furthermore it imposes also the symmetry
  !! q -> -q+G if present.
  !! February 2020: Update (A. Urru) to include the symmetry operations
  !! that require the time reversal operator (meaning that TS is a
  !! symmetry of the crystal). For more information please see:
  !! Phys. Rev. B 100, 045115 (2019).
  !
  USE symm_base, ONLY: t_rev
  !
  implicit none
  !
  integer, intent(in) :: nat
  !! input: the number of atoms
  integer, intent(in) :: s(3, 3, 48)
  !! input: the symmetry matrices
  integer, intent(in) :: irt(48, nat)
  !! input: the rotated of each vector
  integer, intent(in) :: invs(48)
  !! input: the inverse of each matrix
  integer, intent(in) :: nsymq
  !! input: the order of the small group
  integer, intent(in) :: irotmq
  !! input: the rotation sending q ->-q+G
  double precision, intent(in) :: xq(3), at(3, 3), bg(3, 3)
  !! input: the q point, the direct and reciprocal lattice vectors
  double precision, intent(in) :: rtau(3, 48, nat)
  !! input: the R associated at each t
  logical, intent(in) :: minus_q
  !! input: true if a symmetry q->-q+G
  double complex, intent(in) :: fcq(3*nat, 3*nat)
  !! input: the matrix to symmetrize
  double complex, intent(out) :: fcqsymm(3, 3, nat, nat)
  !! output: the matrix to symmetrize
  !
  ! ... local variables
  !
  double precision, parameter :: tpi = 6.283185307179586
  ! 2*pi
  integer :: isymq, sna, snb, irot, na, nb, ipol, jpol, lpol, kpol, icar, jcar, i, j, &
             iflb(nat, nat)
  ! counters, indices, work space

  double precision :: arg
  ! the argument of the phase

  double complex :: phi(3,3,nat,nat), phip(3, 3, nat, nat), work(3, 3), fase, faseq(48)
  ! work space, phase factors

  !CALL scompact_dyn(nat, fcq, phi)
  do i = 1, 3 * nat
    na = (i - 1) / 3 + 1
    icar = i - 3 * (na - 1)
    do j = 1, 3 * nat
       nb = (j - 1) / 3 + 1
       jcar = j - 3 * (nb - 1)
       phi (icar, jcar, na, nb) = fcq (i, j)
       print *, "PHI1:", icar, jcar, na,nb, phi(icar, jcar, na, nb)
    enddo
  enddo

  ! Convert to cartesian coordinates
  do na = 1, nat
    do nb = 1, nat
      call trntnsc(phi(1, 1, na, nb), at, bg, -1)
    end do
  end do

  !
  !    We start by imposing hermiticity
  !
  do na = 1, nat
    do nb = 1, nat
      do ipol = 1, 3
        do jpol = 1, 3
          phi(ipol, jpol, na, nb) = 0.5d0*(phi(ipol, jpol, na, nb) &
                                           + CONJG(phi(jpol, ipol, nb, na)))
          phi(jpol, ipol, nb, na) = CONJG(phi(ipol, jpol, na, nb))
        end do
      end do
    end do
  end do
  !
  !    If no other symmetry is present we quit here
  !
  if ((nsymq == 1) .and. (.not. minus_q)) return
  !
  !    Then we impose the symmetry q -> -q+G if present
  !
  if (minus_q) then
    do na = 1, nat
      do nb = 1, nat
        do ipol = 1, 3
          do jpol = 1, 3
            work(:, :) = (0.d0, 0.d0)
            sna = irt(irotmq, na)
            snb = irt(irotmq, nb)
            arg = 0.d0
            do kpol = 1, 3
              arg = arg + (xq(kpol)*(rtau(kpol, irotmq, na) - &
                                     rtau(kpol, irotmq, nb)))
            end do
            arg = arg*tpi
            fase = CMPLX(cos(arg), sin(arg), kind(1.0D0))
            do kpol = 1, 3
              do lpol = 1, 3
                work(ipol, jpol) = work(ipol, jpol) + &
                                   s(ipol, kpol, irotmq)*s(jpol, lpol, irotmq) &
                                   *phi(kpol, lpol, sna, snb)*fase
              end do
            end do
            phip(ipol, jpol, na, nb) = (phi(ipol, jpol, na, nb) + &
                                        CONJG(work(ipol, jpol)))*0.5d0
          end do
        end do
      end do
    end do
    phi = phip
  end if

  !
  !    Here we symmetrize with respect to the small group of q
  !
  if (nsymq == 1) return

  iflb(:, :) = 0
  do na = 1, nat
    do nb = 1, nat
      if (iflb(na, nb) == 0) then
        work(:, :) = (0.d0, 0.d0)
        do isymq = 1, nsymq
          irot = isymq
          sna = irt(irot, na)
          snb = irt(irot, nb)
          arg = 0.d0
          do ipol = 1, 3
            arg = arg + (xq(ipol)*(rtau(ipol, irot, na) - &
                                   rtau(ipol, irot, nb)))
          end do
          arg = arg*tpi
          faseq(isymq) = CMPLX(cos(arg), sin(arg), kind(1.0D0))
          do ipol = 1, 3
            do jpol = 1, 3
              do kpol = 1, 3
                do lpol = 1, 3
                  IF (t_rev(isymq) == 1) THEN
                    work(ipol, jpol) = work(ipol, jpol) + &
                                       s(ipol, kpol, irot)*s(jpol, lpol, irot) &
                                       *CONJG(phi(kpol, lpol, sna, snb)*faseq(isymq))
                  ELSE
                    work(ipol, jpol) = work(ipol, jpol) + &
                                       s(ipol, kpol, irot)*s(jpol, lpol, irot) &
                                       *phi(kpol, lpol, sna, snb)*faseq(isymq)
                  END IF
                end do
              end do
            end do
          end do
        end do
        do isymq = 1, nsymq
          irot = isymq
          sna = irt(irot, na)
          snb = irt(irot, nb)
          do ipol = 1, 3
            do jpol = 1, 3
              phi(ipol, jpol, sna, snb) = (0.d0, 0.d0)
              do kpol = 1, 3
                do lpol = 1, 3
                  IF (t_rev(isymq) == 1) THEN
                    phi(ipol, jpol, sna, snb) = phi(ipol, jpol, sna, snb) &
                                                + s(ipol, kpol, invs(irot))*s(jpol, lpol, invs(irot)) &
                                                *CONJG(work(kpol, lpol)*faseq(isymq))
                  ELSE
                    phi(ipol, jpol, sna, snb) = phi(ipol, jpol, sna, snb) &
                                                + s(ipol, kpol, invs(irot))*s(jpol, lpol, invs(irot)) &
                                                *work(kpol, lpol)*CONJG(faseq(isymq))
                  END IF
                end do
              end do
            end do
          end do
          iflb(sna, snb) = 1
        end do
      end if
    end do
  end do
  phi(:, :, :, :) = phi(:, :, :, :)/DBLE(nsymq)

  ! Convert back to fractional coordinates
  do na = 1, nat
    do nb = 1, nat
      call trntnsc(phi(1, 1, na, nb), at, bg, +1)
    end do
  end do
  !CALL compact_dyn(nat, fcqsymm, phi)
  fcqsymm = phi

  return
end subroutine symdynph_gq_new

