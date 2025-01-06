!
! Copyright (C) 2001 - 2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine smallg_q(xq, modenum, at, nrot, s, sym, minus_q)
  !-----------------------------------------------------------------------
  !
  ! This routine selects, among the symmetry matrices of the point group
  ! of a crystal, the symmetry operations which leave q unchanged.
  ! Furthermore it checks if one of the above matrices send q --> -q+G.
  ! In this case minus_q is set true.
  !
  !  input-output variables
  !
  use symm_base, only: t_rev

  implicit none

  double precision, parameter :: accep = 1.d-5

  double precision, intent(in) :: at(3, 3), xq(3)
  ! input: the reciprocal lattice vectors
  ! input: the direct lattice vectors
  ! input: the q point of the crystal

  integer, intent(in) :: s(3, 3, 48), nrot, modenum
  ! input: the symmetry matrices
  ! input: number of symmetry operations
  ! input: main switch of the program, used for
  !        q<>0 to restrict the small group of q
  !        to operation such that Sq=q (exactly,
  !        without G vectors) when iswitch = -3.
  logical, intent(inout) :: sym(48), minus_q
  ! input-output: .true. if symm. op. S q = q + G
  ! output: .true. if there is an op. sym.: S q = - q + G
  !
  !  local variables
  !

  double precision :: aq(3), raq(3), zero(3)
  ! q vector in crystal basis
  ! the rotated of the q vector
  ! the zero vector

  integer :: irot, ipol, jpol
  ! counter on symmetry op.
  ! counter on polarizations
  ! counter on polarizations

  logical :: eqvect
  ! logical function, check if two vectors are equa
  !
  ! return immediately (with minus_q=.true.) if xq=(0,0,0)
  !
  minus_q = .true.
  if ((xq(1) == 0.d0) .and. (xq(2) == 0.d0) .and. (xq(3) == 0.d0)) &
    return
  !
  !   Set to zero some variables
  !
  minus_q = .false.
  zero(:) = 0.d0
  !
  !   Transform xq to the crystal basis
  !
  aq = xq
  call cryst_to_cart(1, aq, at, -1)
  !
  !   Test all symmetries to see if this operation send Sq in q+G or in -q+G
  !
  do irot = 1, nrot
    if (.not. sym(irot)) goto 100
    raq(:) = 0.d0
    do ipol = 1, 3
      do jpol = 1, 3
        raq(ipol) = raq(ipol) + dble(s(ipol, jpol, irot))*aq(jpol)
      end do
    end do
    if (t_rev(irot) == 1) raq = -raq
    sym(irot) = eqvect(raq, aq, zero, accep)
    !
    !  if "iswitch.le.-3" (modenum.ne.0) S must be such that Sq=q exactly !
    !
    if (modenum .ne. 0 .and. sym(irot)) then
      do ipol = 1, 3
        sym(irot) = sym(irot) .and. (abs(raq(ipol) - aq(ipol)) < 1.0d-5)
      end do
    end if
    !     if (.not.minus_q) then
    if (sym(irot) .and. .not. minus_q) then
      raq = -raq
      minus_q = eqvect(raq, aq, zero, accep)
    end if
100 continue
  end do
  !
  !  if "iswitch.le.-3" (modenum.ne.0) time reversal symmetry is not included !
  !
  if (modenum .ne. 0) minus_q = .false.
  !
  return
  !
end subroutine smallg_q
!-----------------------------------------------------------------------
subroutine set_giq(xq, lgamma, bg, at, s, nsymq, nsym, irotmq, minus_q, gi, gimq)
  !-----------------------------------------------------------------------
  !
  ! This routine calculates the possible vectors G associated
  ! to the symmetries of the small group of q: Sq -> q + G
  ! Furthermore if minus_q and irotmq are set it finds the G for Sq -> -q+G.
  !
  !USE kinds, ONLY : DP
  !USE cell_base, ONLY : bg, at
  !USE control_lr, ONLY : lgamma
  use symm_base, only: t_rev

  implicit none

  double precision, parameter :: accep = 1.d-5

  double precision, intent(IN) :: xq(3), at(3, 3), bg(3, 3)
  ! input: the q point, the lattice vectors, the reciprocal lattice vectors
  logical, intent(IN) :: lgamma
  double precision, intent(OUT) ::gi(3, 48), gimq(3)
  ! output: the G associated to a symmetry:[S(irotq)*q - q]
  ! output: the G associated to:  [S(irotmq)*q + q]

  logical, intent(IN) :: minus_q
  ! input: .t. if there is sym.ops. such that Sq=-q+G
  integer, intent(IN) :: s(3, 3, 48), nsymq, nsym
  ! input: the symmetry matrices
  ! input: dimension of the small group of q

  integer, intent(OUT) :: irotmq
  ! input: op. symmetry: s_irotmq(q)=-q+G

  double precision :: wrk(3), aq(3), raq(3), zero(3)
  ! additional space to compute gi and gimq
  ! q vector in crystal basis
  ! the rotated of the q vector
  ! the zero vector

  integer :: isym, ipol, jpol
  ! counter on symmetry operations
  ! counter on polarizations
  ! counter on polarizations

  logical :: eqvect
  ! logical function, check if two vectors are equal
  !
  !  Set to zero some variables and transform xq to the crystal basis
  !
  zero = 0.d0
  gi = 0.d0
  gimq = 0.d0
  irotmq = 0
  if (lgamma) then
    irotmq = 1
    return
  end if
  aq = xq
  call cryst_to_cart(1, aq, at, -1)
  !
  !   test all symmetries to see if the operation S sends q in q+G ...
  !
  do isym = 1, nsymq
    raq = 0.d0
    do ipol = 1, 3
      do jpol = 1, 3
        raq(ipol) = raq(ipol) + dble(s(ipol, jpol, isym))* &
                    aq(jpol)
      end do
    end do
    if (t_rev(isym) == 1) raq = -raq
    if (.not. eqvect(raq, aq, zero, accep)) call errore('set_giq', &
                                                        'problems with the input group', 1)
    do ipol = 1, 3
      if (t_rev(isym) == 1) then
        wrk(ipol) = aq(ipol) - raq(ipol)
      else
        wrk(ipol) = raq(ipol) - aq(ipol)
      end if
    end do
    call cryst_to_cart(1, wrk, bg, 1)
    gi(:, isym) = wrk(:)
    if (irotmq == 0) then
      raq = -raq
      if (eqvect(raq, aq, zero, accep)) then
        irotmq = isym
        wrk = aq - raq
        call cryst_to_cart(1, wrk, bg, 1)
        gimq = wrk
      end if
    end if
  end do
  !
  !   ... and in -q+G
  !
  if (minus_q .and. irotmq == 0) then
    do isym = nsymq + 1, nsym
      raq = 0.d0
      do ipol = 1, 3
        do jpol = 1, 3
          raq(ipol) = raq(ipol) + dble(s(ipol, jpol, isym))* &
                      aq(jpol)
        end do
      end do
      raq = -raq
      if (eqvect(raq, aq, zero, accep)) then
        wrk = aq - raq
        call cryst_to_cart(1, wrk, bg, 1)
        gimq(:) = wrk(:)
        irotmq = isym
      end if
      if (irotmq /= 0) exit
    end do
  end if
  if (minus_q .and. irotmq == 0) &
    call errore('set_giq', 'problem with minus_q', 1)
  !
  return
end subroutine set_giq

!-----------------------------------------------------------------------
subroutine sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the vector rtau which contains for each
  !     atom and each rotation the vector S\tau_a - \tau_b, where
  !     b is the rotated a atom, given by the array irt. These rtau are
  !     non zero only if fractional translations are present.
  !
  implicit none
  !
  !     first the dummy variables
  !
  integer, intent(in) :: nsym, s(3, 3, 48), nat, irt(48, nat)
  ! nsym: number of symmetries of the point group
  ! s:    matrices of symmetry operations
  ! nat : number of atoms in the unit cell
  ! irt(n,m) = transformed of atom m for symmetry n
  double precision, intent(in) :: at(3, 3), bg(3, 3), tau(3, nat)
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  ! tau: coordinates of the atoms
  double precision, intent(out):: rtau(3, 48, nat)
  ! rtau: the direct translations
  !
  !    here the local variables
  !
  integer :: na, nb, isym, ipol
  ! counters on: atoms, symmetry operations, polarization
  double precision, allocatable :: xau(:, :)
  double precision :: ft(3)
  !
  allocate (xau(3, nat))
  !
  !   compute the atomic coordinates in crystal axis, xau
  !
  do na = 1, nat
    do ipol = 1, 3
      xau(ipol, na) = bg(1, ipol)*tau(1, na) + &
                      bg(2, ipol)*tau(2, na) + &
                      bg(3, ipol)*tau(3, na)
    end do
  end do
  !
  !    for each symmetry operation, compute the atomic coordinates
  !    of the rotated atom, ft, and calculate rtau = Stau'-tau
  !
  rtau(:, :, :) = 0.0d0
  do isym = 1, nsym
    do na = 1, nat
      nb = irt(isym, na)
      do ipol = 1, 3
        ft(ipol) = s(1, ipol, isym)*xau(1, na) + &
                   s(2, ipol, isym)*xau(2, na) + &
                   s(3, ipol, isym)*xau(3, na) - xau(ipol, nb)
      end do
      do ipol = 1, 3
        rtau(ipol, isym, na) = at(ipol, 1)*ft(1) + &
                               at(ipol, 2)*ft(2) + &
                               at(ipol, 3)*ft(3)
      end do
    end do
  end do
  !
  !    deallocate workspace
  !
  deallocate (xau)
  return
end subroutine sgam_lr
!