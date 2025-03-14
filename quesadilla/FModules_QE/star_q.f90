!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine star_q(xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, verbosity)
  !-----------------------------------------------------------------------
  ! generate the star of q vectors that are equivalent to the input one
  ! NB: input s(:,:,1:nsym) must contain all crystal symmetries,
  ! i.e. not those of the small-qroup of q only
  !
  implicit none
  !
  integer, intent(in) :: nsym, s(3, 3, 48), invs(48)
  ! nsym matrices of symmetry operations
  ! invs: list of inverse operation indices
  double precision, intent(in) :: xq(3), at(3, 3), bg(3, 3)
  ! xq: q vector
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  !
  integer, intent(out) :: nq, isq(48), imq
  ! nq  : degeneracy of the star of q
  ! isq : index of q in the star for a given sym
  ! imq : index of -q in the star (0 if not present)

  double precision, intent(out) :: sxq(3, 48)
  ! list of vectors in the star of q
  logical, intent(in) :: verbosity
  ! if true prints several messages.
  !
  integer :: nsq(48), isym, ism1, iq, i
  ! number of symmetry ops. of bravais lattice
  ! counters on symmetry ops.
  ! index of inverse of isym
  ! counters
  double precision :: saq(3, 48), aq(3), raq(3), zero(3)
  ! auxiliary list of q (crystal coordinates)
  ! input q in crystal coordinates
  ! rotated q in crystal coordinates
  ! coordinates of fractionary translations
  ! a zero vector: used in eqvect

  logical, external :: eqvect
  ! function used to compare two vectors
  !
  zero(:) = 0.d0
  !
  ! go to  crystal coordinates
  !
  do i = 1, 3
    aq(i) = xq(1)*at(1, i) + xq(2)*at(2, i) + xq(3)*at(3, i)
  end do
  !
  ! create the list of rotated q
  !
  do i = 1, 48
    nsq(i) = 0
    isq(i) = 0
  end do
  nq = 0
  do isym = 1, nsym
    ism1 = invs(isym)
    do i = 1, 3
      raq(i) = s(i, 1, ism1)*aq(1) &
               + s(i, 2, ism1)*aq(2) &
               + s(i, 3, ism1)*aq(3)
    end do
    do i = 1, 3
      sxq(i, 48) = bg(i, 1)*raq(1) &
                   + bg(i, 2)*raq(2) &
                   + bg(i, 3)*raq(3)
    end do
    do iq = 1, nq
      if (eqvect(raq, saq(1, iq), zero)) then
        isq(isym) = iq
        nsq(iq) = nsq(iq) + 1
      end if
    end do
    if (isq(isym) == 0) then
      nq = nq + 1
      nsq(nq) = 1
      isq(isym) = nq
      saq(:, nq) = raq(:)
      do i = 1, 3
        sxq(i, nq) = bg(i, 1)*saq(1, nq) &
                     + bg(i, 2)*saq(2, nq) &
                     + bg(i, 3)*saq(3, nq)
      end do
    end if
  end do
  !
  ! set imq index if needed and check star degeneracy
  !
  raq(:) = -aq(:)
  imq = 0
  do iq = 1, nq
    if (eqvect(raq, saq(1, iq), zero)) imq = iq
    if (nsq(iq)*nq /= nsym) call errore('star_q', 'too many symmetries! (is this a supercell?)', iq)
  end do
  !
  ! writes star of q
  !
  if (verbosity) then
    write (*, *)
    write (*, '(5x,a,i4)') 'Number of q in the star = ', nq
    write (*, '(5x,a)') 'List of q in the star:'
    write (*, '(7x,i4,3f14.9)') (iq, (sxq(i, iq), i=1, 3), iq=1, nq)
    call flush ()
    if (imq == 0) then
      write (*, '(5x,a)') 'In addition there is the -q list: '
      write (*, '(7x,i4,3f14.9)') (iq, (-sxq(i, iq), i=1, 3), iq=1, nq)
      call flush ()
    end if
  end if
  return
end subroutine star_q
