!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
subroutine q2qstar_out (phi, xq, at, bg, nat, ntyp, ityp, tau, m_loc, nqs,lgamma, dynqstar)
    !----------------------------------------------------------------------------
    !! A small utility that reads the first q from a dynamical matrix file (either
    !! xml or plain text), recomputes the system symmetry (starting from the lattice)
    !! and generates the star of q.
    !
    !! Useful for debugging and for producing the star of the wannier-phonon code output.
    !
    !! Syntax:  
    !!   \(\texttt{q2qstar.x}\) filein [fileout]
    !
    !! fileout default: rot_filein (old format) or rot_filein.xml (new format) 
    !
    !USE io_global,          ONLY : *!, ionode_id, ionode, 
    ! symmetry
    USE symm_base,          ONLY : s, invs, nsym, find_sym, set_sym_bl, irt, copy_sym, nrot, inverse_s
    ! for reading the dyn.mat.
    ! Provided by the user instead
    !USE cell_base,          ONLY : at, bg
    !USE ions_base,          ONLY : nat, ityp, ntyp => nsp, atm, tau, amass
    ! small group symmetry
    USE lr_symm_base,       ONLY : rtau, nsymq, minus_q, irotmq, gi, gimq
    !
    implicit none
    

    ! -------- INPUT VARIABLES ------------------------------------
    integer, intent(in) :: nat, ntyp, ityp(nat)
    ! number of atoms in the unit cell, and the type of each atom
    double complex, intent(in) :: phi(3, 3, nat, nat)
    ! the input dynamical matrix
    double precision, intent(in)      :: xq(3)
    double precision, intent(in)      :: tau(3, nat)
    ! the atomic positions in cartesian coordinates

    double precision, intent(in) :: at (3, 3), bg (3, 3), m_loc(3, nat)
    ! direct lattice vectors
    ! reciprocal lattice vectors
    ! Magnetic moment of each atom
    logical, intent(in) :: lgamma
    ! -------- OUTPUT VARIABLES -----------------------------------
    double complex, intent(out) :: dynqstar(48,3,3,nat,nat) 
    ! the output dynamical matrix at each q in the star
    integer, intent(out) :: nqs
    ! number of q in the star



    ! -------- LOCAL VARIABLES ------------------------------------
    !integer :: ierr, nargs
    !
    !integer               :: gi(3), gimq(3), irt(3, 48)
    !double precision :: s(3,3,48)
    double precision :: sxq(3, 48)
    double precision :: xqs(3,48)
    integer               :: isq (48)
    integer :: imq
    !integer :: nsym, nsymq

    logical :: sym(48)
    double complex :: d2( 3 * nat, 3 * nat)
    integer :: i,j, icar,jcar, na,nb
    ! ######################### reading ######################### 
    ! ntyp = ntypx
    !  CALL read_dyn_from_file (nqs, xqs, epsilon, lrigid,  &
    !      ntyp, nat, ibrav, celldm, at, atm, amass)
    !  !
    !  IF (ionode) CLOSE(unit=1)
    !  !
    !  xq = xqs(:,1)
    !  ALLOCATE(phi(3,3,nat,nat))
    !  ALLOCATE(tau(3,nat))
    !  ALLOCATE(ityp(nat))
    !  phi  = dq_phiq(:,:,:,:,1)
    !  tau =  dq_tau
    !  ityp = dq_ityp
    !  !zeu =  dq_zeu ! note: zeu from dynamicalq is a real(dp) array, zeu from control_ph is a flag (logical)
    !  amass = amass/amu_ry
    !  !
    !ENDIF XML_FORMAT_READ
    !
    ! regenerate the lattice
    ! CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
    ! at = at / celldm(1)  !  bring at in units of alat
    ! CALL volume(celldm(1),at(1,1),at(1,2),at(1,3),omega)
    ! CALL recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))


    WRITE(*,'(//,5x,a,3f14.9/)') "Dynamical matrix at q =", xq

    ! ~~~~~~~~ setup bravais lattice symmetry ~~~~~~~~ 
    CALL set_sym_bl (at, bg) 
    WRITE(*, '(5x,a,i3)') "Symmetries of bravais lattice: ", nrot
    
    ! ~~~~~~~~ find the symmetries of the crystal ~~~~~~~~
    ! TODO: provide this routine
    CALL find_sym (at, bg, nat, tau, ityp, .false., m_loc) 
    WRITE(*, '(5x,a,i3)') "Symmetries of crystal:         ", nsym
    !
    ! ~~~~~~~~ setup small group of q symmetry ~~~~~~~~ 
    ! part 1: call smallg_q and the copy_sym, 
    minus_q = .true.
    sym = .false.
    sym(1:nsym) = .true.
    ! TODO: provide this routine
    CALL smallg_q(xq, 0, at, bg, nsym, s, sym, minus_q) 
    nsymq = copy_sym(nsym, sym) ! TODO: provide this routine

    ! recompute the inverses as the order of sym.ops. has changed
    ! TODO: provide this routine
    CALL inverse_s ( )  

    ! part 2: this computes gi, gimq
    call set_giq (xq,s,nsymq,nsym,minus_q,gi,gimq,lgamma) ! TODO: provide this routine
    WRITE(*, '(5x,a,i3)') "Symmetries of small group of q:", nsymq
    IF(minus_q) WRITE(*, '(10x,a)') "in addition sym. q -> -q+G:"
    !
    ! finally this does some of the above again and also computes rtau...
    ! TODO: provide this routine
    ALLOCATE(rtau( 3, 48, nat))
    CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
    !
    ! ######################### star of q #########################
    do na = 1, nat
       do nb = 1, nat
          call trntnsc (phi (1, 1, na, nb), at, bg, - 1) ! TODO: provide this routine
       enddo
    enddo
    CALL symdynph_gq_new (xq, phi, s, invs, rtau, irt, nsymq, nat, &
         irotmq, minus_q) ! TODO: provide this routine
    do na = 1, nat
       do nb = 1, nat
          call trntnsc (phi (1, 1, na, nb), at, bg, + 1) ! TODO: provide this routine
       enddo
    enddo
    !
    CALL star_q(xq, at, bg, nsym, s, invs, nqs, sxq, isq, imq, .true. ) ! TODO: provide this routine
    !
    !XML_FORMAT_WRITE : &
    !IF (xmldyn) THEN
    !   nqq=nqs
    !   IF (imq==0) nqq=2*nqs
  ! !     IF (lgamma.AND.done_epsil.AND.done_zeu) THEN
  ! !        CALL write_dyn_mat_header( fildyn, ntyp, nat, ibrav, nspin_mag, &
  ! !             celldm, at, bg, omega, atm, amass, tau, ityp, m_loc, &
  ! !             nqq, epsilon, zstareu, lraman, ramtns)
  ! !     ELSE
    !      CALL write_dyn_mat_header( filout, ntyp, nat, ibrav, nspin_mag, &
    !           celldm, at, bg, omega, atm, amass, tau,ityp,m_loc,nqq)
  ! !     ENDIF
    !ELSE XML_FORMAT_WRITE
    !    OPEN (unit=1, file=filout,status='unknown',form='formatted',iostat=ierr)
    !    IF (ierr /= 0) CALL errore(CODE,'opening output file',1)
    !    CALL write_old_dyn_mat_head(1)
    !ENDIF XML_FORMAT_WRITE
    !
    ! repack phi to 3*nat,3*nat so that it can be repacked and then rerepacked again in q2qstar_ph
    DO i = 1, 3 * nat
      na = (i - 1) / 3 + 1
      icar = i - 3 * (na - 1)
      DO j = 1, 3 * nat
          nb = (j - 1) / 3 + 1
          jcar = j - 3 * (nb - 1)
          d2 (i, j) = phi(icar, jcar, na, nb)
      ENDDO
    ENDDO
    !
    CALL q2qstar_ph (d2, at, bg, nat, nsym, s, invs, irt, rtau, &
                     nqs, sxq, isq, imq, 1) ! TODO: provide this routine
  
    DEALLOCATE(rtau)
    DEALLOCATE(irt) ! from symm_base
end subroutine q2qstar_out