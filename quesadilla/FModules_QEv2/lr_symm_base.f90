!
! Copyright (C) 2001-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for LR_Modules routines
!
MODULE lr_symm_base
  !
  !USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the modes and the small group of q
  !
  SAVE
  !
  INTEGER :: irgq(48), nsymq=0, irotmq
  ! selects the operations of the small group
  ! the number of symmetry of the small group
  ! selects the symmetry sending q <-> -q+G
  double precision, ALLOCATABLE :: rtau(:,:,:) !3, 48, nat)
  ! coordinates of direct translations
  double precision :: gi(3,48), gimq(3)
  ! the possible G associated to each symmetry
  ! the G associated to the symmetry q<->-q+G
  LOGICAL :: minus_q, & ! if .TRUE. there is the symmetry sending q<->-q
             invsymq    ! if .TRUE. the small group of q has inversion
  !
  ! Symmetry representation of the perturbations
  !
  !INTEGER :: lr_npert
  !! Number of perturbations considered at the same time.
  !! e.g., for phonons: dimension of the irreducible representation
  !! e.g., for electric fields: 3
  !COMPLEX(DP), ALLOCATABLE :: upert(:, :, :)
  !! Representation of the symmetry in the perturbation basis. Size (lr_npert, lr_npert, nsymq)
  !! e.g., for phonons: transformation matrix of the patterns
  !! e.g., for electric fields: transformation matrix of Cartesian vectors
  !COMPLEX(DP), ALLOCATABLE :: upert_mq(:, :)
  !! Representation of the symmetry that transforms q to -q. Size (lr_npert, lr_npert)
  !
END MODULE lr_symm_base
!