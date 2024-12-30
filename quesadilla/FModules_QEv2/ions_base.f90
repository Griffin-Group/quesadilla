!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
  MODULE ions_base
!------------------------------------------------------------------------------!
      !! Ions configuration management.
      !
      USE parameters, ONLY : ntypx
      !USE uspp_param, ONLY : nsp
      !
      IMPLICIT NONE
      SAVE

      !     nsp       = number of species
      !     na(is)    = number of atoms of species is
      !     nax       = max number of atoms of a given species
      !     nat       = total number of atoms of all species
 
      INTEGER :: nax      = 0
      INTEGER :: nat      = 0
      INTEGER :: na(ntypx)= 0

      !     ityp( i ) = the type of i-th atom in stdin

      INTEGER,  ALLOCATABLE :: ityp(:)

      !     zv(is)    = (pseudo-)atomic charge
      !     amass(is) = mass of ions, in atomic mass units
      !     rcmax(is) = Ewald radius (for ion-ion interactions)

      double precision :: zv(ntypx)    = 0.0D0
      double precision :: amass(ntypx) = 0.0D0
      double precision :: rcmax(ntypx) = 0.0D0

      !     ityp( i ) = the type of i-th atom in stdin
      !     atm( j )  = name of the type of the j-th atomic specie
      !     tau( 1:3, i ) = position of the i-th atom

      double precision, ALLOCATABLE :: tau(:,:)     !  initial positions read from stdin (in bohr)
      double precision, ALLOCATABLE :: vel(:,:)     !  initial velocities read from stdin (in bohr)
      CHARACTER(LEN=6)      :: atm( ntypx )
      CHARACTER(LEN=80)     :: tau_format   ! format of input atomic positions:
                                            ! 'alat','crystal','bohr','angstrom'

      ! if_pos( x, i ) = 0 : x coordinate of i-th atom will be kept fixed
      INTEGER, ALLOCATABLE :: if_pos(:,:)  ! allowed values: 0 or 1 only
      INTEGER, ALLOCATABLE :: iforce(:,:)  ! if_pos sorted by specie 
      INTEGER :: fixatom   = 0            ! number of frozen atoms
      INTEGER :: ndofp     =-1            ! ionic degree of freedom
      INTEGER :: ndfrz     = 0            ! frozen degrees of freedom

      double precision :: fricp   ! friction parameter for damped dynamics
      double precision :: greasp  ! friction parameter for damped dynamics

      ! ... taui = real ionic positions in the center of mass reference
      ! ... system at istep = 0
      ! ... this array is used to compute mean square displacements,
      ! ... it is initialized when NBEG = -1, NBEG = 0 and TAURDR = .TRUE.
      ! ... first index: x,y,z, second index: atom sorted by specie with respect input
      ! ... this array is saved in the restart file

      double precision, ALLOCATABLE :: taui(:,:)

      ! ... cdmi = center of mass reference system (related to the taui)
      ! ... this vector is computed when NBEG = -1, NBEG = 0 and TAURDR = .TRUE.
      ! ... this array is saved in the restart file

      double precision :: cdmi(3), cdm(3)

      ! ... cdms = center of mass computed for scaled positions (taus)

      double precision :: cdms(3)
      !
      double precision, ALLOCATABLE :: extfor(:,:)     !  external forces on atoms


      
!------------------------------------------------------------------------------!
  END MODULE ions_base
!------------------------------------------------------------------------------!
