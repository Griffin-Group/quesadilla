!
! Copyright (C) 2002-2023 Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
subroutine errore(calling_routine, message, ierr)
  !----------------------------------------------------------------------------
  implicit none
  !
  character(LEN=*), intent(IN) :: calling_routine, message
  ! the name of the calling calling_routine
  ! the output message
  integer, intent(IN) :: ierr
  ! the error flag
  integer, external :: find_free_unit
  character(LEN=6) :: cerr
  !
  if (ierr <= 0) return
  !
  ! ... the error message is written on the "*" unit
  !
  write (cerr, FMT='(I6)') ierr
  write (UNIT=*, FMT='(/,1X,78("%"))')
  write (UNIT=*, FMT='(5X,"Error in routine ",A," (",A,"):")') &
    trim(calling_routine), trim(adjustl(cerr))
  write (UNIT=*, FMT='(5X,A)') trim(message)
  write (UNIT=*, FMT='(1X,78("%"),/)')
  !
  write (*, '("     stopping ...")')
  !
  call flush ()
  !
  stop 1
  !
  return
  !
end subroutine errore
!
!----------------------------------------------------------------------
subroutine infomsg(routine, message)
  !----------------------------------------------------------------------
  !
  ! ... This is a simple routine which writes an info message
  ! ... from a given routine to output.
  !
  implicit none
  !
  character(LEN=*) :: routine, message
  ! the name of the calling routine
  ! the output message
  !
  !
  write (*, '(5X,"Message from routine ",A,":")') routine
  write (*, '(5X,A)') message
  !
  !
  return
  !
end subroutine infomsg
!-----------------------------------------------------------------------
function int_to_char(i)
  !-----------------------------------------------------------------------
  !! Converts an integer number of up to 6 figures into a left-justifed
  !! character variable.
  !
  implicit none
  !
  integer, intent(IN) :: i
  character(LEN=6)   :: int_to_char
  character :: c
  integer   :: n, j, nc
  logical   :: neg
  !
  nc = 6
  !
  if (i < 0) then
    nc = nc - 1
    n = -i
    neg = .true.
  else
    n = i
    neg = .false.
  end if
  !
  j = 1
  do while (j <= nc)
    int_to_char(j:j) = char(mod(n, 10) + ichar('0'))
    n = n/10
    if (n == 0) exit
    j = j + 1
  end do
  !
  if (j <= nc) then
    do n = 1, j/2
      c = int_to_char(n:n)
      int_to_char(n:n) = int_to_char(j - n + 1:j - n + 1)
      int_to_char(j - n + 1:j - n + 1) = c
    end do
    if (j < nc) int_to_char(j + 1:nc) = ' '
  else
    int_to_char(:) = '*'
  end if
  !
  if (neg) then
    do n = nc + 1, 2, -1
      int_to_char(n:n) = int_to_char(n - 1:n - 1)
    end do
    int_to_char(1:1) = '-'
  end if
  !
  return
  !
end function int_to_char
