module version
  implicit none

  contains

subroutine print_version

  write(*,*)
  write(*,'(6x,''* * * * * * * * * * * * * * * * * *'')')
  write(*,'(6x,''*           P l o t M S           *'')')
  write(*,'(6x,''* * * * * * * * * * * * * * * * * *'')')
  write(*,'(6x,''*    QCxMS spectra plotting tool  *'')')
  write(*,'(6x,''*         -  v. 6.1.0  -          *'')')
  write(*,'(6x,''*          25. Jul 2022           *'')')
  write(*,'(6x,''*           *   *   *             *'')')
  write(*,'(6x,''*           S. Grimme             *'')')
  write(*,'(6x,''*           J. Koopman            *'')')
  write(*,'(6x,''* * * * * * * * * * * * * * * * * *'')')
  write(*,'(6x,''   Contributor:  C.Bauer, T.Kind   '')')
  write(*,*)

end subroutine print_version
end module version
