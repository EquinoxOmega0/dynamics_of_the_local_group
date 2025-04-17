PROGRAM zerhacker
IMPLICIT NONE

double precision :: x,y,z
character(len=16) :: filename
integer :: i,n


  READ(*,*) n


filename='F2.dat'

!WRITE(*,*) 'Dateiname:'
!READ(*,*) filename
!WRITE(*,*) 'Zeitschritt:'
!READ(*,*) time
OPEN(51,file=filename)
OPEN(52,file='Ftborigin.dat')
OPEN(53,file='Fmorigin.dat')
OPEN(54,file='Fmmorigin.dat')


DO i=1,n
  
READ(51,*) x,y,z

IF (i<3) THEN
WRITE(54,*) x,y,z
ELSE IF (i<5) THEN
WRITE(53,*) x,y,z
ELSE
WRITE(52,*) x,y,z
END IF



END DO

CLOSE(54)
CLOSE(53)
CLOSE(52)
CLOSE(51)


END PROGRAM 
