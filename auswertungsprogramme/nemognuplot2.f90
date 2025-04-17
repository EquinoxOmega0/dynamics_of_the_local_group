PROGRAM nemognuplot
IMPLICIT NONE

double precision :: x,y,z
character(len=16) :: filename,zahl
integer :: time,i,ii,n,u

time=0

DO u=0,5

write(zahl,*) u
zahl=adjustl(zahl)

filename='rotated'//trim(zahl)//'.dat'

if (u==0) then
filename='vglnemo.dat'
end if

!WRITE(*,*) 'Dateiname:'
!READ(*,*) filename
!WRITE(*,*) 'Zeitschritt:'
!READ(*,*) time
OPEN(51,file=filename)
OPEN(52,file='gnu'//trim(zahl)//'.dat')

DO i=0,time
READ(51,*) n
READ(51,*) 
READ(51,*) 

DO ii=1,n
READ(51,*) 
END DO

DO ii=1,n
READ(51,*) x,y,z
IF (i==time) THEN
WRITE(52,*) x,y,z
END IF
END DO

DO ii=1,n
READ(51,*) 
END DO


END DO

CLOSE(52)
CLOSE(51)


END DO

END PROGRAM 
