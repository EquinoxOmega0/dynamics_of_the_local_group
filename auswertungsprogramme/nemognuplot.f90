PROGRAM nemognuplot
IMPLICIT NONE

double precision :: x,y,z
character(len=16) :: filename
integer :: time,i,ii,n

WRITE(*,*) 'Dateiname:'
READ(*,*) filename
WRITE(*,*) 'Zeitschritt:'
READ(*,*) time
OPEN(51,file=filename)
OPEN(52,file='gnu.dat')

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


END PROGRAM 
