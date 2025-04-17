PROGRAM convergence
IMPLICIT NONE
integer :: n,i,ii,iii,best,gal,ng,nb
character(len=16) :: zahl
double precision :: wert

OPEN(49,file='p.dat')
READ(49,*) n
READ(49,*) nb
READ(49,*) ng
READ(49,*) best
READ(49,*) gal
CLOSE(49)

OPEN(50,file='diagram.dat')

DO i=1,(n+1)

IF (i<(n+1)) THEN
write(zahl,*) i
zahl=adjustl(zahl)
ELSE
zahl=''
END IF

OPEN(51,file='convergence'//trim(zahl)//'.dat')

DO ii=1,nb

DO iii=0,ng

READ(51,*) wert

IF ((ii==best).AND.(iii==gal)) THEN
WRITE(50,*) wert
END IF

END DO

READ(51,*)

END DO

CLOSE(51)
END DO

CLOSE(50)

END PROGRAM 
