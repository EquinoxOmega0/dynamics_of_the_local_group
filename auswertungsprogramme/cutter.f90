 PROGRAM cutter
IMPLICIT NONE
integer :: n,i,np,ii,iii,nt
character(len=16) :: zahl
double precision, allocatable :: pos(:,:),vel(:,:)

OPEN(50,file='c.dat')
READ(50,*) n
READ(50,*) nt
CLOSE(50)

DO i=1,n

write(zahl,*) i
zahl=adjustl(zahl)
OPEN(51,file='fmod'//trim(zahl)//'.dat')

DO ii=0,nt
READ(51,*) np
READ(51,*)
READ(51,*) 

If ((ii==0).AND.(i==1)) THEN
allocate(pos(1:3,1:np))
allocate(vel(1:3,1:np))
END IF


DO iii=1,np
READ (51,*) 
END DO

DO iii=1,np
READ (51,*) pos(1:3,iii)
END DO

DO iii=1,np
READ (51,*) vel(1:3,iii)
END DO


If (ii==0) THEN
OPEN(52,file='mod'//trim(zahl)//'-0.dat')

WRITE(52,*) np
WRITE(52,*) 3
WRITE(52,*) 0.D0

DO iii=1,np
WRITE (52,*) 0
END DO

DO iii=1,np
WRITE (52,*) pos(1:3,iii)
END DO

DO iii=1,np
WRITE (52,*) vel(1:3,iii)
END DO

CLOSE(52)
END IF

If (ii==nt) THEN
OPEN(52,file='mod'//trim(zahl)//'-F.dat')

WRITE(52,*) np
WRITE(52,*) 3
WRITE(52,*) DBLE(nt)

DO iii=1,np
WRITE (52,*) 0
END DO

DO iii=1,np
WRITE (52,*) pos(1:3,iii)
END DO

DO iii=1,np
WRITE (52,*) vel(1:3,iii)
END DO


CLOSE(52)
END IF



END DO

CLOSE(51)

END DO


END PROGRAM