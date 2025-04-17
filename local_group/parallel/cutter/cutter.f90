 PROGRAM cutter
IMPLICIT NONE
include 'mpif.h'
integer :: n,in,np,i,ii,iii,nt,nsplit,rfile,wfile
character(len=16) :: zahl
double precision, allocatable :: m(:),pos(:,:),vel(:,:)
integer :: dd
double precision :: tt
integer ierr, rankrank, nbnodes, namelen
character (len=MPI_MAX_PROCESSOR_NAME) :: name


call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rankrank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nbnodes, ierr)
call MPI_GET_PROCESSOR_NAME(name, namelen, ierr)	


OPEN(50,file='c.dat')
READ(50,*) n
READ(50,*) nt
CLOSE(50)


nsplit=INT(DBLE(n)/DBLE(nbnodes))+1

DO in=1,nsplit

i=(in-1)*nbnodes+rankrank+1

IF (i<=n) THEN

rfile=100+rankrank
wfile=200+rankrank
write(zahl,*) i
zahl=adjustl(zahl)
OPEN(rfile,file='fmod'//trim(zahl)//'.dat')

DO ii=0,nt
READ(rfile,*) np
READ(rfile,*) dd
READ(rfile,*) tt

If ((ii==0).AND.(in==1)) THEN
allocate(m(1:np))
allocate(pos(1:3,1:np))
allocate(vel(1:3,1:np))
END IF


DO iii=1,np
READ(rfile,*) m(iii)
END DO

DO iii=1,np
READ(rfile,*) pos(1:3,iii)
END DO

DO iii=1,np
READ(rfile,*) vel(1:3,iii)
END DO


If (ii==0) THEN
OPEN(wfile,file='mod'//trim(zahl)//'-0.dat')

WRITE(wfile,*) np
WRITE(wfile,*) 3
WRITE(wfile,*) 0.D0

DO iii=1,np
WRITE (wfile,*) m(iii)
END DO

DO iii=1,np
WRITE (wfile,*) pos(1:3,iii)
END DO

DO iii=1,np
WRITE (wfile,*) vel(1:3,iii)
END DO

CLOSE(wfile)
END IF

If (ii==nt) THEN
OPEN(wfile,file='mod'//trim(zahl)//'-F.dat')

WRITE(wfile,*) np
WRITE(wfile,*) 3
WRITE(wfile,*) DBLE(nt)

DO iii=1,np
WRITE(wfile,*) m(iii)
END DO

DO iii=1,np
WRITE(wfile,*) pos(1:3,iii)
END DO

DO iii=1,np
WRITE(wfile,*) vel(1:3,iii)
END DO


CLOSE(wfile)
END IF



END DO

CLOSE(rfile)

END IF

END DO


call MPI_FINALIZE(ierr)


END PROGRAM