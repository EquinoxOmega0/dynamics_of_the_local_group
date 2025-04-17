Program mkinput
! declaration of variables
IMPLICIT NONE
integer :: i,ii,iii,nori,nm,nmodels,nset
double precision :: dbvar,dbx,dby,dbz,PI,tol,L,rmin,vstart,mass!,xcm,ycm,zcm
double precision :: mindist,distvar,streu,polya,polyb1,polyb2,anglediv,orimass,opdist,ovdist
double precision :: oxsize,oysize,ozsize,ovxsize,ovysize,ovzsize
double precision :: mxsize,mysize,mzsize,mvxsize,mvysize,mvzsize
character(len=16) :: filename,modname,zahl
double precision, allocatable :: masses(:)!,pos(:,:)
!intialisation of parameters
PI=ACOS(-1.D0)


OPEN(55,file='origin.dat')
READ(55,*) nm
READ(55,*) nori
READ(55,*) nmodels
READ(55,*) nset
allocate(masses(1:nm))
DO iii=1,nm
READ(55,*) masses(iii)
END DO
READ(55,*) mindist
READ(55,*) distvar
READ(55,*) polya,polyb1,polyb2
READ(55,*) streu
READ(55,*) anglediv
READ(55,*) mxsize,mysize,mzsize
READ(55,*) mvxsize,mvysize,mvzsize
READ(55,*) orimass
READ(55,*) oxsize,oysize,ozsize
READ(55,*) ovxsize,ovysize,ovzsize
READ(55,*) opdist
READ(55,*) ovdist

tol=4.D-2
!allocate(pos(1:3,1:(nm+nori)))


call random_seed()
DO i=1,nmodels
WRITE(zahl,*) i
zahl=adjustl(zahl)
filename='old'//trim(zahl)//'.dat'

OPEN(50,file=filename)
!number of particles
WRITE(50,*) nm
WRITE(50,*) nori
WRITE(50,*) nset
! name for model
modname='mod'//trim(zahl)//'.dat'

WRITE(50,*) modname
! masses

DO iii=1,nm
WRITE(50,*) masses(iii)
END DO



! initial distance (x*3kpc)
call random_number(dbvar)
rmin=dbvar*distvar+mindist
WRITE(50,*) rmin



call random_number(dbvar) !initial velocity

rmin=rmin/(1.D0+masses(1)/masses(2))
vstart=polya+polyb1*rmin+polyb2*rmin**2
vstart=vstart*(1.D0+masses(1)/masses(2))
vstart=vstart+dbvar*streu*2.D0-streu

WRITE(50,*) vstart



! angle of velocity
call random_number(dbvar)
dbvar=dbvar*anglediv*2.D0+90.D0-anglediv
WRITE(50,*) dbvar



DO ii=1,(nm-2)



call random_number(dbvar)
call random_number(dbx)
dbx=dbx*mxsize
IF (dbvar>0.5D0) THEN
dbx=-dbx
ENDIF


call random_number(dbvar)
dby=dby*mysize
call random_number(dby)
IF (dbvar>0.5D0) THEN
dby=-dby
ENDIF

call random_number(dbvar)
dbz=dbz*mzsize
call random_number(dbz)
IF (dbvar>0.5D0) THEN
dbz=-dbz
ENDIF



WRITE(50,*) dbx,dby,dbz



call random_number(dbvar)
call random_number(dbx)
dbx=dbx*mvxsize
IF (dbvar>0.5D0) THEN
dbx=-dbx
ENDIF


call random_number(dbvar)
call random_number(dby)
dby=dby*mvysize
IF (dbvar>0.5D0) THEN
dby=-dby
ENDIF

call random_number(dbvar)
call random_number(dbz)
dbz=dbz*mvzsize
IF (dbvar>0.5D0) THEN
dbz=-dbz
ENDIF

WRITE(50,*) dbx,dby,dbz


END DO


WRITE(50,*) orimass


DO ii=1,nori

!102 CONTINUE

call random_number(dbvar)
call random_number(dbx)
dbx=dbx*oxsize
IF (dbvar>0.5D0) THEN
dbx=-dbx
ENDIF


call random_number(dbvar)
call random_number(dby)
dby=dby*oysize
IF (dbvar>0.5D0) THEN
dby=-dby
ENDIF

call random_number(dbvar)
dby=dby*ozsize
call random_number(dbz)
IF (dbvar>0.5D0) THEN
dbz=-dbz
ENDIF

WRITE(50,*) dbx,dby,dbz
!pos(1,ii+nm)=dbx
!pos(2,ii+nm)=dby
!pos(3,ii+nm)=dbz

!DO iii=1,(ii+nm)

!IF (((pos(1,ii+nm)-pos(1,iii))**2+(pos(2,ii+nm)-pos(2,iii))**2+(pos(3,ii+nm)-pos(3,iii))**2)<tol) THEN

!IF ((ii+nm).NE.iii) THEN
!GOTO 102
!ENDIF
!ENDIF
!END DO

call random_number(dbvar)
call random_number(dbx)
dbx=dbx*ovxsize
IF (dbvar>0.5D0) THEN
dbx=-dbx
ENDIF


call random_number(dbvar)
call random_number(dby)
dby=dby*ovysize
IF (dbvar>0.5D0) THEN
dby=-dby
ENDIF

call random_number(dbvar)
call random_number(dbz)
dbz=dbz*ovzsize
IF (dbvar>0.5D0) THEN
dbz=-dbz
ENDIF

WRITE(50,*) dbx,dby,dbz

WRITE(50,*) opdist

WRITE(50,*) ovdist


END DO






CLOSE(50)


END DO


CLOSE(55)


END PROGRAM 
