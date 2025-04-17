!Modul Para(meter) for "global" variables
MODULE Para  
IMPLICIT NONE
DOUBLE PRECISION, ALLOCATABLE :: pos(:,:),vel(:,:),masses(:),tbpos(:,:),tbvel(:,:)
DOUBLE PRECISION, ALLOCATABLE :: m(:),sqrtm(:)
INTEGER :: n,nm,ntb,nnn,hubbleon
DOUBLE PRECISION :: tstop,massless,tstep,a0,sqrta0,soft,toffset,goldenratio,hadjust
END MODULE


MODULE Forces
CONTAINS 

SUBROUTINE derivs(x,y,dydx)
 USE nrtype
 USE Para
 IMPLICIT NONE
 REAL(SP), INTENT(IN) :: x
 !Array von Variablen die abgleitet werden
 REAL(SP), DIMENSION(:), INTENT(IN) :: y
 !Ableitungen von y-Array nach x
 REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
 INTEGER :: u,uu
 DOUBLE PRECISION :: dist2,hsoft



DO u=1,(nm+1)
dydx((u-1)*6+1)=y((u-1)*6+4)
dydx((u-1)*6+2)=y((u-1)*6+5)
dydx((u-1)*6+3)=y((u-1)*6+6)
END DO


Do u=1,(nm+1)
dydx((u-1)*6+4)=0.D0
dydx((u-1)*6+5)=0.D0
dydx((u-1)*6+6)=0.D0
END DO

! MOND gravitation
DO u=1,nm
DO uu=1,nm
IF (u.NE.uu) THEN

hsoft=soft*(sqrtm(u)/sqrta0)**2

dist2=((y((u-1)*6+1)-y((uu-1)*6+1))**2+(y((u-1)*6+2)-y((uu-1)*6+2))**2+(y((u-1)*6+3)-y((uu-1)*6+3))**2+hsoft)


dydx((uu-1)*6+4)=dydx((uu-1)*6+4)+sqrtm(u)*(y((u-1)*6+1)-y((uu-1)*6+1))/dist2
dydx((uu-1)*6+5)=dydx((uu-1)*6+5)+sqrtm(u)*(y((u-1)*6+2)-y((uu-1)*6+2))/dist2
dydx((uu-1)*6+6)=dydx((uu-1)*6+6)+sqrtm(u)*(y((u-1)*6+3)-y((uu-1)*6+3))/dist2
ENDIF
END DO
END DO


Do u=1,nm
dydx((u-1)*6+4)=sqrta0*dydx((u-1)*6+4)
dydx((u-1)*6+5)=sqrta0*dydx((u-1)*6+5)
dydx((u-1)*6+6)=sqrta0*dydx((u-1)*6+6)
END DO

!Hubble Expansion
IF (hubbleon==1) THEN

DO u=1,nm
dydx((u-1)*6+4)=dydx((u-1)*6+4)+y((u-1)*6+1)*2.D0/(3.D0*(x**2))/goldenratio*hadjust
dydx((u-1)*6+5)=dydx((u-1)*6+5)+y((u-1)*6+2)*2.D0/(3.D0*(x**2))/goldenratio*hadjust
dydx((u-1)*6+6)=dydx((u-1)*6+6)+y((u-1)*6+3)*2.D0/(3.D0*(x**2))/goldenratio*hadjust
END DO
  
ENDIF



! MOND gravitation
DO u=1,nm

hsoft=soft*(sqrtm(u)/sqrta0)**2

dist2=((y((u-1)*6+1)-y(nm*6+1))**2+(y((u-1)*6+2)-y(nm*6+2))**2+(y((u-1)*6+3)-y(nm*6+3))**2+hsoft)


dydx(nm*6+4)=dydx(nm*6+4)+sqrtm(u)*(y((u-1)*6+1)-y(nm*6+1))/dist2
dydx(nm*6+5)=dydx(nm*6+5)+sqrtm(u)*(y((u-1)*6+2)-y(nm*6+2))/dist2
dydx(nm*6+6)=dydx(nm*6+6)+sqrtm(u)*(y((u-1)*6+3)-y(nm*6+3))/dist2
END DO


dydx(nm*6+4)=sqrta0*dydx(nm*6+4)
dydx(nm*6+5)=sqrta0*dydx(nm*6+5)
dydx(nm*6+6)=sqrta0*dydx(nm*6+6)



!Hubble Expansion
IF (hubbleon==1) THEN

dydx(nm*6+4)=dydx(nm*6+4)+y(nm*6+1)*2.D0/(3.D0*(x**2))/goldenratio*hadjust
dydx(nm*6+5)=dydx(nm*6+5)+y(nm*6+2)*2.D0/(3.D0*(x**2))/goldenratio*hadjust
dydx(nm*6+6)=dydx(nm*6+6)+y(nm*6+3)*2.D0/(3.D0*(x**2))/goldenratio*hadjust

ENDIF




END SUBROUTINE derivs

END MODULE




PROGRAM deep_MOND_integrator 
! DeMONi
! use numerical recipies
 USE nr
 USE nrtype
 USE nrutil
 USE Para
 USE Forces
 USE ode_path
!Declaration of Variables
IMPLICIT NONE
include 'mpif.h'
Integer :: i,ii
Double precision :: t0,t1,h1,hmin,genau,hnext
Double precision, allocatable :: hopti(:)
character(len=16) :: filein,fileout
Double precision :: comx,comy,comz,comvx,comvy,comvz,summass
REAL(SP), allocatable :: ystart(:)
integer ierr, rankrank, nbnodes, namelen,npara
Double precision, allocatable :: parellelpv(:),hilfsverteiler(:)
character (len=MPI_MAX_PROCESSOR_NAME) :: name
!EXTERNAL  derivs



call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rankrank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nbnodes, ierr)
call MPI_GET_PROCESSOR_NAME(name, namelen, ierr)	


!Defination of golden ratio
goldenratio=(1.D0+SQRT(5.D0))/2.D0


! Read parameters for integrator and general settings
OPEN(52,file='para.dat')
READ(52,*) filein
READ(52,*) fileout
READ(52,*) tstop
READ(52,*) tstep
READ(52,*) a0            
READ(52,*) soft
READ(52,*) massless
READ(52,*) toffset
READ(52,*) hubbleon
READ(52,*) hadjust
CLOSE(52)

!WRITE(*,*) 'read in'
 !unitmass = 10^11 Msol
 !unitspeed = 378.66876 km/s
 !unitlength = 3kpc
 !unittime = 7.74655044869528 Myr
 !unita0 = 1.54444 10^9 m/s²






sqrta0=SQRT(a0)


! Read data of particle masses, position and velocity from snapshot
OPEN(51,file=filein)

READ(51,*) n
READ(51,*)
READ(51,*)

nm=0

ALLOCATE(m(1:n))

DO i=1,n
READ(51,*) m(i)
IF (m(i)>massless) THEN
nm=nm+1
END IF
END DO


ALLOCATE(pos(1:3,1:nm))
ALLOCATE(vel(1:3,1:nm))
ALLOCATE(masses(1:nm))
ALLOCATE(sqrtm(1:nm))

ALLOCATE(ystart(1:(nm*6+6)))

DO i=1,nm
masses(i)=m(i)
sqrtm(i)=SQRT(masses(i))
END DO

ntb=n-nm

!WRITE(*,*) 
!WRITE(*,*) 'DEMONI'
!WRITE(*,*) 'is a deep-MOND-integrator'
!WRITE(*,*) 
!WRITE(*,*) 'number of massiv particles: ',nm
!WRITE(*,*) 'number of testparticles: ',ntb
!WRITE(*,*) 
!WRITE(*,*) 'integrating at time:'

ALLOCATE(tbpos(1:3,1:ntb))
ALLOCATE(tbvel(1:3,1:ntb))


DO i=1,nm
READ(51,*) pos(1:3,i)
END DO
DO i=1,ntb
READ(51,*) tbpos(1:3,i)
END DO

DO i=1,nm
READ(51,*) vel(1:3,i)
END DO
DO i=1,ntb
READ(51,*) tbvel(1:3,i)
END DO

CLOSE(51)


t1=toffset

ALLOCATE(hopti(1:ntb))
DO i=1,ntb
hopti(i)=1.D-2
END DO

hmin=1.D-6
genau=1.D-3





!calculate total mass
summass=0.D0
DO i=1,nm 
summass=summass+masses(i)
END DO

!calculate center of mass

comx=0.D0
comy=0.D0
comz=0.D0
comvx=0.D0
comvy=0.D0
comvz=0.D0

DO i=1,nm
comx=comx+pos(1,i)*masses(i)
comy=comy+pos(2,i)*masses(i)
comz=comz+pos(3,i)*masses(i)
comvx=comvx+vel(1,i)*masses(i)
comvy=comvy+vel(2,i)*masses(i)
comvz=comvz+vel(3,i)*masses(i)
END DO

comx=comx/summass
comy=comy/summass
comz=comz/summass
comvx=comvx/summass
comvy=comvy/summass
comvz=comvz/summass


!shift to center of mass

DO i=1,nm
pos(1,i)=pos(1,i)-comx
pos(2,i)=pos(2,i)-comy
pos(3,i)=pos(3,i)-comz
vel(1,i)=vel(1,i)-comvx
vel(2,i)=vel(2,i)-comvy
vel(3,i)=vel(3,i)-comvz
END DO

DO i=1,ntb
tbpos(1,i)=tbpos(1,i)-comx
tbpos(2,i)=tbpos(2,i)-comy
tbpos(3,i)=tbpos(3,i)-comz
tbvel(1,i)=tbvel(1,i)-comvx
tbvel(2,i)=tbvel(2,i)-comvy
tbvel(3,i)=tbvel(3,i)-comvz
END DO





if (rankrank==0) THEN
OPEN(52,file=fileout)

WRITE(52,*) n
WRITE(52,*) 3
WRITE(52,*) 0.D0

DO i=1,nm
WRITE(52,*) masses(i)
END DO
DO i=1,ntb
WRITE(52,*) massless
END DO

DO i=1,nm
WRITE(52,*) pos(1:3,i)
END DO

DO i=1,ntb
WRITE(52,*) tbpos(1:3,i)
END DO

DO i=1,nm
WRITE(52,*) vel(1:3,i)
END DO

DO i=1,ntb
WRITE(52,*) tbvel(1:3,i)
END DO
END IF


npara=INT(DBLE(ntb)/DBLE(nbnodes))+1

ALLOCATE(parellelpv(1:(6*npara)))
ALLOCATE(hilfsverteiler(1:(6*npara*nbnodes)))



DO i=1,(6*npara*nbnodes)
hilfsverteiler(i)=0.D0
END DO



101 CONTINUE
t0=t1
t1=t0+tstep

!WRITE(*,*) (t0-toffset)

DO ii=1,ntb
hilfsverteiler((ii-1)*6+1)=tbpos(1,ii)
hilfsverteiler((ii-1)*6+2)=tbpos(2,ii)
hilfsverteiler((ii-1)*6+3)=tbpos(3,ii)
hilfsverteiler((ii-1)*6+4)=tbvel(1,ii)
hilfsverteiler((ii-1)*6+5)=tbvel(2,ii)
hilfsverteiler((ii-1)*6+6)=tbvel(3,ii)
END DO




call MPI_Scatter(hilfsverteiler,(6*npara),MPI_DOUBLE_PRECISION,parellelpv,(6*npara),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)



 DO i=1,npara




DO ii=1,nm
ystart((ii-1)*6+1)=pos(1,ii)
ystart((ii-1)*6+2)=pos(2,ii)
ystart((ii-1)*6+3)=pos(3,ii)
ystart((ii-1)*6+4)=vel(1,ii)
ystart((ii-1)*6+5)=vel(2,ii)
ystart((ii-1)*6+6)=vel(3,ii)
END DO



 h1=hopti(i)
IF (h1<hmin) THEN
h1=hmin*10.D0
END IF
! nnn=i+nm

DO ii=1,6
ystart(nm*6+ii)=parellelpv((i-1)*6+ii)
END DO


CALL odeint(ystart,t0,t1,genau,h1,hmin,derivs,rkqs,hnext)


DO ii=1,6
parellelpv((i-1)*6+ii)=ystart(nm*6+ii)
END DO


	  hopti(i)=hnext
END DO

!end loop for all test particles



call MPI_AllGather(parellelpv,(6*npara),MPI_DOUBLE_PRECISION,hilfsverteiler,(6*npara),MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)



DO ii=1,ntb
tbpos(1,ii)=hilfsverteiler((ii-1)*6+1)
tbpos(2,ii)=hilfsverteiler((ii-1)*6+2)
tbpos(3,ii)=hilfsverteiler((ii-1)*6+3)
tbvel(1,ii)=hilfsverteiler((ii-1)*6+4)
tbvel(2,ii)=hilfsverteiler((ii-1)*6+5)
tbvel(3,ii)=hilfsverteiler((ii-1)*6+6)
END DO




DO ii=1,nm
pos(1,ii)=ystart((ii-1)*6+1)
pos(2,ii)=ystart((ii-1)*6+2)
pos(3,ii)=ystart((ii-1)*6+3)
vel(1,ii)=ystart((ii-1)*6+4)
vel(2,ii)=ystart((ii-1)*6+5)
vel(3,ii)=ystart((ii-1)*6+6)
END DO



!calculate center of mass

comx=0.D0
comy=0.D0
comz=0.D0
comvx=0.D0
comvy=0.D0
comvz=0.D0

DO i=1,nm
comx=comx+pos(1,i)*masses(i)
comy=comy+pos(2,i)*masses(i)
comz=comz+pos(3,i)*masses(i)
comvx=comvx+vel(1,i)*masses(i)
comvy=comvy+vel(2,i)*masses(i)
comvz=comvz+vel(3,i)*masses(i)
END DO

comx=comx/summass
comy=comy/summass
comz=comz/summass
comvx=comvx/summass
comvy=comvy/summass
comvz=comvz/summass


!shift to center of mass

DO i=1,nm
pos(1,i)=pos(1,i)-comx
pos(2,i)=pos(2,i)-comy
pos(3,i)=pos(3,i)-comz
vel(1,i)=vel(1,i)-comvx
vel(2,i)=vel(2,i)-comvy
vel(3,i)=vel(3,i)-comvz
END DO

DO i=1,ntb
tbpos(1,i)=tbpos(1,i)-comx
tbpos(2,i)=tbpos(2,i)-comy
tbpos(3,i)=tbpos(3,i)-comz
tbvel(1,i)=tbvel(1,i)-comvx
tbvel(2,i)=tbvel(2,i)-comvy
tbvel(3,i)=tbvel(3,i)-comvz
END DO




IF (rankrank==0) THEN

WRITE(52,*) n
WRITE(52,*) 3
WRITE(52,*) (t1-toffset)

DO i=1,n
WRITE(52,*) masses(i)
END DO

DO i=1,nm
WRITE(52,*) pos(1:3,i)
END DO

DO i=1,ntb
WRITE(52,*) tbpos(1:3,i)
END DO

DO i=1,nm
WRITE(52,*) vel(1:3,i)
END DO

DO i=1,ntb
WRITE(52,*) tbvel(1:3,i)
END DO

END IF



IF (t1<(tstop+toffset)) THEN
  GOTO 101
ENDIF

IF (rankrank==0) THEN
CLOSE(52)
END IF

!WRITE(*,*) (t1-toffset)
!WRITE(*,*) 
!WRITE(*,*) 'program finished'
call MPI_FINALIZE(ierr)





END PROGRAM