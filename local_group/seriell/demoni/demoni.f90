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
Integer :: i,ii
Double precision :: t0,t1,h1,hmin,genau,hnext
Double precision, allocatable :: hopti(:)
character(len=16) :: filein,fileout
REAL(SP), allocatable :: ystart(:)
Double precision :: compos(1:3),comvel(1:3),mtot
!EXTERNAL  derivs

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

WRITE(*,*) 
WRITE(*,*) 'DEMONI'
WRITE(*,*) 'is a deep-MOND-integrator'
WRITE(*,*) 
WRITE(*,*) 'number of massiv particles: ',nm
WRITE(*,*) 'number of testparticles: ',ntb
WRITE(*,*) 
WRITE(*,*) 'integrating at time:'

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


!shift to center of mass

mtot=0.D0
compos(1)=0.D0
compos(2)=0.D0
compos(3)=0.D0
comvel(1)=0.D0
comvel(2)=0.D0
comvel(3)=0.D0

DO i=1,nm
mtot=mtot+masses(i)
compos(1)=compos(1)+pos(1,i)*masses(i)
compos(2)=compos(2)+pos(2,i)*masses(i)
compos(3)=compos(3)+pos(3,i)*masses(i)
comvel(1)=comvel(1)+vel(1,i)*masses(i)
comvel(2)=comvel(2)+vel(2,i)*masses(i)
comvel(3)=comvel(3)+vel(3,i)*masses(i)
END DO

compos(1)=compos(1)/mtot
compos(2)=compos(2)/mtot
compos(3)=compos(3)/mtot
comvel(1)=comvel(1)/mtot
comvel(2)=comvel(2)/mtot
comvel(3)=comvel(3)/mtot

DO i=1,n
pos(1,i)=pos(1,i)-compos(1)
pos(2,i)=pos(2,i)-compos(2)
pos(3,i)=pos(3,i)-compos(3)
vel(1,i)=vel(1,i)-comvel(1)
vel(2,i)=vel(2,i)-comvel(2)
vel(3,i)=vel(3,i)-comvel(3)
END DO


OPEN(52,file=fileout)

WRITE(52,*) n
WRITE(52,*) 3
WRITE(52,*) 0

DO i=1,nm
WRITE (52,*) masses(i)
END DO
DO i=1,ntb
WRITE (52,*) massless
END DO

DO i=1,nm
WRITE (52,*) pos(1:3,i)
END DO

DO i=1,ntb
WRITE (52,*) tbpos(1:3,i)
END DO

DO i=1,nm
WRITE (52,*) vel(1:3,i)
END DO

DO i=1,ntb
WRITE (52,*) tbvel(1:3,i)
END DO




101 CONTINUE
t0=t1
t1=t0+tstep

WRITE(*,*) (t0-toffset)


 DO i=1,ntb
 h1=hopti(i)

IF (h1<hmin) THEN
h1=hmin*10.D0
END IF
 nnn=i

DO ii=1,nm
ystart((ii-1)*6+1)=pos(1,ii)
ystart((ii-1)*6+2)=pos(2,ii)
ystart((ii-1)*6+3)=pos(3,ii)
ystart((ii-1)*6+4)=vel(1,ii)
ystart((ii-1)*6+5)=vel(2,ii)
ystart((ii-1)*6+6)=vel(3,ii)
END DO

ystart(nm*6+1)=tbpos(1,i)
ystart(nm*6+2)=tbpos(2,i)
ystart(nm*6+3)=tbpos(3,i)
ystart(nm*6+4)=tbvel(1,i)
ystart(nm*6+5)=tbvel(2,i)
ystart(nm*6+6)=tbvel(3,i)

CALL odeint(ystart,t0,t1,genau,h1,hmin,derivs,rkqs,hnext)

tbpos(1,i)=ystart(nm*6+1)
tbpos(2,i)=ystart(nm*6+2)
tbpos(3,i)=ystart(nm*6+3)
tbvel(1,i)=ystart(nm*6+4)
tbvel(2,i)=ystart(nm*6+5)
tbvel(3,i)=ystart(nm*6+6)

	  hopti(i)=hnext
END DO


DO ii=1,nm
pos(1,ii)=ystart((ii-1)*6+1)
pos(2,ii)=ystart((ii-1)*6+2)
pos(3,ii)=ystart((ii-1)*6+3)
vel(1,ii)=ystart((ii-1)*6+4)
vel(2,ii)=ystart((ii-1)*6+5)
vel(3,ii)=ystart((ii-1)*6+6)
END DO





!shift to center of mass

compos(1)=0.D0
compos(2)=0.D0
compos(3)=0.D0
comvel(1)=0.D0
comvel(2)=0.D0
comvel(3)=0.D0

DO i=1,nm
compos(1)=compos(1)+pos(1,i)*masses(i)
compos(2)=compos(2)+pos(2,i)*masses(i)
compos(3)=compos(3)+pos(3,i)*masses(i)
comvel(1)=comvel(1)+vel(1,i)*masses(i)
comvel(2)=comvel(2)+vel(2,i)*masses(i)
comvel(3)=comvel(3)+vel(3,i)*masses(i)
END DO

compos(1)=compos(1)/mtot
compos(2)=compos(2)/mtot
compos(3)=compos(3)/mtot
comvel(1)=comvel(1)/mtot
comvel(2)=comvel(2)/mtot
comvel(3)=comvel(3)/mtot

DO i=1,n
pos(1,i)=pos(1,i)-compos(1)
pos(2,i)=pos(2,i)-compos(2)
pos(3,i)=pos(3,i)-compos(3)
vel(1,i)=vel(1,i)-comvel(1)
vel(2,i)=vel(2,i)-comvel(2)
vel(3,i)=vel(3,i)-comvel(3)
END DO

! output results





WRITE(52,*) n
WRITE(52,*) 3
WRITE(52,*) (t1-toffset)


DO i=1,n
WRITE (52,*) masses(i)
END DO

DO i=1,nm
WRITE (52,*) pos(1:3,i)
END DO

DO i=1,ntb
WRITE (52,*) tbpos(1:3,i)
END DO

DO i=1,nm
WRITE (52,*) vel(1:3,i)
END DO

DO i=1,ntb
WRITE (52,*) tbvel(1:3,i)
END DO




IF (t1<(tstop+toffset)) THEN
  GOTO 101
ENDIF


CLOSE(52)

WRITE(*,*) (t1-toffset)
WRITE(*,*) 
WRITE(*,*) 'program finished'



END PROGRAM