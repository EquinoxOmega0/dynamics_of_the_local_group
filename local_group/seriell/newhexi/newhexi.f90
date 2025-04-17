!Modul Para(meter) for "global" variables
MODULE Para  
IMPLICIT NONE
DOUBLE PRECISION, ALLOCATABLE :: pos(:,:),vel(:,:),masses(:),tbpos(:,:),tbvel(:,:)
DOUBLE PRECISION, ALLOCATABLE :: m(:),rhalo(:)
INTEGER :: n,nm,ntb,nnn,haloon,dynfricton,hubbleon
DOUBLE PRECISION :: tstop,massless,tstep,soft,toffset,hPI,lncoulomb,scalefrict,goldenratio,hadjust
END MODULE

MODULE mathfunction
CONTAINS 

!Faktorielle
FUNCTION fact(bb)
IMPLICIT NONE
INTEGER , INTENT(IN) :: bb
INTEGER :: fact
INTEGER :: a,cc

a=1

DO cc=1,bb 
a=a*cc
END DO

fact=a

END FUNCTION


! Errorfunction
FUNCTION errorf(z)
USE Para
IMPLICIT NONE
DOUBLE PRECISION , INTENT(IN):: z
DOUBLE PRECISION :: errorf
DOUBLE PRECISION :: zwischen
INTEGER :: k



zwischen=0.D0

if ((z>2.D0).or.(z<-2.D0)) then

if (z>2.D0) then
zwischen=1.D0
endif

if (z<-2.D0) then
zwischen=-1.D0
endif


else

DO k=0,8

zwischen=zwischen+(((-1.D0)**k)*(z**(2*k+1))/(DBLE(2*k+1)*DBLE(fact(k))))

END DO

zwischen=2.D0/SQRT(hPI)*zwischen

if (zwischen>1.D0) then
zwischen=1.D0
endif

if (zwischen<-1.D0) then
zwischen=-1.D0
endif

end if

errorf=zwischen

END FUNCTION



END MODULE






MODULE Forces
CONTAINS 

SUBROUTINE derivs(x,y,dydx)
 USE nrtype
 USE Para
 USE mathfunction
 IMPLICIT NONE
 REAL(SP), INTENT(IN) :: x
 !Array von Variablen die abgleitet werden
 REAL(SP), DIMENSION(:), INTENT(IN) :: y
 !Ableitungen von y-Array nach x
 REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
 INTEGER :: u,uu
 DOUBLE PRECISION :: dist3,dist2,dist1,vrel,xpar,scalarp,beta
 DOUBLE PRECISION :: r(1:3),vr(1:3),ffrict(1:3)
 

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


! Newtonian Gravitation
DO u=1,nm
DO uu=1,nm
IF (u.NE.uu) THEN

r(1)=y((uu-1)*6+1)-y((u-1)*6+1)
r(2)=y((uu-1)*6+2)-y((u-1)*6+2)
r(3)=y((uu-1)*6+3)-y((u-1)*6+3)

IF (haloon==1) THEN

dist2=((r(1)**2)+(r(2)**2)+(r(3)**2)+(soft**2))


IF (dist2<(rhalo(uu)**2)) THEN

dydx((u-1)*6+4)=dydx((u-1)*6+4)+(r(1)*masses(uu)/(dist2*rhalo(uu)))
dydx((u-1)*6+5)=dydx((u-1)*6+5)+(r(2)*masses(uu)/(dist2*rhalo(uu)))
dydx((u-1)*6+6)=dydx((u-1)*6+6)+(r(3)*masses(uu)/(dist2*rhalo(uu)))

ELSE
  
dist1=SQRT(dist2)
dist3=dist1**3

dydx((u-1)*6+4)=dydx((u-1)*6+4)+(r(1)*masses(uu)/dist3)
dydx((u-1)*6+5)=dydx((u-1)*6+5)+(r(2)*masses(uu)/dist3)
dydx((u-1)*6+6)=dydx((u-1)*6+6)+(r(3)*masses(uu)/dist3)

END IF


ELSE
dist1=SQRT(((r(1)**2)+(r(2)**2)+(r(3)**2)+(soft**2)))
dist3=dist1**3

dydx((u-1)*6+4)=dydx((u-1)*6+4)+(r(1)*masses(uu)/dist3)
dydx((u-1)*6+5)=dydx((u-1)*6+5)+(r(2)*masses(uu)/dist3)
dydx((u-1)*6+6)=dydx((u-1)*6+6)+(r(3)*masses(uu)/dist3)
END IF

ENDIF

END DO
END DO

!Hubble Expansion
IF (hubbleon==1) THEN
DO u=1,nm
dydx((u-1)*6+4)=dydx((u-1)*6+4)+y((u-1)*6+1)*2.D0/(3.D0*(x**2))/goldenratio*hadjust
dydx((u-1)*6+5)=dydx((u-1)*6+5)+y((u-1)*6+2)*2.D0/(3.D0*(x**2))/goldenratio*hadjust
dydx((u-1)*6+6)=dydx((u-1)*6+6)+y((u-1)*6+3)*2.D0/(3.D0*(x**2))/goldenratio*hadjust
END DO
END IF



!Dynamical Friction
if (dynfricton==1) THEN
DO u=1,nm
DO uu=1,nm
IF (u.NE.uu) THEN
IF (masses(uu)>masses(u)) THEN

r(1)=y((uu-1)*6+1)-y((u-1)*6+1)
r(2)=y((uu-1)*6+2)-y((u-1)*6+2)
r(3)=y((uu-1)*6+3)-y((u-1)*6+3)

dist2=((r(1)**2)+(r(2)**2)+(r(3)**2)+(soft**2))


IF ((dist2<(rhalo(uu)**2)).AND.(dist2>(2.D0*(soft**2)))) THEN

vr(1)=y((u-1)*6+4)-y((uu-1)*6+4)
vr(2)=y((u-1)*6+5)-y((uu-1)*6+5)
vr(3)=y((u-1)*6+6)-y((uu-1)*6+6)

vrel=SQRT((vr(1)**2)+(vr(2)**2)+(vr(3)**2)+1.D-4)

xpar=vrel/SQRT(masses(uu)/rhalo(uu))

xpar=(errorf(xpar)-2.D0*xpar/SQRT(hPI)*EXP(-(xpar**2)))

dist1=SQRT(dist2)
lncoulomb=LOG(1+masses(uu)*dist1/(masses(u)*rhalo(uu)))

! change direction of friction (minus wegen meiner anderen Konvektion)
scalarp=vr(1)*r(1)+vr(2)*r(2)+vr(3)*r(3)

!beta=ACOS(scalarp/(vrel*dist1))
!beta=0.75D0*COS(beta-hPI/2.D0)   simplified by mathematica
beta=0.75D0*SQRT(1.D0-(scalarp/(vrel*dist1))**2)
if (beta<0) then
beta=-beta
end if

ffrict(1)=-(r(1)/dist1-vr(1)*scalarp/(dist1*vrel*vrel))
ffrict(2)=-(r(2)/dist1-vr(2)*scalarp/(dist1*vrel*vrel))
ffrict(3)=-(r(3)/dist1-vr(3)*scalarp/(dist1*vrel*vrel))

scalarp=SQRT(ffrict(1)**2+ffrict(2)**2+ffrict(3)**2)

ffrict(1)=ffrict(1)/scalarp
ffrict(2)=ffrict(2)/scalarp
ffrict(3)=ffrict(3)/scalarp

ffrict(1)=COS(beta)*vr(1)/vrel+ffrict(1)*SIN(beta)
ffrict(2)=COS(beta)*vr(2)/vrel+ffrict(2)*SIN(beta)
ffrict(3)=COS(beta)*vr(3)/vrel+ffrict(3)*SIN(beta)

ffrict(1)=-masses(uu)/(rhalo(uu)*dist2)*lncoulomb/(vrel**2)*masses(u)*xpar*ffrict(1)*scalefrict
ffrict(2)=-masses(uu)/(rhalo(uu)*dist2)*lncoulomb/(vrel**2)*masses(u)*xpar*ffrict(2)*scalefrict
ffrict(3)=-masses(uu)/(rhalo(uu)*dist2)*lncoulomb/(vrel**2)*masses(u)*xpar*ffrict(3)*scalefrict

dydx((u-1)*6+4)=dydx((u-1)*6+4)+ffrict(1)
dydx((u-1)*6+5)=dydx((u-1)*6+5)+ffrict(2)
dydx((u-1)*6+6)=dydx((u-1)*6+6)+ffrict(3)
! actio est reactio
dydx((uu-1)*6+4)=dydx((uu-1)*6+4)-ffrict(1)*masses(u)/masses(uu)
dydx((uu-1)*6+5)=dydx((uu-1)*6+5)-ffrict(2)*masses(u)/masses(uu)
dydx((uu-1)*6+6)=dydx((uu-1)*6+6)-ffrict(3)*masses(u)/masses(uu)


END IF

END IF
END IF
END DO
END DO
END IF

!Dynamical Friction 2 - less acuarte
if (dynfricton==2) THEN
DO u=1,nm
DO uu=1,nm
IF (u.NE.uu) THEN
IF (masses(uu)>masses(u)) THEN

r(1)=y((uu-1)*6+1)-y((u-1)*6+1)
r(2)=y((uu-1)*6+2)-y((u-1)*6+2)
r(3)=y((uu-1)*6+3)-y((u-1)*6+3)

dist2=((r(1)**2)+(r(2)**2)+(r(3)**2)+(soft**2))


IF ((dist2<(rhalo(uu)**2)).AND.(dist2>(2.D0*(soft**2)))) THEN

vr(1)=y((u-1)*6+4)-y((uu-1)*6+4)
vr(2)=y((u-1)*6+5)-y((uu-1)*6+5)
vr(3)=y((u-1)*6+6)-y((uu-1)*6+6)

vrel=SQRT((vr(1)**2)+(vr(2)**2)+(vr(3)**2)+1.D-4)

xpar=vrel/SQRT(masses(uu)/rhalo(uu))

xpar=(errorf(xpar)-2.D0*xpar/SQRT(hPI)*EXP(-(xpar**2)))


lncoulomb=LOG(1+masses(uu)*dist1/(masses(u)*rhalo(uu)))

ffrict(1)=-masses(uu)/(rhalo(uu)*dist2)*lncoulomb/(vrel**3)*masses(u)*xpar*vr(1)*scalefrict
ffrict(2)=-masses(uu)/(rhalo(uu)*dist2)*lncoulomb/(vrel**3)*masses(u)*xpar*vr(2)*scalefrict
ffrict(3)=-masses(uu)/(rhalo(uu)*dist2)*lncoulomb/(vrel**3)*masses(u)*xpar*vr(3)*scalefrict

dydx((u-1)*6+4)=dydx((u-1)*6+4)+ffrict(1)
dydx((u-1)*6+5)=dydx((u-1)*6+5)+ffrict(2)
dydx((u-1)*6+6)=dydx((u-1)*6+6)+ffrict(3)
! actio est reactio
dydx((uu-1)*6+4)=dydx((uu-1)*6+4)-ffrict(1)*masses(u)/masses(uu)
dydx((uu-1)*6+5)=dydx((uu-1)*6+5)-ffrict(2)*masses(u)/masses(uu)
dydx((uu-1)*6+6)=dydx((uu-1)*6+6)-ffrict(3)*masses(u)/masses(uu)


END IF

END IF
END IF
END DO
END DO
END IF



! Newtonian Gravitation
DO u=1,nm

r(1)=y((u-1)*6+1)-y(nm*6+1)
r(2)=y((u-1)*6+2)-y(nm*6+2)
r(3)=y((u-1)*6+3)-y(nm*6+3)

IF (haloon==1) THEN

dist2=((r(1)**2)+(r(2)**2)+(r(3)**2)+(soft**2))

IF (dist2<(rhalo(u)**2)) THEN

dydx(nm*6+4)=dydx(nm*6+4)+(r(1)*masses(u)/(dist2*rhalo(u)))
dydx(nm*6+5)=dydx(nm*6+5)+(r(2)*masses(u)/(dist2*rhalo(u)))
dydx(nm*6+6)=dydx(nm*6+6)+(r(3)*masses(u)/(dist2*rhalo(u)))

ELSE
  
dist1=SQRT(dist2)
dist3=dist1**3

dydx(nm*6+4)=dydx(nm*6+4)+(r(1)*masses(u)/dist3)
dydx(nm*6+5)=dydx(nm*6+5)+(r(2)*masses(u)/dist3)
dydx(nm*6+6)=dydx(nm*6+6)+(r(3)*masses(u)/dist3)

END IF



ELSE
dist1=SQRT(((r(1)**2)+(r(2)**2)+(r(3)**2)+(soft**2)))
dist3=dist1**3

dydx(nm*6+4)=dydx(nm*6+4)+(r(1)*masses(u)/dist3)
dydx(nm*6+5)=dydx(nm*6+5)+(r(2)*masses(u)/dist3)
dydx(nm*6+6)=dydx(nm*6+6)+(r(3)*masses(u)/dist3)
END IF

END DO

!Hubble Expansions
IF (hubbleon==1) THEN
dydx(nm*6+4)=dydx(nm*6+4)+y(nm*6+1)*2.D0/(3.D0*(x**2))/goldenratio*hadjust
dydx(nm*6+5)=dydx(nm*6+5)+y(nm*6+2)*2.D0/(3.D0*(x**2))/goldenratio*hadjust
dydx(nm*6+6)=dydx(nm*6+6)+y(nm*6+3)*2.D0/(3.D0*(x**2))/goldenratio*hadjust
END IF


!Dynamical Friction

if (dynfricton==1) THEN
DO u=1,nm

r(1)=y((u-1)*6+1)-y(nm*6+1)
r(2)=y((u-1)*6+2)-y(nm*6+2)
r(3)=y((u-1)*6+3)-y(nm*6+3)

dist2=((r(1)**2)+(r(2)**2)+(r(3)**2)+(soft**2))

IF ((dist2<(rhalo(u)**2)).AND.(dist2>(2.D0*(soft**2)))) THEN

vr(1)=y(nm*6+4)-y((u-1)*6+4)
vr(2)=y(nm*6+5)-y((u-1)*6+5)
vr(3)=y(nm*6+6)-y((u-1)*6+6)

vrel=SQRT((vr(1)**2)+(vr(2)**2)+(vr(3)**2)+1.D-4)

xpar=vrel/SQRT(masses(u)/rhalo(u))

xpar=(errorf(xpar)-2.D0*xpar/SQRT(hPI)*EXP(-(xpar**2)))

dist1=SQRT(dist2)

lncoulomb=LOG(1+masses(u)*dist1/(massless*rhalo(u)))


! change direction of friction (minus wegen meiner anderen Konvektion)
scalarp=vr(1)*r(1)+vr(2)*r(2)+vr(3)*r(3)

!beta=ACOS(scalarp/(vrel*dist1))
!beta=0.75D0*COS(beta-hPI/2.D0)   simplified by mathematica
beta=0.75D0*SQRT(1.D0-(scalarp/(vrel*dist1))**2)
if (beta<0) then
beta=-beta
end if

ffrict(1)=-(r(1)/dist1-vr(1)*scalarp/(dist1*vrel*vrel))
ffrict(2)=-(r(2)/dist1-vr(2)*scalarp/(dist1*vrel*vrel))
ffrict(3)=-(r(3)/dist1-vr(3)*scalarp/(dist1*vrel*vrel))

scalarp=SQRT(ffrict(1)**2+ffrict(2)**2+ffrict(3)**2)

ffrict(1)=ffrict(1)/scalarp
ffrict(2)=ffrict(2)/scalarp
ffrict(3)=ffrict(3)/scalarp

ffrict(1)=COS(beta)*vr(1)/vrel+ffrict(1)*SIN(beta)
ffrict(2)=COS(beta)*vr(2)/vrel+ffrict(2)*SIN(beta)
ffrict(3)=COS(beta)*vr(3)/vrel+ffrict(3)*SIN(beta)

ffrict(1)=-masses(u)/(rhalo(u)*dist2)*lncoulomb/(vrel**2)*massless*xpar*ffrict(1)*scalefrict
ffrict(2)=-masses(u)/(rhalo(u)*dist2)*lncoulomb/(vrel**2)*massless*xpar*ffrict(2)*scalefrict
ffrict(3)=-masses(u)/(rhalo(u)*dist2)*lncoulomb/(vrel**2)*massless*xpar*ffrict(3)*scalefrict

dydx(nm*6+4)=dydx(nm*6+4)+ffrict(1)
dydx(nm*6+5)=dydx(nm*6+5)+ffrict(2)
dydx(nm*6+6)=dydx(nm*6+6)+ffrict(3)

END IF

END DO
END IF


!Dynamical Friction 2 - less acuarte

if (dynfricton==1) THEN
DO u=1,nm

r(1)=y((u-1)*6+1)-y(nm*6+1)
r(2)=y((u-1)*6+2)-y(nm*6+2)
r(3)=y((u-1)*6+3)-y(nm*6+3)

dist2=((r(1)**2)+(r(2)**2)+(r(3)**2)+(soft**2))

IF ((dist2<(rhalo(u)**2)).AND.(dist2>(2.D0*(soft**2)))) THEN

vr(1)=y(nm*6+4)-y((u-1)*6+4)
vr(2)=y(nm*6+5)-y((u-1)*6+5)
vr(3)=y(nm*6+6)-y((u-1)*6+6)

vrel=SQRT((vr(1)**2)+(vr(2)**2)+(vr(3)**2)+1.D-4)

xpar=vrel/SQRT(masses(u)/rhalo(u))

xpar=(errorf(xpar)-2.D0*xpar/SQRT(hPI)*EXP(-(xpar**2)))

lncoulomb=LOG(1+masses(u)*dist1/(massless*rhalo(u)))


ffrict(1)=-masses(u)/(rhalo(u)*dist2)*lncoulomb/(vrel**3)*massless*xpar*vr(1)*scalefrict
ffrict(2)=-masses(u)/(rhalo(u)*dist2)*lncoulomb/(vrel**3)*massless*xpar*vr(2)*scalefrict
ffrict(3)=-masses(u)/(rhalo(u)*dist2)*lncoulomb/(vrel**3)*massless*xpar*vr(3)*scalefrict

dydx(nm*6+4)=dydx(nm*6+4)+ffrict(1)
dydx(nm*6+5)=dydx(nm*6+5)+ffrict(2)
dydx(nm*6+6)=dydx(nm*6+6)+ffrict(3)

END IF

END DO
END IF



END SUBROUTINE derivs

END MODULE


PROGRAM newtonian_hubble_expansion_integrator 
! NewHExI
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

!Defination of PI
hPI=ACOS(-1.D0)

!Defination of golden ratio
goldenratio=(1.D0+SQRT(5.D0))/2.D0

! Read parameters for integrator and general settings
OPEN(52,file='para.dat')
READ(52,*) filein
READ(52,*) fileout
READ(52,*) tstop
READ(52,*) tstep
READ(52,*) soft
READ(52,*) massless
READ(52,*) toffset
READ(52,*) haloon
READ(52,*) dynfricton
READ(52,*) hubbleon
READ(52,*) scalefrict
READ(52,*) hadjust
CLOSE(52)




 !unitmass = 10^11 Msol
 !unitspeed = 378.66876 km/s
 !unitlength = 3kpc
 !unittime = 7.74655044869528 Myr




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
ALLOCATE(ystart(1:((nm+1)*6)))

DO i=1,nm
masses(i)=m(i)
END DO

ntb=n-nm



!get halo extension
ALLOCATE(rhalo(1:nm))
  
IF ((dynfricton.NE.0).OR.(haloon==1)) THEN
  
OPEN(53,file='halo.dat')

DO i=1,nm
READ(53,*) rhalo(i)
END DO

CLOSE(53)
  
END IF



WRITE(*,*) 
WRITE(*,*) 'NewHExI'
WRITE(*,*) 'is a NEWtonian Hubble EXpansion Integrator'
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

DO i=1,nm
pos(1,i)=pos(1,i)-compos(1)
pos(2,i)=pos(2,i)-compos(2)
pos(3,i)=pos(3,i)-compos(3)
vel(1,i)=vel(1,i)-comvel(1)
vel(2,i)=vel(2,i)-comvel(2)
vel(3,i)=vel(3,i)-comvel(3)
END DO

DO i=1,ntb
tbpos(1,i)=tbpos(1,i)-compos(1)
tbpos(2,i)=tbpos(2,i)-compos(2)
tbpos(3,i)=tbpos(3,i)-compos(3)
tbvel(1,i)=tbvel(1,i)-comvel(1)
tbvel(2,i)=tbvel(2,i)-comvel(2)
tbvel(3,i)=tbvel(3,i)-comvel(3)
END DO



! output initial coordinates
OPEN(52,file=fileout)

WRITE(52,*) n
WRITE(52,*) 3
WRITE(52,*) 0

DO i=1,n
WRITE (52,*) m(i)
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
 nnn=i+nm

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

DO i=1,nm
pos(1,i)=pos(1,i)-compos(1)
pos(2,i)=pos(2,i)-compos(2)
pos(3,i)=pos(3,i)-compos(3)
vel(1,i)=vel(1,i)-comvel(1)
vel(2,i)=vel(2,i)-comvel(2)
vel(3,i)=vel(3,i)-comvel(3)
END DO

DO i=1,ntb
tbpos(1,i)=tbpos(1,i)-compos(1)
tbpos(2,i)=tbpos(2,i)-compos(2)
tbpos(3,i)=tbpos(3,i)-compos(3)
tbvel(1,i)=tbvel(1,i)-comvel(1)
tbvel(2,i)=tbvel(2,i)-comvel(2)
tbvel(3,i)=tbvel(3,i)-comvel(3)
END DO

! output results
WRITE(52,*) n
WRITE(52,*) 3
WRITE(52,*) (t1-toffset)

DO i=1,n
WRITE (52,*) m(i)
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