PROGRAM makemodel
IMPLICIT NONE
! declartion of variables
! masses for masssive particles
Double Precision,allocatable  :: masses(:)
!all phase space coordinates for massive particles
Double Precision,allocatable :: x(:),y(:),z(:),vx(:),vy(:),vz(:)
!all phase space coordiantes for test particles
Double Precision,allocatable :: tbx(:),tby(:),tbz(:),tbvx(:),tbvy(:),tbvz(:)
! origin of testparticles
Double Precision,allocatable :: orix(:),oriy(:),oriz(:),orivx(:),orivy(:),orivz(:),oriradius(:),orivdisp(:)
! variables for initial condiations, centre of mass and random generator
!Double Precision :: E,L,theta,startdist
Double Precision :: startdist,vabs,vangle,tol
!Double Precision :: mu,km,startr,startx,starty,startvx,startvy,drdt,dthetadt,scalef
Double Precision :: xcm,ycm,zcm,vxcm,vycm,vzcm,summass,rr,rtheta,rphi,rv,rvphi,rvtheta,PI,orimasse
! particles counters
Integer :: nm,ntb,nori,nset,i,ii,ntot,iii,schleife
!filename for output
character(len=16) :: filename

!Defination of PI
PI=ACOS(-1.D0)



OPEN(50,file='input.dat')
!Read in needed lengthes
READ(50,*) nm
READ(50,*) nori
READ(50,*) nset

ntb=nori*nset

! allocate all arrays
allocate(masses(1:nm))
allocate(x(1:nm))
allocate(y(1:nm))
allocate(z(1:nm))
allocate(vx(1:nm))
allocate(vy(1:nm))
allocate(vz(1:nm))
allocate(tbx(1:ntb))
allocate(tby(1:ntb))
allocate(tbz(1:ntb))
allocate(tbvx(1:ntb))
allocate(tbvy(1:ntb))
allocate(tbvz(1:ntb))
allocate(orix(1:nori))
allocate(oriy(1:nori))
allocate(oriz(1:nori))
allocate(orivx(1:nori))
allocate(orivy(1:nori))
allocate(orivz(1:nori))
allocate(oriradius(1:nori))
allocate(orivdisp(1:nori))

! filename for output
READ(50,*) filename

! read in masses
DO i=1,nm 
READ(50,*) masses(i)
END DO

!intial condiation for 2 dominat galaxies
READ(50,*) startdist
READ(50,*) vabs
READ(50,*) vangle
vangle=vangle*PI/180.D0

!read in intial condiations of other massive particles
IF (nm>2) THEN
DO i=3,nm
READ(50,*) x(i),y(i),z(i)
READ(50,*) vx(i),vy(i),vz(i)
END DO
END IF


!read in intial conditions for test particles
READ(50,*) orimasse

DO i=1,nori
READ(50,*) orix(i),oriy(i),oriz(i)
READ(50,*) orivx(i),orivy(i),orivz(i)
READ(50,*) oriradius(i)
READ(50,*) orivdisp(i)
END DO


! Positions of main galaxies in system of galaxy 1
x(1)=0.D0
y(1)=0.D0
z(1)=0.D0
vx(1)=0.D0
vy(1)=0.D0
vz(1)=0.D0

x(2)=startdist
y(2)=0.D0
z(2)=0.D0
vx(2)=vabs*COS(vangle)
vy(2)=vabs*SIN(vangle)
vz(2)=0.D0


! move to center of mass frame of galaxy 1 and 2

summass=masses(1)+masses(2)

xcm=x(1)*masses(1)+x(2)*masses(2)
ycm=y(1)*masses(1)+y(2)*masses(2)
zcm=z(1)*masses(1)+z(2)*masses(2)

xcm=xcm/summass
ycm=ycm/summass
zcm=zcm/summass

vxcm=vx(1)*masses(1)+vx(2)*masses(2)
vycm=vy(1)*masses(1)+vy(2)*masses(2)
vzcm=vz(1)*masses(1)+vz(2)*masses(2)

vxcm=vxcm/summass
vycm=vycm/summass
vzcm=vzcm/summass

DO i=1,2

x(i)=x(i)-xcm
y(i)=y(i)-ycm
z(i)=z(i)-zcm

vx(i)=vx(i)-vxcm
vy(i)=vy(i)-vycm
vz(i)=vz(i)-vzcm

END DO



! shift to centre of mass of whole system

xcm=0.D0
ycm=0.D0
zcm=0.D0
vxcm=0.D0
vycm=0.D0
vzcm=0.D0
summass=0.D0


DO i=1,nm 

xcm=xcm+masses(i)*x(i)
ycm=ycm+masses(i)*y(i)
zcm=zcm+masses(i)*z(i)

vxcm=vxcm+masses(i)*vx(i)
vycm=vycm+masses(i)*vy(i)
vycm=vycm+masses(i)*vz(i)

summass=summass+masses(i)

END DO

xcm=xcm/summass
ycm=ycm/summass
zcm=zcm/summass


vxcm=vxcm/summass
vycm=vycm/summass
vzcm=vzcm/summass


DO i=1,nm
x(i)=x(i)-xcm
y(i)=y(i)-ycm
z(i)=z(i)-zcm
vx(i)=vx(i)-vxcm
vy(i)=vy(i)-vycm
vz(i)=vz(i)-vzcm
END DO



tol=(oriradius(i))/((FLOAT(nset))**(1/3))
tol=(tol/100.D0)**2

!create randomized test particles

call random_seed()

DO i=1,nori

!one testparticle in origin of tporigin
tbx((i-1)*nset+1)=orix(i)
tby((i-1)*nset+1)=oriy(i)
tbz((i-1)*nset+1)=oriz(i)
tbvx((i-1)*nset+1)=orivx(i)
tbvy((i-1)*nset+1)=orivy(i)
tbvz((i-1)*nset+1)=orivz(i)
  

If (nset>1) then
DO ii=2,nset

schleife=0

103 CONTINUE 

schleife=schleife+1

! all others randomized around origin
call random_number(rr)
call random_number(rphi)
call random_number(rtheta)
call random_number(rv)
call random_number(rvphi)
call random_number(rvtheta)

rr=rr*oriradius(i)
rphi=rphi*2.D0*PI
rtheta=rtheta*PI
rv=rv*orivdisp(i)
rvphi=rvphi*2.D0*PI
rvtheta=rvtheta*PI

tbx((i-1)*nset+ii)=orix(i)+rr*SIN(rvtheta)*COS(rphi)
tby((i-1)*nset+ii)=oriy(i)+rr*SIN(rvtheta)*SIN(rphi)
tbz((i-1)*nset+ii)=oriz(i)+rr*COS(rvtheta)

DO iii=1,((i-1)*nset+ii)
rr=(tbx((i-1)*nset+ii)-tbx(iii))**2
rr=rr+(tby((i-1)*nset+ii)-tby(iii))**2
rr=rr+(tbz((i-1)*nset+ii)-tbz(iii))**2
IF (rr<tol) THEN
IF ((iii.NE.((i-1)*nset+ii)).AND.(schleife<1000)) THEN
GOTO 103
END IF
END IF
END DO

tbvx((i-1)*nset+ii)=orivx(i)+rv*SIN(rvtheta)*COS(rvphi)
tbvy((i-1)*nset+ii)=orivy(i)+rv*SIN(rvtheta)*SIN(rvphi)
tbvz((i-1)*nset+ii)=orivz(i)+rv*COS(rvtheta)

END DO
end if
END DO




! output model

OPEN(60,file=filename)

! header of file
ntot=nm+ntb
WRITE(60,*) ntot
WRITE(60,*) 3
WRITE(60,*) 0

! output masses
DO i=1,nm
WRITE(60,*) masses(i)
END DO


! output coordinates
DO i=1,ntb
WRITE(60,*) orimasse
END DO

DO i=1,nm
WRITE(60,*) x(i),y(i),z(i)
END DO

DO i=1,ntb
WRITE(60,*) tbx(i),tby(i),tbz(i)
END DO


!output velocities
DO i=1,nm
WRITE(60,*) vx(i),vy(i),vz(i)
END DO

DO i=1,ntb
WRITE(60,*) tbvx(i),tbvy(i),tbvz(i)
END DO

!finish
CLOSE(50)
CLOSE(60)

WRITE(*,*) 'Modell erfolgreich erstellt'

END 
