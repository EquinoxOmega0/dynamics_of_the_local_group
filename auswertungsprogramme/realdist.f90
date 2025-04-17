PROGRAM Real_distribution
IMPLICIT NONE
! Counter variables
integer :: n,i,nm,modus
! Spherical (galactic) coordinates variables
double precision, allocatable :: r(:),l(:),b(:),vr(:),vt(:),alpha(:),m(:)
! carthesian coordinates
double precision, allocatable :: x(:),y(:),z(:),vx(:),vy(:),vz(:)
!help variables
double precision :: vtx,vty,vtz,rho,norm,PI,w1,w2,w3,mm
double precision :: xcross,ycross,zcross,rcross,xh,yh,zh
!Defination of PI
PI=ACOS(-1.D0)

!Read Input (spherical data)
OPEN(10,file='vglinput.dat')

!Read modus: 1 for normal vgl-file output; 0: for nemo-file output
! 2: only koordinates, keep radial velocity
READ(10,*) modus

!get number of galaxies
READ(10,*) n
READ(10,*) nm


! allocate arrays
ALLOCATE(r(1:(n+1)))
ALLOCATE(l(1:(n+1)))
ALLOCATE(b(1:(n+1)))
ALLOCATE(vr(1:(n+1)))
ALLOCATE(vt(1:(n+1)))
ALLOCATE(alpha(1:(n+1)))

ALLOCATE(m(1:nm))

ALLOCATE(x(1:(n+1)))
ALLOCATE(y(1:(n+1)))
ALLOCATE(z(1:(n+1)))
ALLOCATE(vx(1:(n+1)))
ALLOCATE(vy(1:(n+1)))
ALLOCATE(vz(1:(n+1)))


! Read data from file
DO i=1,nm
READ(10,*) m(i)
END DO


DO i=1,n 

READ(10,*) l(i),b(i),r(i),vr(i),vt(i),alpha(i)

! transfrom ° to rad
b(i)=b(i)*PI/180
l(i)=l(i)*PI/180
alpha(i)=alpha(i)*PI/180

END DO

CLOSE(10)

!save origin (=Earth)
x(n+1)=0.D0
y(n+1)=0.D0
z(n+1)=0.D0
vx(n+1)=0.D0
vy(n+1)=220.D0  !eventuell vorzeichen ändern, nochmal durchdenken!!!
vz(n+1)=0.D0



DO i=1,n

!Calculates carthesian coordiantes
x(i)=r(i)*COS(b(i))*COS(l(i))
y(i)=r(i)*COS(b(i))*SIN(l(i))
z(i)=r(i)*SIN(b(i))


! Calculate radial velocity part in carthesian coordiantes
vx(i)=(x(i)/r(i))*vr(i)
vy(i)=(y(i)/r(i))*vr(i)
vz(i)=(z(i)/r(i))*vr(i)

! check if there is also a tangential velocity
IF (vt(i).NE.0.D0) THEN

! transform tangential velocity into carthesian coordiantes
rho=SQRT(x(i)**2+y(i)**2)
norm=SQRT(2.D0*(y(i)**2)*(z(i)**2)+rho**2)

vty=(SIN(alpha(i))*norm+z(i)*rho*COS(alpha(i))-x(i)/y(i)*COS(alpha(i))*(rho**3))
vty=vty/(-z(i)*x(i)+y(i)*z(i)+(x(i)**2)*(rho**2)/(y(i)*z(i))+y(i)*(rho**2)/z(i))

vtx=(COS(alpha(i))*rho+vty*x(i))/y(i)

vtz=-x(i)*COS(alpha(i))*rho/(z(i)*y(i))-(x(i)**2)*vty/(y(i)*z(i))-y(i)*vty/z(i)

! add tangential velocity to total velocity
vx(i)=vx(i)+vt(i)*vtx
vy(i)=vy(i)+vt(i)*vty
vz(i)=vz(i)+vt(i)*vtz

ENDIF

END DO





! Shift system to Milkyway-centered

DO i=2,(n+1)
x(i)=x(i)-x(1)
y(i)=y(i)-y(1)
z(i)=z(i)-z(1)
END DO
x(1)=0.D0
y(1)=0.D0
z(1)=0.D0



! rotated to MW-M31-M33 plane
xcross=y(2)*z(3)-y(3)*z(2)
ycross=z(2)*x(3)-z(3)*x(2)
zcross=x(2)*y(3)-x(3)*y(2)

rcross=xcross**2+ycross**2+zcross**2
rcross=SQRT(rcross)

xcross=xcross/rcross
ycross=ycross/rcross
zcross=zcross/rcross

w1=-atan2(ycross,xcross)
w2=-acos(zcross)



DO i=1,(n+1)

xh=x(i)*cos(w1)-y(i)*sin(w1)
yh=x(i)*sin(w1)+y(i)*cos(w1)
zh=z(i)

x(i)=xh
y(i)=yh
z(i)=zh

xh=vx(i)*cos(w1)-vy(i)*sin(w1)
yh=vx(i)*sin(w1)+vy(i)*cos(w1)
zh=vz(i)

vx(i)=xh
vy(i)=yh
vz(i)=zh

END DO



DO i=1,(n+1)

xh=cos(w2)*x(i)+sin(w2)*z(i)
yh=y(i)
zh=-sin(w2)*x(i)+cos(w2)*z(i)

x(i)=xh
y(i)=yh
z(i)=zh

xh=cos(w2)*vx(i)+sin(w2)*vz(i)
yh=vy(i)
zh=-sin(w2)*vx(i)+cos(w2)*vz(i)

vx(i)=xh
vy(i)=yh
vz(i)=zh

END DO


w3=-atan2(y(2),x(2))


! Rotate M31 on x-axis
DO i=1,(n+1)

xh=x(i)*cos(w3)-y(i)*sin(w3)
yh=x(i)*sin(w3)+y(i)*cos(w3)
zh=z(i)

x(i)=xh
y(i)=yh
z(i)=zh

xh=vx(i)*cos(w3)-vy(i)*sin(w3)
yh=vx(i)*sin(w3)+vy(i)*cos(w3)
zh=vz(i)

vx(i)=xh
vy(i)=yh
vz(i)=zh

END DO


                            !unitmass = 10^11 Msol
rho=3.D0                    !unitlength = 3kpc
norm=378.668761569285D0     !unitspeed = 378.66876 km/s
                            !unittime = 7.74655044869528 Myr


! Rescale System
DO i=1,(n+1)
x(i)=x(i)/rho
y(i)=y(i)/rho
z(i)=z(i)/rho

vx(i)=vx(i)/norm
vy(i)=vy(i)/norm
vz(i)=vz(i)/norm

END DO

! save radial velocity instead of other velocity if modus is 2

IF (modus==2) THEN


DO i=1,n
vx(i)=vr(i)/norm
vy(i)=vt(i)/norm
vz(i)=alpha(i)
END DO




END IF




! normal vgl-file output
IF (modus.ne.0) THEN


! Output carthesian data
OPEN(11,file='vgl.dat')


! number of massiv galaxies
WRITE(11,*) nm
! number of massless galaxies
WRITE(11,*) (n-nm)

DO i=1,(n+1)

!write carthesian positions
WRITE(11,*) x(i),y(i),z(i)
!write carthesian velocities
WRITE(11,*) vx(i),vy(i),vz(i)

END DO

CLOSE(11)


ELSE


OPEN(11,file='vglnemo.dat')

! Write header
WRITE(11,*) n
WRITE(11,*) 3
WRITE(11,*) 0.D0

! write massiv particles
DO i=1,nm 
WRITE(11,*) m(i)
END DO

! output mass for massless particles
mm=1.D-8

DO i=nm+1,n
WRITE(11,*) mm
END DO

! Write particle positions
DO i=1,n
WRITE(11,*) x(i),y(i),z(i)
END DO

! Write particle velocites
DO i=1,n
WRITE(11,*) vx(i),vy(i),vz(i)
END DO

CLOSE(11)

END IF


END PROGRAM 