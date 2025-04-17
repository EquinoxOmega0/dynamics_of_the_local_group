
!main program
PROGRAM auswertung
IMPLICIT NONE

! declaration of variables
character(len=16) :: filename,vglfile,newpar,oldpar,fname,zahl,protocol
double precision, allocatable :: vglx(:,:),vglv(:,:),vgltbx(:,:),vgltbv(:,:),pos(:,:),vel(:,:),p0(:,:),v0(:,:)
integer :: i,ii,vglnm,vglntb,nm,nori,nset,n,counter,nmodel,nnewparents,nmut
integer :: onlyrad
Double Precision,allocatable :: orix(:,:),oriy(:,:),oriz(:,:),orivx(:,:),orivy(:,:),orivz(:,:),oriradius(:,:)
Double Precision,allocatable :: orivdisp(:,:)
Double Precision,allocatable :: x(:,:),y(:,:),z(:,:),vx(:,:),vy(:,:),vz(:,:),masses(:,:)
Double Precision :: distm,vdistm,disttotb,vdisttotb
integer,allocatable:: bestindex(:),bestmodel(:)
Double Precision,allocatable :: orimasse(:),gmodel(:)
Double Precision,allocatable :: startdist(:),vabsolut(:),vangle(:)
Double Precision :: gmpos,gmvel,gtbpos,gtbvel,massivg,tbgtot
Double Precision :: xobs,yobs,zobs,vxobs,vyobs,vzobs,vradmod,vradvgl!,PI
Double Precision, dimension(:) :: va(1:3),vb(1:3)
double precision :: xcross,ycross,zcross,rcross,xh,yh,zh,w1,w2,w3
double precision,allocatable :: gewichte(:,:)

! open file for settings of genetic algoritm
OPEN(51,file='in.dat')

READ(51,*) vglfile
READ(51,*) filename
READ(51,*) oldpar
READ(51,*) newpar
READ(51,*) protocol
READ(51,*) nmodel
READ(51,*) gmpos
READ(51,*) gmvel
READ(51,*) gtbpos
READ(51,*) gtbvel
READ(51,*) nnewparents
READ(51,*) nmut
READ(51,*) onlyrad

CLOSE(51)

OPEN(60,file='ergebnis.dat')


! Read in data from compare-file
OPEN(52,file=vglfile)

READ(52,*) vglnm
READ(52,*) vglntb

ALLOCATE(vglx(1:3,1:vglnm))
ALLOCATE(vglv(1:3,1:vglnm))

ALLOCATE(vgltbx(1:3,1:vglntb))
ALLOCATE(vgltbv(1:3,1:vglntb))

DO i=1,vglnm
READ(52,*) vglx(1:3,i)
READ(52,*) vglv(1:3,i)
END DO

DO i=1,vglntb
READ(52,*) vgltbx(1:3,i)
READ(52,*) vgltbv(1:3,i)
END DO

READ(52,*) xobs,yobs,zobs
READ(52,*) vxobs,vyobs,vzobs

CLOSE(52)



! Read in array-sizes form one old parameter file
OPEN(53,file=trim(oldpar)//'1.dat')

READ(53,*) nm
READ(53,*) nori
READ(53,*) nset

CLOSE(53)



! set size of arrays
ALLOCATE(masses(1:nm,1:nmodel))
ALLOCATE(x(1:nm,1:nmodel))
ALLOCATE(y(1:nm,1:nmodel))
ALLOCATE(z(1:nm,1:nmodel))
ALLOCATE(vx(1:nm,1:nmodel))
ALLOCATE(vy(1:nm,1:nmodel))
ALLOCATE(vz(1:nm,1:nmodel))

allocate(bestmodel(1:nnewparents))

ALLOCATE(startdist(1:nmodel))
ALLOCATE(vabsolut(1:nmodel))
ALLOCATE(vangle(1:nmodel))


ALLOCATE(orimasse(1:nmodel))

if (nori<vglntb) then
nori=vglntb
end if

ALLOCATE(orix(1:nori,1:nmodel))
ALLOCATE(oriy(1:nori,1:nmodel))
ALLOCATE(oriz(1:nori,1:nmodel))
ALLOCATE(orivx(1:nori,1:nmodel))
ALLOCATE(orivy(1:nori,1:nmodel))
ALLOCATE(orivz(1:nori,1:nmodel))
ALLOCATE(oriradius(1:nori,1:nmodel))
ALLOCATE(orivdisp(1:nori,1:nmodel))

ALLOCATE(gmodel(1:nmodel))

n=nm+nori*nset

ALLOCATE(pos(1:3,1:n))
ALLOCATE(vel(1:3,1:n))

ALLOCATE(p0(1:3,1:n))
ALLOCATE(v0(1:3,1:n))

ALLOCATE(bestindex(1:vglntb))

ALLOCATE(gewichte(0:(nm+nori),1:nmodel))



! Loop for all different models
DO counter=1,nmodel

! Open old parameter file for model
write(zahl,*) counter
zahl=adjustl(zahl)
OPEN(53,file=trim(oldpar)//trim(zahl)//'.dat')
  Write(*,*) 'OK'//trim(zahl)


! Read in header
READ(53,*) nm
READ(53,*) nori
READ(53,*) nset
READ(53,*) fname

! Read in initial conditions of for old model
DO i=1,nm 
READ(53,*) masses(i,counter)
END DO


READ(53,*) startdist(counter)
READ(53,*) vabsolut(counter)
READ(53,*) vangle(counter)


IF (nm>2) THEN
DO i=3,nm
READ(53,*) x(i,counter),y(i,counter),z(i,counter)
READ(53,*) vx(i,counter),vy(i,counter),vz(i,counter)
END DO
END IF

READ(53,*) orimasse(counter)

DO i=1,nori
READ(53,*) orix(i,counter),oriy(i,counter),oriz(i,counter)
READ(53,*) orivx(i,counter),orivy(i,counter),orivz(i,counter)
READ(53,*) oriradius(i,counter)
READ(53,*) orivdisp(i,counter)
END DO

CLOSE(53)



! Read in particle data for time 0
OPEN(50,file=trim(filename)//trim(zahl)//'-0.dat')

READ(50,*) n
READ(50,*) 
READ(50,*) 

DO i=1,n 
READ(50,*) 
END DO

DO i=1,n 
READ(50,*) p0(1:3,i)
END DO

DO i=1,n 
READ(50,*) v0(1:3,i)
END DO

CLOSE(50)

 


! Read in particle data
OPEN(50,file=trim(filename)//trim(zahl)//'-F.dat')

READ(50,*) n
READ(50,*) 
READ(50,*) 

DO i=1,n 
READ(50,*) 
END DO

DO i=1,n 
READ(50,*) pos(1:3,i)
END DO

DO i=1,n 
READ(50,*) vel(1:3,i)
END DO

CLOSE(50)

! Shift to Milky Way 

DO i=2,n

pos(1,i)=pos(1,i)-pos(1,1)
pos(2,i)=pos(2,i)-pos(2,1)
pos(3,i)=pos(3,i)-pos(3,1)
vel(1,i)=vel(1,i)-vel(1,1)
vel(2,i)=vel(2,i)-vel(2,1)
vel(3,i)=vel(3,i)-vel(3,1)

END DO

pos(1,1)=0.D0
pos(2,1)=0.D0
pos(3,1)=0.D0
vel(1,1)=0.D0
vel(2,1)=0.D0
vel(3,1)=0.D0


! rotated to MW-M31-M33 plane
xcross=pos(2,2)*pos(3,3)-pos(3,2)*pos(2,3)
ycross=pos(3,2)*pos(1,3)-pos(1,2)*pos(3,3)
zcross=pos(1,2)*pos(2,3)-pos(2,2)*pos(1,3)

rcross=xcross**2+ycross**2+zcross**2
rcross=SQRT(rcross)

xcross=xcross/rcross
ycross=ycross/rcross
zcross=zcross/rcross

w1=-atan2(ycross,xcross)
w2=-acos(zcross)


DO i=1,n

xh=pos(1,i)*cos(w1)-pos(2,i)*sin(w1)
yh=pos(1,i)*sin(w1)+pos(2,i)*cos(w1)
zh=pos(3,i)

pos(1,i)=xh
pos(2,i)=yh
pos(3,i)=zh

xh=vel(1,i)*cos(w1)-vel(2,i)*sin(w1)
yh=vel(1,i)*sin(w1)+vel(2,i)*cos(w1)
zh=vel(3,i)

vel(1,i)=xh
vel(2,i)=yh
vel(3,i)=zh

END DO



DO i=1,n

xh=cos(w2)*pos(1,i)+sin(w2)*pos(3,i)
yh=pos(2,i)
zh=-sin(w2)*pos(1,i)+cos(w2)*pos(3,i)

pos(1,i)=xh
pos(2,i)=yh
pos(3,i)=zh

xh=cos(w2)*vel(1,i)+sin(w2)*vel(3,i)
yh=vel(2,i)
zh=-sin(w2)*vel(1,i)+cos(w2)*vel(3,i)

vel(1,i)=xh
vel(2,i)=yh
vel(3,i)=zh

END DO


w3=-atan2(pos(2,2),pos(1,2))


! Rotate M31 on x-axis
DO i=1,n

xh=pos(1,i)*cos(w3)-pos(2,i)*sin(w3)
yh=pos(1,i)*sin(w3)+pos(2,i)*cos(w3)
zh=pos(3,i)

pos(1,i)=xh
pos(2,i)=yh
pos(3,i)=-zh

xh=vel(1,i)*cos(w3)-vel(2,i)*sin(w3)
yh=vel(1,i)*sin(w3)+vel(2,i)*cos(w3)
zh=vel(3,i)

vel(1,i)=xh
vel(2,i)=yh
vel(3,i)=-zh

END DO





! get fitness value for massiv particles

DO i=1,vglnm

! fitness for distance of massiv particles
distm=(vglx(1,i)-pos(1,i))**2+(vglx(2,i)-pos(2,i))**2+(vglx(3,i)-pos(3,i))**2

! fitness for velocity of massiv particles
IF (onlyrad==0) THEN
  !full velocity
vdistm=(vglv(1,i)-vel(1,i))**2+(vglv(2,i)-vel(2,i))**2+(vglv(3,i)-vel(3,i))**2

ELSE
!only radial velocity to observer
va(1)=vel(1,i)-vxobs
va(2)=vel(2,i)-vyobs
va(3)=vel(3,i)-vzobs

vb(1)=xobs-pos(1,i)
vb(2)=yobs-pos(2,i)
vb(3)=zobs-pos(3,i)

vradmod=SQRT(vb(1)**2+vb(2)**2+vb(3)**2)

vradmod=(va(1)*vb(1)+va(2)*vb(2)+va(3)*vb(3))/vradmod

vradvgl=vglv(1,i)

vdistm=(vradmod-vradvgl)**2

ENDIF

distm=SQRT(distm)
vdistm=SQRT(vdistm)

WRITE(60,*) distm
WRITE(60,*) vdistm
  
END DO


! get fitness value for massless particles


DO ii=1,vglntb

  
! fitness for distance of massiv particles
disttotb=(vgltbx(1,ii)-pos(1,ii+vglnm))**2+(vgltbx(2,ii)-pos(2,ii+vglnm))**2+(vgltbx(3,ii)-pos(3,ii+vglnm))**2

 Write(*,*) 'OK'
! fitness for velocity of massiv particles
IF (onlyrad==0) then
    !full velocity
vdisttotb=(vgltbv(1,ii)-vel(1,ii+vglnm))**2+(vgltbv(2,ii)-vel(2,ii+vglnm))**2+(vgltbv(3,ii)-vel(3,ii+vglnm))**2

ELSE
!only radial velocity to observer

va(1)=vel(1,ii+vglnm)-vxobs
va(2)=vel(2,ii+vglnm)-vyobs
va(3)=vel(3,ii+vglnm)-vzobs

vb(1)=xobs-pos(1,ii+vglnm)
vb(2)=yobs-pos(2,ii+vglnm)
vb(3)=zobs-pos(3,ii+vglnm)

vradmod=SQRT(vb(1)**2+vb(2)**2+vb(3)**2)

vradmod=(va(1)*vb(1)+va(2)*vb(2)+va(3)*vb(3))/vradmod

vradvgl=vgltbv(1,ii)

vdisttotb=(vradmod-vradvgl)**2

ENDIF

disttotb=SQRT(disttotb)
vdisttotb=SQRT(vdisttotb)

WRITE(60,*) disttotb
WRITE(60,*) vdisttotb
  

  
END DO




WRITE(60,*) '------------------------------------------------'







OPEN(70,file='rotated'//trim(zahl)//'.dat')

WRITE(70,*) n
WRITE(70,*) 3
WRITE(70,*) 1640

DO i=1,n 
WRITE(70,*) '1'
END DO

DO i=1,n 
WRITE(70,*) pos(1:3,i)
END DO

DO i=1,n 
WRITE(70,*) vel(1:3,i)
END DO

CLOSE(70)









END DO



CLOSE(60)






END PROGRAM