PROGRAM Creator
use struct
IMPLICIT NONE

real(kind=4) :: rpar(5)
integer :: ipar(5)
type(particle),allocatable :: pd(:)
character(len=10) :: filename
integer :: k,nbd_tot,modus,m_tot

200 FORMAT(1X,E16.4,1X,E16.4,1X,E16.4)


write(*,*) 'Binäre Datei erzeugen = 1'
write(*,*) 'Model-Datei lesen =2'
write(*,*) 'Model-Datei lesen nur Koordinaten =3'
read(*,*) modus

if (modus==1) THEN

write(*,*) 'Dateiname eingeben'
read(*,*) filename 

open(10,file=filename,form='formatted')
open(11,file='mout00.bin',form='unformatted')
read(10,*) nbd_tot
read(10,*) m_tot


ipar(1)=nbd_tot
ipar(2)=0
ipar(3)=1
ipar(4)=0
ipar(5)=0

rpar(1)=m_tot
rpar(2)=0
rpar(3)=0
rpar(4)=0
rpar(5)=0

write(11),ipar(1:5)
write(11),rpar(1:5)

allocate(pd(nbd_tot))

do k=1,nbd_tot
write(*,*) k
read(10,*) 
read(10,200) pd(k)%posvel(1),pd(k)%posvel(2),pd(k)%posvel(3)
read(10,200) pd(k)%posvel(4),pd(k)%posvel(5),pd(k)%posvel(6)

write(11),pd(k)%posvel(1:6)

end do


ELSE IF (modus==2) THEN



write(*,*) 'Dateiname eingeben'
read(*,*) filename

open(10,file=filename,form='unformatted')
open(11,file='lesbar.dat',form='formatted')

read(10),ipar(1:5)
read(10),rpar(1:5)

m_tot=rpar(1)
nbd_tot=ipar(1)
allocate(pd(nbd_tot))

write(11,*),nbd_tot
write(11,*),m_tot
write(*,*) nbd_tot


do k=1,nbd_tot

read(10),pd(k)%posvel(1:6)
write(11,*) ' '
write(11,200) pd(k)%posvel(1),pd(k)%posvel(2),pd(k)%posvel(3)
write(11,200) pd(k)%posvel(4),pd(k)%posvel(5),pd(k)%posvel(6)

end do

ELSE



write(*,*) 'Dateiname eingeben'
read(*,*) filename

open(10,file=filename,form='unformatted')
open(11,file='lesbar.dat',form='formatted')

read(10),ipar(1:5)
read(10),rpar(1:5)

m_tot=rpar(1)
nbd_tot=ipar(1)
allocate(pd(nbd_tot))



do k=1,nbd_tot

read(10),pd(k)%posvel(1:6)

write(11,*) pd(k)%posvel(1),pd(k)%posvel(2),pd(k)%posvel(3)

end do



ENDIF

close(10)
close(11)



write(*,*) 'finished'


END
