!####################################################################################
!                          NUMERICAL RECIPES ROUTINES
!####################################################################################

!####################################################################################
!                          NUMERICAL INTEGRATION
!####################################################################################

!================================================================
!                          ODEINT 
!================================================================

MODULE ode_path
	USE nrtype
	INTEGER(I4B) :: nok,nbad,kount
	LOGICAL(LGT), SAVE :: save_steps=.false.
	REAL(SP) :: dxsav
	REAL(SP), DIMENSION(:), POINTER :: xp
	REAL(SP), DIMENSION(:,:), POINTER :: yp
END MODULE ode_path

	SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs,hnext)
	USE nrtype; USE nrutil, ONLY : nrerror,reallocate
	USE ode_path
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: ystart
	REAL(SP), INTENT(IN) :: x1,x2,eps,h1,hmin
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
!BL
		SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
		REAL(SP), INTENT(INOUT) :: x
		REAL(SP), INTENT(IN) :: htry,eps
		REAL(SP), INTENT(OUT) :: hdid,hnext
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			IMPLICIT NONE
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE rkqs
	END INTERFACE
	REAL(SP), PARAMETER :: TINY=1.0e-30_sp
	INTEGER(I4B), PARAMETER :: MAXSTP=100000000
	INTEGER(I4B) :: nstp
	REAL(SP) :: h,hdid,hnext,x,xsav
	REAL(SP), DIMENSION(size(ystart)) :: dydx,y,yscal
	x=x1
	h=sign(h1,x2-x1)
	nok=0
	nbad=0
	kount=0
	y(:)=ystart(:)
	if (save_steps) then
		xsav=x-2.0_sp*dxsav
		nullify(xp,yp)
		allocate(xp(256))
		allocate(yp(size(ystart),size(xp)))
	end if
	do nstp=1,MAXSTP
		call derivs(x,y,dydx)
		yscal(:)=abs(y(:))+abs(h*dydx(:))+TINY
		if (save_steps .and. (abs(x-xsav) > abs(dxsav))) &
			call save_a_step
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x
		call rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)
		if (hdid == h) then
			nok=nok+1
		else
			nbad=nbad+1
		end if
		if ((x-x2)*(x2-x1) >= 0.0) then
			ystart(:)=y(:)
			if (save_steps) call save_a_step
			RETURN
		end if
		if (abs(hnext) < hmin)&
			hnext=hmin!call nrerror('stepsize smaller than minimum in odeint')
		h=hnext
	end do
	!call nrerror('too many steps in odeint')
	CONTAINS
!BL
	SUBROUTINE save_a_step
	kount=kount+1
	if (kount > size(xp)) then
		xp=>reallocate(xp,2*size(xp))
		yp=>reallocate(yp,size(yp,1),size(xp))
	end if
	xp(kount)=x
	yp(:,kount)=y(:)
	xsav=x
	END SUBROUTINE save_a_step

	END SUBROUTINE odeint
!================================================================
!                          rkck
!================================================================

SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
  USE nrtype; USE nrutil, ONLY : assert_eq
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
  REAL(SP), INTENT(IN) :: x,h
  REAL(SP), DIMENSION(:), INTENT(OUT) :: yout,yerr
  INTERFACE
     SUBROUTINE derivs(x,y,dydx)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: x
       REAL(SP), DIMENSION(:), INTENT(IN) :: y
       REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
     END SUBROUTINE derivs
  END INTERFACE
  INTEGER(I4B) :: ndum
  REAL(SP), DIMENSION(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
  REAL(SP), PARAMETER :: A2=0.2_sp,A3=0.3_sp,A4=0.6_sp,A5=1.0_sp,&
       A6=0.875_sp,B21=0.2_sp,B31=3.0_sp/40.0_sp,B32=9.0_sp/40.0_sp,&
       B41=0.3_sp,B42=-0.9_sp,B43=1.2_sp,B51=-11.0_sp/54.0_sp,&
       B52=2.5_sp,B53=-70.0_sp/27.0_sp,B54=35.0_sp/27.0_sp,&
       B61=1631.0_sp/55296.0_sp,B62=175.0_sp/512.0_sp,&
       B63=575.0_sp/13824.0_sp,B64=44275.0_sp/110592.0_sp,&
       B65=253.0_sp/4096.0_sp,C1=37.0_sp/378.0_sp,&
       C3=250.0_sp/621.0_sp,C4=125.0_sp/594.0_sp,&
       C6=512.0_sp/1771.0_sp,DC1=C1-2825.0_sp/27648.0_sp,&
       DC3=C3-18575.0_sp/48384.0_sp,DC4=C4-13525.0_sp/55296.0_sp,&
       DC5=-277.0_sp/14336.0_sp,DC6=C6-0.25_sp
  ndum=assert_eq(size(y),size(dydx),size(yout),size(yerr),'rkck')
  ytemp=y+B21*h*dydx
  call derivs(x+A2*h,ytemp,ak2)
  ytemp=y+h*(B31*dydx+B32*ak2)
  call derivs(x+A3*h,ytemp,ak3)
  ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
  call derivs(x+A4*h,ytemp,ak4)
  ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
  call derivs(x+A5*h,ytemp,ak5)
  ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
  call derivs(x+A6*h,ytemp,ak6)
  yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
  yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
END SUBROUTINE rkck

!================================================================
!                          rkqs
!================================================================

SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  USE nr, ONLY : rkck
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
  REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
  REAL(SP), INTENT(INOUT) :: x
  REAL(SP), INTENT(IN) :: htry,eps
  REAL(SP), INTENT(OUT) :: hdid,hnext
  INTERFACE
     SUBROUTINE derivs(x,y,dydx)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: x
       REAL(SP), DIMENSION(:), INTENT(IN) :: y
       REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
     END SUBROUTINE derivs
  END INTERFACE
  INTEGER(I4B) :: ndum
  REAL(SP) :: errmax,h,htemp,xnew
  REAL(SP), DIMENSION(size(y)) :: yerr,ytemp
  REAL(SP), PARAMETER :: SAFETY=0.9_sp,PGROW=-0.2_sp,PSHRNK=-0.25_sp,&
       ERRCON=1.89e-4
  ndum=assert_eq(size(y),size(dydx),size(yscal),'rkqs')
  h=htry
  do
     call rkck(y,dydx,x,h,ytemp,yerr,derivs)
     errmax=maxval(abs(yerr(:)/yscal(:)))/eps
     if (errmax <= 1.0) exit
     htemp=SAFETY*h*(errmax**PSHRNK)
     h=sign(max(abs(htemp),0.1_sp*abs(h)),h)
     xnew=x+h
!        modified
!     if (xnew == x) call nrerror('stepsize underflow in rkqs')
     IF (xnew == x) THEN
        CALL nrerror('stepsize underflow in rkqs')
        RETURN
     ENDIF
  end do
  if (errmax > ERRCON) then
     hnext=SAFETY*h*(errmax**PGROW)
  else
     hnext=5.0_sp*h
  end if
  hdid=h
  x=x+h
  y(:)=ytemp(:)
END SUBROUTINE rkqs


!####################################################################################
!                          RANDOM NUMBER GENERATOR
!####################################################################################

!================================================================
!                          ran_state
!================================================================

MODULE ran_state
  USE nrtype
  IMPLICIT NONE
  INTEGER, PARAMETER :: K4B=selected_int_kind(9)
  INTEGER(K4B), PARAMETER :: hg=huge(1_K4B), hgm=-hg, hgng=hgm-1
  INTEGER(K4B), SAVE :: lenran=0, seq=0
  INTEGER(K4B), SAVE :: iran0,jran0,kran0,nran0,mran0,rans
  INTEGER(K4B), DIMENSION(:,:), POINTER, SAVE :: ranseeds
  INTEGER(K4B), DIMENSION(:), POINTER, SAVE :: iran,jran,kran, &
       nran,mran,ranv
  REAL(SP), SAVE :: amm
  INTERFACE ran_hash
     MODULE PROCEDURE ran_hash_s, ran_hash_v
  END INTERFACE
CONTAINS
  !BL
  SUBROUTINE ran_init(length)
    USE nrtype; USE nrutil, ONLY : arth,nrerror,reallocate
    IMPLICIT NONE
    INTEGER(K4B), INTENT(IN) :: length
    INTEGER(K4B) :: new,j,hgt
    if (length < lenran) RETURN
    hgt=hg
    if (hg /= 2147483647) call nrerror('ran_init: arith assump 1 fails')
    if (hgng >= 0)        call nrerror('ran_init: arith assump 2 fails')
    if (hgt+1 /= hgng)    call nrerror('ran_init: arith assump 3 fails')
    if (not(hg) >= 0)     call nrerror('ran_init: arith assump 4 fails')
    if (not(hgng) < 0)    call nrerror('ran_init: arith assump 5 fails')
    if (hg+hgng >= 0)     call nrerror('ran_init: arith assump 6 fails')
    if (not(-1_k4b) < 0)  call nrerror('ran_init: arith assump 7 fails')
    if (not(0_k4b) >= 0)  call nrerror('ran_init: arith assump 8 fails')
    if (not(1_k4b) >= 0)  call nrerror('ran_init: arith assump 9 fails')
    if (lenran > 0) then
       ranseeds=>reallocate(ranseeds,length,5)
       ranv=>reallocate(ranv,length-1)
       new=lenran+1
    else
       allocate(ranseeds(length,5))
       allocate(ranv(length-1))
       new=1
       amm=nearest(1.0_sp,-1.0_sp)/hgng
       if (amm*hgng >= 1.0 .or. amm*hgng <= 0.0) &
            call nrerror('ran_init: arth assump 10 fails')
    end if
    ranseeds(new:,1)=seq
    ranseeds(new:,2:5)=spread(arth(new,1,size(ranseeds(new:,1))),2,4)
    do j=1,4
       call ran_hash(ranseeds(new:,j),ranseeds(new:,j+1))
    end do
    where (ranseeds(new:,1:3) < 0) &
         ranseeds(new:,1:3)=not(ranseeds(new:,1:3))
    where (ranseeds(new:,4:5) == 0) ranseeds(new:,4:5)=1
    if (new == 1) then
       iran0=ranseeds(1,1)
       jran0=ranseeds(1,2)
       kran0=ranseeds(1,3)
       mran0=ranseeds(1,4)
       nran0=ranseeds(1,5)
       rans=nran0
    end if
    if (length > 1) then
       iran => ranseeds(2:,1)
       jran => ranseeds(2:,2)
       kran => ranseeds(2:,3)
       mran => ranseeds(2:,4)
       nran => ranseeds(2:,5)
       ranv = nran
    end if
    lenran=length
  END SUBROUTINE ran_init
  !BL
  SUBROUTINE ran_deallocate
    if (lenran > 0) then
       deallocate(ranseeds,ranv)
       nullify(ranseeds,ranv,iran,jran,kran,mran,nran)
       lenran = 0
    end if
  END SUBROUTINE ran_deallocate
  !BL
  SUBROUTINE ran_seed(sequence,size,put,get)
    IMPLICIT NONE
    INTEGER, OPTIONAL, INTENT(IN) :: sequence
    INTEGER, OPTIONAL, INTENT(OUT) :: size
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: put
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: get
    if (present(size)) then
       size=5*lenran
    else if (present(put)) then
       if (lenran == 0) RETURN
       ranseeds=reshape(put,shape(ranseeds))
       where (ranseeds(:,1:3) < 0) ranseeds(:,1:3)=not(ranseeds(:,1:3))
       where (ranseeds(:,4:5) == 0) ranseeds(:,4:5)=1
       iran0=ranseeds(1,1)
       jran0=ranseeds(1,2)
       kran0=ranseeds(1,3)
       mran0=ranseeds(1,4)
       nran0=ranseeds(1,5)
    else if (present(get)) then
       if (lenran == 0) RETURN
       ranseeds(1,1:5)=(/ iran0,jran0,kran0,mran0,nran0 /)
       get=reshape(ranseeds,shape(get))
    else if (present(sequence)) then
       call ran_deallocate
       seq=sequence
    end if
  END SUBROUTINE ran_seed
  !BL
  SUBROUTINE ran_hash_s(il,ir)
    IMPLICIT NONE
    INTEGER(K4B), INTENT(INOUT) :: il,ir
    INTEGER(K4B) :: is,j
    do j=1,4
       is=ir
       ir=ieor(ir,ishft(ir,5))+1422217823
       ir=ieor(ir,ishft(ir,-16))+1842055030
       ir=ieor(ir,ishft(ir,9))+80567781
       ir=ieor(il,ir)
       il=is
    end do
  END SUBROUTINE ran_hash_s
  !BL
  SUBROUTINE ran_hash_v(il,ir)
    IMPLICIT NONE
    INTEGER(K4B), DIMENSION(:), INTENT(INOUT) :: il,ir
    INTEGER(K4B), DIMENSION(size(il)) :: is
    INTEGER(K4B) :: j
    do j=1,4
       is=ir
       ir=ieor(ir,ishft(ir,5))+1422217823
       ir=ieor(ir,ishft(ir,-16))+1842055030
       ir=ieor(ir,ishft(ir,9))+80567781
       ir=ieor(il,ir)
       il=is
    end do
  END SUBROUTINE ran_hash_v
END MODULE ran_state

!================================================================
!                          ran1
!================================================================

SUBROUTINE ran1_s(harvest)
  USE nrtype
  USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
       iran0,jran0,kran0,nran0,mran0,rans
  IMPLICIT NONE
  REAL(SP), INTENT(OUT) :: harvest
  if (lenran < 1) call ran_init(1)
  rans=iran0-kran0
  if (rans < 0) rans=rans+2147483579_k4b
  iran0=jran0
  jran0=kran0
  kran0=rans
  nran0=ieor(nran0,ishft(nran0,13))
  nran0=ieor(nran0,ishft(nran0,-17))
  nran0=ieor(nran0,ishft(nran0,5))
  if (nran0 == 1) nran0=270369_k4b
  mran0=ieor(mran0,ishft(mran0,5))
  mran0=ieor(mran0,ishft(mran0,-13))
  mran0=ieor(mran0,ishft(mran0,6))
  rans=ieor(nran0,rans)+mran0
  harvest=amm*merge(rans,not(rans), rans<0 )
END SUBROUTINE ran1_s

SUBROUTINE ran1_v(harvest)
  USE nrtype
  USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
       iran,jran,kran,nran,mran,ranv
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
  INTEGER(K4B) :: n
  n=size(harvest)
  if (lenran < n+1) call ran_init(n+1)
  ranv(1:n)=iran(1:n)-kran(1:n)
  where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
  iran(1:n)=jran(1:n)
  jran(1:n)=kran(1:n)
  kran(1:n)=ranv(1:n)
  nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
  nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
  nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
  where (nran(1:n) == 1) nran(1:n)=270369_k4b
  mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
  mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
  mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
  ranv(1:n)=ieor(nran(1:n),ranv(1:n))+mran(1:n)
  harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
END SUBROUTINE ran1_v

!####################################################################################
!                          ROOT FINDING,
!####################################################################################

!================================================================
!                          RTSAFE
!================================================================

	FUNCTION rtsafe(funcd,x1,x2,xacc)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,xacc
	REAL(SP) :: rtsafe
	INTERFACE
		SUBROUTINE funcd(x,fval,fderiv)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: fval,fderiv
		END SUBROUTINE funcd
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXIT=100
	INTEGER(I4B) :: j
	REAL(SP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
	call funcd(x1,fl,df)
	call funcd(x2,fh,df)
	if ((fl > 0.0 .and. fh > 0.0) .or. &
		(fl < 0.0 .and. fh < 0.0)) &
		call nrerror('root must be bracketed in rtsafe')
	if (fl == 0.0) then
		rtsafe=x1
		RETURN
	else if (fh == 0.0) then
		rtsafe=x2
		RETURN
	else if (fl < 0.0) then
		xl=x1
		xh=x2
	else
		xh=x1
		xl=x2
	end if
	rtsafe=0.5_sp*(x1+x2)
	dxold=abs(x2-x1)
	dx=dxold
	call funcd(rtsafe,f,df)
	do j=1,MAXIT
		if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) >= 0.0 .or. &
			abs(2.0_sp*f) > abs(dxold*df) ) then
			dxold=dx
			dx=0.5_sp*(xh-xl)
			rtsafe=xl+dx
			if (xl == rtsafe) RETURN
		else
			dxold=dx
			dx=f/df
			temp=rtsafe
			rtsafe=rtsafe-dx
			if (temp == rtsafe) RETURN
		end if
		if (abs(dx) < xacc) RETURN
		call funcd(rtsafe,f,df)
		if (f < 0.0) then
			xl=rtsafe
		else
			xh=rtsafe
		end if
	end do
	call nrerror('rtsafe: exceeded maximum iterations')
	END FUNCTION rtsafe
