!************************************************************************************
SUBROUTINE f_RealRef(mode,q,res,Nlayers,LayerMatrix,fRealRef)    ! Incorporation of Instrumental Resolution in the Reflectivity Calculation
IMPLICIT NONE
integer, intent(in) :: mode	                   ! 0 = default (17 points integration), 1 = integration with 2*mode+1 points
real*8, intent(in) :: q                       ! neutron wavelength (A)
real*8, intent(in) :: res                     ! instrumental resolution dq/q
integer, intent(in) :: Nlayers                 ! Number of layers in the model
real*8, dimension(0:Nlayers+1,0:4), intent(in) :: LayerMatrix
                                               ! Matrix containing the real (1st column) and imaginary (2nd column) scattering
                                               ! length density of each layer. The 3rd column corresponds to the thickness (A) of
                                               ! each layer and the 4th column to the roughness (A).
                                               ! row 1 and Nlayers+1 concern the fronting and backing material respectively
real*8, intent(out) :: fRealRef
!---------------------------------------------------------------------
INTEGER i,j                                    ! counter
REAL*8 qq                                      ! wavevector of points around q, = q[1+(a/n)(dq/q)], where dq/q is the resolution
REAL*8 gweight                                 ! Gaussian weight = exp[-2(a/n)^2], n = number of points considered, a = point index
REAL*8 pi
REAL*8 deltaq
COMPLEX*16 ii
COMPLEX*16, dimension(0:Nlayers+1) :: l1
REAL*8 deltaq2,gfact,dx

pi=3.141592653589793238D0
ii=(0.0D0,1.0D0)

fRealRef=0.0D0
deltaq=q*res/2.354820D0 ! FWHM


! Layer Matric calcs
! Force the correct branch cut for sqrt below the critical edge by addition of a tiny number
do j=0,Nlayers+1
	if(LayerMatrix(Nlayers+1,4).EQ.1.0D0.AND.LayerMatrix(0,4).EQ.0.0D0) then
    	l1(j)=(LayerMatrix(j,0)*(1.0D0-LayerMatrix(j,4))+LayerMatrix(Nlayers+1,0)*LayerMatrix(j,4))+ii*(LayerMatrix(j,1)*&
        (1.0D0-LayerMatrix(j,4))+LayerMatrix(Nlayers+1,1)*LayerMatrix(j,4))-LayerMatrix(0,0)+ii*1e-30
    endif
        
    if(LayerMatrix(Nlayers+1,4).EQ.0.0D0.AND.LayerMatrix(0,4).EQ.1.0D0) then
    	l1(j)=(LayerMatrix(j,0)*(1.0D0-LayerMatrix(j,4))+LayerMatrix(0,0)*LayerMatrix(j,4))+ii*(LayerMatrix(j,1)*&
        (1.0D0-LayerMatrix(j,4))+LayerMatrix(0,1)*LayerMatrix(j,4))-LayerMatrix(0,0)+ii*1e-30
    endif
        
    if(LayerMatrix(Nlayers+1,4).EQ.0.0D0.AND.LayerMatrix(0,4).EQ.0.0) then
    	l1(j)=LayerMatrix(j,0)+ii*LayerMatrix(j,1)-LayerMatrix(0,0)+ii*1e-30
    endif
enddo


! Gaussian convolution
! from -4*sigma to 4*sigma, 17 point evaluation
if(deltaq.EQ.0.0D0) then
    qq=q
    fRealRef=fRealRef+Reflectivity()
else
    if(mode.EQ.0) then ! integration with midpoint rule, 17 point evaluation
    	gfact=(1.0D0/(deltaq*sqrt(2.0D0*pi)))
    	deltaq2=2.0D0*(deltaq**2)
        do i=-8,8
        	dx=8.0D0*deltaq/17.0D0
            qq=q+dble(i)*dx
            gweight=gfact*dexp(-((qq-q)**2)/(deltaq2))*dx
            fRealRef=fRealRef+gweight*Reflectivity()
        enddo
					
    else ! integration with midpoint rule, 2*mode+1 points evaluation
    	gfact=(1.0D0/(deltaq*sqrt(2.0D0*pi)))
    	deltaq2=2.0D0*(deltaq**2)
        do i=-mode,mode
        	dx=8.0D0*deltaq/dble(2*mode+1)
            qq=q+dble(i)*dx
            gweight=gfact*dexp(-((qq-q)**2)/(deltaq2))*dx
            fRealRef=fRealRef+gweight*Reflectivity()
        enddo
    endif
endif

contains
    !************************************************************************************
    real*8 FUNCTION Reflectivity()   ! Calculation of Reflectivity based on the Matrix Method
    !IMPLICIT NONE
    COMPLEX*16 kz
    !COMPLEX*16 k1,k2  ! k_{n},k_{n+1}    
	COMPLEX*16, dimension(0:Nlayers+1) :: k
    COMPLEX*16 b

    COMPLEX*8 r
    REAL*8 d
    COMPLEX*16 m11,m12,m21,m22
    COMPLEX*16 k11,k12,k21,k22
    COMPLEX*16 h11,h12,h21,h22
    !COMPLEX*16 ii
    !COMPLEX*16 l1,l2

    !INTEGER j                            ! counter
    COMPLEX*16 kz2
    COMPLEX*16 iixb,cdexpiixb,cdexpmiixb

                                                   !  theory: kn=sqrt[kz^2-4pi(rhon-rho0)]
                                                   !          r_{n,n+1}=[(k_{n}-k_{n+1})/(k_{n}+k_{n+1})] exp(-2k_{n}k_{n+1}sigma_{n,n+1}^2)
                                                   !          b_{0}=0, b_{n}=ik_{n}d_{n}
                                                   !
                                                   !                  |        exp(-bn)           r_{n,n+1} exp(-bn) |
                                                   !          C_{n} = |  r_{n,n+1} exp(-bn)            exp (-bn)     |
                                                   !
                                                   !          M=C0 C1 C2 ... Cn,        R=|M10/M00|^2


    !ii=(0.0D0,1.0D0)
    !kz=CMPLX(qq/2.0D0,0.0D0)
    !kz=qq/2.0D0
	!kz2=kz**2


    m11=(1.0D0,0.0D0)
    m12=(0.0D0,0.0D0)
    m21=(0.0D0,0.0D0)
    m22=(1.0D0,0.0D0)

	do j=0,Nlayers+1
		k(j)=sqrt((qq**2)/4.0D0-4.0D0*pi*l1(j))
	enddo

    do j=0,Nlayers
		
        r=((k(j)-k(j+1))/(k(j)+k(j+1)))*exp(-2.0D0*k(j)*k(j+1)*LayerMatrix(j,3)**2)

        d=LayerMatrix(j,2)

        if(j.EQ.0) then
            b=0.0D0
        else
            b=k(j)*d
      	endif

		iixb=ii*b
		cdexpiixb=exp(iixb)
		cdexpmiixb=1.0D0/cdexpiixb !exp(-iixb)
		
        k11=cdexpiixb
        k12=r*cdexpiixb
        k21=r*cdexpmiixb
        k22=cdexpmiixb

        h11=m11 
        h12=m12 
        h21=m21 
        h22=m22

        m11=h11*k11+h12*k21 
        m12=h11*k12+h12*k22 
        m21=h21*k11+h22*k21 
        m22=h21*k12+h22*k22

    enddo

    Reflectivity=REAL((m21*CONJG(m21))/(m11*CONJG(m11)))
    END FUNCTION
    !************************************************************************************
END SUBROUTINE

!************************************************************************************
SUBROUTINE f_profilecalc(z, LayerMatrix, Nlayers,fProfile)
!z, distance from the surface along the z axis
!LayerMatrix, matrix containing the real (1st column) and imaginary (2nd column) scattering
! length density of each layer. The 3rd column corresponds to the thickness (A) of
! each layer and the 4th column to the roughness (A). 5th column layer hydration
! row 1 and Nlayers+1 concern the fronting and backing material respectively
! Nlayers, Number of layers in the model
IMPLICIT NONE
real*8, intent(in) :: z                       ! distance from the surface along the z axis
real*8, dimension(0:Nlayers+1,0:4), intent(in) :: LayerMatrix
                                               ! Matrix containing the real (1st column) and imaginary (2nd column) scattering
                                               ! length density of each layer. The 3rd column corresponds to the thickness (A) of
                                               ! each layer and the 4th column to the roughness (A).
                                               ! row 1 and Nlayers+1 concern the fronting and backing 
integer, intent(in) :: Nlayers                 ! Number of layers in the model
real*8, intent(out) :: fProfile
!---------------------------------------------------------------------
INTEGER i,j                                    ! counter
real*8 zi,lmi1,lmi2,rho

rho=LayerMatrix(0,0)
do i=2,Nlayers+2
	zi=0
	do j=1,i-1
		zi=zi+LayerMatrix(j-1,2)
	enddo

	if(LayerMatrix(Nlayers+1,4).EQ.1.0D0.AND.LayerMatrix(0,4).EQ.0.0D0) then
		lmi1=LayerMatrix(i-1,0)*(1.0-LayerMatrix(i-1,4))+LayerMatrix(Nlayers+1,0)*LayerMatrix(i-1,4)
		lmi2=LayerMatrix(i-2,0)*(1.0-LayerMatrix(i-2,4))+LayerMatrix(Nlayers+1,0)*LayerMatrix(i-2,4)
	endif
	if(LayerMatrix(Nlayers+1,4).EQ.0.0D0.AND.LayerMatrix(0,4).EQ.1.0D0) then
		lmi1=LayerMatrix(i-1,0)*(1.0-LayerMatrix(i-1,4))+LayerMatrix(0,0)*LayerMatrix(i-1,4)
		lmi2=LayerMatrix(i-2,0)*(1.0-LayerMatrix(i-2,4))+LayerMatrix(0,0)*LayerMatrix(i-2,4)
	endif
	if(LayerMatrix(Nlayers+1,4).EQ.0.0D0.AND.LayerMatrix(0,4).EQ.0.0D0) then
		lmi1=LayerMatrix(i-1,0)
		lmi2=LayerMatrix(i-2,0)
	endif

	if(LayerMatrix(i-2,3).NE.0.0D0) then
		rho=rho+((lmi1-lmi2)/2.0D0)*(1.0D0+erf((z-zi)/(sqrt(2.0D0)*LayerMatrix(i-2,3))))
	else
		if(z.GE.zi) then
			rho=rho+((lmi1-lmi2)/2.0D0)*(2.0D0)
		else
			rho=rho+0.0D0	
		endif
	endif
enddo

fProfile=rho

END SUBROUTINE

!************************************************************************************
SUBROUTINE f_solventcalc(z, LayerMatrix, Nlayers,fsolvent)
!z, distance from the surface along the z axis
!LayerMatrix, matrix containing the real (1st column) and imaginary (2nd column) scattering
! length density of each layer. The 3rd column corresponds to the thickness (A) of
! each layer and the 4th column to the roughness (A). 5th column layer hydration
! row 1 and Nlayers+1 concern the fronting and backing material respectively
! Nlayers, Number of layers in the model
IMPLICIT NONE
real*8, intent(in) :: z                       ! distance from the surface along the z axis
real*8, dimension(0:Nlayers+1,0:4), intent(in) :: LayerMatrix
                                               ! Matrix containing the real (1st column) and imaginary (2nd column) scattering
                                               ! length density of each layer. The 3rd column corresponds to the thickness (A) of
                                               ! each layer and the 4th column to the roughness (A).
                                               ! row 1 and Nlayers+1 concern the fronting and backing 
integer, intent(in) :: Nlayers                 ! Number of layers in the model
real*8, intent(out) :: fsolvent
!---------------------------------------------------------------------
real*8 rho,zi,lmi1,lmi2
integer i,j


if(LayerMatrix(Nlayers+1,4).EQ.0.0D0.AND.LayerMatrix(0,4).EQ.0.0D0) then
	rho=0.0
else
	if(LayerMatrix(Nlayers+1,4).EQ.1.0D0.AND.LayerMatrix(0,4).EQ.0.0D0) rho=0.0
	if(LayerMatrix(Nlayers+1,4).EQ.0.0D0.AND.LayerMatrix(0,4).EQ.1.0D0) rho=1.0

	do i=2,Nlayers+2
		zi=0
		do j=1,i-1
			zi=zi+LayerMatrix(j-1,2)
		enddo

		if(LayerMatrix(Nlayers+1,4).EQ.1.0D0.AND.LayerMatrix(0,4).EQ.0.0D0) then
			lmi1=LayerMatrix(i-1,4)
			lmi2=LayerMatrix(i-2,4)
		endif
		if(LayerMatrix(Nlayers+1,4).EQ.0.0D0.AND.LayerMatrix(0,4).EQ.1.0D0) then
			lmi1=LayerMatrix(i-1,4)
			lmi2=LayerMatrix(i-2,4)
		endif

		if(LayerMatrix(i-2,3).NE.0.0D0) then
			rho=rho+((lmi1-lmi2)/2.0)*(1.0+erf((z-zi)/(sqrt(2.0)*LayerMatrix(i-2,3))))
		else
			if(z.GT.zi) then
				rho=rho+((lmi1-lmi2)/2.0)*(2.0)
			else
				rho=rho+0.0D0
			endif
		endif	
	enddo

endif

fsolvent=rho

END SUBROUTINE