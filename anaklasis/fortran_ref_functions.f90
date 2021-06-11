!************************************************************************************
SUBROUTINE f_RealRef(mode,q,res,Nlayers,LayerMatrix,fRealRef)    ! Incorporation of Instrumental Resolution in the Reflectivity Calculation
IMPLICIT NONE
integer, intent(in) :: mode	                   ! 0 = normal, 1 = slow
real(8), intent(in) :: q                       ! neutron wavelength (A)
real(8), intent(in) :: res                     ! instrumental resolution dq/q
integer, intent(in) :: Nlayers                 ! Number of layers in the model
real(8), dimension(0:Nlayers+1,0:4), intent(in) :: LayerMatrix
                                               ! Matrix containing the real (1st column) and imaginary (2nd column) scattering
                                               ! length density of each layer. The 3rd column corresponds to the thickness (A) of
                                               ! each layer and the 4th column to the roughness (A).
                                               ! row 1 and Nlayers+1 concern the fronting and backing material respectively
real(8), intent(out) :: fRealRef
!---------------------------------------------------------------------
INTEGER i,j                                    ! counter
REAL*8 qq                                      ! wavevector of points around q, = q[1+(a/n)(dq/q)], where dq/q is the resolution
REAL*8 gweight                                 ! Gaussian weight = exp[-2(a/n)^2], n = number of points considered, a = point index
REAL*8 pi
REAL*8 deltaq
COMPLEX*16 ii
COMPLEX(8), dimension(0:Nlayers) :: l1,l2
REAL*8 deltaq2,gfact
real(8), dimension(0:8) :: gfactor

gfactor(0)=0.19741265136584743
gfactor(1)=0.17466632194020804
gfactor(2)=0.12097757871001297
gfactor(3)=0.06559061680303818
gfactor(4)=0.027834684208772387
gfactor(6)=0.009244709419990171
gfactor(7)=0.002402738192663789
gfactor(8)=0.0004886077571899516

pi=3.141592653589793238D0
ii=(0.0D0,1.0D0)

fRealRef=0.0D0
deltaq=q*res/2.0D0 ! the width of the Gaussian should be equal to 2*DeltaQ

! Layer Matric calcs
do j=0,Nlayers
	if(LayerMatrix(Nlayers+1,4).EQ.1.0D0.AND.LayerMatrix(0,4).EQ.0.0D0) then
    	l1(j)=(LayerMatrix(j,0)*(1.0D0-LayerMatrix(j,4))+LayerMatrix(Nlayers+1,0)*LayerMatrix(j,4))+ii*(LayerMatrix(j,1)*&
        (1.0D0-LayerMatrix(j,4))+LayerMatrix(Nlayers+1,1)*LayerMatrix(j,4))-LayerMatrix(0,0)-ii*LayerMatrix(0,1)
        l2(j)=(LayerMatrix(j+1,0)*(1.0D0-LayerMatrix(j+1,4))+LayerMatrix(Nlayers+1,0)*LayerMatrix(j+1,4))+ii*(LayerMatrix(j+1,1)&
        *(1.0D0-LayerMatrix(j+1,4))+LayerMatrix(Nlayers+1,1)*LayerMatrix(j+1,4))-LayerMatrix(0,0)-ii*LayerMatrix(0,1)
    endif
        
    if(LayerMatrix(Nlayers+1,4).EQ.0.0D0.AND.LayerMatrix(0,4).EQ.1.0D0) then
    	l1(j)=(LayerMatrix(j,0)*(1.0D0-LayerMatrix(j,4))+LayerMatrix(0,0)*LayerMatrix(j,4))+ii*(LayerMatrix(j,1)*&
        (1.0D0-LayerMatrix(j,4))+LayerMatrix(0,1)*LayerMatrix(j,4))-LayerMatrix(0,0)-ii*LayerMatrix(0,1)
        l2(j)=(LayerMatrix(j+1,0)*(1.0D0-LayerMatrix(j+1,4))+LayerMatrix(0,0)*LayerMatrix(j+1,4))+ii*(LayerMatrix(j+1,1)&
        *(1.0D0-LayerMatrix(j+1,4))+LayerMatrix(0,1)*LayerMatrix(j+1,4))-LayerMatrix(0,0)-ii*LayerMatrix(0,1)
    endif
        
    if(LayerMatrix(Nlayers+1,4).EQ.0.0D0.AND.LayerMatrix(0,4).EQ.0.0) then
    	l1(j)=(LayerMatrix(j,0)*(1.0D0-LayerMatrix(j,4)))+ii*(LayerMatrix(j,1)*&
        (1.0D0-LayerMatrix(j,4)))-LayerMatrix(0,0)-ii*LayerMatrix(0,1)
        l2(j)=(LayerMatrix(j+1,0)*(1.0D0-LayerMatrix(j+1,4)))+ii*(LayerMatrix(j+1,1)&
        *(1.0D0-LayerMatrix(j+1,4)))-LayerMatrix(0,0)-ii*LayerMatrix(0,1)
    endif
enddo


! Gaussian convolution
if(deltaq.EQ.0.0D0) then
    qq=q
    fRealRef=fRealRef+Reflectivity()
else
    if(mode.EQ.0) then
    	!gfact=(1.0D0/(deltaq*sqrt(2.0D0*pi)))
    	deltaq2=2.0D0*deltaq**2
        do i=-8,8
            qq=q+dble(i)*deltaq/2.0D0
            !gweight=gfact*exp(-((qq-q)**2)/(deltaq2))*(deltaq/2.0D0)
            fRealRef=fRealRef+gfactor(abs(i))*Reflectivity()
        enddo
    else
    	gfact=(1.0D0/(deltaq*sqrt(2.0D0*pi)))
    	deltaq2=2.0D0*deltaq**2
        do i=-27,27
            qq=q+dble(i)*deltaq/9.0D0
            gweight=gfact*exp(-((qq-q)**2)/(deltaq2))*(deltaq/9.0D0)
            fRealRef=fRealRef+gweight*Reflectivity()
        enddo
    endif
endif

contains
    !************************************************************************************
    real(8) FUNCTION Reflectivity()   ! Calculation of Reflectivity based on the Matrix Method
    !IMPLICIT NONE
    COMPLEX*16 kz
    COMPLEX*16 k1,k2  ! k_{n},k_{n+1}                              
    COMPLEX*16 b

    COMPLEX*16 r
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
    kz=CMPLX(qq/2.0D0,0.0D0)
	kz2=kz**2


    m11=(1.0D0,0.0D0)
    m12=(0.0D0,0.0D0)
    m21=(0.0D0,0.0D0)
    m22=(1.0D0,0.0D0)

    do j=0,Nlayers
		
        k1=cdsqrt(kz2-4.0D0*pi*l1(j))
        k2=cdsqrt(kz2-4.0D0*pi*l2(j))

        r=((k1-k2)/(k1+k2))*cdexp(-2.0D0*k1*k2*LayerMatrix(j,3)**2)

        d=LayerMatrix(j,2)

        if(j.EQ.0) then
            b=0.0D0
        else
            b=k1*d
      	endif

		iixb=ii*b
		cdexpiixb=cdexp(iixb)
		cdexpmiixb=cdexp(-iixb)
		
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

    Reflectivity=REAL((m21*DCONJG(m21))/(m11*DCONJG(m11)))
    END FUNCTION
    !************************************************************************************
END SUBROUTINE

