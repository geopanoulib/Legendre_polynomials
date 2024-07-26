! 
! Copyright (c) 2024, Georgios Panou
! 
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your
! option) any later version.
! 
! This program is distributed in the hope that it will be useful, but WITHOUT 
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License 
! along with this program. If not, see <https://www.gnu.org/licenses/>.
! 
! Authors: Iossifidis C., Koci J. and Panou G. <geopanou@survey.ntua.gr>
! 

program geolp

implicit none
integer*8 NX, NMAX, n, p, q, t
real*8 theta, PIHALF, raddeg, cosi_1, cosi_2, coeff, cosix
real*8 A, B, coss, Pi, pnn, pn0, C, Tni

parameter (NMAX=1000)

real*8 cosn(0:NMAX) ! multiple angle cosines
real*8 Pn_0(0:NMAX) ! Legendre Polynomials (m=0)

! Initialize
NX=10 ! The maximum degree, less than NMAX
theta=30.0d0 ! The co-latitude in degrees
PIHALF=atan(1.d0)*2.d0
raddeg=PIHALF/90.d0
theta=raddeg*theta

! Computes the multiple angle cosines by using the Chebyshev's algorithm
! https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Chebyshev_method

cosi_1=dcos(theta)
coeff=2.d0*cosi_1  
cosi_2=1.d0
cosn(0)=1.d0
cosn(1)=cosi_1
do n=2,NX
    cosix=coeff*cosi_1-cosi_2
    cosn(n)=cosix
    cosi_2=cosi_1
    cosi_1=cosix
enddo

write(*,"(a20,i20)") "NX=",NX
write(*,"(a20,1pe20.8)") "theta(deg)=",theta/raddeg
write(*,"(1a10,1a25)") "n", "Pn(theta)"

! Computes the fully normalized Legendre Polynomials (m = 0)

n=0
Pn_0(n)=1.d0

n=1
Pn_0(n)=dsqrt(3.0d0)*cosn(1)

! Computes the odd Legendre polynomials (m = 0)
pnn=1.d0
coss=cosn(1)
do n=3,NX,2
    p=n-1
    ! Computes the Dn_n00 ratio
    A=1.d0/p
    B=1.d0/n
    coeff=1.d0-0.75d0*B
    pnn=pnn*(1.d0-A*coeff)
    Pi=coss
    q=n+2
    t=3
    do while (p.gt.0)
        ! Computes the Tni ratio
        B=1.d0-1.d0/p
        C=1.d0+1.d0/q
        Tni=B*C
        Pi=Pi*Tni
        Pi=Pi+cosn(t)
        t=t+2
        p=p-2
        q=q+2
    enddo
    coeff=dsqrt(2.0d0*n+1)
    Pn_0(n)=Pi*pnn*coeff
enddo

! Computes the even Legendre Polynomials (m = 0)
pn0=1.d0
pnn=2.d0
coss=cosn(2)
do n=2,NX,2
    p=n-1
    
    ! Computes the Dn_n00 ratio
    A=1.d0/p
    B=1.d0/n
    coeff=1.d0-0.75d0*B
    pnn=pnn*(1.d0-A*coeff)
    
    ! Computes the tn ratio
    A=1.d0-B
    pn0=pn0*A*A
    Pi=coss
    
    ! Computation of the pni coefficients and each Polynomial from them */
    t=4
    q=n+3
    p=p-1
    do while (p.gt.0)
        ! Computes the Tni ratio
        B=1.d0-1.d0/p
        C=1.d0+1.d0/q
        Tni=B*C
        Pi=Pi*Tni
        Pi=Pi+cosn(t)
        t=t+2
        p=p-2
        q=q+2
    enddo
    coeff=dsqrt(2.0d0*n+1)
    Pn_0(n)=coeff*(Pi*pnn+pn0)
enddo

! Printing the Legendre Polynomials
do n=0,NX
    write(*,"(1i10,1p1e25.16)") n,Pn_0(n)
enddo

stop
end program
