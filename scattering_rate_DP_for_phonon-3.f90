program main
implicit none

integer, parameter :: dp = kind(1.0D0)
integer imode,nband,nmode,z,nk1, nk2, nk3, nk, add1,delta
integer i, j, k, ii, jj, kk, ikq, jkq, kkq, iband, jband, ibndmax, ibndmin
real(dp) ef,ef0,smear,fk,fkq,kB,kf,temp,x,arg,weight,sqrtpm1,def,pi,elast,volume,rad2ev,ev2thz, fsthick,ev2j
real(dp),allocatable::dk(:,:,:,:),ene(:,:,:,:),wgk(:,:,:),ddk(:,:),dene(:,:),dwgk(:)
integer,allocatable::add2(:,:,:)
real(dp),allocatable::freq(:,:,:,:),dfreq(:,:),tau(:,:)

kB = 1.380649D-23  ! J/K
sqrtpm1 = 1.0D0/1.77245385090551602729D0
pi = 3.14159265358979323846D0
rad2ev   = 6.58551D-4
ev2thz = 241.7988407D0
ev2j = 1.602176621D-19  ! J
open(unit=9,file="out")
open(unit=10,file='inputDP')
read(10,*) nk1,nk2,nk3,nmode,temp,ef,smear 
read(10,*) def,elast,volume,fsthick
read(10,*) delta
close(10)
write(9,*) "Read inputDP completed"
open(unit=11,file='EIGENVAL')
read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*) z, nk, nband
allocate(add2(nk1,nk2,nk3),dk(3,nk1,nk2,nk3),ene(nband,nk1,nk2,nk3),wgk(nk1,nk2,nk3))
allocate(ddk(3,nk),dene(nband,nk),dwgk(nk))
do k = 1, nk3
   do j = 1, nk2
      do i = 1, nk1
         read(11,*)
         read(11,*) dk(1:3,i,j,k), wgk(i,j,k)
         do iband = 1, nband
            read(11,*) z, ene(iband,i,j,k)     ! eV
         enddo
      enddo
   enddo
enddo
close(11)
write(9,*) "Read EIGENVAL completed"
allocate(freq(nmode,nk1,nk2,nk3),dfreq(nmode,nk),tau(nmode,nk))
open(unit=12,file='BTE.omega_full')
do k = 1, nk3
   do j = 1, nk2
      do i = 1, nk1
         read(12,*) freq(:,i,j,k)   ! rad/ps
      enddo
   enddo
enddo
close(12)
freq=freq*rad2ev*ev2thz    ! in THz
write(9,*) "Read BTE.omega_full completed"
add1=0
do k = 1, nk3
   do j = 1, nk2
      do i = 1, nk1
         add2(i,j,k)=add1+i
         ddk(:,add2(i,j,k))=dk(:,i,j,k)
         dene(:,add2(i,j,k))=ene(:,i,j,k)
         dwgk(add2(i,j,k))=wgk(i,j,k)
         dfreq(:,add2(i,j,k))=freq(:,i,j,k)
      enddo
      add1=add1+nk1
   enddo
enddo
write(9,*) "Change 3D parameters to 1D completed"
!ibndmin=100000
!ibndmax=0
!do k = 1, nk3
!   do j = 1, nk2
!      do i = 1, nk1
!         do iband = 1, nband
!            if (abs(ene(iband,i,j,k)-ef0).lt.fsthick) then
!               ibndmin = min(iband,ibndmin)
!               ibndmax = max(iband,ibndmax)
!            endif
!         enddo
!      enddo
!   enddo
!enddo
!write(9,*) 'ibndmin = ', ibndmin, 'ibndmax = ', ibndmax
tau=0.0
do k = 1, nk3
   do j = 1, nk2
      do i = 1, nk1
         do kk = 1, nk3
            do jj = 1, nk2
               do ii = 1, nk1
                  if ((i+ii)>(nk1+1)) then
                     ikq = i+ii-nk1-1
                  else
                     ikq = i+ii-1
                  endif
                  if ((j+jj)>(nk2+1)) then
                     jkq = j+jj-nk2-1
                  else
                     jkq = j+jj-1
                  endif
                  if ((k+kk)>(nk3+1)) then
                     kkq = k+kk-nk3-1
                  else
                     kkq = k+kk-1
                  endif
                  do iband = 1, nband
                     do jband = 1, nband
                        if ((abs(dene(iband,add2(ii,jj,kk))-ef).lt.fsthick).and.&
                           (abs(dene(jband,add2(ikq,jkq,kkq))-ef).lt.fsthick))  then
                           fk = 1.0D0/(1.0D0+exp((dene(iband,add2(ii,jj,kk))-ef)*ev2j/kB/temp))
                           fkq = 1.0D0/(1.0D0+exp((dene(jband,add2(ikq,jkq,kkq))-ef)*ev2j/kB/temp))
                           do imode = 1, nmode
                              if (delta==1) then   ! derivative of Methfessel-Paxton
                                 x = (dene(iband,add2(ii,jj,kk))-dene(jband,add2(ikq,jkq,kkq))-dfreq(imode,add2(i,j,k))/ev2thz)/smear
                                 arg = min(200.0,x**2)
                                 weight = exp(-arg)*sqrtpm1/smear
                              else if (delta==2) then ! Gaussian 
                                 x = (dene(iband,add2(ii,jj,kk))-dene(jband,add2(ikq,jkq,kkq))-dfreq(imode,add2(i,j,k))/ev2thz)/smear/sqrt(2.0D0)
                                 weight = exp(-x**2)/smear/sqrt(2.0D0*pi)
                              else if (delta==3) then ! Lorenz
                                 x = dene(iband,add2(ii,jj,kk))-dene(jband,add2(ikq,jkq,kkq))-dfreq(imode,add2(i,j,k))/ev2thz
                                 weight=smear/(x**2+smear**2)/pi
                              endif
                              tau(imode,add2(i,j,k))=tau(imode,add2(i,j,k))+160.2176621D0*dwgk(add2(ii,jj,kk))*pi*(def**2)*dfreq(imode,add2(i,j,k))*(fkq-fk)*weight/volume/elast   ! THz
                           enddo ! imode
                          endif
                        enddo ! jband
                     enddo ! iband
               enddo  ! ii
            enddo ! jj
         enddo  ! kk
      enddo ! i
   enddo ! j
enddo  ! k
write(9,*) "Calculate tau^-1 completed"
open(unit=17,file="tau-1.txt")
write(17,"(2A20)") "Freq (THz)", "tau^-1 (THz)"
do k = 1, nk3
   do j = 1, nk2
      do i = 1, nk1
         do imode = 1, nmode
            write(17,"(50000ES20.10)") dfreq(imode,add2(i,j,k)),tau(imode,add2(i,j,k))
         enddo
      enddo
   enddo
enddo

open(unit=18,file="linewidth.phself")
write(18,*)
write(18,"(A9,A6,2A22)") "Q-point", "Mode", "Phonon freq (meV) ","Phonon linewidth (meV)"
do k = 1, nk3
   do j = 1, nk2
      do i = 1, nk1
         do imode = 1, nmode
            write(18,"(I9,I6,2ES20.8)") add2(i,j,k),imode,dfreq(imode,add2(i,j,k))*4.1356,tau(imode,add2(i,j,k))*4.1356
         enddo
      enddo
   enddo
enddo

close(9)
close(17)
close(18)
      
end
