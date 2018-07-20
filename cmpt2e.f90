subroutine cmpt2e()
use nrtype; use molprops ; use s_sp_l_terms ; use s_d_l_terms ; use s_p_d_terms
implicit none

integer(i4b) :: ab1, ab2, sh1,sh2,sh3,sh4, ao1, ao2,ao3,ao4, nschw,atc,atd, n1, n2, n3, n4
integer(i4b) :: ii, jj, i, j, k, l, ij, kl
integer(i4b) :: lij,lkl,prevdupl(5000),lijkl

real(dp) :: essss, esssp(3), esspp(3,3), espsp(3,3), esppp(3,3,3), epppp(3,3,3,3)
real(dp) :: esssd(5), essdd(5,5), esdsd(5,5), esddd(5,5,5), edddd(5,5,5,5)
real(dp) :: esspd(3,5), espsd(3,5), esppd(3,3,5), esdpp(5,3,3)
real(dp) :: espdd(3,5,5), esdpd(5,3,5)
real(dp) :: epppd(3,3,3,5), eppdd(3,3,5,5), epdpd(3,5,3,5), epddd(3,5,5,5)
real(dp) :: esssl(4), eslsl(4,4), essll(4,4), eslll(4,4,4), ellll(4,4,4,4)
real(dp) :: epppl(3,3,3,4), eplpl(3,4,3,4), eppll(3,3,4,4), eplll(3,4,4,4)
real(dp) :: esspl(3,4), espsl(3,4), espll(3,4,4), eslpl(4,3,4)
real(dp) :: esppl(3,3,4), eslpp(4,3,3)
real(dp) :: essdl(5,4), esdsl(5,4), esddl(5,5,4), esldd(4,5,5)
real(dp) :: esdll(5,4,4), esldl(4,5,4)
real(dp) :: edddl(5,5,5,4), eddll(5,5,4,4), edldl(5,4,5,4), edlll(5,4,4,4)

! This subroutine manages computation of the two-=electron integrals 
! To do this, we make several classes of shell pairs and pairs of pairs, based on:
!  * their angular momentum s, p, d, l
!  * whether the bra or ket pairs are 'diagonal' i.e. involve the same shell
!  * whether the bra and ket pair and the same
! Each unique integral is then stored based on its ijkl index. Duplicates are avoided as much as possible, but clearly 
!   they will simply overlap each other, so it does not really matter
! Order of funcs, bra/kets: s, p, d, l, |ss), |sp), |sd), |sl), |pp), |pd), |pl), |dd), |dl), |ll)

! First work out how many TEIs of each type there may be, and create arrays to store them
! Except for the fully permutational sort (ab|cd) - we'll estimate those more carefully shortly

allocate(i2eaaaa(nb),aonaaaa(nb))
i2eaaaa=0._dp
aonaaaa=0
naaaa=0
maxaabb=nb*(nb-1)/2
allocate(i2eaabb(maxaabb),aonaabb(maxaabb))
i2eaabb=0._dp
aonaabb=0
naabb=0
maxabab=nb*(nb+1)
allocate(i2eabab(maxabab),aonabab(maxabab))
i2eabab=0._dp
aonabab=0
nabab=0
maxaabc=nb*(nb-1)*(nb-2)
allocate(i2eaabc(maxaabc),aonaabc(2,maxaabc))
i2eaabc=0._dp
aonaabc=0
naabc=0
maxabcd_d=6*nb2
allocate(i2eabcd_d(maxabcd_d),aonabcd_d(2,maxabcd_d))
i2eabcd_d=0._dp
aonabcd_d=0
nabcd_d=0
schw=0._dp

! First calculate the 'diagonal' or 'Schwartz' integrals that have the same bra and ket - but different shells within these, i.e. (ab|ab)
! In due course I'll save these for Schwartz screening, but not in this version

! Counters for number of significant pairs
nsigss=0 ; nsigsp=0 ; nsigsd=0 ; nsigsl=0 ; nsigpp=0 ; nsigpd=0 ; nsigpl=0 ; nsigdd=0 ; nsigdl=0 ; nsigll=0

! By shell pair type.
do ii=1,nss 
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  call i2e_ssss(sh1,sh2,sh1,sh2,essss)
  schw(ab1)=sqrt(abs(essss))
  if (schw(ab1).gt.shthresh) then
    nsigss=nsigss+1
  end if
  if (abs(essss).gt.intthresh) then
    call save_abcd_integral_d(ao1,ao2,ao1,ao2,essss)
  end if
end do

! Now do (sp| pairs, which yield 3 (ab|ab) integrals and 3 unique (ab|ac)
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)   ! First problem. sh2 can be smaller than sh1
  ao2=aosh(sh2)
  call i2e_spsp(sh1,sh2,sh1,sh2,espsp)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(espsp)))
  if (schw(ab1).gt.shthresh) then
    nsigsp=nsigsp+1
  end if
  ! Here by permutation, we know that espsp(i,j) must equal espsp(j,i) so only store one of them
  do i=1,3
    do j=1,i
      if (abs(espsp(i,j)).gt.intthresh) then
        call save_abcd_integral_d(ao1,ao2+i-1,ao1,ao2+j-1,espsp(i,j))
      end if
    end do
  end do
end do

! Now do (sd| pairs
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)   ! First problem. sh2 can be smaller than sh1
  ao2=aosh(sh2)
  call i2e_sdsd(sh1,sh2,sh1,sh2,esdsd)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(esdsd)))
  if (schw(ab1).gt.shthresh) then
    nsigsd=nsigsd+1
  end if
  ! Here by permutation, we know that esdsd(i,j) must equal esdsd(j,i) so only store one of them
  do i=1,5
    do j=1,i
      if (abs(esdsd(i,j)).gt.intthresh) then
        call save_abcd_integral_d(ao1,ao2+i-1,ao1,ao2+j-1,esdsd(i,j))
      end if
    end do
  end do
end do

! Now do (sl| pairs
do ii=1,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)   ! First problem. sh2 can be smaller than sh1
  ao2=aosh(sh2)
  call i2e_slsl(sh1,sh2,sh1,sh2,eslsl)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(eslsl)))
  if (schw(ab1).gt.shthresh) then
    nsigsl=nsigsl+1
  end if
  do i=1,4
    ij=ind2(ao1,ao2+i-1)
    do j=1,i
      if (abs(eslsl(i,j)).gt.intthresh) then
        call save_abcd_integral_d(ao1,ao2+i-1,ao1,ao2+j-1,eslsl(i,j))
      end if
    end do
  end do
end do

! Now do (pp| pairs
do ii=1,npp
  sh1=ppprs(1,ii)
  ao1=aosh(sh1)
  sh2=ppprs(2,ii)
  ao2=aosh(sh2)
  call i2e_pppp(sh1,sh2,sh1,sh2,epppp)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(epppp)))
  if (schw(ab1).gt.shthresh) then
    nsigpp=nsigpp+1
  end if
  prevdupl=0
  do i=1,3
    do j=1,3
      if (sh1.eq.sh2) then
        lij=ind2(i,j)
      else
        lij=ind2(i,j+3)
      end if
      do k=1,i
        do l=1,3
          if (abs(epppp(i,j,k,l)).gt.intthresh) then
            if (sh1.eq.sh2) then
              lkl=ind2(k,l)
            else
              lkl=ind2(k,l+3)
            end if
            if (lij.gt.lkl) then
              lijkl=ioff(lij)+lkl
            else
              lijkl=ioff(lkl)+lij
            end if
            if (prevdupl(lijkl).ne.0) cycle
            prevdupl(lijkl)=1
            call save_abcd_integral_d(ao1+i-1,ao2+j-1,ao1+k-1,ao2+l-1,epppp(i,j,k,l))
          end if
        end do
      end do
    end do
  end do
end do

! Now do (pd| pairs
do ii=1,npd
  sh1=pdprs(1,ii)
  ao1=aosh(sh1)
  sh2=pdprs(2,ii)
  ao2=aosh(sh2)
  call i2e_pdpd(sh1,sh2,sh1,sh2,epdpd)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(epdpd)))
  if (schw(ab1).gt.shthresh) then
    nsigpd=nsigpd+1
  end if
  prevdupl=0
  do i=1,3
    do j=1,5
      lij=ind2(i,j+3)
      do k=1,3
        do l=1,5
          if (abs(epdpd(i,j,k,l)).gt.intthresh) then
            lkl=ind2(k,l+3)
            if (lij.gt.lkl) then
              lijkl=ioff(lij)+lkl
            else
              lijkl=ioff(lkl)+lij
            end if
            if (prevdupl(lijkl).ne.0) cycle
            prevdupl(lijkl)=1
            call save_abcd_integral_d(ao1+i-1,ao2+j-1,ao1+k-1,ao2+l-1,epdpd(i,j,k,l))
          end if
        end do
      end do
    end do
  end do
end do

! Now do (pl| pairs
do ii=1,npl
  sh1=plprs(1,ii)
  ao1=aosh(sh1)
  sh2=plprs(2,ii)   ! First problem. sh2 can be smaller than sh1
  ao2=aosh(sh2)
  call i2e_plpl(sh1,sh2,sh1,sh2,eplpl)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(eplpl)))
  if (schw(ab1).gt.shthresh) then
    nsigpl=nsigpl+1
  end if
  prevdupl=0
  do i=1,3
    do j=1,4
      lij=ind2(i,j+3)
      do k=1,i
        do l=1,4
          if (abs(eplpl(i,j,k,l)).gt.intthresh) then
            lkl=ind2(k,l+3)
            if (lij.gt.lkl) then
              lijkl=ioff(lij)+lkl
            else
              lijkl=ioff(lkl)+lij
            end if
            if (prevdupl(lijkl).ne.0) cycle
            prevdupl(lijkl)=1
            call save_abcd_integral_d(ao1+i-1,ao2+j-1,ao1+k-1,ao2+l-1,eplpl(i,j,k,l))
          end if
        end do
      end do
    end do
  end do
end do

! Now do (dd| pairs
do ii=1,ndd
  sh1=ddprs(1,ii)
  ao1=aosh(sh1)
  sh2=ddprs(2,ii)
  ao2=aosh(sh2)
  call i2e_dddd(sh1,sh2,sh1,sh2,edddd)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(edddd)))
  if (schw(ab1).gt.shthresh) then
    nsigdd=nsigdd+1
  end if
  prevdupl=0
  do i=1,5
    do j=1,5
      if (sh1.eq.sh2) then
        lij=ind2(i,j)
      else
        lij=ind2(i,j+5)
      end if
      do k=1,i
        do l=1,5
          if (abs(edddd(i,j,k,l)).gt.intthresh) then
            if (sh1.eq.sh2) then
              lkl=ind2(k,l)
            else
              lkl=ind2(k,l+5)
            end if
            if (lij.gt.lkl) then
              lijkl=ioff(lij)+lkl
            else
              lijkl=ioff(lkl)+lij
            end if
            if (prevdupl(lijkl).ne.0) cycle
            prevdupl(lijkl)=1
            call save_abcd_integral_d(ao1+i-1,ao2+j-1,ao1+k-1,ao2+l-1,edddd(i,j,k,l))
          end if
        end do
      end do
    end do
  end do
end do

! Now do (dl| pairs
do ii=1,ndl
  sh1=dlprs(1,ii)
  ao1=aosh(sh1)
  sh2=dlprs(2,ii)
  ao2=aosh(sh2)
  call i2e_dldl(sh1,sh2,sh1,sh2,edldl)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(edldl)))
  if (schw(ab1).gt.shthresh) then
    nsigdl=nsigdl+1
  end if
  prevdupl=0
  do i=1,5
    do j=1,4
      lij=ind2(i,j+5)
      do k=1,5
        do l=1,4
          if (abs(edldl(i,j,k,l)).gt.intthresh) then
            lkl=ind2(k,l+5)
            if (lij.gt.lkl) then
              lijkl=ioff(lij)+lkl
            else
              lijkl=ioff(lkl)+lij
            end if
            if (prevdupl(lijkl).ne.0) cycle
            prevdupl(lijkl)=1
            call save_abcd_integral_d(ao1+i-1,ao2+j-1,ao1+k-1,ao2+l-1,edldl(i,j,k,l))
          end if
        end do
      end do
    end do
  end do
end do

! Now do (ll'| pairs
do ii=1,nll
  sh1=llprs(1,ii)
  ao1=aosh(sh1)
  sh2=llprs(2,ii)
  ao2=aosh(sh2)
  call i2e_llll(sh1,sh2,sh1,sh2,ellll)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(ellll)))
  if (schw(ab1).gt.shthresh) then
    nsigll=nsigll+1
  end if
  prevdupl=0
  do i=1,4
    do j=1,4
      if (sh1.eq.sh2) then
        lij=ind2(i,j)
      else
        lij=ind2(i,j+4)
      end if
      do k=1,i
        do l=1,4
          if (abs(ellll(i,j,k,l)).gt.intthresh) then
            if (sh1.eq.sh2) then
              lkl=ind2(k,l)
            else
              lkl=ind2(k,l+4)
            end if
            if(lij.gt.lkl) then
              lijkl=ioff(lij)+lkl
            else
              lijkl=ioff(lkl)+lij
            end if
            if (prevdupl(lijkl).ne.0) cycle
            prevdupl(lijkl)=1
            call save_abcd_integral_d(ao1+i-1,ao2+j-1,ao1+k-1,ao2+l-1,ellll(i,j,k,l))
          end if
        end do
      end do
    end do
  end do
end do

call time_checker(-1,"Two electron integrals - diagonal terms, timing:")
write (9,*) "Significant / Total  number of each type of shell pairs:"
write (9,'(A,I6,"/",I6,A,I6,"/",I6)') "(ss|",nsigss,nss,"     (sp|",nsigsp,nsp
write (9,'(A,I6,"/",I6,A,I6,"/",I6)') "(sd|",nsigsd,nsd,"     (sl|",nsigsl,nsl
write (9,'(A,I6,"/",I6,A,I6,"/",I6)') "(pp|",nsigpp,npp,"     (pd|",nsigpd,npd
write (9,'(A,I6,"/",I6,A,I6,"/",I6)') "(pl|",nsigpl,npl,"     (dd|",nsigdd,ndd
write (9,'(A,I6,"/",I6,A,I6,"/",I6)') "(dl|",nsigdl,ndl,"     (ll|",nsigll,nll
flush(9)


! ********************************************
! Now loop over off-diagonal shell pair pairs.
! First work out how many fully permutationally relevant i2es there may be
! ********************************************
maxabcd=nsigss*(nsigss-1)/2+nsigpp*(nsigpp-1)/2*81+nsigdd*(nsigdd-1)/2*625
maxabcd=maxabcd+nsigll*(nsigll-1)/2*256
maxabcd=maxabcd+nsigsp*(nsigsp-1)/2*9+nsigsd*(nsigsd-1)/2*25+nsigsl*(nsigsl-1)/2*16
maxabcd=maxabcd+nsigpd*(nsigpd-1)/2*225+nsigpl*(nsigpl-1)/2*144+nsigdl*(nsigdl-1)/2*400

maxabcd=maxabcd+nsigss*nsigsp*3+nsigss*nsigsd*5+nsigss*nsigsl*4+nsigss*nsigpp*9
maxabcd=maxabcd+nsigss*nsigpd*15+nsigss*nsigpl*12+nsigss*nsigdd*25+nsigss*nsigdl*20
maxabcd=maxabcd+nsigss*nsigll*16+nsigsp*nsigsd*15+nsigsp*nsigsl*12+nsigsp*nsigpp*27
maxabcd=maxabcd+nsigsp*nsigpd*45+nsigsp*nsigpl*36+nsigsp*nsigdd*75+nsigsp*nsigdl*60
maxabcd=maxabcd+nsigsp*nsigll*48+nsigsd*nsigsl*20+nsigsd*nsigpp*45+nsigsd*nsigpd*75
maxabcd=maxabcd+nsigsd*nsigpl*60+nsigsd*nsigdd*125+nsigsd*nsigdl*100+nsigsd*nsigll*80
maxabcd=maxabcd+nsigsl*nsigpp*36+nsigsl*nsigpd*60+nsigsl*nsigpl*48+nsigsl*nsigdd*100
maxabcd=maxabcd+nsigsl*nsigdl*80+nsigsl*nsigll*64+nsigpp*nsigpd*135+nsigpp*nsigpl*108
maxabcd=maxabcd+nsigpp*nsigdd*225+nsigpp*nsigdl*180+nsigpp*nsigll*144+nsigpd*nsigpl*180
maxabcd=maxabcd+nsigpd*nsigdd*375+nsigpd*nsigdl*300+nsigpd*nsigll*240+nsigdd*nsigdl*500
maxabcd=maxabcd+nsigdd*nsigll*400+nsigdl*nsigll*320



write (9,*) "Based on Schwartz screen values for each shell, the number of integrals has been estimated:"
write (9,'(I12,A,I12)') maxabcd,"  vs. a theoretical maximum of ",nb4
allocate(i2eabcd(maxabcd),aonabcd(2,maxabcd))
i2eabcd=0._dp
aonabcd=0
nabcd=0
write (9,*)
write (9,*) "Type of shell quadruples, Total number, n. based on sig, number screened"

! Now calculate them.
! Start with (ss|ss)
nschw=0
do ii=2,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=ssprs(1,jj)
    ao3=aosh(sh3)
    sh4=ssprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ssss(sh1,sh2,sh3,sh4,essss)
    if (abs(essss).gt.intthresh) then
      call save_abcd_integral(ao1,ao2,ao3,ao4,essss)
    end if
  end do
end do
if (nss.gt.0) then
  write (9,'(A,3I12)') "(ss|ss)",nss*(nss-1)/2,nsigss*(nsigss-1)/2,nschw
end if

nschw=0
! Then (ss|sp)
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nsp
    sh3=spprs(1,jj)
    ao3=aosh(sh3)
    sh4=spprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sssp(sh1,sh2,sh3,sh4,esssp)
    do i=1,3
      if (abs(esssp(i)).gt.intthresh) then
        call save_abcd_integral(ao1,ao2,ao3,ao4+i-1,esssp(i))
      end if
    end do
  end do
end do
if (nss*nsp.gt.0) then
  write (9,'(A,3I12)') "(ss|sp)",nss*nsp,nsigss*nsigsp,nschw
end if

nschw=0
! Then (ss|sd)
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nsd
    sh3=sdprs(1,jj)
    ao3=aosh(sh3)
    sh4=sdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sssd(sh1,sh2,sh3,sh4,esssd)
    do i=1,5
      if (abs(esssd(i)).gt.intthresh) then
        call save_abcd_integral(ao1,ao2,ao3,ao4+i-1,esssd(i))
      end if
    end do
  end do
end do
if(nss*nsd.gt.0) then
  write (9,'(A,3I12)') "(ss|sd)",nss*nsd,nsigss*nsigsd,nschw
end if

! Then (ss|sl)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nsl
    sh3=slprs(1,jj)
    ao3=aosh(sh3)
    sh4=slprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sssl(sh1,sh2,sh3,sh4,esssl)
    do i=1,4
      if (abs(esssl(i)).gt.intthresh) then
        call save_abcd_integral(ao1,ao2,ao3,ao4+i-1,esssl(i))
      end if
    end do
  end do
end do
if(nss*nsl.gt.0) then
  write (9,'(A,3I12)') "(ss|sl)",nss*nsl,nsigss*nsigsl,nschw
end if

! Then (ss|pp), which can be (ab|cd) or (aa|bc) or (ab|cc) or (aa|cc)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npp
    sh3=ppprs(1,jj)
    ao3=aosh(sh3)
    sh4=ppprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sspp(sh1,sh2,sh3,sh4,esspp)
    do i=1,3
      do j=1,3
        if ((sh3.eq.sh4).and.(j.gt.i)) cycle
        if (abs(esspp(i,j)).gt.intthresh) then
          call save_abcd_integral(ao1,ao2,ao3+i-1,ao4+j-1,esspp(i,j))
        end if
      end do
    end do
  end do
end do
if(nss*npp.gt.0) then
  write (9,'(A,3I12)') "(ss|pp)",nss*npp,nsigss*nsigpp,nschw
end if

! Then (ss|pd)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npd
    sh3=pdprs(1,jj)
    ao3=aosh(sh3)
    sh4=pdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sspd(sh1,sh2,sh3,sh4,esspd)
    do i=1,3
      do j=1,5
        if (abs(esspd(i,j)).gt.intthresh) then
          call save_abcd_integral(ao1,ao2,ao3+i-1,ao4+j-1,esspd(i,j))
        end if
      end do
    end do
  end do
end do
if(nss*npd.gt.0) then
  write (9,'(A,3I12)') "(ss|pd)",nss*npd,nsigss*nsigpd,nschw
end if

! Then (ss|pl)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npl
    sh3=plprs(1,jj)
    ao3=aosh(sh3)
    sh4=plprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sspl(sh1,sh2,sh3,sh4,esspl)
    do i=1,3
      do j=1,4
        if (abs(esspl(i,j)).gt.intthresh) then
          call save_abcd_integral(ao1,ao2,ao3+i-1,ao4+j-1,esspl(i,j))
        end if
      end do
    end do
  end do
end do
if(nss*npl.gt.0) then
  write (9,'(A,3I12)') "(ss|pl)",nss*npl,nsigss*nsigpl,nschw
end if

! Then (ss|dd)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndd
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ssdd(sh1,sh2,sh3,sh4,essdd)
    if (sh3.eq.sh4) then
      do i=1,5
        do j=1,i
          if (abs(essdd(i,j)).gt.intthresh) then
            call save_abcd_integral(ao1,ao2,ao3+i-1,ao4+j-1,essdd(i,j))
          end if
        end do
      end do
    else
      do i=1,5
        do j=1,5
          if (abs(essdd(i,j)).gt.intthresh) then
            call save_abcd_integral(ao1,ao2,ao3+i-1,ao4+j-1,essdd(i,j))
          end if
        end do
      end do
    end if
  end do
end do
if (nss*ndd.gt.0) then
  write (9,'(A,3I12)') "(ss|dd)",nss*ndd,nsigss*nsigdd,nschw
end if

! Then (ss|dl)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndl
    sh3=dlprs(1,jj)
    ao3=aosh(sh3)
    sh4=dlprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ssdl(sh1,sh2,sh3,sh4,essdl)
    do i=1,5
      do j=1,4
        if (abs(essdl(i,j)).gt.intthresh) then
          call save_abcd_integral(ao1,ao2,ao3+i-1,ao4+j-1,essdl(i,j))
        end if
      end do
    end do
  end do
end do
if (nss*ndl.gt.0) then
  write (9,'(A,3I12)') "(ss|dl)",nss*ndl,nsigss*nsigdl,nschw
end if

! Then (ss|ll)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ssll(sh1,sh2,sh3,sh4,essll)
    if (sh3.eq.sh4) then
      do i=1,4
        do j=1,i
          if (abs(essll(i,j)).gt.intthresh) then
            call save_abcd_integral(ao1,ao2,ao3+i-1,ao4+j-1,essll(i,j))
          end if
        end do
      end do
    else
      do i=1,4
        do j=1,4
          if (abs(essll(i,j)).gt.intthresh) then
            call save_abcd_integral(ao1,ao2,ao3+i-1,ao4+j-1,essll(i,j))
          end if
        end do
      end do
    end if
  end do
end do
if (nss*nll.gt.0) then
  write (9,'(A,3I12)') "(ss|ll)",nss*nll,nsigss*nsigll,nschw
end if

! Then (sp|sp)
nschw=0
do ii=2,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=spprs(1,jj)
    ao3=aosh(sh3)
    sh4=spprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    atc=atsh(sh3) ; atd=atsh(sh4)
    call i2e_spsp(sh1,sh2,sh3,sh4,espsp)
    do i=1,3
      do j=1,3
        if (abs(espsp(i,j)).gt.intthresh) then
          nabcd=nabcd+1
          i2eabcd(nabcd)=espsp(i,j)
          aonabcd(:,nabcd)=(/ind2(ao1,ao2+i-1),ind2(ao3,ao4+j-1)/)
        end if
      end do
    end do
  end do
end do
if (nsp.gt.0) then
  write (9,'(A,3I12)') "(sp|sp)",nsp*(nsp-1)/2,nsigsp*(nsigsp-1)/2,nschw
end if

! Then (sp|sd)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nsd
    sh3=sdprs(1,jj)
    ao3=aosh(sh3)
    sh4=sdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_spsd(sh1,sh2,sh3,sh4,espsd)
    do i=1,3
      do j=1,5
        if (abs(espsd(i,j)).gt.intthresh) then
          call save_abcd_integral(ao1,ao2+i-1,ao3,ao4+j-1,espsd(i,j))
        end if
      end do
    end do
  end do
end do
if (nsp*nsd.gt.0) then
  write (9,'(A,3I12)') "(sp|sd)",nsp*nsd,nsigsp*nsigsd,nschw
end if

! Then (sp|sl)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nsl
    sh3=slprs(1,jj)
    ao3=aosh(sh3)
    sh4=slprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_spsl(sh1,sh2,sh3,sh4,espsl)
    do i=1,3
      do j=1,4
        if (abs(espsl(i,j)).gt.intthresh) then
          call save_abcd_integral(ao1,ao2+i-1,ao3,ao4+j-1,espsl(i,j))
        end if
      end do
    end do
  end do
end do
if (nsp*nsl.gt.0) then
  write (9,'(A,3I12)') "(sp|sl)",nsp*nsl,nsigsp*nsigsl,nschw
end if

! Then (sp|pp)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npp
    sh3=ppprs(1,jj)
    ao3=aosh(sh3)
    sh4=ppprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sppp(sh1,sh2,sh3,sh4,esppp)
    do l=1,3
      do k=1,3
        if ((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,3
          if (abs(esppp(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,esppp(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsp*npp.gt.0) then
  write (9,'(A,3I12)') "(sp|pp)",nsp*npp,nsigsp*nsigpp,nschw
end if

! Then (sp|pd)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npd
    sh3=pdprs(1,jj)
    ao3=aosh(sh3)
    sh4=pdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sppd(sh1,sh2,sh3,sh4,esppd)
    do l=1,5
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,3
          if (abs(esppd(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,esppd(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsp*npd.gt.0) then
  write (9,'(A,3I12)') "(sp|pd)",nsp*npd,nsigsp*nsigpd,nschw
end if

! Then (sp|pl)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npl
    sh3=plprs(1,jj)
    ao3=aosh(sh3)
    sh4=plprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sppl(sh1,sh2,sh3,sh4,esppl)
    do l=1,4
      do k=1,3
        do j=1,3
          if (abs(esppl(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,esppl(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsp*npl.gt.0) then
  write (9,'(A,3I12)') "(sp|pl)",nsp*npl,nsigsp*nsigpl,nschw
end if

! Then (sp|dd)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndd
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_spdd(sh1,sh2,sh3,sh4,espdd)
    do l=1,5
      do k=1,5
        if ((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,3
          if (abs(espdd(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,espdd(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsp*ndd.gt.0) then
  write (9,'(A,3I12)') "(sp|dd)",nsp*ndd,nsigsp*nsigdd,nschw
end if

! Then (sp|dl)
!nschw=0
!do ii=1,nsp
!  sh1=spprs(1,ii)
!  ao1=aosh(sh1)
!  sh2=spprs(2,ii)
!  ao2=aosh(sh2)
!  ab1=ind2(sh1,sh2)
!  do jj=1,ndl
!    sh3=dlprs(1,jj)
!    ao3=aosh(sh3)
!    sh4=dlprs(2,jj)
!    ao4=aosh(sh4)
!    ab2=ind2(sh3,sh4)
!    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
!      nschw=nschw+1
!      cycle
!    end if
!    call i2e_spdl(sh1,sh2,sh3,sh4,espdl)
!    do l=1,4
!      do k=1,5
!        kl=ind2(ao3+k-1,ao4+l-1)
!        do j=1,3
!          if (abs(espdl(j,k,l)).gt.intthresh) then 
!            ij=ind2(ao1,ao2+j-1)
!            ijkl=ind2(ij,kl)
!            i2e(ijkl)=espdl(j,k,l)
!          end if
!        end do
!      end do
!    end do
!  end do
!end do
!write (9,'(A,I8)') "Number of screened (sp|dl) integrals:",nschw

! Then (sp|ll)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_spll(sh1,sh2,sh3,sh4,espll)
    do l=1,4
      do k=1,4
        if ((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,3
          if (abs(espll(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,espll(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsp*nll.gt.0) then
  write (9,'(A,3I12)') "(sp|ll)",nsp*nll,nsigsp*nsigll,nschw
end if

! Then (sd|sd)
nschw=0
do ii=2,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=sdprs(1,jj)
    ao3=aosh(sh3)
    sh4=sdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sdsd(sh1,sh2,sh3,sh4,esdsd)
    do l=1,5
      do j=1,5
        if (abs(esdsd(j,l)).gt.intthresh) then 
          call save_abcd_integral(ao1,ao2+j-1,ao3,ao4+l-1,esdsd(j,l))
        end if
      end do
    end do
  end do
end do
if (nsd.gt.0) then
  write (9,'(A,3I12)') "(sd|sd)",nsd*(nsd-1)/2,nsigsd*(nsigsd-1)/2,nschw
end if

! Then (sd|sl)
nschw=0
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nsl
    sh3=slprs(1,jj)
    ao3=aosh(sh3)
    sh4=slprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sdsl(sh1,sh2,sh3,sh4,esdsl)
    do l=1,4
      do j=1,5
        if (abs(esdsl(j,l)).gt.intthresh) then 
          call save_abcd_integral(ao1,ao2+j-1,ao3,ao4+l-1,esdsl(j,l))
        end if
      end do
    end do
  end do
end do
if (nsd*nsl.gt.0) then
  write (9,'(A,3I12)') "(sd|sl)",nsd*nsl,nsigsd*nsigsl,nschw
end if
!call time_checker(-1,"Two electron integrals - (sd|sl) terms, timing:")

! Then (sd|pp)
nschw=0
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npp
    sh3=ppprs(1,jj)
    ao3=aosh(sh3)
    sh4=ppprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sdpp(sh1,sh2,sh3,sh4,esdpp)
    do l=1,3
      do k=1,3
        if ((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,5
          if (abs(esdpp(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,esdpp(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsd*npp.gt.0) then
  write (9,'(A,3I12)') "(sd|pp)",nsd*npp,nsigsd*nsigpp,nschw
end if

! Then (sd|pd)
nschw=0
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npd
    sh3=pdprs(1,jj)
    ao3=aosh(sh3)
    sh4=pdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sdpd(sh1,sh2,sh3,sh4,esdpd)
    do l=1,5
      do k=1,3
        do j=1,5
          if (abs(esdpd(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,esdpd(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsd*npd.gt.0) then
  write (9,'(A,3I12)') "(sd|pd)",nsd*npd,nsigsd*nsigpd,nschw
end if

! Then (sd|pl)
!nschw=0
!do ii=1,nsd
!  sh1=sdprs(1,ii)
!  ao1=aosh(sh1)
!  sh2=sdprs(2,ii)
!  ao2=aosh(sh2)
!  ab1=ind2(sh1,sh2)
!  do jj=1,npl
!    sh3=plprs(1,jj)
!    ao3=aosh(sh3)
!    sh4=plprs(2,jj)
!    ao4=aosh(sh4)
!    ab2=ind2(sh3,sh4)
!    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
!      nschw=nschw+1
!      cycle
!    end if
!    call i2e_sdpl(sh1,sh2,sh3,sh4,esdpl)
!    do l=1,4
!      do k=1,3
!        kl=ind2(ao3+k-1,ao4+l-1)
!        do j=1,5
!          if (abs(esdpl(j,k,l)).gt.intthresh) then 
!            ij=ind2(ao1,ao2+j-1)
!            ijkl=ind2(ij,kl)
!            i2e(ijkl)=esdpl(j,k,l)
!          end if
!        end do
!      end do
!    end do
!  end do
!end do
!write (9,'(A,I8)') "Number of screened (sd|pl) integrals:",nschw

! Then (sd|dd)
nschw=0
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndd
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sddd(sh1,sh2,sh3,sh4,esddd)
    do l=1,5
      do k=1,5
        if ((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,5
          if (abs(esddd(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,esddd(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsd*ndd.gt.0) then
  write (9,'(A,3I12)') "(sd|dd)",nsd*ndd,nsigsd*nsigdd,nschw
end if

! Then (sd|dl)
nschw=0
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndl
    sh3=dlprs(1,jj)
    ao3=aosh(sh3)
    sh4=dlprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sddl(sh1,sh2,sh3,sh4,esddl)
    do l=1,4
      do k=1,5
        do j=1,5
          if (abs(esddl(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,esddl(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsd*ndl.gt.0) then
  write (9,'(A,3I12)') "(sd|dl)",nsd*ndl,nsigsd*nsigdl,nschw
end if

! Then (sd|ll)
nschw=0
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sdll(sh1,sh2,sh3,sh4,esdll)
    do l=1,4
      do k=1,4
        if ((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,5
          if (abs(esdll(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,esdll(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsd*nll.gt.0) then
  write (9,'(A,3I12)') "(sd|ll)",nsd*nll,nsigsd*nsigll,nschw
end if

! Then (sl|sl)
nschw=0
n1=0;n2=0;n3=0;n4=0
do ii=2,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=slprs(1,jj)
    ao3=aosh(sh3)
    sh4=slprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_slsl(sh1,sh2,sh3,sh4,eslsl)
    do l=1,4
      do j=1,4
        if (abs(eslsl(j,l)).gt.intthresh) then 
          call save_abcd_integral(ao1,ao2+j-1,ao3,ao4+l-1,eslsl(j,l))
        end if
      end do
    end do
  end do
end do
if (nsl.gt.0) then
  write (9,'(A,3I12)') "(sl|sl)",nsl*(nsl-1)/2,nsigsl*(nsigsl-1)/2,nschw
end if

! Then (sl|pp)
nschw=0
do ii=1,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npp
    sh3=ppprs(1,jj)
    ao3=aosh(sh3)
    sh4=ppprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_slpp(sh1,sh2,sh3,sh4,eslpp)
    do l=1,3
      do k=1,3
        if ((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,4
          if (abs(eslpp(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,eslpp(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsl*npp.gt.0) then
  write (9,'(A,3I12)') "(sl|pp)",nsl*npp,nsigsl*nsigpp,nschw
end if

! Then (sl|pd)
!nschw=0
!do ii=1,nsl
!  sh1=slprs(1,ii)
!  ao1=aosh(sh1)
!  sh2=slprs(2,ii)
!  ao2=aosh(sh2)
!  ab1=ind2(sh1,sh2)
!  do jj=1,npd
!    sh3=pdprs(1,jj)
!    ao3=aosh(sh3)
!    sh4=pdprs(2,jj)
!    ao4=aosh(sh4)
!    ab2=ind2(sh3,sh4)
!    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
!      nschw=nschw+1
!      cycle
!    end if
!    call i2e_slpd(sh1,sh2,sh3,sh4,eslpd)
!    do l=1,5
!      do k=1,3
!        kl=ind2(ao3+k-1,ao4+l-1)
!        do j=1,4
!          if (abs(eslpd(j,k,l)).gt.intthresh) then 
!            ij=ind2(ao1,ao2+j-1)
!            ijkl=ind2(ij,kl)
!            i2e(ijkl)=eslpd(j,k,l)
!          end if
!        end do
!      end do
!    end do
!  end do
!end do

! Then (sl|pl)
nschw=0
do ii=1,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npl
    sh3=plprs(1,jj)
    ao3=aosh(sh3)
    sh4=plprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_slpl(sh1,sh2,sh3,sh4,eslpl)
    do l=1,4
      do k=1,3
        do j=1,4
          if (abs(eslpl(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,eslpl(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsl*npl.gt.0) then
  write (9,'(A,3I12)') "(sl|pl)",nsl*npl,nsigsl*nsigpl,nschw
end if

! Then (sl|dd)
nschw=0
do ii=1,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndd
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sldd(sh1,sh2,sh3,sh4,esldd)
    do l=1,5
      do k=1,5
        if ((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,4
          if (abs(esldd(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,esldd(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsl*ndd.gt.0) then
  write (9,'(A,3I12)') "(sl|dd)",nsl*ndd,nsigsl*nsigdd,nschw
end if

! Then, (sl|dl)
nschw=0
do ii=1,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndl
    sh3=dlprs(1,jj)
    ao3=aosh(sh3)
    sh4=dlprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sldl(sh1,sh2,sh3,sh4,esldl)
    do l=1,4
      do k=1,5
        do j=1,4
          if (abs(esldl(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,esldl(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (ndl*nsl.gt.0) then
  write (9,'(A,3I12)') "(sl|dl)",nsl*ndl,nsigsl*nsigdl,nschw
end if

! Then (sl|ll)
nschw=0
do ii=1,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_slll(sh1,sh2,sh3,sh4,eslll)
    do l=1,4
      do k=1,4
        if ((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,4
          if (abs(eslll(j,k,l)).gt.intthresh) then 
            call save_abcd_integral(ao1,ao2+j-1,ao3+k-1,ao4+l-1,eslll(j,k,l))
          end if
        end do
      end do
    end do
  end do
end do
if (nsl*nll.gt.0) then
  write (9,'(A,3I12)') "(sl|ll)",nsl*nll,nsigsl*nsigll,nschw
end if

! Then (pp|pp)
nschw=0
do ii=2,npp
  sh1=ppprs(1,ii)
  ao1=aosh(sh1)
  sh2=ppprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=ppprs(1,jj)
    ao3=aosh(sh3)
    sh4=ppprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_pppp(sh1,sh2,sh3,sh4,epppp)
    do l=1,3
      do k=1,3
        if ((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,3
          do i=1,3
            if ((sh1.eq.sh2).and.(i.gt.j)) cycle
            if (abs(epppp(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,epppp(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (npp.gt.0) then
  write (9,'(A,3I12)') "(pp|pp)",npp*(npp-1)/2,nsigpp*(nsigpp-1)/2,nschw
end if

! Then (pp|pd)
nschw=0
do ii=1,npp
  sh1=ppprs(1,ii)
  ao1=aosh(sh1)
  sh2=ppprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npd  
    sh3=pdprs(1,jj)
    ao3=aosh(sh3)
    sh4=pdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_pppd(sh1,sh2,sh3,sh4,epppd)
    do l=1,5
      do k=1,3
        do j=1,3
          do i=1,3
            if ((sh1.eq.sh2).and.(j.gt.i)) cycle
            if (abs(epppd(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,epppd(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (npp*npd.gt.0) then
  write (9,'(A,3I12)') "(pp|pd)",npp*npd,nsigpp*nsigpd,nschw
end if

! Then (pp|pl)
nschw=0
do ii=1,npp
  sh1=ppprs(1,ii)
  ao1=aosh(sh1)
  sh2=ppprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npl  
    sh3=plprs(1,jj)
    ao3=aosh(sh3)
    sh4=plprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_pppl(sh1,sh2,sh3,sh4,epppl)
    do l=1,4
      do k=1,3
        do j=1,3
          do i=1,3
            if ((sh1.eq.sh2).and.(j.gt.i)) cycle
            if (abs(epppl(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,epppl(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (npp*npl.gt.0) then
write (9,'(A,I8,A,I8)') "Number of screened (pp|pl) integral sets:",nschw," out of:",npp*npl
end if

! Then (pp|dd)
nschw=0
do ii=1,npp
  sh1=ppprs(1,ii)
  ao1=aosh(sh1)
  sh2=ppprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndd  
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ppdd(sh1,sh2,sh3,sh4,eppdd)
    do l=1,5
      do k=1,5
        if ((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,3
          do i=1,3
            if ((sh1.eq.sh2).and.(i.gt.j)) cycle
            if (abs(eppdd(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,eppdd(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (npp*ndd.gt.0) then
  write (9,'(A,3I12)') "(pp|dd)",npp*ndd,nsigpp*nsigdd,nschw
end if

! Then (pp|dl) which i'll skip

! Then (pp|ll)
nschw=0
do ii=1,npp
  sh1=ppprs(1,ii)
  ao1=aosh(sh1)
  sh2=ppprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll  
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ppll(sh1,sh2,sh3,sh4,eppll)
    do l=1,4
      do k=1,4
        if((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,3
          do i=1,3
            if((sh1.eq.sh2).and.(i.gt.j)) cycle
            if (abs(eppll(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,eppll(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (npp*nll.gt.0) then
  write (9,'(A,3I12)') "(pp|ll)",npp*nll,nsigpp*nsigll,nschw
end if

! Then (pd|pd)
nschw=0
do ii=2,npd
  sh1=pdprs(1,ii)
  ao1=aosh(sh1)
  sh2=pdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=pdprs(1,jj)
    ao3=aosh(sh3)
    sh4=pdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_pdpd(sh1,sh2,sh3,sh4,epdpd)
    do l=1,5
      do k=1,3
        do j=1,5
          do i=1,3
            if (abs(epdpd(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,epdpd(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (npd.gt.0) then
  write (9,'(A,3I12)') "(pd|pd)",npd*(npd-1)/2,nsigpd*(nsigpd-1)/2,nschw
end if

! Gonna skip (pd|pl)

! Then (pd|dd)
nschw=0
do ii=1,npd
  sh1=pdprs(1,ii)
  ao1=aosh(sh1)
  sh2=pdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndd
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_pddd(sh1,sh2,sh3,sh4,epddd)
    do l=1,5
      do k=1,5
        if((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,5
          do i=1,3
            if (abs(epddd(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,epddd(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (npd*ndd.gt.0) then
  write (9,'(A,3I12)') "(pd|dd)",npd*ndd,nsigpd*nsigdd,nschw
end if


! I will skip (pd|dl) and (pd|ll)

! Then (pl|pl)
nschw=0
do ii=2,npl
  sh1=plprs(1,ii)
  ao1=aosh(sh1)
  sh2=plprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=plprs(1,jj)
    ao3=aosh(sh3)
    sh4=plprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_plpl(sh1,sh2,sh3,sh4,eplpl)
    do l=1,4
      do k=1,3
        do j=1,4
          do i=1,3
            if (abs(eplpl(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,eplpl(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (npl.gt.0) then
  write (9,'(A,3I12)') "(pl|pl)",npl*(npl-1)/2,nsigpl*(nsigpl-1)/2,nschw
end if

! Also skipping (pl|dd) and (pl|dl)

! Then (pl|ll)
nschw=0
do ii=1,npl
  sh1=plprs(1,ii)
  ao1=aosh(sh1)
  sh2=plprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll  
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_plll(sh1,sh2,sh3,sh4,eplll)
    do l=1,4
      do k=1,4
        if((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,4
          do i=1,3
            if (abs(eplll(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,eplll(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (npl*nll.gt.0) then
write (9,'(A,I8,A,I8)') "Number of screened (pl|ll) integral sets:",nschw," out of:",npl*nll
end if

! Then (dd|dd)
nschw=0
do ii=2,ndd
  sh1=ddprs(1,ii)
  ao1=aosh(sh1)
  sh2=ddprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_dddd(sh1,sh2,sh3,sh4,edddd)
    do l=1,5
      do k=1,5
        if((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,5
          do i=1,5
            if((sh1.eq.sh2).and.(i.gt.j)) cycle
            if (abs(edddd(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,edddd(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (ndd.gt.0) then
  write (9,'(A,3I12)') "(dd|dd)",ndd*(ndd-1)/2,nsigdd*(nsigdd-1)/2,nschw
end if
!call time_checker(-1,"Two electron integrals - (dd|dd) terms, timing:")

! Then (dd|dl)
nschw=0
do ii=1,ndd
  sh1=ddprs(1,ii)
  ao1=aosh(sh1)
  sh2=ddprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndl
    sh3=dlprs(1,jj)
    ao3=aosh(sh3)
    sh4=dlprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_dddl(sh1,sh2,sh3,sh4,edddl)
    do l=1,4
      do k=1,5
        do j=1,5
          do i=1,5
            if((sh1.eq.sh2).and.(j.gt.i)) cycle
            if (abs(edddl(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,edddl(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (ndl*ndd.gt.0) then
  write (9,'(A,3I12)') "(dl|dd)",ndl*ndd,nsigdl*nsigdd,nschw
end if
!call time_checker(-1,"Two electron integrals - (dd|dl) terms, timing:")

! Then (dd|ll)
nschw=0
do ii=1,ndd
  sh1=ddprs(1,ii)
  ao1=aosh(sh1)
  sh2=ddprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ddll(sh1,sh2,sh3,sh4,eddll)
    do l=1,4
      do k=1,4
        if ((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,5
          do i=1,5
            if((sh1.eq.sh2).and.(j.gt.i)) cycle
            if (abs(eddll(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,eddll(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (ndd*nll.gt.0) then
  write (9,'(A,3I12)') "(dd|ll)",ndd*nll,nsigdd*nsigll,nschw
end if
!call time_checker(-1,"Two electron integrals - (dd|ll) terms, timing:")

! Then (dl|dl)
nschw=0
do ii=2,ndl
  sh1=dlprs(1,ii)
  ao1=aosh(sh1)
  sh2=dlprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=dlprs(1,jj)
    ao3=aosh(sh3)
    sh4=dlprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_dldl(sh1,sh2,sh3,sh4,edldl)
    do l=1,4
      do k=1,5
        do j=1,4
          do i=1,5
            if (abs(edldl(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,edldl(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (ndl.gt.0) then
  write (9,'(A,3I12)') "(dl|dl)",ndl*(ndl-1)/2,nsigdl*(nsigdl-1)/2,nschw
end if
!call time_checker(-1,"Two electron integrals - (dl|dl) terms, timing:")

! Then (dl|ll)
nschw=0
do ii=1,ndl
  sh1=dlprs(1,ii)
  ao1=aosh(sh1)
  sh2=dlprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_dlll(sh1,sh2,sh3,sh4,edlll)
    do l=1,4
      do k=1,4
        if((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,4
          do i=1,5
            if (abs(edlll(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,edlll(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (ndl*nll.gt.0) then
  write (9,'(A,3I12)') "(dl|ll)",ndl*nll,nsigdl*nsigll,nschw
end if
!call time_checker(-1,"Two electron integrals - (dl|ll) terms, timing:")

! Then (ll|ll)
nschw=0
do ii=2,nll
  sh1=llprs(1,ii)
  ao1=aosh(sh1)
  sh2=llprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_llll(sh1,sh2,sh3,sh4,ellll)
    do l=1,4
      do k=1,4
        if((sh3.eq.sh4).and.(k.gt.l)) cycle
        do j=1,4
          do i=1,4
            if((sh1.eq.sh2).and.(j.gt.i)) cycle
            if (abs(ellll(i,j,k,l)).gt.intthresh) then
              call save_abcd_integral(ao1+i-1,ao2+j-1,ao3+k-1,ao4+l-1,ellll(i,j,k,l))
            end if
          end do
        end do
      end do
    end do
  end do
end do
if (nll.gt.0) then
  write (9,'(A,3I12)') "(ll|ll)",nll*(nll-1)/2,nsigll*(nsigll-1)/2,nschw
end if
!call time_checker(-1,"Two electron integrals - (ll|ll) terms, timing:")
write (9,*) "Final number of integrals used:"
write (9,'(A,2I12)') "(aa|aa)",naaaa,nb
write (9,'(A,2I12)') "(aa|bb)",naabb,maxaabb
write (9,'(A,2I12)') "(ab|ab)",nabab,maxabab
write (9,'(A,2I12)') "(aa|bc)",naabc,maxaabc
write (9,'(A,2I12,A)') "(ab|cd)",nabcd_d,maxabcd_d,"coming from the diag part"
write (9,'(A,2I12)') "(ab|cd)",nabcd,maxabcd


!write (902,*) "(aa|aa) class"
!do i=1,naaaa
!  ia=aonaaaa(i)
!  ij=ind2(ia,ia)
!  ijkl=ioff(ij)+ij
!  write (902,'(I6,3x,4I3,F12.7)') ijkl,ia,ia,ia,ia,i2eaaaa(i)
!  if (ni2e(ijkl).eq.1) then
!    write (902,*) "++++ present two times in new integrals"
!  end if
!  ni2e(ijkl)=1
!  if (oi2e(ijkl).ne.1) then
!    write (902,*) "#### missing in old integrals"
!  end if
!end do
!write (902,*) "(aa|bb) class"
!do i=1,naabb
!  ia=revind2(1,aonaabb(i)) ; ka=revind2(2,aonaabb(i))
!  ij=ind2(ia,ia) ; kl=ind2(ka,ka)
!  if (ij.gt.kl) then
!    ijkl=ioff(ij)+kl
!  else
!    ijkl=ioff(kl)+ij
!  end if
!  write (902,'(I6,3x,4I3,F12.7)') ijkl,ia,ia,ka,ka,i2eaabb(i)
!  if (ni2e(ijkl).eq.1) then
!    write (902,*) "++++ present two times in new integrals"
!  end if
!  ni2e(ijkl)=1
!  if (oi2e(ijkl).ne.1) then
!    write (902,*) "#### missing in old integrals"
!  end if
!end do
!write (902,*) "(aa|bc) class"
!do i=1,naabc
!  ia=aonaabc(1,i)
!  ka=revind2(1,aonaabc(2,i)) ; la=revind2(2,aonaabc(2,i))
!  ij=ind2(ia,ia) ; kl=ind2(ka,la)
!  if (ij.gt.kl) then
!    ijkl=ioff(ij)+kl
!  else
!    ijkl=ioff(kl)+ij
!  end if
!  write (902,'(I6,3x,4I3,F12.7)') ijkl,ia,ia,ka,la,i2eaabc(i)
!  if (ni2e(ijkl).eq.1) then
!    write (902,*) "++++ present two times in new integrals"
!  end if
!  ni2e(ijkl)=1
!  if (oi2e(ijkl).ne.1) then
!    write (902,*) "#### missing in old integrals"
!  end if
!end do
!write (902,*) "(ab|ab) class"
!do i=1,nabab
!  ia=revind2(1,aonabab(i)) ; ja=revind2(2,aonabab(i))
!  ij=ind2(ia,ja)
!  ijkl=ioff(ij)+ij
!  write (902,'(I6,3x,4I3,F12.7)') ijkl,ia,ja,ia,ja,i2eabab(i)
!  if (ni2e(ijkl).eq.1) then
!    write (902,*) "++++ present two times in new integrals"
!  end if
!  ni2e(ijkl)=1
!  if (oi2e(ijkl).ne.1) then
!    write (902,*) "#### missing in old integrals"
!  end if
!end do
!write (902,*) "(ab|cd) class"
!do i=1,nabcd1
!  ia=revind2(1,aonabcd1(1,i)) ; ja=revind2(2,aonabcd1(1,i))
!  ka=revind2(1,aonabcd1(2,i)) ; la=revind2(2,aonabcd1(2,i))
!  ij=ind2(ia,ja) ; kl=ind2(ka,la)
!  if (ij.gt.kl) then
!    ijkl=ioff(ij)+kl
!  else
!    ijkl=ioff(kl)+ij
!  end if
!  write (902,'(I6,3x,4I3,F12.7)') ijkl,ia,ja,ka,la,i2eabcd1(i)
!  if (ni2e(ijkl).eq.1) then
!    write (902,*) "++++ present two times in new integrals"
!  end if
!  ni2e(ijkl)=1
!  if (oi2e(ijkl).ne.1) then
!    write (902,*) "#### missing in old integrals"
!  end if
!end do
!if (nabcdblocks.gt.1) then
!write (902,*) "(ab|cd) class"
!do i=1,nabcd2
!  ia=revind2(1,aonabcd2(1,i)) ; ja=revind2(2,aonabcd2(1,i))
!  ka=revind2(1,aonabcd2(2,i)) ; la=revind2(2,aonabcd2(2,i))
!  ij=ind2(ia,ja) ; kl=ind2(ka,la)
!  if (ij.gt.kl) then
!    ijkl=ioff(ij)+kl
!  else
!    ijkl=ioff(kl)+ij
!  end if
!  write (902,'(I6,3x,4I3,F12.7)') ijkl,ia,ja,ka,la,i2eabcd2(i)
!  if (ni2e(ijkl).eq.1) then
!    write (902,*) "++++ present two times in new integrals"
!  end if
!  ni2e(ijkl)=1
!  if (oi2e(ijkl).ne.1) then
!    write (902,*) "#### missing in old integrals"
!  end if
!end do
!end if
!write (902,*) "Total number non-zero:"
!write (902,*) naaaa+naabb+nabab+naabc+nabcd1+nabcd2+nabcd3


end subroutine cmpt2e

