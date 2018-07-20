Module s_sp_l_terms
use nrtype
implicit none

contains
pure subroutine i2e_ssss(sh1,sh2,sh3,sh4,essss)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: essss

integer(i4b) :: ab, cd, zzab, zzcd
real(dp) :: ppq, qpq, rpq2
real(dp) :: bfac

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
essss=0._dp
do zzab=1,nksqu(ab)
    if (negab(zzab,ab)) cycle
    do zzcd=1,nksqu(cd)
        if (negab(zzcd,cd)) cycle
        ppq=pab(zzab,ab)+pab(zzcd,cd)
        qpq=pab(zzab,ab)*pab(zzcd,cd)/ppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        call boys_0(qpq*rpq2,bfac)
        ! Build (ss|ss)^m, m=0,1,2
        essss=essss+kab(zzab,ab)*kab(zzcd,cd)*bfac/sqrt(ppq)
    end do
end do
return

end subroutine i2e_ssss


pure subroutine i2e_sssp(sh1,sh2,sh3,sh4,esssp)
use nrtype ; use molprops ; use boys
implicit none

! This version does simple O-S recursion to get these (ss|sp) integrals

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esssp(3)

integer(i4b) :: ab, cd, zzab, zzcd
real(dp) :: ooppq, qpq, qoppq, rpq2
real(dp) :: bfac(0:1), os(0:1), rdq(3), rqw(3), wx(3)

esssp=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzab=1,nksqu(ab)
    if (negab(zzab,ab)) cycle
    do zzcd=1,nksqu(cd)
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)/(pab(zzab,ab)+pab(zzcd,cd))
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        call boys_1(qpq*rpq2,bfac)
        ! Build (ss|ss)^m, m=0,1
        os(:)=kab(zzab,ab)*kab(zzcd,cd)*bfac(:)*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)
        esssp(:)=esssp(:)+rdq(:)*os(0)+rqw(:)*os(1)
    end do
end do

return

end subroutine i2e_sssp


pure subroutine i2e_sspp(sh1,sh2,sh3,sh4,esspp)
use nrtype ; use molprops ; use boys
implicit none

! This version exploits the horizontal relation, eq. 18 of HG & P 1988
! Need to compute the uncontr. and contr. (ss|ss), 3 (ss|sp) and 6 (ss|sd)
! integrals

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esspp(3,3)

integer(i4b) :: ab, cd, zzab, zzcd, m, l, lp1, lm1, jj, k
real(dp) :: ooppq, qpq, poppq, qoppq, rpq2, fac1
real(dp) :: bfac(0:2), os(10,0:2), osc(10), rdq(3), rqw(3), rcd(3), wx(3)

! First assemble the [ss|sx] (x=s,p,d) integrals and contract them

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzab=1,nksqu(ab)
    if (negab(zzab,ab)) cycle
    do zzcd=1,nksqu(cd)
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        call boys_2(qpq*rpq2,bfac)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2
        os(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac(:)*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)
        do m=0,1
          os(2:4,m)=os(1,m)*rdq(:)+rqw(:)*os(1,m+1)
        end do
        fac1=(os(1,0)-poppq*os(1,1))*hpm1(zzcd,cd)
        do lp1=5,10
          do jj=1,3
            l=minus1s(jj,lp1)
            if (l.eq.0) cycle
            os(lp1,0)=os(l,0)*rdq(jj)+os(l,1)*rqw(jj)
            lm1=minus1s(jj,l)
            if (lm1.eq.0) cycle
            os(lp1,0)=os(lp1,0)+fac1
          end do
        end do
        osc=osc+os(:,0)
    end do
end do

rcd=xsh(:,sh4)-xsh(:,sh3)
do lp1=5,10
  do k=2,4
    l=minus1s(k-1,lp1)
    if (l.eq.0) cycle
    esspp(k-1,l-1)=osc(lp1)+rcd(k-1)*osc(l)
  end do
end do

end subroutine i2e_sspp



pure subroutine i2e_spsp(sh1,sh2,sh3,sh4,espsp)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: espsp(3,3)

integer(i4b) :: ab, cd, zzab, zzcd, m, j
real(dp) :: ooppq, qpq, qoppq, rpq2, fac1
real(dp) :: bfac(0:2), os(0:3,0:3,0:2), rbp(3), rpw(3), rdq(3), rqw(3), wx(3)

! Simply assemble the [sx|sy] (x,y=s,p) integrals and contract them

espsp=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzab=1,nksqu(ab)
      if (negab(zzab,ab)) cycle
      rbp=Px(:,zzab,ab)-xsh(:,sh2)
      do zzcd=1,nksqu(cd)
            if (negab(zzcd,cd)) cycle
            ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
            qoppq=pab(zzcd,cd)*ooppq
            qpq=pab(zzab,ab)*qoppq
            rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
            Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
            rdq=Px(:,zzcd,cd)-xsh(:,sh4)
            rpw=Wx-Px(:,zzab,ab)
            rqw=Wx-Px(:,zzcd,cd)
            ! Build (ss|ss)^m, m=0,1,2
            call boys_2(qpq*rpq2,bfac)
            os(0,0,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
            ! OS recursion: Build (ss|sp)^m, m=0,1
            do m=0,1
              os(0,1:3,m)=rdq(:)*os(0,0,m)+rqw(:)*os(0,0,m+1)
            end do
            ! OS recursion: Build (sp|sp)
            fac1=.5_dp*os(0,0,1)*ooppq
            do j=1,3
              os(j,:,0)=rbp(j)*os(0,:,0)+rpw(j)*os(0,:,1)
              os(j,j,0)=os(j,j,0)+fac1
            end do
            espsp=espsp+os(1:3,1:3,0)
      end do
end do

end subroutine i2e_spsp



pure subroutine i2e_sppp(sh1,sh2,sh3,sh4,esppp)
use nrtype ; use molprops ; use boys
implicit none

! This version exploits the horizontal relation, eq. 18 of HG & P 1988
! Need to compute the uncontr. [ss|ss], [ss|sp], [sp|sp] and [sp|sd]
! integrals

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esppp(3,3,3)

integer(i4b) :: ab, cd, zzab, zzcd, m, l, lm1, lp1, jj
integer(i4b) :: kp1, jrec
real(dp) :: ooppq, qpq, poppq, qoppq, rpq2
real(dp) :: bfac(0:3), os(4,10,0:3), osc(4,10), rpw(3), rbp(3), rdq(3), rqw(3), rcd(3)
real(dp) :: wx(3), afac(0:1)

! First assemble the [sx|sy] (x=s,p; y=s,p,d) integrals and contract them
! We use the matrix os(j,l1,l2,m) depending on the one angular momentum at b, the two at d, and the Boys order

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzab=1,nksqu(ab)
    if (negab(zzab,ab)) cycle
    rbp=Px(:,zzab,ab)-xsh(:,sh2)
    do zzcd=1,nksqu(cd)
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rbp=Px(:,zzab,ab)-xsh(:,sh2)
        rpw=Wx-Px(:,zzab,ab)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3
        call boys_n(3,qpq*rpq2,bfac)
        os(1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,3
        do m=0,2
          os(1,2:4,m)=rdq(:)*os(1,1,m)+rqw(:)*os(1,1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0,2
        afac(0:1)=(os(1,1,0:1)-poppq*os(1,1,1:2))*hpm1(zzcd,cd)
        do lp1=5,10
          do jj=1,3
            l=minus1s(jj,lp1)
            if(l.eq.0) cycle
            jrec=jj
            exit
          end do
          os(1,lp1,0:1)=os(1,l,0:1)*rdq(jrec)+os(1,l,1:2)*rqw(jrec)
          lm1=minus1s(jrec,l)
          if (lm1.eq.0) cycle
          os(1,lp1,0:1)=os(1,lp1,0:1)+afac(0:1)
        end do
        do jj=1,3
          kp1=jj+1
          os(kp1,2:10,0)=rbp(jj)*os(1,2:10,0)+rpw(jj)*os(1,2:10,1)
          do l=2,10
            lm1=minus1s(jj,l)
            if(lm1.eq.0) cycle
            os(kp1,l,0)=os(kp1,l,0)+.5_dp*ooppq*os(1,lm1,1)*real(indexps(jj,l),dp)
            !os(kp1,l,0)=os(kp1,l,0)+.5_dp*ooppq*os(1,lm1,1)
          end do
        end do
        osc=osc+os(:,:,0)
    end do
end do

rcd=xsh(:,sh4)-xsh(:,sh3)
do jj=1,3
  do l=5,10
    lm1=minus1s(jj,l)
    if (lm1.eq.0) cycle
    esppp(:,jj,lm1-1)=osc(2:4,l)+rcd(jj)*osc(2:4,lm1)
  end do
end do

end subroutine i2e_sppp



pure subroutine i2e_pppp(sh1,sh2,sh3,sh4,epppp)
use nrtype ; use molprops ; use boys
implicit none

! This version exploits the electron transfer eq. of Giese & the HRR,  eq. 18 of HG & P 1988
! Start with [ss|ss] up to [gs|ss] using the VRR, then use ETR to [sd|sd] finally use HRR

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: epppp(3,3,3,3)

integer(i4b) :: ab, cd, zb, zd, zzab, zzcd
integer(i4b) :: m, i, ip1, im1, k, kp1, km1, l, lp1, ijk
real(dp) :: ooppq, qpq, poq, qoppq, rpq2, afac(0:2)
real(dp) :: bfac(0:4), os(35,0:4), oset(35,10), osc(10,10), rpw(3), rap(3)
real(dp) :: rba(3), rdc(3), delv(3), ketmat(10,3,3), wx(3)

osc=0._dp
ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rdc=xsh(:,sh3)-xsh(:,sh4)
rba=xsh(:,sh1)-xsh(:,sh2)
do zzab=1,nksqu(ab)
    if (negab(zzab,ab)) cycle
    rap=Px(:,zzab,ab)-xsh(:,sh1)
    zb=mod((zzab-1),nzet(sh2))+1
    do zzcd=1,nksqu(cd)
        if (negab(zzcd,cd)) cycle
        zd=mod((zzcd-1),nzet(sh4))+1
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rpw=Wx-Px(:,zzab,ab)
        ! Build (ss|ss)^m, m=0,1,2,3,4
        call boys_n(4,qpq*rpq2,bfac)
        os(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ps|ss)^m, m=0,3
        do m=0,3
          os(2:4,m)=rap(:)*os(1,m)+rpw(:)*os(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0,2
        afac(0:2)=(os(1,0:2)-qoppq*os(1,1:3))*hpm1(zzab,ab)
        do ip1=4,10
          do ijk=1,3
            i=minus1s(ijk,ip1)
            if(i.eq.0) cycle
            os(ip1,0:2)=os(i,0:2)*rap(ijk)+os(i,1:3)*rpw(ijk)
            im1=minus1s(ijk,i)
            if (im1.eq.0) cycle
            os(ip1,0:2)=os(ip1,0:2)+afac(0:2)*real(indexps(ijk,i),dp)
          end do
        end do
        do ip1=11,20
          do ijk=1,3
            i=minus1(ijk,ip1)
            if(i.eq.0) cycle
            os(ip1,0:1)=os(i,0:1)*rap(ijk)+os(i,1:2)*rpw(ijk)
            im1=minus1(ijk,i)
            if (im1.eq.0) cycle
            os(ip1,0:1)=os(ip1,0:1)+real(indexps(ijk,i),dp)*(os(im1,0:1)-qoppq*os(im1,1:2))*hpm1(zzab,ab)
          end do
        end do
        do ip1=21,35
          do ijk=1,3
            i=minus1(ijk,ip1)
            if(i.eq.0) cycle
            os(ip1,0)=os(i,0)*rap(ijk)+os(i,1)*rpw(ijk)
            im1=minus1(ijk,i)
            if (im1.eq.0) cycle
            os(ip1,0)=os(ip1,0)+real(indexps(ijk,i),dp)*(os(im1,0)-qoppq*os(im1,1))*hpm1(zzab,ab)
          end do
        end do
        
        poq=pab(zzab,ab)*pm1(zzcd,cd)
        delv=-pm1(zzcd,cd)*(zeta(zb,sh2)*rba(:)+zeta(zd,sh4)*rdc(:))
        oset(:,1)=os(:,0)
        do i=2,20
          do ijk=1,3
            ip1=plus1(ijk,i)
            kp1=ijk+1
            oset(i,kp1)=delv(ijk)*oset(i,1)-poq*oset(ip1,1)
            im1=minus1(ijk,i)
            if (im1.eq.0) cycle
            oset(i,kp1)=oset(i,kp1)+hpm1(zzcd,cd)*real(indexps(ijk,i),dp)*oset(im1,1)
          end do
        end do
         
        ! Then we assemble [ps|ds], which requires [ps|ps], [ds|ps], [ss|ps], [ss|ss] and [ps|ss]
        do i=1,4
          do k=1,4
            do ijk=1,3
              ip1=plus1s(ijk,i)
              kp1=plus1s(ijk,k)
              oset(i,kp1)=delv(ijk)*oset(i,k)-poq*oset(ip1,k)
              km1=minus1s(ijk,k)
              im1=minus1s(ijk,i)
              if (km1.ne.0) then
                oset(i,kp1)=oset(i,kp1)+real(indexps(ijk,k),dp)*hpm1(zzcd,cd)*oset(i,km1)
              end if
              if (im1.ne.0) then
                oset(i,kp1)=oset(i,kp1)+real(indexps(ijk,i),dp)*hpm1(zzcd,cd)*oset(im1,k)
              end if
            end do
          end do
        end do
        ! Then we assemble [ds|ds], which requires [ds|ps], [fs|ps], [ps|ps] and [ds|ss]
        do i=2,10
          do k=2,4 
            do ijk=1,3
              ip1=plus1s(ijk,i)
              kp1=plus1s(ijk,k)
              oset(i,kp1)=delv(ijk)*oset(i,k)-poq*oset(ip1,k)
              km1=minus1s(ijk,k)
              im1=minus1s(ijk,i)
              if (km1.ne.0) then
                oset(i,kp1)=oset(i,kp1)+real(indexps(ijk,k),dp)*hpm1(zzcd,cd)*oset(i,km1)
              end if
              if (im1.ne.0) then
                oset(i,kp1)=oset(i,kp1)+real(indexps(ijk,i),dp)*hpm1(zzcd,cd)*oset(im1,k)
              end if
            end do
          end do
        end do
        osc=osc+oset(1:10,:)
    end do
end do

! Then the HRR. Start with the ket
! Here i Use as input the (ds|ds), (ds|ps), (ps|ds), (ps|ps) i.e. osc(2:10,2:10)
do ijk=1,3
  do k=2,4
    kp1=plus1s(ijk,k)
    l=1
    lp1=plus1s(ijk,l)
    ketmat(:,k-1,lp1-1)=osc(:,kp1)+rdc(ijk)*osc(:,k)
  end do
end do

! Now the bra
! Here i Use as input the (ds|pp), (ps|pp), i.e. ketmat(2:10,1:3,1:3)
do ijk=1,3
  do i=2,4
    ip1=plus1s(ijk,i)
    epppp(i-1,ijk,:,:)=ketmat(ip1,:,:)+rba(ijk)*ketmat(i,:,:)
  end do
end do
return

end subroutine i2e_pppp


pure subroutine i2e_sspl(sh1,sh2,sh3,sh4,esspl)
use nrtype ; use molprops ; use boys
implicit none

! This version uses standard OS VRR

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esspl(3,4)

integer(i4b) :: ab, cd, zzab, zzcd, m
integer(i4b) :: k, jj
real(dp) :: ooppq,poppq, rpq2, fac, qpq
real(dp) :: bfac(0:2), os(4,4,0:2), osc(4,4), rcq(3), rdq(3), rqw(3), wx(3)

! First assemble the [ss|sx] (x=s,p,d) integrals and contract them

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzab=1,nksqu(ab)
    if (negab(zzab,ab)) cycle
    do zzcd=1,nksqu(cd)
        if (negab(zzcd,cd)) cycle
        rcq=Px(:,zzcd,cd)-xsh(:,sh3)
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        call boys_2(qpq*rpq2,bfac)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2
        os(1,1,0:2)=kab(zzab,ab)*kab(zzcd,cd)*bfac(0:2)*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,1
        do m=0,1
          os(1,2:4,m)=(rdq(:)*os(1,1,m)+rqw(:)*os(1,1,m+1))*ddabl(zzcd,cd)
          os(2:4,1,m)=rcq(:)*os(1,1,m)+rqw(:)*os(1,1,m+1)
        end do
        ! OS recursion: Build (ss|pp) parts
        fac=(os(1,1,0)-poppq*os(1,1,1))*hpm1(zzcd,cd)*ddabl(zzcd,cd)
        do jj=1,3
          k=1+jj
          os(k,2:4,0)=rcq(jj)*os(1,2:4,0)+rqw(jj)*os(1,2:4,1)
          os(k,k,0)=os(k,k,0)+fac
        end do
        osc=osc+os(:,:,0)
    end do
end do

esspl=osc(2:4,:)

end subroutine i2e_sspl



pure subroutine i2e_spsl(sh1,sh2,sh3,sh4,espsl)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: espsl(3,4)

integer(i4b) :: ab, cd, zzab, zzcd, m
integer(i4b) :: j
real(dp) :: ooppq, qpq, qoppq, rpq2, fac1, wx(3)
real(dp) :: bfac(0:2), os(4,4,0:2), rbp(3), rpw(3), rdq(3), rqw(3)

! Simply assemble the [sx|sy] (x,y=s,p) integrals and contract them

espsl=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzab=1,nksqu(ab)
    if (negab(zzab,ab)) cycle
    rbp=Px(:,zzab,ab)-xsh(:,sh2)
    do zzcd=1,nksqu(cd)
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rpw=Wx-Px(:,zzab,ab)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2
        call boys_2(qpq*rpq2,bfac)
        os(1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,1
        do m=0,1
          os(1,2:4,m)=(rdq(:)*os(1,1,m)+rqw(:)*os(1,1,m+1))*ddabl(zzcd,cd)
        end do
        os(2:4,1,0)=rbp(:)*os(1,1,0)+rpw(:)*os(1,1,1)
        ! OS recursion: Build (sp|sp)
        fac1=.5_dp*os(1,1,1)*ooppq*ddabl(zzcd,cd)
        do j=1,3
          os(1+j,2:4,0)=rbp(j)*os(1,2:4,0)+rpw(j)*os(1,2:4,1)
          os(j+1,j+1,0)=os(j+1,j+1,0)+fac1
        end do
        espsl=espsl+os(2:4,:,0)
    end do
end do

end subroutine i2e_spsl



pure subroutine i2e_sppl(sh1,sh2,sh3,sh4,esppl)
use nrtype ; use molprops ; use boys
implicit none

! This version simply uses O-S recursion

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esppl(3,3,4)

integer(i4b) :: ab, cd, zzab, zzcd, m, jj, k
real(dp) :: ooppq, qpq, poppq, qoppq, rpq2
real(dp) :: fackl(0:1), facjl(0:1), facjk(0:1)
real(dp) :: bfac(0:3), os(4,4,4,0:3), rpw(3), rbp(3), rcq(3),rdq(3), rqw(3), wx(3)

esppl=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzcd=1,nksqu(cd)
    if (negab(zzcd,cd)) cycle
    rcq=Px(:,zzcd,cd)-xsh(:,sh3)
    rdq=Px(:,zzcd,cd)-xsh(:,sh4)
    do zzab=1,nksqu(ab)
        if (negab(zzab,ab)) cycle
        rbp=Px(:,zzab,ab)-xsh(:,sh2)
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rpw=Wx-Px(:,zzab,ab)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3
        call boys_n(3,qpq*rpq2,bfac)
        os(1,1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)  !   1 element out of 64 done
        ! OS recursion: Build (ss|sp)^m, m=0,2
        do m=0,2
          os(1,1,2:4,m)=(rdq(:)*os(1,1,1,m)+rqw(:)*os(1,1,1,m+1))*ddabl(zzcd,cd)
          os(1,2:4,1,m)=(rcq(:)*os(1,1,1,m)+rqw(:)*os(1,1,1,m+1)) 
          os(2:4,1,1,m)=(rbp(:)*os(1,1,1,m)+rpw(:)*os(1,1,1,m+1)) 
        end do
        ! OS recursion: Build (ss|pp), (sp|sp) and (sp|ps)
        fackl=(os(1,1,1,0:1)-poppq*os(1,1,1,1:2))*hpm1(zzcd,cd)*ddabl(zzcd,cd)
        facjl=os(1,1,1,1:2)*.5_dp*ooppq*ddabl(zzcd,cd)
        facjk=os(1,1,1,1:2)*.5_dp*ooppq
        do jj=1,3
          k=jj+1
          os(1,k,2:4,0:1)=(rcq(jj)*os(1,1,2:4,0:1)+rqw(jj)*os(1,1,2:4,1:2))
          os(1,k,k,0:1)=os(1,k,k,0:1)+fackl
          os(k,1,2:4,0:1)=(rbp(jj)*os(1,1,2:4,0:1)+rpw(jj)*os(1,1,2:4,1:2))
          os(k,1,k,0:1)=os(k,1,k,0:1)+facjl
          os(k,2:4,1,0:1)=(rbp(jj)*os(1,2:4,1,0:1)+rpw(jj)*os(1,2:4,1,1:2))
          os(k,k,1,0:1)=os(k,k,1,0:1)+facjk
        end do
        ! Last step: build the (sp|pp) bit
        do jj=1,3
          k=jj+1
          os(k,2:4,2:4,0)=(os(1,2:4,2:4,0)*rbp(jj)+os(1,2:4,2:4,1)*rpw(jj))
          os(k,2:4,k,0)=os(k,2:4,k,0)+.5_dp*ooppq*os(1,2:4,1,1)*ddabl(zzcd,cd)
          os(k,k,2:4,0)=os(k,k,2:4,0)+.5_dp*ooppq*os(1,1,2:4,1)
        end do
        esppl=esppl+os(2:4,2:4,:,0)
    end do
end do

end subroutine i2e_sppl



pure subroutine i2e_spll(sh1,sh2,sh3,sh4,espll)
use nrtype ; use molprops ; use boys
implicit none

! This version simply uses O-S recursion

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: espll(3,4,4)

integer(i4b) :: ab, cd, zc, zd, zzab, zzcd, m, jj, k
real(dp) :: ooppq, qpq, poppq, qoppq, rpq2
real(dp) :: fackl(0:1), facjl(0:1), facjk(0:1)
real(dp) :: bfac(0:3), os(4,4,4,0:3), rpw(3), rbp(3), rcq(3),rdq(3), rqw(3), wx(3)

espll=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzcd=1,nksqu(cd)
    if (negab(zzcd,cd)) cycle
    zc=int((zzcd-1)/nzet(sh4))+1
    zd=mod((zzcd-1),nzet(sh4))+1
    rcq=Px(:,zzcd,cd)-xsh(:,sh3)
    rdq=Px(:,zzcd,cd)-xsh(:,sh4)
    do zzab=1,nksqu(ab)
        if (negab(zzab,ab)) cycle
        rbp=Px(:,zzab,ab)-xsh(:,sh2)
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rpw=Wx-Px(:,zzab,ab)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3
        call boys_n(3,qpq*rpq2,bfac)
        os(1,1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)  !   1 element out of 64 done
        ! OS recursion: Build (ss|sp)^m, m=0,2
        do m=0,2
          os(1,1,2:4,m)=(rdq(:)*os(1,1,1,m)+rqw(:)*os(1,1,1,m+1))*djkl(zd,sh4)  
          os(1,2:4,1,m)=(rcq(:)*os(1,1,1,m)+rqw(:)*os(1,1,1,m+1))*djkl(zc,sh3)  
          os(2:4,1,1,m)=(rbp(:)*os(1,1,1,m)+rpw(:)*os(1,1,1,m+1))
        end do
        ! OS recursion: Build (ss|pp), (sp|sp) and (sp|ps)
        fackl=(os(1,1,1,0:1)-poppq*os(1,1,1,1:2))*hpm1(zzcd,cd)*ddabl(zzcd,cd)
        facjl=os(1,1,1,1:2)*.5_dp*ooppq*djkl(zd,sh4)
        facjk=os(1,1,1,1:2)*.5_dp*ooppq*djkl(zc,sh3)
        do jj=1,3
          k=jj+1
          os(1,k,2:4,0:1)=(rcq(jj)*os(1,1,2:4,0:1)+rqw(jj)*os(1,1,2:4,1:2))*djkl(zc,sh3)  ! 19
          os(1,k,k,0:1)=os(1,k,k,0:1)+fackl
          os(k,1,2:4,0:1)=(rbp(jj)*os(1,1,2:4,0:1)+rpw(jj)*os(1,1,2:4,1:2))
          os(k,1,k,0:1)=os(k,1,k,0:1)+facjl
          os(k,2:4,1,0:1)=(rbp(jj)*os(1,2:4,1,0:1)+rpw(jj)*os(1,2:4,1,1:2))
          os(k,k,1,0:1)=os(k,k,1,0:1)+facjk
        end do
        ! Last step: build the (sp|pp) bit
        do jj=1,3
          k=jj+1
          os(k,2:4,2:4,0)=(os(1,2:4,2:4,0)*rbp(jj)+os(1,2:4,2:4,1)*rpw(jj))
          os(k,2:4,k,0)=os(k,2:4,k,0)+.5_dp*ooppq*os(1,2:4,1,1)*djkl(zd,sh4)
          os(k,k,2:4,0)=os(k,k,2:4,0)+.5_dp*ooppq*os(1,1,2:4,1)*djkl(zc,sh3)
        end do
        espll=espll+os(2:4,:,:,0)
    end do
end do

end subroutine i2e_spll


pure subroutine i2e_slpp(sh1,sh2,sh3,sh4,eslpp)
use nrtype ; use molprops ; use boys
implicit none

! This version simply uses O-S recursion

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: eslpp(4,3,3)

integer(i4b) :: ab, cd, zzab, zzcd, m, jj, k
real(dp) :: ooppq, qpq, poppq, qoppq, rpq2
real(dp) :: fackl(0:1), facjl(0:1), facjk(0:1)
real(dp) :: bfac(0:3), os(4,4,4,0:3), rpw(3), rbp(3), rcq(3),rdq(3), rqw(3), wx(3)

eslpp=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzcd=1,nksqu(cd)
    if (negab(zzcd,cd)) cycle
    rcq=Px(:,zzcd,cd)-xsh(:,sh3)
    rdq=Px(:,zzcd,cd)-xsh(:,sh4)
    do zzab=1,nksqu(ab)
        if (negab(zzab,ab)) cycle
        rbp=Px(:,zzab,ab)-xsh(:,sh2)
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rpw=Wx-Px(:,zzab,ab)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3
        call boys_n(3,qpq*rpq2,bfac)
        os(1,1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)  !   1 element out of 64 done
        ! OS recursion: Build (ss|sp)^m, m=0,2
        do m=0,2
          os(1,1,2:4,m)=(rdq(:)*os(1,1,1,m)+rqw(:)*os(1,1,1,m+1))
          os(1,2:4,1,m)=(rcq(:)*os(1,1,1,m)+rqw(:)*os(1,1,1,m+1))
          os(2:4,1,1,m)=(rbp(:)*os(1,1,1,m)+rpw(:)*os(1,1,1,m+1))*ddabl(zzab,ab)
        end do
        ! OS recursion: Build (ss|pp), (sp|sp) and (sp|ps)
        fackl=(os(1,1,1,0:1)-poppq*os(1,1,1,1:2))*hpm1(zzcd,cd)
        facjl=os(1,1,1,1:2)*.5_dp*ooppq*ddabl(zzab,ab)
        facjk=os(1,1,1,1:2)*.5_dp*ooppq*ddabl(zzab,ab)
        do jj=1,3
          k=jj+1
          os(1,k,2:4,0:1)=(rcq(jj)*os(1,1,2:4,0:1)+rqw(jj)*os(1,1,2:4,1:2))
          os(1,k,k,0:1)=os(1,k,k,0:1)+fackl
          os(k,1,2:4,0:1)=(rbp(jj)*os(1,1,2:4,0:1)+rpw(jj)*os(1,1,2:4,1:2))*ddabl(zzab,ab)! 28
          os(k,1,k,0:1)=os(k,1,k,0:1)+facjl
          os(k,2:4,1,0:1)=(rbp(jj)*os(1,2:4,1,0:1)+rpw(jj)*os(1,2:4,1,1:2))*ddabl(zzab,ab)! 37
          os(k,k,1,0:1)=os(k,k,1,0:1)+facjk
        end do
        ! Last step: build the (sp|pp) bit
        do jj=1,3
          k=jj+1
          os(k,2:4,2:4,0)=(os(1,2:4,2:4,0)*rbp(jj)+os(1,2:4,2:4,1)*rpw(jj))*ddabl(zzab,ab)! 64
          os(k,2:4,k,0)=os(k,2:4,k,0)+.5_dp*ooppq*os(1,2:4,1,1)*ddabl(zzab,ab)
          os(k,k,2:4,0)=os(k,k,2:4,0)+.5_dp*ooppq*os(1,1,2:4,1)*ddabl(zzab,ab)
        end do
        eslpp=eslpp+os(:,2:4,2:4,0)
    end do
end do

end subroutine i2e_slpp


pure subroutine i2e_slpl(sh1,sh2,sh3,sh4,eslpl)
use nrtype ; use molprops ; use boys
implicit none

! This version simply uses O-S recursion

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: eslpl(4,3,4)

integer(i4b) :: ab, cd, zzab, zzcd, m, jj, k
real(dp) :: ooppq, qpq, poppq, qoppq, rpq2
real(dp) :: ggjl, fackl(0:1), facjl(0:1), facjk(0:1)
real(dp) :: bfac(0:3), os(4,4,4,0:3), rpw(3), rbp(3), rcq(3),rdq(3), rqw(3), wx(3)

eslpl=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzcd=1,nksqu(cd)
    if (negab(zzcd,cd)) cycle
    rcq=Px(:,zzcd,cd)-xsh(:,sh3)
    rdq=Px(:,zzcd,cd)-xsh(:,sh4)
    do zzab=1,nksqu(ab)
        if (negab(zzab,ab)) cycle
        ggjl=ddabl(zzab,ab)*ddabl(zzcd,cd)
        rbp=Px(:,zzab,ab)-xsh(:,sh2)
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rpw=Wx-Px(:,zzab,ab)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3
        call boys_n(3,qpq*rpq2,bfac)
        os(1,1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)  !   1 element out of 64 done
        ! OS recursion: Build (ss|sp)^m, m=0,2
        do m=0,2
          os(1,1,2:4,m)=(rdq(:)*os(1,1,1,m)+rqw(:)*os(1,1,1,m+1))*ddabl(zzcd,cd) !  4 / 64
          os(1,2:4,1,m)=(rcq(:)*os(1,1,1,m)+rqw(:)*os(1,1,1,m+1))
          os(2:4,1,1,m)=(rbp(:)*os(1,1,1,m)+rpw(:)*os(1,1,1,m+1))*ddabl(zzab,ab) ! 10 / 64
        end do
        ! OS recursion: Build (ss|pp), (sp|sp) and (sp|ps)
        fackl=(os(1,1,1,0:1)-poppq*os(1,1,1,1:2))*hpm1(zzcd,cd)*ddabl(zzcd,cd)
        facjl=os(1,1,1,1:2)*.5_dp*ooppq*ggjl
        facjk=os(1,1,1,1:2)*.5_dp*ooppq*ddabl(zzab,ab)
        do jj=1,3
          k=jj+1
          os(1,k,2:4,0:1)=(rcq(jj)*os(1,1,2:4,0:1)+rqw(jj)*os(1,1,2:4,1:2))
          os(1,k,k,0:1)=os(1,k,k,0:1)+fackl
          os(k,1,2:4,0:1)=(rbp(jj)*os(1,1,2:4,0:1)+rpw(jj)*os(1,1,2:4,1:2))*ddabl(zzab,ab)! 28
          os(k,1,k,0:1)=os(k,1,k,0:1)+facjl
          os(k,2:4,1,0:1)=(rbp(jj)*os(1,2:4,1,0:1)+rpw(jj)*os(1,2:4,1,1:2))*ddabl(zzab,ab)! 37
          os(k,k,1,0:1)=os(k,k,1,0:1)+facjk
        end do
        ! Last step: build the (sp|pp) bit
        do jj=1,3
          k=jj+1
          os(k,2:4,2:4,0)=(os(1,2:4,2:4,0)*rbp(jj)+os(1,2:4,2:4,1)*rpw(jj))*ddabl(zzab,ab)! 64
          os(k,2:4,k,0)=os(k,2:4,k,0)+.5_dp*ooppq*os(1,2:4,1,1)*ggjl
          os(k,k,2:4,0)=os(k,k,2:4,0)+.5_dp*ooppq*os(1,1,2:4,1)*ddabl(zzab,ab)
        end do
        eslpl=eslpl+os(:,2:4,:,0)
    end do
end do

end subroutine i2e_slpl



pure subroutine i2e_plpl(sh1,sh2,sh3,sh4,eplpl)
use nrtype ; use molprops ; use boys
implicit none

! This version uses only O-S recursion

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: eplpl(3,4,3,4)

integer(i4b) :: ab, cd, zzab, zzcd, m, jj, k
real(dp) :: ooppq, qpq, poppq, qoppq, rpq2, fac1(3)
real(dp) :: bfac(0:4), os(4,4,4,4,0:4), rap(3), rpw(3), rbp(3), rcq(3), rdq(3), rqw(3)
real(dp) :: wx(3), eplplc(3,4,3,4)

! First assemble the [sx|sy] (x=s,p; y=s,p,d) integrals and contract them
! We use the matrix os(j1,j2,l1,l2,m) depending on the two angular momenta directions (0 = none,
! 1 = x, 2 = y, 3 = z) at b, the two at d, and the Boys order m.

eplpl=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzab=1,nksqu(ab)
    if (negab(zzab,ab)) cycle
    rap=Px(:,zzab,ab)-xsh(:,sh1)
    rbp=Px(:,zzab,ab)-xsh(:,sh2)
    eplplc=0._dp
    do zzcd=1,nksqu(cd)
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rcq=Px(:,zzcd,cd)-xsh(:,sh3)
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rpw=Wx-Px(:,zzab,ab)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4
        call boys_n(4,qpq*rpq2,bfac)
        os(1,1,1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,3
        do m=0,3
          os(1,1,1,2:4,m)=rdq(:)*os(1,1,1,1,m)+rqw(:)*os(1,1,1,1,m+1)
          os(1,1,2:4,1,m)=(rcq(:)*os(1,1,1,1,m)+rqw(:)*os(1,1,1,1,m+1))
          os(1,2:4,1,1,m)=rbp(:)*os(1,1,1,1,m)+rpw(:)*os(1,1,1,1,m+1)
          os(2:4,1,1,1,m)=(rap(:)*os(1,1,1,1,m)+rpw(:)*os(1,1,1,1,m+1))
        end do
        ! OS recursion: Build (sp|sp)^m, m=0,1,2 and other s,s,p,p terms
        fac1=.5_dp*ooppq*os(1,1,1,1,1:3)
        do jj=1,3
          k=jj+1
          os(1,k,1,2:4,0:2)=(rbp(jj)*os(1,1,1,2:4,0:2)+rpw(jj)*os(1,1,1,2:4,1:3))
          os(1,k,1,k,0:2)=os(1,k,1,k,0:2)+fac1
          os(k,1,1,2:4,0:2)=(rap(jj)*os(1,1,1,2:4,0:2)+rpw(jj)*os(1,1,1,2:4,1:3))
          os(k,1,1,k,0:2)=os(k,1,1,k,0:2)+fac1
          os(1,k,2:4,1,0:2)=(rbp(jj)*os(1,1,2:4,1,0:2)+rpw(jj)*os(1,1,2:4,1,1:3))
          os(1,k,k,1,0:2)=os(1,k,k,1,0:2)+fac1
          os(k,1,2:4,1,0:2)=(rap(jj)*os(1,1,2:4,1,0:2)+rpw(jj)*os(1,1,2:4,1,1:3))
          os(k,1,k,1,0:2)=os(k,1,k,1,0:2)+fac1
          os(1,1,k,2:4,0:2)=(rcq(jj)*os(1,1,1,2:4,0:2)+rqw(jj)*os(1,1,1,2:4,1:3))
          os(1,1,k,k,0:2)=os(1,1,k,k,0:2)+(os(1,1,1,1,0:2)-poppq*os(1,1,1,1,1:3))*hpm1(zzcd,cd)
          os(k,2:4,1,1,0:2)=(rap(jj)*os(1,2:4,1,1,0:2)+rpw(jj)*os(1,2:4,1,1,1:3))
          os(k,k,1,1,0:2)=os(k,k,1,1,0:2)+(os(1,1,1,1,0:2)-qoppq*os(1,1,1,1,1:3))*hpm1(zzab,ab)
        end do
        ! OS recursion: Build 4 different (sp|pp)^m, m=0,1 parts. Only need m=1 for pne of them
        do jj=1,3
          k=jj+1
          os(1,k,2:4,2:4,0:1)=(os(1,1,2:4,2:4,0:1)*rbp(jj)+os(1,1,2:4,2:4,1:2)*rpw(jj))
          os(1,k,2:4,k,0:1)=os(1,k,2:4,k,0:1)+.5_dp*ooppq*os(1,1,2:4,1,1:2)
          os(1,k,k,2:4,0:1)=os(1,k,k,2:4,0:1)+.5_dp*ooppq*os(1,1,1,2:4,1:2)
          os(k,1,2:4,2:4,0)=(os(1,1,2:4,2:4,0)*rap(jj)+os(1,1,2:4,2:4,1)*rpw(jj))
          os(k,1,2:4,k,0)=os(k,1,2:4,k,0)+.5_dp*ooppq*os(1,1,2:4,1,1)
          os(k,1,k,2:4,0)=os(k,1,k,2:4,0)+.5_dp*ooppq*os(1,1,1,2:4,1)
          os(2:4,2:4,k,1,0)=(os(2:4,2:4,1,1,0)*rcq(jj)+os(2:4,2:4,1,1,1)*rqw(jj))
          os(k,2:4,k,1,0)=os(k,2:4,k,1,0)+.5_dp*ooppq*os(1,2:4,1,1,1)
          os(2:4,k,k,1,0)=os(2:4,k,k,1,0)+.5_dp*ooppq*os(2:4,1,1,1,1)
          os(2:4,2:4,1,k,0)=(os(2:4,2:4,1,1,0)*rdq(jj)+os(2:4,2:4,1,1,1)*rqw(jj))
          os(k,2:4,1,k,0)=os(k,2:4,1,k,0)+.5_dp*ooppq*os(1,2:4,1,1,1)
          os(2:4,k,1,k,0)=os(2:4,k,1,k,0)+.5_dp*ooppq*os(2:4,1,1,1,1)
        end do
        ! Final recursion: The (pp|pp) part
        do jj=1,3
          k=jj+1
          os(k,2:4,2:4,2:4,0)=(os(1,2:4,2:4,2:4,0)*rap(jj)+os(1,2:4,2:4,2:4,1)*rpw(jj))
          os(k,k,2:4,2:4,0)=os(k,k,2:4,2:4,0)+(os(1,1,2:4,2:4,0)-qoppq*os(1,1,2:4,2:4,1))*hpm1(zzab,ab)
          os(k,2:4,k,2:4,0)=os(k,2:4,k,2:4,0)+.5_dp*ooppq*os(1,2:4,1,2:4,1)
          os(k,2:4,2:4,k,0)=os(k,2:4,2:4,k,0)+.5_dp*ooppq*os(1,2:4,2:4,1,1)
        end do
        eplplc(:,:,:,1)=eplplc(:,:,:,1)+os(2:4,:,2:4,1,0)
        eplplc(:,:,:,2:4)=eplplc(:,:,:,2:4)+os(2:4,:,2:4,2:4,0)*ddabl(zzcd,cd)
    end do
    eplpl(:,1,:,:)=eplpl(:,1,:,:)+eplplc(:,1,:,:)
    eplpl(:,2:4,:,:)=eplpl(:,2:4,:,:)+eplplc(:,2:4,:,:)*ddabl(zzab,ab)
end do

end subroutine i2e_plpl


pure subroutine i2e_pppl(sh1,sh2,sh3,sh4,epppl)
use nrtype ; use molprops ; use boys
implicit none

! This version uses only O-S recursion

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: epppl(3,3,3,4)

integer(i4b) :: ab, cd, zzab, zzcd, m, jj, k
real(dp) :: ooppq, qpq, poppq, qoppq, rpq2, fac1(3)
real(dp) :: bfac(0:4), os(4,4,4,4,0:4), rap(3), rpw(3), rbp(3), rcq(3), rdq(3), rqw(3)
real(dp) :: wx(3)

! First assemble the [sx|sy] (x=s,p; y=s,p,d) integrals and contract them
! We use the matrix os(j1,j2,l1,l2,m) depending on the two angular momenta directions (0 = none,
! 1 = x, 2 = y, 3 = z) at b, the two at d, and the Boys order m.

epppl=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzab=1,nksqu(ab)
    if (negab(zzab,ab)) cycle
    rap=Px(:,zzab,ab)-xsh(:,sh1)
    rbp=Px(:,zzab,ab)-xsh(:,sh2)
    do zzcd=1,nksqu(cd)
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rcq=Px(:,zzcd,cd)-xsh(:,sh3)
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rpw=Wx-Px(:,zzab,ab)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4
        call boys_n(4,qpq*rpq2,bfac)
        os(1,1,1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,3
        do m=0,3
          os(1,1,1,2:4,m)=ddabl(zzcd,cd)*(rdq(:)*os(1,1,1,1,m)+rqw(:)*os(1,1,1,1,m+1))
          os(1,1,2:4,1,m)=(rcq(:)*os(1,1,1,1,m)+rqw(:)*os(1,1,1,1,m+1))
          os(1,2:4,1,1,m)=(rbp(:)*os(1,1,1,1,m)+rpw(:)*os(1,1,1,1,m+1))
          os(2:4,1,1,1,m)=(rap(:)*os(1,1,1,1,m)+rpw(:)*os(1,1,1,1,m+1))
        end do
        ! OS recursion: Build (sp|sp)^m, m=0,1,2 and other s,s,p,p terms
        fac1=.5_dp*ooppq*os(1,1,1,1,1:3)
        do jj=1,3
          k=jj+1
          os(1,k,1,2:4,0:2)=(rbp(jj)*os(1,1,1,2:4,0:2)+rpw(jj)*os(1,1,1,2:4,1:3))
          os(1,k,1,k,0:2)=os(1,k,1,k,0:2)+fac1*ddabl(zzcd,cd)
          os(k,1,1,2:4,0:2)=(rap(jj)*os(1,1,1,2:4,0:2)+rpw(jj)*os(1,1,1,2:4,1:3))
          os(k,1,1,k,0:2)=os(k,1,1,k,0:2)+fac1*ddabl(zzcd,cd)
          os(1,k,2:4,1,0:2)=(rbp(jj)*os(1,1,2:4,1,0:2)+rpw(jj)*os(1,1,2:4,1,1:3))
          os(1,k,k,1,0:2)=os(1,k,k,1,0:2)+fac1
          os(k,1,2:4,1,0:2)=(rap(jj)*os(1,1,2:4,1,0:2)+rpw(jj)*os(1,1,2:4,1,1:3))
          os(k,1,k,1,0:2)=os(k,1,k,1,0:2)+fac1

          os(1,1,k,2:4,0:2)=(rcq(jj)*os(1,1,1,2:4,0:2)+rqw(jj)*os(1,1,1,2:4,1:3))
          os(1,1,k,k,0:2)=os(1,1,k,k,0:2)+.5_dp*(os(1,1,1,1,0:2)-poppq*os(1,1,1,1,1:3))*pm1(zzcd,cd)*ddabl(zzcd,cd)

          os(k,2:4,1,1,0:2)=(rap(jj)*os(1,2:4,1,1,0:2)+rpw(jj)*os(1,2:4,1,1,1:3))
          os(k,k,1,1,0:2)=os(k,k,1,1,0:2)+.5_dp*(os(1,1,1,1,0:2)-qoppq*os(1,1,1,1,1:3))*pm1(zzab,ab)
        end do
        ! OS recursion: Build 4 different (sp|pp)^m, m=0,1 parts. Only need m=1 for pne of them
        do jj=1,3
          k=jj+1
          os(1,k,2:4,2:4,0:1)=(os(1,1,2:4,2:4,0:1)*rbp(jj)+os(1,1,2:4,2:4,1:2)*rpw(jj))
          os(1,k,2:4,k,0:1)=os(1,k,2:4,k,0:1)+.5_dp*ooppq*os(1,1,2:4,1,1:2)*ddabl(zzcd,cd)
          os(1,k,k,2:4,0:1)=os(1,k,k,2:4,0:1)+.5_dp*ooppq*os(1,1,1,2:4,1:2)
          os(k,1,2:4,2:4,0)=(os(1,1,2:4,2:4,0)*rap(jj)+os(1,1,2:4,2:4,1)*rpw(jj))
          os(k,1,2:4,k,0)=os(k,1,2:4,k,0)+.5_dp*ooppq*os(1,1,2:4,1,1)*ddabl(zzcd,cd)
          os(k,1,k,2:4,0)=os(k,1,k,2:4,0)+.5_dp*ooppq*os(1,1,1,2:4,1)
          os(2:4,2:4,k,1,0)=(os(2:4,2:4,1,1,0)*rcq(jj)+os(2:4,2:4,1,1,1)*rqw(jj))
          os(k,2:4,k,1,0)=os(k,2:4,k,1,0)+.5_dp*ooppq*os(1,2:4,1,1,1)
          os(2:4,k,k,1,0)=os(2:4,k,k,1,0)+.5_dp*ooppq*os(2:4,1,1,1,1)
          os(2:4,2:4,1,k,0)=(os(2:4,2:4,1,1,0)*rdq(jj)+os(2:4,2:4,1,1,1)*rqw(jj))*ddabl(zzcd,cd)
          os(k,2:4,1,k,0)=os(k,2:4,1,k,0)+.5_dp*ooppq*os(1,2:4,1,1,1)*ddabl(zzcd,cd)
          os(2:4,k,1,k,0)=os(2:4,k,1,k,0)+.5_dp*ooppq*os(2:4,1,1,1,1)*ddabl(zzcd,cd)
        end do
        ! Final recursion: The (pp|pp) part
        do jj=1,3
          k=jj+1
          os(k,2:4,2:4,2:4,0)=(os(1,2:4,2:4,2:4,0)*rap(jj)+os(1,2:4,2:4,2:4,1)*rpw(jj))
          os(k,k,2:4,2:4,0)=os(k,k,2:4,2:4,0)+.5_dp*(os(1,1,2:4,2:4,0)-qoppq*os(1,1,2:4,2:4,1))*pm1(zzab,ab)
          os(k,2:4,k,2:4,0)=os(k,2:4,k,2:4,0)+.5_dp*ooppq*os(1,2:4,1,2:4,1)
          os(k,2:4,2:4,k,0)=os(k,2:4,2:4,k,0)+.5_dp*ooppq*os(1,2:4,2:4,1,1)*ddabl(zzcd,cd)
        end do
        epppl=epppl+os(2:4,2:4,2:4,:,0)
    end do
end do

end subroutine i2e_pppl


pure subroutine i2e_ppll(sh1,sh2,sh3,sh4,eppll)
use nrtype ; use molprops ; use boys
implicit none

! This version uses only O-S recursion

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: eppll(3,3,4,4)

integer(i4b) :: ab, cd, zc, zd, zzab, zzcd, m, jj, k
real(dp) :: ooppq, qpq, poppq, qoppq, rpq2, fac1(3)
real(dp) :: bfac(0:4), os(4,4,4,4,0:4), rap(3), rpw(3), rbp(3), rcq(3), rdq(3), rqw(3)
real(dp) :: wx(3), eppllc(3,3,4,4)

! First assemble the [sx|sy] (x=s,p; y=s,p,d) integrals and contract them
! We use the matrix os(j1,j2,l1,l2,m) depending on the two angular momenta directions (0 = none,
! 1 = x, 2 = y, 3 = z) at b, the two at d, and the Boys order m.

eppll=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzcd=1,nksqu(cd)
    if (negab(zzcd,cd)) cycle
    zc=int((zzcd-1)/nzet(sh4))+1
    zd=mod((zzcd-1),nzet(sh4))+1
    rcq=Px(:,zzcd,cd)-xsh(:,sh3)
    rdq=Px(:,zzcd,cd)-xsh(:,sh4)
    eppllc=0._dp
    do zzab=1,nksqu(ab)
        if (negab(zzab,ab)) cycle
        rap=Px(:,zzab,ab)-xsh(:,sh1)
        rbp=Px(:,zzab,ab)-xsh(:,sh2)
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rpw=Wx-Px(:,zzab,ab)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4
        call boys_n(4,qpq*rpq2,bfac)
        os(1,1,1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,3
        do m=0,3
          os(1,1,1,2:4,m)=djkl(zd,sh4)*(rdq(:)*os(1,1,1,1,m)+rqw(:)*os(1,1,1,1,m+1))
          os(1,1,2:4,1,m)=djkl(zc,sh3)*(rcq(:)*os(1,1,1,1,m)+rqw(:)*os(1,1,1,1,m+1))
          os(1,2:4,1,1,m)=(rbp(:)*os(1,1,1,1,m)+rpw(:)*os(1,1,1,1,m+1))
          os(2:4,1,1,1,m)=(rap(:)*os(1,1,1,1,m)+rpw(:)*os(1,1,1,1,m+1))
        end do
        ! OS recursion: Build (sp|sp)^m, m=0,1,2 and other s,s,p,p terms
        fac1=.5_dp*ooppq*os(1,1,1,1,1:3)
        do jj=1,3
          k=jj+1
          os(1,k,1,2:4,0:2)=(rbp(jj)*os(1,1,1,2:4,0:2)+rpw(jj)*os(1,1,1,2:4,1:3))
          os(1,k,1,k,0:2)=os(1,k,1,k,0:2)+fac1
          os(k,1,1,2:4,0:2)=(rap(jj)*os(1,1,1,2:4,0:2)+rpw(jj)*os(1,1,1,2:4,1:3))
          os(k,1,1,k,0:2)=os(k,1,1,k,0:2)+fac1
          os(1,k,2:4,1,0:2)=(rbp(jj)*os(1,1,2:4,1,0:2)+rpw(jj)*os(1,1,2:4,1,1:3))
          os(1,k,k,1,0:2)=os(1,k,k,1,0:2)+fac1
          os(k,1,2:4,1,0:2)=(rap(jj)*os(1,1,2:4,1,0:2)+rpw(jj)*os(1,1,2:4,1,1:3))
          os(k,1,k,1,0:2)=os(k,1,k,1,0:2)+fac1

          os(1,1,k,2:4,0:2)=(rcq(jj)*os(1,1,1,2:4,0:2)+rqw(jj)*os(1,1,1,2:4,1:3))
          os(1,1,k,k,0:2)=os(1,1,k,k,0:2)+.5_dp*(os(1,1,1,1,0:2)-poppq*os(1,1,1,1,1:3))*pm1(zzcd,cd)

          os(k,2:4,1,1,0:2)=(rap(jj)*os(1,2:4,1,1,0:2)+rpw(jj)*os(1,2:4,1,1,1:3))
          os(k,k,1,1,0:2)=os(k,k,1,1,0:2)+.5_dp*(os(1,1,1,1,0:2)-qoppq*os(1,1,1,1,1:3))*pm1(zzab,ab)
        end do
        ! OS recursion: Build 4 different (sp|pp)^m, m=0,1 parts. Only need m=1 for pne of them
        do jj=1,3
          k=jj+1
          os(1,k,2:4,2:4,0:1)=(os(1,1,2:4,2:4,0:1)*rbp(jj)+os(1,1,2:4,2:4,1:2)*rpw(jj))
          os(1,k,2:4,k,0:1)=os(1,k,2:4,k,0:1)+.5_dp*ooppq*os(1,1,2:4,1,1:2)
          os(1,k,k,2:4,0:1)=os(1,k,k,2:4,0:1)+.5_dp*ooppq*os(1,1,1,2:4,1:2)
          os(k,1,2:4,2:4,0)=(os(1,1,2:4,2:4,0)*rap(jj)+os(1,1,2:4,2:4,1)*rpw(jj))
          os(k,1,2:4,k,0)=os(k,1,2:4,k,0)+.5_dp*ooppq*os(1,1,2:4,1,1)
          os(k,1,k,2:4,0)=os(k,1,k,2:4,0)+.5_dp*ooppq*os(1,1,1,2:4,1)
          os(2:4,2:4,k,1,0)=(os(2:4,2:4,1,1,0)*rcq(jj)+os(2:4,2:4,1,1,1)*rqw(jj))
          os(k,2:4,k,1,0)=os(k,2:4,k,1,0)+.5_dp*ooppq*os(1,2:4,1,1,1)
          os(2:4,k,k,1,0)=os(2:4,k,k,1,0)+.5_dp*ooppq*os(2:4,1,1,1,1)
          os(2:4,2:4,1,k,0)=(os(2:4,2:4,1,1,0)*rdq(jj)+os(2:4,2:4,1,1,1)*rqw(jj))
          os(k,2:4,1,k,0)=os(k,2:4,1,k,0)+.5_dp*ooppq*os(1,2:4,1,1,1)
          os(2:4,k,1,k,0)=os(2:4,k,1,k,0)+.5_dp*ooppq*os(2:4,1,1,1,1)
        end do
        ! Final recursion: The (pp|pp) part
        do jj=1,3
          k=jj+1
          os(k,2:4,2:4,2:4,0)=(os(1,2:4,2:4,2:4,0)*rap(jj)+os(1,2:4,2:4,2:4,1)*rpw(jj))
          os(k,k,2:4,2:4,0)=os(k,k,2:4,2:4,0)+.5_dp*(os(1,1,2:4,2:4,0)-qoppq*os(1,1,2:4,2:4,1))*pm1(zzab,ab)
          os(k,2:4,k,2:4,0)=os(k,2:4,k,2:4,0)+.5_dp*ooppq*os(1,2:4,1,2:4,1)
          os(k,2:4,2:4,k,0)=os(k,2:4,2:4,k,0)+.5_dp*ooppq*os(1,2:4,2:4,1,1)
        end do
        eppllc=eppllc+os(2:4,2:4,:,:,0)
    end do
    eppll(:,:,1,1)=eppll(:,:,1,1)+eppllc(:,:,1,1)
    eppll(:,:,1,2:4)=eppll(:,:,1,2:4)+eppllc(:,:,1,2:4)*djkl(zd,sh4)
    eppll(:,:,2:4,1)=eppll(:,:,2:4,1)+eppllc(:,:,2:4,1)*djkl(zc,sh3)
    eppll(:,:,2:4,2:4)=eppll(:,:,2:4,2:4)+eppllc(:,:,2:4,2:4)*ddabl(zzcd,cd)
end do

end subroutine i2e_ppll



pure subroutine i2e_plll(sh1,sh2,sh3,sh4,eplll)
use nrtype ; use molprops ; use boys
implicit none

! This version uses only O-S recursion

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: eplll(3,4,4,4)

integer(i4b) :: ab, cd, zc, zd, zzab, zzcd, m, jj, k
real(dp) :: ooppq, qpq, poppq, qoppq, rpq2, fac1(3)
real(dp) :: bfac(0:4), os(4,4,4,4,0:4), rap(3), rpw(3), rbp(3), rcq(3), rdq(3), rqw(3)
real(dp) :: wx(3), eplllc(3,4,4,4)

! First assemble the [sx|sy] (x=s,p; y=s,p,d) integrals and contract them
! We use the matrix os(j1,j2,l1,l2,m) depending on the two angular momenta directions (0 = none,
! 1 = x, 2 = y, 3 = z) at b, the two at d, and the Boys order m.

eplll=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zzcd=1,nksqu(cd)
    if (negab(zzcd,cd)) cycle
    zc=int((zzcd-1)/nzet(sh4))+1
    zd=mod((zzcd-1),nzet(sh4))+1
    rcq=Px(:,zzcd,cd)-xsh(:,sh3)
    rdq=Px(:,zzcd,cd)-xsh(:,sh4)
    eplllc=0._dp
    do zzab=1,nksqu(ab)
        if (negab(zzab,ab)) cycle
        rap=Px(:,zzab,ab)-xsh(:,sh1)
        rbp=Px(:,zzab,ab)-xsh(:,sh2)
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rpw=Wx-Px(:,zzab,ab)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4
        call boys_n(4,qpq*rpq2,bfac)
        os(1,1,1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,3
        do m=0,3
          os(1,1,1,2:4,m)=(rdq(:)*os(1,1,1,1,m)+rqw(:)*os(1,1,1,1,m+1))
          os(1,1,2:4,1,m)=(rcq(:)*os(1,1,1,1,m)+rqw(:)*os(1,1,1,1,m+1))
          os(1,2:4,1,1,m)=(rbp(:)*os(1,1,1,1,m)+rpw(:)*os(1,1,1,1,m+1))
          os(2:4,1,1,1,m)=(rap(:)*os(1,1,1,1,m)+rpw(:)*os(1,1,1,1,m+1))
        end do
        ! OS recursion: Build (sp|sp)^m, m=0,1,2 and other s,s,p,p terms
        fac1=.5_dp*ooppq*os(1,1,1,1,1:3)
        do jj=1,3
          k=jj+1
          os(1,k,1,2:4,0:2)=(rbp(jj)*os(1,1,1,2:4,0:2)+rpw(jj)*os(1,1,1,2:4,1:3))
          os(1,k,1,k,0:2)=os(1,k,1,k,0:2)+fac1
          os(k,1,1,2:4,0:2)=(rap(jj)*os(1,1,1,2:4,0:2)+rpw(jj)*os(1,1,1,2:4,1:3))
          os(k,1,1,k,0:2)=os(k,1,1,k,0:2)+fac1
          os(1,k,2:4,1,0:2)=(rbp(jj)*os(1,1,2:4,1,0:2)+rpw(jj)*os(1,1,2:4,1,1:3))
          os(1,k,k,1,0:2)=os(1,k,k,1,0:2)+fac1
          os(k,1,2:4,1,0:2)=(rap(jj)*os(1,1,2:4,1,0:2)+rpw(jj)*os(1,1,2:4,1,1:3))
          os(k,1,k,1,0:2)=os(k,1,k,1,0:2)+fac1

          os(1,1,k,2:4,0:2)=(rcq(jj)*os(1,1,1,2:4,0:2)+rqw(jj)*os(1,1,1,2:4,1:3))
          os(1,1,k,k,0:2)=os(1,1,k,k,0:2)+.5_dp*(os(1,1,1,1,0:2)-poppq*os(1,1,1,1,1:3))*pm1(zzcd,cd)

          os(k,2:4,1,1,0:2)=(rap(jj)*os(1,2:4,1,1,0:2)+rpw(jj)*os(1,2:4,1,1,1:3))
          os(k,k,1,1,0:2)=os(k,k,1,1,0:2)+.5_dp*(os(1,1,1,1,0:2)-qoppq*os(1,1,1,1,1:3))*pm1(zzab,ab)
        end do
        ! OS recursion: Build 4 different (sp|pp)^m, m=0,1 parts. Only need m=1 for pne of them
        do jj=1,3
          k=jj+1
          os(1,k,2:4,2:4,0:1)=(os(1,1,2:4,2:4,0:1)*rbp(jj)+os(1,1,2:4,2:4,1:2)*rpw(jj))
          os(1,k,2:4,k,0:1)=os(1,k,2:4,k,0:1)+.5_dp*ooppq*os(1,1,2:4,1,1:2)
          os(1,k,k,2:4,0:1)=os(1,k,k,2:4,0:1)+.5_dp*ooppq*os(1,1,1,2:4,1:2)
          os(k,1,2:4,2:4,0)=(os(1,1,2:4,2:4,0)*rap(jj)+os(1,1,2:4,2:4,1)*rpw(jj))
          os(k,1,2:4,k,0)=os(k,1,2:4,k,0)+.5_dp*ooppq*os(1,1,2:4,1,1)
          os(k,1,k,2:4,0)=os(k,1,k,2:4,0)+.5_dp*ooppq*os(1,1,1,2:4,1)
          os(2:4,2:4,k,1,0)=(os(2:4,2:4,1,1,0)*rcq(jj)+os(2:4,2:4,1,1,1)*rqw(jj))
          os(k,2:4,k,1,0)=os(k,2:4,k,1,0)+.5_dp*ooppq*os(1,2:4,1,1,1)
          os(2:4,k,k,1,0)=os(2:4,k,k,1,0)+.5_dp*ooppq*os(2:4,1,1,1,1)
          os(2:4,2:4,1,k,0)=(os(2:4,2:4,1,1,0)*rdq(jj)+os(2:4,2:4,1,1,1)*rqw(jj))
          os(k,2:4,1,k,0)=os(k,2:4,1,k,0)+.5_dp*ooppq*os(1,2:4,1,1,1)
          os(2:4,k,1,k,0)=os(2:4,k,1,k,0)+.5_dp*ooppq*os(2:4,1,1,1,1)
        end do
        ! Final recursion: The (pp|pp) part
        do jj=1,3
          k=jj+1
          os(k,2:4,2:4,2:4,0)=(os(1,2:4,2:4,2:4,0)*rap(jj)+os(1,2:4,2:4,2:4,1)*rpw(jj))
          os(k,k,2:4,2:4,0)=os(k,k,2:4,2:4,0)+.5_dp*(os(1,1,2:4,2:4,0)-qoppq*os(1,1,2:4,2:4,1))*pm1(zzab,ab)
          os(k,2:4,k,2:4,0)=os(k,2:4,k,2:4,0)+.5_dp*ooppq*os(1,2:4,1,2:4,1)
          os(k,2:4,2:4,k,0)=os(k,2:4,2:4,k,0)+.5_dp*ooppq*os(1,2:4,2:4,1,1)
        end do
        eplllc(:,1,:,:)=eplllc(:,1,:,:)+os(2:4,1,:,:,0)
        eplllc(:,2:4,:,:)=eplllc(:,2:4,:,:)+os(2:4,2:4,:,:,0)*ddabl(zzab,ab)
    end do
    eplll(:,:,1,1)=eplll(:,:,1,1)+eplllc(:,:,1,1)
    eplll(:,:,1,2:4)=eplll(:,:,1,2:4)+eplllc(:,:,1,2:4)*djkl(zd,sh4)
    eplll(:,:,2:4,1)=eplll(:,:,2:4,1)+eplllc(:,:,2:4,1)*djkl(zc,sh3)
    eplll(:,:,2:4,2:4)=eplll(:,:,2:4,2:4)+eplllc(:,:,2:4,2:4)*ddabl(zzcd,cd)
end do

end subroutine i2e_plll





end module s_sp_l_terms









