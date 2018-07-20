subroutine fock_build_incore(f,p)
use nrtype; use molprops
implicit none

real(dp), intent(in) :: p(nb,nb)
real(dp), intent(inout) :: f(nb,nb)

integer(i4b) :: ij,kl,ia,ja,ka,la,i,j
integer(i8b) :: ii
real(dp) :: fac

kl=1; ij=1
do i=1,naaaa
  ia=aonaaaa(i)
  f(ia,ia)=f(ia,ia)+.5_dp*p(ia,ia)*i2eaaaa(i)
end do
!write (*,*) "Done",naaaa,"(aa|aa) ints"
do i=1,naabb
  ia=revind2(1,aonaabb(i)) ; ka=revind2(2,aonaabb(i))
  f(ia,ia)=f(ia,ia)+p(ka,ka)*i2eaabb(i)
  fac=.5_dp*p(ia,ka)*i2eaabb(i)
  f(ia,ka)=f(ia,ka)-fac
!  f(ka,ia)=f(ka,ia)-fac
  f(ka,ka)=f(ka,ka)+p(ia,ia)*i2eaabb(i)
end do
!!write (*,*) "Done",naabb,"(aa|bb) ints"
!write (*,*) "Starting",naabc,"(aa|bc) ints",size(i2eaabc)
do i=1,naabc
  ia=aonaabc(1,i)
  ka=revind2(1,aonaabc(2,i)) ; la=revind2(2,aonaabc(2,i))
  f(ia,ia)=f(ia,ia)+2._dp*p(ka,la)*i2eaabc(i)
  fac=.5_dp*p(ia,ka)*i2eaabc(i)
  if (ia.gt.la) then
    f(ia,la)=f(ia,la)-fac
  else if (ia.eq.la) then
    f(ia,la)=f(ia,la)-fac*2._dp
  else
    f(la,ia)=f(la,ia)-fac
  end if
  fac=.5_dp*p(ia,la)*i2eaabc(i)
  if (ia.gt.ka) then
    f(ia,ka)=f(ia,ka)-fac
  else if (ia.eq.ka) then
    f(ia,ka)=f(ia,ka)-fac*2._dp
  else
    f(ka,ia)=f(ka,ia)-fac
  end if
  fac=p(ia,ia)*i2eaabc(i)
  f(ka,la)=f(ka,la)+fac
!  f(la,ka)=f(la,ka)+fac
end do
!write (*,*) "Done",naabc,"(aa|bc) ints"
do i=1,nabab
  ia=revind2(1,aonabab(i)) ; ja=revind2(2,aonabab(i))
  fac=1.5_dp*p(ja,ia)*i2eabab(i)
  F(ia,ja) = F(ia,ja) + fac
!  F(ja,ia) = F(ja,ia) + fac
  F(ia,ia) = F(ia,ia) - .5_dp*p(ja,ja) * i2eabab(i)
  F(ja,ja) = F(ja,ja) - .5_dp*p(ia,ia) * i2eabab(i)
end do
!write (*,*) "Done",nabab,"(ab|ab) ints"
do i=1,nabcd_d
  ia=revind2(1,aonabcd_d(1,i)) ; ja=revind2(2,aonabcd_d(1,i))
  ka=revind2(1,aonabcd_d(2,i)) ; la=revind2(2,aonabcd_d(2,i))
  fac=2._dp*  p(la,ka) * i2eabcd_d(i)
  F(ia,ja) = F(ia,ja) +  fac
!  F(ja,ia) = F(ja,ia) +  fac
  fac=.5_dp*p(ja,ka)*i2eabcd_d(i)
  if (ia.gt.la) then
    F(ia,la) = F(ia,la) - fac
  else if (ia.eq.la) then
    F(ia,la) = F(ia,la) - 2._dp*fac
  else
    F(la,ia) = F(la,ia) - fac
  end if
  fac=.5_dp*p(ja,la)*i2eabcd_d(i)
  if (ia.gt.ka) then
    F(ia,ka) = F(ia,ka) - fac
  else if (ia.eq.ka) then
    F(ia,ka) = F(ia,ka) - fac*2._dp
  else
    F(ka,ia) = F(ka,ia) - fac
  end if
  fac=.5_dp*p(ia,ka)*i2eabcd_d(i)
  if (ja.gt.la) then
    F(ja,la) = F(ja,la) - fac
  else if (ja.eq.la) then
    F(ja,la) = F(ja,la) - fac*2._dp
  else
    F(la,ja) = F(la,ja) - fac
  end if
  fac=.5_dp*p(ia,la)*i2eabcd_d(i)
  if (ja.gt.ka) then
    F(ja,ka) = F(ja,ka) - fac
  else if (ja.eq.ka) then
    F(ja,ka) = F(ja,ka) - fac*2._dp
  else
    F(ka,ja) = F(ka,ja) - fac
  end if
  fac=2._dp*p(ja,ia)*i2eabcd_d(i)
  F(ka,la) = F(ka,la) +  fac
!  F(la,ka) = F(la,ka) +   fac
end do
do ii=1,nabcd
  ia=revind2(1,aonabcd(1,ii)) ; ja=revind2(2,aonabcd(1,ii))
  ka=revind2(1,aonabcd(2,ii)) ; la=revind2(2,aonabcd(2,ii))
  fac=2._dp*  p(la,ka) * i2eabcd(ii)
  F(ia,ja) = F(ia,ja) +  fac
!  F(ja,ia) = F(ja,ia) +  fac
  fac=.5_dp*p(ja,ka)*i2eabcd(ii)
  if (ia.gt.la) then
    F(ia,la) = F(ia,la) - fac
  else if (ia.eq.la) then
    F(ia,la) = F(ia,la) - 2._dp*fac
  else
    F(la,ia) = F(la,ia) - fac
  end if
  fac=.5_dp*p(ja,la)*i2eabcd(ii)
  if (ia.gt.ka) then
    F(ia,ka) = F(ia,ka) - fac
  else if (ia.eq.ka) then
    F(ia,ka) = F(ia,ka) - fac*2._dp
  else
    F(ka,ia) = F(ka,ia) - fac
  end if
  fac=.5_dp*p(ia,ka)*i2eabcd(ii)
  if (ja.gt.la) then
    F(ja,la) = F(ja,la) - fac
  else if (ja.eq.la) then
    F(ja,la) = F(ja,la) - fac*2._dp
  else
    F(la,ja) = F(la,ja) - fac
  end if
  fac=.5_dp*p(ia,la)*i2eabcd(ii)
  if (ja.gt.ka) then
    F(ja,ka) = F(ja,ka) - fac
  else if (ja.eq.ka) then
    F(ja,ka) = F(ja,ka) - fac*2._dp
  else
    F(ka,ja) = F(ka,ja) - fac
  end if
  fac=2._dp*p(ja,ia)*i2eabcd(ii)
  F(ka,la) = F(ka,la) +  fac
!  F(la,ka) = F(la,ka) +   fac
end do
!write (*,*) "Done",nabcd,"(ab|cd) ints"
do i=2,nb
  do j=1,i-1
    f(j,i)=f(i,j)
  end do
end do

end subroutine fock_build_incore

