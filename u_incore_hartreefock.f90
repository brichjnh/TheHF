subroutine u_incore_hartreefock()
use nrtype; use molprops
implicit none

integer(i4b), parameter :: maxiter=50
real(dp) :: fa(nb,nb),fprima(nb,nb), fb(nb,nb), fprimb(nb,nb)
real(dp) :: cmata(nb,nb), pa(nb,nb), polda(nb,nb), cmatb(nb,nb), pb(nb,nb), poldb(nb,nb), pj(nb,nb)
real(dp) :: mprod,invseigs(nb,nb), trace
real(dp) :: p1(nb,nb),p2(nb,nb),p3(nb,nb),ptot(nb,nb)
real(dp) :: folda(nb,nb,7),fdiisa(nb,nb),foldb(nb,nb,7),fdiisb(nb,nb)
real(dp) :: hcore(nb,nb), sminhalf(nb,nb)
real(dp) :: resmaxa,resmaxb,resida(nb,nb),residoa(nb,nb,7),residb(nb,nb),residob(nb,nb,7)
real(dp) :: bmat(8,8),vecb(8),solb(8)
real(dp) :: eigsa(nb),eigsb(nb), eeold, deltae, deltapa, deltapb
real(dp), parameter :: conv=1.d-6, conve=1.d-9
integer(i4b) :: i, j,ij,k,jk, mu, nu,oo(14)
logical(lgt) :: diis

! Build the guess for the Fock matrix - either core hamiltonian or Huckel
if (densup) then
  write (*,*) "Density matrix sum initial guess not yet implemented for UHF"
  stop
else
  hcore=t1e+v1e
end if

! Then obtain the matrix S^-1/2 that orthogonalizes the basis
call wrapper_dsyev(nb,s1e,eigsa,pa)
invseigs=0._dp
ij=0._dp
do i=1,nb
    invseigs(i,i)=1._dp/sqrt(eigsa(i))
end do
sminhalf=matmul(pa,matmul(invseigs,transpose(pa)))

! Now transform the core matrix and diagonalize for initial guess

fa=matmul(sminhalf,matmul(hcore,transpose(sminhalf)))
call wrapper_dsyev(nb,fa,eigsa,cmata)
cmata=matmul(sminhalf,cmata)
cmatb=cmata
if (prtlevl.ge.2) then
  write (9,*) ""
  write (9,*) "Eigenvalues for core or Huckel Hamiltonian:"
  write (9,'(5F12.5)') eigsa
end if
if (prtlevl.ge.3) then
  write (9,*) "Eigenvectors for core or Huckel Hamiltonian:"
  do i=1,nb
    write (9,'(7F12.6)') cmata(i,:)
  end do
end if
eelec=0._dp
pa=matmul(cmata(:,1:nocca),transpose(cmata(:,1:nocca)))
pb=matmul(cmata(:,1:noccb),transpose(cmata(:,1:noccb)))
pj=pa+pb
! Finally assign hcore 'properly':
hcore=t1e+v1e

! Set up some things for DIIS
residoa=1000._dp
residob=1000._dp
folda=0._dp
foldb=0._dp
diis=.false.

write (9,*) ""
write (9,*) "Starting HF cycles: cycle, etot, deltae, deltap"
do i=1,maxiter
    fa=hcore
    fb=hcore
    eeold=eelec
    polda=pa
    poldb=pb
    pa=matmul(cmata(:,1:nocca),transpose(cmata(:,1:nocca)))
    pb=matmul(cmatb(:,1:noccb),transpose(cmatb(:,1:noccb)))
    if (i.le.8) then
        pa=.5_dp*pa+.5_dp*polda
        pb=.5_dp*pb+.5_dp*poldb
    end if
    pj=pa+pb
    fa=hcore
    fb=hcore
    call u_fock_build_incore(fa,fb,pa,pb,pj)
    resida=matmul(fa,matmul(pa,s1e))-matmul(s1e,matmul(pa,fa))
    residb=matmul(fb,matmul(pb,s1e))-matmul(s1e,matmul(pb,fb))
    do j=1,6
        residoa(:,:,8-j)=residoa(:,:,7-j)
        residob(:,:,8-j)=residob(:,:,7-j)
        folda(:,:,8-j)=folda(:,:,7-j)
        foldb(:,:,8-j)=foldb(:,:,7-j)
    end do
    residoa(:,:,1)=resida
    residob(:,:,1)=residb
    resmaxa=maxval(abs(resida))
    resmaxb=maxval(abs(residb))
    folda(:,:,1)=fa
    foldb(:,:,1)=fb
    if ((i.ge.6).and.(resmaxa.lt.0.1).and.(resmaxb.lt.0.1).and.(.not.diis)) then
       diis=.true.
       write (9,*) "Starting DIIS"
    end if
    if (diis) then
!      build b matrix
       bmat=0._dp
       bmat(8,:)=-1._dp; bmat(:,8)=-1._dp; bmat(8,8)=0._dp
       vecb=0._dp ; vecb(8)=-1._dp
       do j=1,7
          do k=1,j
              mprod=0._dp
              do jk=1,nb
                  mprod=mprod+sum(residoa(:,jk,j)*residoa(:,jk,k))+sum(residob(:,jk,j)*residob(:,jk,k))
              end do
              bmat(j,k)=mprod
              bmat(k,j)=bmat(j,k)
          end do
       end do
       call wrapper_dsysv(8,bmat,vecb,solb)
       fdiisa=0._dp
       fdiisb=0._dp
       do j=1,7
           fdiisa=fdiisa+solb(j)*folda(:,:,j)
           fdiisb=fdiisb+solb(j)*foldb(:,:,j)
       end do
       fprima=matmul(transpose(sminhalf),matmul(fdiisa,sminhalf))
       fprimb=matmul(transpose(sminhalf),matmul(fdiisb,sminhalf))
    else
       fprima=matmul(transpose(sminhalf),matmul(fa,sminhalf))
       fprimb=matmul(transpose(sminhalf),matmul(fb,sminhalf))
    end if
    call wrapper_dsyev(nb,fprima,eigsa,cmata)
    call wrapper_dsyev(nb,fprimb,eigsb,cmatb)
    cmata=matmul(sminhalf,cmata)
    cmatb=matmul(sminhalf,cmatb)

    eelec=0._dp
    do nu = 1, nb
        do mu = 1, nb
            eelec=eelec+pa(nu,mu)*(hcore(mu,nu)+Fa(mu,nu))
            eelec=eelec+pb(nu,mu)*(hcore(mu,nu)+Fb(mu,nu))
        end do
    end do
    eelec=eelec/2._dp
    deltae=eelec-eeold
    deltapa=maxval(abs(pa-polda))
    deltapb=maxval(abs(pb-poldb))
    etot=eelec+unuc
    write (9,'(I3,4F20.10)') i,etot,deltae,deltapa,deltapb
    if (abs(deltae).lt.conve) then
        if ((deltapa.lt.conv).and.(deltapb.lt.conv)) then
             write (9,*) "SCF Converged"
             exit
        end if
    end if
end do
write (9,*) "Final Alpha Eigenvalues:"
write (9,'(5F12.5)') eigsa
write (9,*) "Final Beta  Eigenvalues:"
write (9,'(5F12.5)') eigsb
if (prtlevl.ge.3) then
  write (9,*) "Final Alpha Eigenvectors:"
  do i=1,nb
    write (9,'(7F12.6)') cmata(i,:)
  end do
  write (9,*) "Final Beta  Eigenvectors:"
  do i=1,nb
    write (9,'(7F12.6)') cmatb(i,:)
  end do
end if
if (prtlevl.ge.4) then
  write (9,*) "Final Alpha Density Matrix:"
  trace=0._dp
  do i=1,nb
    write (9,'(7F12.6)') pa(i,:)
    trace=trace+pa(i,i)
  end do
  write (9,*) "With trace:",trace
  write (9,*) "Final Beta  Density Matrix:"
  trace=0._dp
  do i=1,nb
    write (9,'(7F12.6)') pb(i,:)
    trace=trace+pb(i,i)
  end do
  write (9,*) "With trace:",trace
  write (9,*) "Final Total Density Matrix:"
  p1=pa+pb
  trace=0._dp
  do i=1,nb
    write (9,'(7F12.6)') p1(i,:)
    trace=trace+p1(i,i)
  end do
  write (9,*) "With trace:",trace
! for re-use, make sym d-mat
  oo=(/1,2,4,5,3,6,8,9,7,10,11,12,13,14/)
  do i=1,9
    do j=1,9
      p2(i,j)=p1(oo(i),oo(j))
    end do
  end do
  do i=1,9
    do j=1,9
      p3(i,j)=p2(oo(i),oo(j))
    end do
  end do
  ptot=(p1+p2+p3)/3._dp
  write (9,*) "Symmetrised Final Total Density Matrix:"
  trace=0._dp
  do i=1,nb
    write (9,'(7F12.6)') ptot(i,:)
    trace=trace+ptot(i,i)
  end do
  write (9,*) "With trace:",trace
end if
allocate(cij(nb,nb))
cij=cmata
return

end subroutine u_incore_hartreefock

