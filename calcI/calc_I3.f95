module constants
        integer,parameter::n=129
        real(kind=8),parameter::pi=3.1415926
        real(kind=8),parameter::mu=4*pi*0.0000001
end module constants

program calcI
use constants
implicit none
real(kind=8)::start,finish1,finish11,finish2,finish3,finish,inte1,BR,BZ,RRR,ZZZ
integer::i,j,k,numm,flag2,l,rr,u,d
real(kind=8),dimension(n,n)::R,Z,psi,q,Ic,B_R,B_Z
integer,dimension(n,n)::num,order
real(kind=8),dimension(n,n,0:10000)::Rm,Zm,psim,dll,nablapsi
!Initialization
dll(:,:,:)=0
nablapsi(:,:,:)=0
call cpu_time(start)
open(11,file='q.plt')
do j=1,n
        do i=1,n
        read(11,*)R(i,j),Z(i,j),psi(i,j),q(i,j)
        enddo
enddo
close(11)
psim(:,:,:)=0;Rm(:,:,:)=0;Zm(:,:,:)=0
do i=1,n
do j=1,n
Rm(i,j,0)=R(i,j);Zm(i,j,0)=Z(i,j);psim(i,j,0)=psi(i,j)
enddo
enddo

!calculate BR,BZ

do i=1,n
   do j=1,n
      if(i>1.and.i<n)then
         B_Z(i,j)=(psi(i+1,j)-psi(i-1,j))/(2*pi*R(i,j)*(R(i+1,j)-R(i-1,j)))
      elseif(i==1)then
         B_Z(i,j)=(psi(i+1,j)-psi(i,j))/(2*pi*R(i,j)*(R(i+1,j)-R(i,j)))
      elseif(i==n)then
         B_Z(i,j)=(psi(i,j)-psi(i-1,j))/(2*pi*R(i,j)*(R(i,j)-R(i-1,j)))
      endif
      if(j>1.and.j<n)then
         B_R(i,j)=-(psi(i,j+1)-psi(i,j-1))/(2*pi*R(i,j)*(Z(i,j+1)-Z(i,j-1)))
      elseif(j==1)then
         B_R(i,j)=-(psi(i,j+1)-psi(i,j))/(2*pi*R(i,j)*(Z(i,j+1)-Z(i,j)))
      elseif(j==n)then
         B_R(i,j)=-(psi(i,j)-psi(i,j-1))/(2*pi*R(i,j)*(Z(i,j)-Z(i,j-1)))
      endif
   enddo
enddo

call findflux(R,Z,B_R,B_Z,psi,Rm,Zm,psim,order,numm,num)
call cpu_time(finish1)

do i=1,n
do j=1,n
        if(order(i,j)<=numm.and.0<order(i,j))then
        do k=0,num(i,j)
           RRR=Rm(i,j,k)
           ZZZ=Zm(i,j,k)
           call calc_BRBZ(RRR,ZZZ,R,Z,B_R,B_Z,psi,Rm,Zm,Br,Bz,l,rr,u,d)
           nablapsi(i,j,k)=2*pi*Rm(i,j,k)*sqrt(BR**2+BZ**2)

           !          call calc_nablapsi(i,j,k,R,Z,psi,Rm,Zm,psim,nablapsi,order,numm,num)
        enddo
        endif
enddo
enddo
call cpu_time(finish2)
do i=1,n
do j=1,n
        if(order(i,j)<=numm.and.0<order(i,j))then
        do k=0,num(i,j)
                call calc_dl(i,j,k,R,Z,num,Rm,Zm,dll)
        enddo
        endif
enddo
enddo
print*,'calc I'
Ic(:,:)=0
do i=1,n
do j=1,n
if(0<order(i,j).and.order(i,j)<=numm)then
inte1=0
        do k=0,num(i,j)
        inte1=inte1+dll(i,j,k)/(nablapsi(i,j,k)*Rm(i,j,k))
        end do
        !print*,'order=',order(i0,j0),psi(i0,j0),'inte=',inte1,i0,j0
        Ic(i,j)=2*pi*q(i,j)/(inte1*mu)
endif
enddo
enddo
call cpu_time(finish3)
open(66,file='I.plt')
do k=1,numm-1
flag2=0
do i=1,n
do j=1,n
        if(order(i,j)==k.and.Ic(i,j)<20000000)then
                write(66,*)i,j,psi(i,j),q(i,j),Ic(i,j)
                flag2=1
                exit
        endif
enddo
if(flag2==1)then
exit
endif
enddo
enddo
close(66)
call cpu_time(finish)
print*,'Findflux Time=',finish1-start,'seconds'
print*,'calc_nabla time',finish2-finish1,'seconds'
print*,'Time=',finish-start,'seconds'
end program calcI

subroutine findflux(R,Z,B_R,B_Z,psi,Rm,Zm,psim,order,numm,num)
use constants
implicit none
integer::i0,j0,i,j,k,ii,jj,kk,iii,jjj,kkk,flag1,numm,flag2
real(kind=8)::min11,psi_min,psi_max,f1,f2,f3,f4
real(kind=8),dimension(n,n)::R,Z,psi,B_R,B_Z
integer,dimension(n,n)::num,m,order
character(len=3)::filename1,filename2
real(kind=8),dimension(n,n,0:10000)::Rm,Zm,psim
real(kind=8)::psi_quick(10001)
!initialization
numm=0
psi_quick(:)=0
num(:,:)=0
m(:,:)=0
order(:,:)=0
psi_min=psi(1,1);psi_max=psi(1,1)
do i=1,n
do j=1,n
if(psi(i,j)<psi_min)then
psi_min=psi(i,j)
endif
if(psi(i,j)>psi_max)then
psi_max=psi(i,j)
endif
enddo
enddo
!do i=50,70
do i=1,n
!        do j=50,70
        do j=1,n
        flag1=1
!        if(i==40.and.j==60)then
        call findfluxI(i,j,R,Z,B_R,B_Z,psi,Rm,Zm,psim,num,flag1,psi_min,psi_max)
!        endif
        if(flag1==1)then
        numm=numm+1
        m(i,j)=numm
        psi_quick(numm)=psi(i,j)
        endif
        enddo
enddo
print*,'numm=',numm
i=40;j=60
open(1111,file='psim24060.plt')
do k=0,num(40,60)
write(1111,*)Rm(i,j,k),Zm(i,j,k),psim(i,j,k)
enddo
close(1111)
!ci mian paixu
i=1
call quick_sort(psi_quick,i,numm)
!print*,psi_quick(20:40)
do k=1,numm
        flag2=0
        do i=1,n
        do j=1,n
                if(psi_quick(k)==psi(i,j))then
                        order(i,j)=k;flag2=1
                        exit
                endif
        enddo
                if(flag2==1)then
                exit
                endif
        enddo
enddo
!!ci mian dian shu
open(22,file='num.dat')
!do k=1,numm
!flag2=0
do i=1,n
do j=1,n
if(0<order(i,j).and.order(i,j)<=k)then
write(22,*)'i=',i,',j=',j,',num=',num(i,j),order(i,j),psi(i,j)
                write(filename1,'(i3)')k 
                open(23,file='./psim/psim'//trim(adjustl(filename1))//'.plt')
                do kk=0,num(i,j)
                        write(23,*)Rm(i,j,kk),Zm(i,j,kk),psim(i,j,kk)
                enddo
                close(23)
!flag2=1
!exit
endif
enddo
!if(flag2==1)then
!exit
!endif
enddo
!enddo
close(22)
end subroutine findflux

recursive subroutine  quick_sort(a,l,r)
use constants
implicit none
integer::i,j,l,r
real(kind=8),dimension(n,n)::psi
real(kind=8)::a(10001),x
if(l<r)then
        i=l;j=r;x=a(l)
        do while(i<j)
                do while(i<j.and.a(j)>=x)
                        j=j-1
                enddo
                if(i<j)then
                        a(i)=a(j)
                        i=i+1
                endif
                do while(i<j.and.a(i)<x)
                        i=i+1
                enddo
                if(i<j)then
                        a(j)=a(i)
                j=j-1
                endif
        enddo
        a(i)=x
        call quick_sort(a,l,i-1)
        call quick_sort(a,i+1,r)
endif
end subroutine quick_sort

subroutine findfluxI(i0,j0,R,Z,B_R,B_Z,psi,Rm,Zm,psim,num,flag1,psi_min,psi_max)
use constants
implicit none
integer::i0,j0,i,j,k,l,rr,u,d,flag1
real(kind=8)::ddl,Br,Bz,dl,Br1,Br2,Br3,Br4,Bz1,Bz2,Bz3,Bz4,rrr,zzz,k1,k2,k3,k4,temp1,temp2,psi_min,psi_max,dl_min,dl_max
real(kind=8),dimension(n,n)::R,Z,psi,B_R,B_Z
integer,dimension(n,n)::num
real(kind=8),dimension(n,n,0:10000)::Rm,Zm,psim
dl_max=0.001
dl_min=0.01
ddl=dl_max+(psi(i0,j0)-psi_max)/(psi_min-psi_max)*(dl_min-dl_max)
temp1=0;temp2=0;
k=0
do while(dl(i0,j0,k,R,Z,Rm,Zm)>=ddl.or.num(i0,j0)<10)
        !calculate k1
        if(num(i0,j0)==0)then
        RRR=R(i0,j0);zzz=Z(i0,j0);
        else if(num(i0,j0)>0)then
        RRR=Rm(i0,j0,num(i0,j0));ZZZ=Zm(i0,j0,num(i0,j0))
        endif
        call calc_BRBZ(RRR,ZZZ,R,Z,B_R,B_Z,psi,Rm,Zm,Br,Bz,l,rr,u,d)
        BR1=BR;BZ1=Bz;k1=Bz1/Br1;
        !calculate k2
        if(l<rr.and.d<u)then
                if(num(i0,j0)==0)then
                RRR=R(i0,j0)+ddl*Br1/(2*sqrt(Br1**2+Bz1**2))
                ZZZ=Z(i0,j0)+ddl*Bz1/(2*sqrt(Br1**2+Bz1**2))
                elseif(num(i0,j0)>0)then
                RRR=Rm(i0,j0,num(i0,j0))+ddl*Br1/(2*sqrt(Br1**2+Bz1**2))
                ZZZ=Zm(i0,j0,num(i0,j0))+ddl*Bz1/(2*sqrt(Br1**2+Bz1**2))
                endif
        else
        flag1=0
        endif
        call calc_BRBZ(RRR,ZZZ,R,Z,B_R,B_Z,psi,Rm,Zm,Br,Bz,l,rr,u,d)
        BR2=BR;BZ2=Bz;k2=Bz2/Br2;
        !calculate k3
        if(l<rr.and.d<u)then
                if(num(i0,j0)==0)then
                RRR=R(i0,j0)+ddl*Br2/(2*sqrt(Br2**2+Bz2**2))
                ZZZ=Z(i0,j0)+ddl*Bz2/(2*sqrt(Br2**2+Bz2**2))
                elseif(num(i0,j0)>0)then
                RRR=Rm(i0,j0,num(i0,j0))+ddl*Br2/(2*sqrt(Br2**2+Bz2**2))
                ZZZ=Zm(i0,j0,num(i0,j0))+ddl*Bz2/(2*sqrt(Br2**2+Bz2**2))
                endif
        else
        flag1=0
        endif
        call calc_BRBZ(RRR,ZZZ,R,Z,B_R,B_Z,psi,Rm,Zm,Br,Bz,l,rr,u,d)
        BR3=BR;BZ3=Bz;k3=Bz3/Br3;
        !calculate k4
        if(l<rr.and.d<u)then
                if(num(i0,j0)==0)then
                RRR=R(i0,j0)+ddl*Br3/(sqrt(Br3**2+Bz3**2))
                ZZZ=Z(i0,j0)+ddl*Bz3/(sqrt(Br3**2+Bz3**2))
                elseif(num(i0,j0)>0)then
                RRR=Rm(i0,j0,num(i0,j0))+ddl*Br3/(sqrt(Br3**2+Bz3**2))
                ZZZ=Zm(i0,j0,num(i0,j0))+ddl*Bz3/(sqrt(Br3**2+Bz3**2))
                endif
        else
        flag1=0
        endif
        call calc_BRBZ(RRR,ZZZ,R,Z,B_R,B_Z,psi,Rm,Zm,Br,Bz,l,rr,u,d)
        BR4=BR;BZ4=Bz;k4=Bz4/Br4;
        if(l<rr.and.u>d)then
        else
        flag1=0
        endif

        if(flag1==1)then
                num(i0,j0)=num(i0,j0)+1
                k=num(i0,j0)
                temp1=(Bz1/sqrt(Br1**2+Bz1**2)+2*Bz2/(sqrt(Br2**2+Bz2**2))+2*Bz3/(sqrt(Br3**2+Bz3**2))+Bz4/sqrt(Br4**2+Bz4**2))/6.0
                temp2=(Br1/sqrt(Br1**2+Bz1**2)+2*Br2/(sqrt(Br2**2+Bz2**2))+2*Br3/(sqrt(Br3**2+Bz3**2))+Br4/sqrt(Br4**2+Bz4**2))/6.0
                if(num(i0,j0)>3000)then
                exit
                else if(num(i0,j0)>1.and.num(i0,j0)<=3000)then
                Zm(i0,j0,num(i0,j0))=Zm(i0,j0,num(i0,j0)-1)+ddl*temp1
                Rm(i0,j0,num(i0,j0))=Rm(i0,j0,num(i0,j0)-1)+ddl*temp2
                elseif(num(i0,j0)==1)then
                Zm(i0,j0,num(i0,j0))=Z(i0,j0)+ddl*temp1
                Rm(i0,j0,num(i0,j0))=R(i0,j0)+ddl*temp2
                endif
                psim(i0,j0,num(i0,j0))=psi(i0,j0)
        elseif(flag1==0)then
        exit
        endif
        Rm(i0,j0,num(i0,j0))=int(Rm(i0,j0,num(i0,j0))*10**9)
        Rm(i0,j0,num(i0,j0))=Rm(i0,j0,num(i0,j0))/(1000000000.0)
        Zm(i0,j0,num(i0,j0))=int(Zm(i0,j0,num(i0,j0))*10**9)
        Zm(i0,j0,num(i0,j0))=Zm(i0,j0,num(i0,j0))/(1000000000.0)
        psim(i0,j0,num(i0,j0))=int(psim(i0,j0,num(i0,j0))*10**9)
        psim(i0,j0,num(i0,j0))=psim(i0,j0,num(i0,j0))/(1000000000.0)
enddo
if(num(i0,j0)>5000.or.num(i0,j0)<50)then
flag1=0
endif
end subroutine findfluxI

subroutine calc_BRBZ(RRR,ZZZ,R,Z,B_R,B_Z,psi,Rm,Zm,Br,Bz,l,rr,u,d)
use constants
implicit none
integer::i,j,ii,jj,l,rr,u,d
real(kind=8)::Br,Bz,brld,brrd,brru,brlu,bzld,bzlu,bzrd,bzru,ds1,ds2,ds3,ds4,RRR,ZZZ
real(kind=8),dimension(n,n)::R,Z,psi,B_R,B_Z
real(kind=8),dimension(n,n,0:10000)::Rm,Zm
integer,dimension(n,n)::num
l=0;rr=0;u=0;d=0;Br=0;Bz=0;
!find lr
do i=1,n-1
        if(RRR>=R(i,1).and.RRR<R(i+1,1))then
        l=i;rr=i+1;
        endif
enddo
!find ud
do j=1,n-1
        if(ZZZ>=Z(1,j).and.ZZZ<Z(1,j+1))then
        d=j;u=j+1;
        endif
enddo
if(l<rr.and.u>d)then
ds1=RRR-R(l,d)
ds2=R(rr,d)-RRR
ds3=ZZZ-Z(rr,d)
ds4=Z(rr,u)-ZZZ
!Br=(ds1*ds4*Brrd+ds4*ds2*Brld+ds3*ds1*Brru+ds3*ds2*Brlu)/((ds3+ds4)*(ds1+ds2))
!Bz=(ds1*ds4*Bzrd+ds4*ds2*Bzld+ds3*ds1*Bzru+ds3*ds2*Bzlu)/((ds3+ds4)*(ds1+ds2))
Br=(ds1*ds4*B_R(rr,d)+ds4*ds2*B_R(l,d)+ds3*ds1*B_R(rr,u)+ds3*ds2*B_R(l,u))/((ds3+ds4)*(ds1+ds2))
Bz=(ds1*ds4*B_Z(rr,d)+ds4*ds2*B_Z(l,d)+ds3*ds1*B_Z(rr,u)+ds3*ds2*B_Z(l,u))/((ds3+ds4)*(ds1+ds2))

endif
end subroutine


function dl(i0,j0,k0,R,Z,Rm,Zm)
use constants
implicit none
integer::i0,j0,k0
real(kind=8)::dl,l1,l2
real(kind=8),dimension(n,n,0:10000)::Rm,Zm
real(kind=8),dimension(n,n)::R,Z
if(k0==0)then
l1=0;l2=0
else
l1=Rm(i0,j0,k0)-R(i0,j0)
l2=Zm(i0,j0,k0)-Z(i0,j0)
endif
dl=sqrt(l1**2+l2**2)
end function dl

subroutine calc_nablapsi(i0,j0,k0,R,Z,psi,Rm,Zm,psim,nablapsi,order,numm,num)
use constants
implicit none
integer::i0,j0,k0,numm
real(kind=8),dimension(n,n)::R,Z,psi
integer,dimension(n,n)::order,num
real(kind=8),dimension(n,n,0:10000)::Rm,Zm,psim,nablapsi,dpsidr,dpsidz
call dpsi_dr(i0,j0,k0,R,Z,psi,Rm,Zm,psim,dpsidr,order,numm,num)
call dpsi_dz(i0,j0,k0,R,Z,psi,Rm,Zm,psim,dpsidz,order,numm,num)
!if(order(i0,j0)==113.or.order(i0,j0)==114)then
!print*,order(i0,j0),k0,dpsidr(i0,j0,k0),dpsidz(i0,j0,k0)
!endif
nablapsi(i0,j0,k0)=sqrt(dpsidr(i0,j0,k0)**2+dpsidz(i0,j0,k0)**2)
end subroutine calc_nablapsi

subroutine dpsi_dr(i0,j0,k0,R,Z,psi,Rm,Zm,psim,dpsidr,order,numm,num)
use constants
implicit none
real(kind=8)::dpsi,dr,temp
integer::i0,j0,k0,i,j,k,l,h,nexti,nextj,numm,number1,flag2,flag3
integer,dimension(2)::ll,hh
real(kind=8),dimension(n,n)::R,Z,psi
integer,dimension(n,n)::order,num
real(kind=8),dimension(n,n,0:10000)::Rm,Zm,psim,dpsidr
number1=0;nexti=0;nextj=0;
flag2=0
flag3=0
l=0;h=0
if(100<order(i0,j0))then
!xun zhao xiang lin ci mian
flag2=0
        do i=1,n
        do j=1,n
        if(order(i,j)==order(i0,j0)-100)then
        nexti=i;nextj=j;
        flag2=1
        exit
        endif
        enddo
        if(flag2==1)then
        exit
        endif
        enddo
flag2=0
        !gu ding z zhi,
        do k=0,num(nexti,nextj)
                flag2=0;flag3=0;
                if(0<=k.and.k<num(nexti,nextj))then
                        if(Zm(i0,j0,k0)>=Zm(nexti,nextj,k).and.Zm(i0,j0,k0)<Zm(nexti,nextj,k+1))then
                        flag2=1
                        endif
                        if(Zm(i0,j0,k0)>Zm(nexti,nextj,k+1).and.Zm(i0,j0,k0)<=Zm(nexti,nextj,k))then
                        flag3=1
                        endif
                        if(flag2==1.or.flag3==1)then
                        number1=number1+1
                        ll(number1)=k;hh(number1)=k+1
                        l=k;h=k+1
                        if(number1==2)then
                                if(abs(Rm(nexti,nextj,ll(1))-Rm(i0,j0,k0))<abs(Rm(nexti,nextj,ll(2))-Rm(i0,j0,k0)))then
                                l=ll(1);h=hh(1)
                                else
                                ll(1)=ll(2);hh(1)=hh(2)
                                l=ll(1);h=hh(1)
                                endif
                                number1=number1-1
                        endif
                        endif
                else if(k==num(nexti,nextj))then
                        if(Zm(i0,j0,k0)>=Zm(nexti,nextj,k).and.Zm(i0,j0,k0)<Z(nexti,nextj))then
                        flag2=1
                        endif
                        if(Zm(i0,j0,k0)>Z(nexti,nextj).and.Zm(i0,j0,k0)<=Zm(nexti,nextj,k))then
                        flag3=1
                        endif
                        if(flag2==1.or.flag3==1)then
                        number1=number1+1
                        ll(number1)=k;hh(number1)=0
                        l=k;h=0;
                        if(number1==2)then
                                if(abs(Rm(nexti,nextj,ll(1))-Rm(i0,j0,k0))<abs(Rm(nexti,nextj,ll(2))-Rm(i0,j0,k0)))then
                                l=ll(1);h=hh(1)
                                else
                                ll(1)=ll(2);hh(1)=hh(2)
                                l=ll(1);h=hh(1)
                                endif
                                number1=number1-1
                        endif
                        endif
                        
                endif
        enddo
endif
!线性差分，需要提高精度
temp=(Rm(nexti,nextj,l)-Rm(nexti,nextj,h))*(Zm(nexti,nextj,l)-Zm(i0,j0,k0))/(Zm(nexti,nextj,h)-Zm(nexti,nextj,l))
dr=Rm(i0,j0,k0)-(Rm(nexti,nextj,l)+temp)
dpsi=psi(i0,j0)-psi(nexti,nextj)
dpsidr(i0,j0,k0)=dpsi/dr
!if(order(i0,j0)==42.and.k0==117)then
!print*,k0,dpsi,dr,dpsidr(i0,j0,k0)
!print*,Zm(nexti,nextj,l),Zm(nexti,nextj,h),l,h
!endif
end subroutine dpsi_dr

subroutine dpsi_dz(i0,j0,k0,R,Z,psi,Rm,Zm,psim,dpsidz,order,numm,num)
use constants
implicit none
real(kind=8)::dpsi,dz,RR,ZZ,temp
integer::i0,j0,k0,i,j,k,l,h,nexti,nextj,numm,number1,flag2,flag3
integer,dimension(2)::ll,hh
real(kind=8),dimension(n,n)::R,Z,psi
integer,dimension(n,n)::order,num
real(kind=8),dimension(n,n,0:10000)::Rm,Zm,psim,dpsidz
number1=0
if(50<order(i0,j0))then
!xun zhao xiang lin ci mian
flag2=0
        do i=1,n
        do j=1,n
        if(order(i,j)==order(i0,j0)-50)then
        nexti=i;nextj=j;
        flag2=1
        exit
        endif
        enddo
        if(flag2==1)then
        exit
        endif
        enddo
        !gu ding r zhi
        do k=0,num(nexti,nextj)
                flag2=0;flag3=0;
                if(0<=k.and.k<num(nexti,nextj))then
                        if(Rm(i0,j0,k0)>=Rm(nexti,nextj,k).and.Rm(i0,j0,k0)<Rm(nexti,nextj,k+1))then
                        flag2=1
                        endif
                        if(Rm(i0,j0,k0)>=Rm(nexti,nextj,k+1).and.Rm(i0,j0,k0)<Rm(nexti,nextj,k))then
                        flag3=1
                        endif
                        if(flag2==1.or.flag3==1)then
                        number1=number1+1
                        ll(number1)=k;hh(number1)=k+1
                        l=k;h=k+1
                        if(number1==2)then
                                if(abs(Zm(nexti,nextj,ll(1))-Zm(i0,j0,k0))<abs(Zm(nexti,nextj,ll(2))-Zm(i0,j0,k0)))then
                                l=ll(1);h=hh(1)
                                else
                                ll(1)=ll(2);hh(1)=hh(2)
                                l=ll(1);h=hh(1)
                                endif
                                number1=number1-1
                        endif
                        endif
                else if(k==num(nexti,nextj))then
                        if(Rm(i0,j0,k0)>=Rm(nexti,nextj,k).and.Rm(i0,j0,k0)<R(nexti,nextj))then
                        flag2=1
                        endif
                        if(Rm(i0,j0,k0)>=R(nexti,nextj).and.Rm(i0,j0,k0)<Rm(nexti,nextj,k))then
                        flag3=1
                        endif
                        if(flag2==1.or.flag3==1)then
                        number1=number1+1
                        ll(number1)=k;hh(number1)=0
                        l=k;h=0
                        if(number1==2)then
                                if(abs(Zm(nexti,nextj,ll(1))-Zm(i0,j0,k0))<abs(Zm(nexti,nextj,ll(2))-Zm(i0,j0,k0)))then
                                l=ll(1);h=hh(1)
                                else
                                ll(1)=ll(2);hh(1)=hh(2)
                                l=ll(1);h=hh(1)
                                endif
                                number1=number1-1
                        endif
                        endif
                endif
        enddo
endif
!线性差分，需要提高精度
temp=(Zm(nexti,nextj,l)-Zm(nexti,nextj,h))*(Rm(nexti,nextj,l)-Rm(i0,j0,k0))/(Rm(nexti,nextj,h)-Rm(nexti,nextj,l))
dz=Zm(i0,j0,k0)-(Zm(nexti,nextj,l)+temp)
dpsi=psi(i0,j0)-psi(nexti,nextj)
dpsidz(i0,j0,k0)=dpsi/dz
end subroutine dpsi_dz

subroutine calc_dl(i0,j0,k0,R,Z,num,Rm,Zm,dll)
use constants
implicit none
integer::i0,j0,k0
real(kind=8)::trans(2,2),cos1,sin1,ds,m1,m2
real(kind=8),dimension(n,n,0:10000)::Rm,Zm
real(kind=8),dimension(n,n)::R,Z
integer,dimension(n,n)::num
real(kind=8),dimension(n,n,0:10000)::dll
real(kind=8)::k1(2,1),k2(2,1),k3(2,1),k4(2,1)
real(kind=8)::k10(2,1),k20(2,1),k30(2,1),k40(2,1)
if(k0<num(i0,j0)-1.and.k0>=1)then
!计算m1,m2
ds=sqrt((Rm(i0,j0,k0)-Rm(i0,j0,k0+1))**2+(Zm(i0,j0,k0)-Zm(i0,j0,k0+1))**2)
cos1=(Rm(i0,j0,k0+1)-Rm(i0,j0,k0))/ds
sin1=(Zm(i0,j0,k0+1)-Zm(i0,j0,k0))/ds
trans(1,1)=cos1
trans(1,2)=sin1
trans(2,1)=-sin1
trans(2,2)=cos1
k10(1,1)=Rm(i0,j0,k0-1)-Rm(i0,j0,k0)
k10(2,1)=Zm(i0,j0,k0-1)-Zm(i0,j0,k0)
k20(1,1)=Rm(i0,j0,k0)-Rm(i0,j0,k0)
k20(2,1)=Zm(i0,j0,k0)-Zm(i0,j0,k0)
k30(1,1)=Rm(i0,j0,k0+1)-Rm(i0,j0,k0)
k40(2,1)=Zm(i0,j0,k0+2)-Zm(i0,j0,k0)
k40(1,1)=Rm(i0,j0,k0+2)-Rm(i0,j0,k0)
k30(2,1)=Zm(i0,j0,k0+1)-Zm(i0,j0,k0)
k1=matmul(trans,k10)
k2=matmul(trans,k20)
k3=matmul(trans,k30)
k4=matmul(trans,k40)
m1=(k3(2,1)-k1(2,1))/(k3(1,1)-k1(1,1))
m2=(k4(2,1)-k2(2,1))/(k4(1,1)-k2(1,1))
dll(i0,j0,k0)=ds*(1+(2*m1**2+2*m2**2-m1*m2)/30.0)
else if(k0==num(i0,j0))then
ds=sqrt((Rm(i0,j0,k0)-R(i0,j0))**2+(Zm(i0,j0,k0)-Z(i0,j0))**2)
cos1=(R(i0,j0)-Rm(i0,j0,k0))/ds
sin1=(Z(i0,j0)-Zm(i0,j0,k0))/ds
trans(1,1)=cos1
trans(1,2)=sin1
trans(2,1)=-sin1
trans(2,2)=cos1
k10(1,1)=Rm(i0,j0,k0-1)-Rm(i0,j0,k0)
k10(2,1)=Zm(i0,j0,k0-1)-Zm(i0,j0,k0)
k20(1,1)=Rm(i0,j0,k0)-Rm(i0,j0,k0)
k20(2,1)=Zm(i0,j0,k0)-Zm(i0,j0,k0)
k30(1,1)=R(i0,j0)-Rm(i0,j0,k0)
k40(2,1)=Zm(i0,j0,1)-Zm(i0,j0,k0)
k40(1,1)=Rm(i0,j0,1)-Rm(i0,j0,k0)
k30(2,1)=Z(i0,j0)-Zm(i0,j0,k0)
k1=matmul(trans,k10)
k2=matmul(trans,k20)
k3=matmul(trans,k30)
k4=matmul(trans,k40)
m1=(k3(2,1)-k1(2,1))/(k3(1,1)-k1(1,1))
m2=(k4(2,1)-k2(2,1))/(k4(1,1)-k2(1,1))
dll(i0,j0,k0)=ds*(1+(2*m1**2+2*m2**2-m1*m2)/30.0)
else if(k0==num(i0,j0)-1)then
ds=sqrt((Rm(i0,j0,k0)-Rm(i0,j0,k0+1))**2+(Zm(i0,j0,k0)-Zm(i0,j0,k0+1))**2)
cos1=(Rm(i0,j0,k0+1)-Rm(i0,j0,k0))/ds
sin1=(Zm(i0,j0,k0+1)-Zm(i0,j0,k0))/ds
trans(1,1)=cos1
trans(1,2)=sin1
trans(2,1)=-sin1
trans(2,2)=cos1
k10(1,1)=Rm(i0,j0,k0-1)-Rm(i0,j0,k0)
k10(2,1)=Zm(i0,j0,k0-1)-Zm(i0,j0,k0)
k20(1,1)=Rm(i0,j0,k0)-Rm(i0,j0,k0)
k20(2,1)=Zm(i0,j0,k0)-Zm(i0,j0,k0)
k30(1,1)=Rm(i0,j0,k0+1)-Rm(i0,j0,k0)
k40(2,1)=Z(i0,j0)-Zm(i0,j0,k0)
k40(1,1)=R(i0,j0)-Rm(i0,j0,k0)
k30(2,1)=Zm(i0,j0,k0+1)-Zm(i0,j0,k0)
k1=matmul(trans,k10)
k2=matmul(trans,k20)
k3=matmul(trans,k30)
k4=matmul(trans,k40)
m1=(k3(2,1)-k1(2,1))/(k3(1,1)-k1(1,1))
m2=(k4(2,1)-k2(2,1))/(k4(1,1)-k2(1,1))
dll(i0,j0,k0)=ds*(1+(2*m1**2+2*m2**2-m1*m2)/30.0)
else if(k0==0)then
ds=sqrt((R(i0,j0)-Rm(i0,j0,k0+1))**2+(Z(i0,j0)-Zm(i0,j0,k0+1))**2)
cos1=(Rm(i0,j0,k0+1)-R(i0,j0))/ds
sin1=(Zm(i0,j0,k0+1)-Z(i0,j0))/ds
trans(1,1)=cos1
trans(1,2)=sin1
trans(2,1)=-sin1
trans(2,2)=cos1
k10(1,1)=Rm(i0,j0,num(i0,j0))-R(i0,j0)
k10(2,1)=Zm(i0,j0,num(i0,j0))-Z(i0,j0)
k20(1,1)=R(i0,j0)-R(i0,j0)
k20(2,1)=Z(i0,j0)-Z(i0,j0)
k30(1,1)=Rm(i0,j0,k0+1)-R(i0,j0)
k40(2,1)=Zm(i0,j0,k0+2)-Z(i0,j0)
k40(1,1)=Rm(i0,j0,k0+2)-R(i0,j0)
k30(2,1)=Zm(i0,j0,k0+1)-Z(i0,j0)
k1=matmul(trans,k10)
k2=matmul(trans,k20)
k3=matmul(trans,k30)
k4=matmul(trans,k40)
m1=(k3(2,1)-k1(2,1))/(k3(1,1)-k1(1,1))
m2=(k4(2,1)-k2(2,1))/(k4(1,1)-k2(1,1))
dll(i0,j0,k0)=ds*(1+(2*m1**2+2*m2**2-m1*m2)/30.0)
endif
end subroutine calc_dl
