program rainflow
!  Downing and Socie, Int J Fatigue 1982

implicit none

integer, parameter :: np = 14000 ! 4096
integer, parameter :: ns = 400
double precision, parameter :: dsa = 0.005d0

integer :: ip,i,j,istart,iread
integer, dimension(ns) :: n,ncum
double precision :: junk,slope,x,y,p,Peq,m,xmean,xrange
double precision, dimension(ns) :: srange,ds
double precision, dimension(2*np) :: e

open (UNIT=48,FILE='cycles.txt',STATUS='REPLACE')
open (UNIT=49,FILE='data.txt',STATUS='OLD',ACTION='READ')

ip = 2
iread = 0
j = 0
istart = 1

read (49,*) junk,e(1)
iread = iread + 1

do                                                !  100
   read (49,*) junk,e(2)
   iread = iread + 1
   if (e(1) .ne. e(2)) then
      exit
   end if
end do                                            !  go to 100

slope = 1.d0
if (e(1) .gt. e(2)) then
   slope = -1.d0
end if

do                                                !  1
   read (49,*) junk,p
   iread = iread + 1
   if (iread .eq. np) then
      exit                                        !  go to 6
   end if
   ip = ip + 1
   slope = slope*(-1.d0)
   e(ip) = p
   do                                             !  2
      if (ip .lt. istart+1) then
         exit                                     !  go to 1
      end if
      x = slope*(e(ip) - e(ip-1))
      if (x .le. 0.d0) then                       !  go to 200
         ip = ip - 1                              !  200
         e(ip) = e(ip+1)
         slope = slope*(-1.d0)                    !  go to 2
      else if (ip .lt. istart+2) then
         exit                                     !  go to 1
      else
         y = slope*(e(ip-2) - e(ip-1))
         if (x .lt. y) then
            exit                                  !  go to 1
         else if (x .eq. y .and. &
                  istart .eq. ip-2) then
            exit                                  !  go to 1
         else if (x .gt. y .and. &
                 istart .eq. ip-2) then
            istart = istart + 1                   !  4
            exit                                  !  go to 1
         else if (x .ge. y .and. &
                  istart .ne. ip-2) then
            xrange = y                             !  5
            xmean = 0.5d0*(e(ip-1) + e(ip-2))
write (48,'(1X,2ES13.4)') xrange,xmean
            ip = ip - 2
            e(ip) = e(ip+2)
         end if                                   !  go to 2
      end if                                      !  go to 2
   end do                                         !  go to 2
end do                                            !  go to 1

do                                                !  6
   j = j + 1
   if (j .gt. istart) then
      exit                                        !  stop
   end if
   ip = ip + 1
   slope = slope*(-1.d0)
   e(ip) = e(j)
   do                                             !  7
      if (ip .lt. istart+1) then
         exit                                     !  go to 6
      end if
      x = slope*(e(ip) - e(ip-1))
      if (x .le. 0.d0) then                       !  go to 300
         ip = ip - 1                              !  300
         e(ip) = e(ip+1)
         slope = slope*(-1.d0)
      else if (ip .lt. istart+2) then
         exit                                     !  go to 6
      else
         y = slope*(e(ip-2) - e(ip-1))
         if (x .lt. y) then                       !  8
            exit                                  !  go to 6
         end if
         xrange = y
         xmean = 0.5d0*(e(ip-1) + e(ip-2))
write (48,'(1X,2ES13.4)') xrange,xmean
         ip = ip - 2
         e(ip) = e(ip+2)
      end if                                      !  go to 7
   end do                                         !  go to 7
end do                                            !  go to 6

close(49)  !  Close input file.
close(48)  !  Close output file.
!  Reopen old output file for reading and sorting.
open (UNIT=49,FILE='cycles.txt',STATUS='OLD',ACTION='READ')
open (UNIT=48,FILE='ncum.txt',STATUS='REPLACE')
open (UNIT=47,FILE='ncyc.txt',STATUS='REPLACE')

!  Set levels.
do i = 1,Ns
   srange(i) = dsa*dble(i)
   ds(i) = dsa
end do

!  Read cycles.
ncum = 0
n = 0
do
   read (49,*,IOSTAT=j) xrange,xmean
   if (j .ne. 0) then
      exit
   end if
   do i = 1,Ns
      if (abs(xrange) .gt. srange(i)) then
         ncum(i) = ncum(i) + 1
      end if
      if ((abs(xrange) .gt. srange(i) - 0.5d0*ds(i)) .and. &
          (abs(xrange) .le. srange(i) + 0.5d0*ds(i))) then
         n(i) = n(i) + 1
      end if
   end do   
end do

do i = Ns,1,-1
   write(48,'(1X,ES12.3,I6)') srange(i),ncum(i)
   write(47,'(1X,ES12.3,I6)') srange(i),n(i)
end do

Peq = 0.d0
m = 10.d0
do i = 1,Ns
  Peq = Peq + (srange(i)**m)*n(i)
end do
Peq = (Peq/300.d0)**(1.d0/m)

print *,Peq

close(47)
close(48)
close(49)

end program