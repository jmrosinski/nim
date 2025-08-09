program DiffStats
implicit none
 integer,parameter :: Nfiles=6,MaxVars=1000
 integer           :: Nfile,Nvar,i
 integer           :: TotPts(MaxVars)
 integer           :: lun1=28,lun2=29,lunNML=10
 integer           :: nz1,glvl1,its1,time1,dim11,nip1
 integer           :: nz2,glvl2,its2,time2,dim12,nip2
 logical           :: PrintHeader
 real              :: digits(MaxVars) = 0.
 real              :: CutoffFactor = .01 !Fraction of abs(var1+var2) below which BinDigits is not counted.
 integer           :: BinDigits(0:8,MaxVars) = 0. !Histogram of matching digits
 real,allocatable  :: var1(:,:),var2(:,:)
 character(len=130):: path1,path2,file1,file2
 character(80)     :: header(10)
 character(13)     :: s1,s2,s3,s4,s5
 character(2)      :: timeunit1, timeunit2
 character(12)     :: yyyymmddhhmm1,yyyymmddhhmm2
 character(6)      :: suffix
 character(4)      :: VarName1,VarName2,VarName(MaxVars)
 character(4)      :: file(Nfiles) = (/'uZZZ','vZZZ','tZZZ','trpZ','pZZZ','wZZZ'/)
 integer :: ios

 namelist /FileInfo/ path1,path2,suffix,CutoffFactor

 open (lunNML,file="DIFFnamelist")
 READ (lunNML,NML=FileInfo)
 close(lunNML)

 PrintHeader = .true.
 print*,'For ',suffix,' comparing ',trim(path1),' with ',trim(path2)
 do Nfile = 1,Nfiles
   file1 = trim(path1) // 'out_' // trim(file(Nfile)) // suffix
   file2 = trim(path2) // 'out_' // trim(file(Nfile)) // suffix
   open (lun1,file=file1,form="unformatted",status='old', iostat=ios)
   if (ios == 0) then
!    write(6,*)'Successfully opened input file=',trim(file1)
   else
     write(6,*)'cannot open file1=', trim(file1)
     stop 999
   end if
   open (lun2,file=file2,form="unformatted",status='old', iostat=ios)
   if (ios == 0) then
!    write(6,*)'Successfully opened input file=',trim(file2)
   else
     write(6,*)'cannot open file1=', trim(file2)
     stop 999
   end if
   do Nvar = 1,MaxVars
      read (lun1,END=10) header
      call parse_hdr_1 (header(1), VarName1, yyyymmddhhmm1)
      call parse_hdr_2 (header(2), nz1, glvl1, timeunit1, its1, time1)
      call parse_hdr_3 (header(3), dim11, nip1)
      allocate( var1(dim11,nip1) )
      read (lun1) var1

      read (lun2) header
      call parse_hdr_1 (header(1), VarName2, yyyymmddhhmm2)
      call parse_hdr_2 (header(2), nz2, glvl2, timeunit2, its2, time2)
      call parse_hdr_3 (header(3), dim12, nip2)

 !    write(6,*)'varnames=', varname1, varname2
 !    write(6,*)'yyyymmddhhmm=', yyyymmddhhmm1, yyyymmddhhmm2
 !    write(6,*)'nz=', nz1, nz2
 !    write(6,*)'glvl=', glvl1, glvl2
 !    write(6,*)'its=', its1, its2
 !    write(6,*)'time=', time1, time2
      call check_header(VarName1,yyyymmddhhmm1,nz1,glvl1,its1,time1,dim11,nip1,&
                        VarName2,yyyymmddhhmm2,nz2,glvl2,its2,time2,dim12,nip2)
      allocate( var2(dim12,nip2) )
      read (lun2) var2
      TotPts(Nfile) = size(var1)
      call diff_vectors(VarName1,var1,var2,size(var1),PrintHeader,CutoffFactor,digits(Nfile),BinDigits(0,Nfile))
      VarName(Nfile) = VarName1
      deallocate(var1,var2)
   enddo
   print*,'Error in DiffStats: Should not be here'
10 continue
   close(lun1)
   close(lun2)
   PrintHeader = .false.
enddo
print*
print"('Histogram of the number of digits that match. CutoffFactor =',f10.7)",CutoffFactor
print"(' Digits        0         1         2         3         4         5         6         7         8   Skipped     Total')"
do Nvar=1,Nfiles
  print"(A6,11i10)",VarName(Nvar),BinDigits(:,Nvar),TotPts(Nvar)-sum(BinDigits(:,Nvar)),TotPts(Nvar)
enddo
print*
print"('Field:  ', 14a6)",(VarName(i),i=1,Nfiles),'  Mean'
print"('Digits: ', 14f6.1)",(digits(i),i=1,Nfiles),sum(digits)/Nfiles

end program DiffStats
