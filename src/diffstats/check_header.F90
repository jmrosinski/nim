subroutine check_header(VarName1,yyyymmddhhmm1,nz1,glvl1,its1,time1,dim11,nip1,&
                        VarName2,yyyymmddhhmm2,nz2,glvl2,its2,time2,dim12,nip2)
 character(4) ,intent(IN) :: VarName1,VarName2
 character(12),intent(IN) :: yyyymmddhhmm1,yyyymmddhhmm2
 integer      ,intent(IN) :: nz1,glvl1,its1,time1,dim11,nip1
 integer      ,intent(IN) :: nz2,glvl2,its2,time2,dim12,nip2

 if(VarName1 /= Varname2) then
   print*,'Variable names are not the same',VarName1,VarName2
 endif
 if(yyyymmddhhmm1 /= yyyymmddhhmm2) then
   print*,'Date-Time is not the same',yyyymmddhhmm1,yyyymmddhhmm2
 endif
 if(nz1 /= nz2) then
   print*,'nz is not the same',nz1,nz2
 endif
 if(glvl1 /= glvl2) then
   print*,'glvl is not the same',glvl1,glvl2
 endif
 if(its1 /= its2) then
   print*,'its is not the same',its1,its2
 endif
 if(time1 /= time2) then
   print*,'time is not the same',time1,time2
 endif
 if(dim11 /= dim12) then
   print*,'dim1 is not the same',dim11,dim12
 endif
 if(nip1 /= nip2) then
   print*,'nip is not the same',nip1,nip2
 endif

end subroutine check_header
