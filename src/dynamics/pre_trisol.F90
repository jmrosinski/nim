! Copyright 2016
! United States Government as represented by the National Oceanic and
! Atmospheric Administration.
! No copyright is claimed in the United States under Title 17, U.S.Code.
! All Other Rights Reserved.
!
! Copyright 2016 Regents of the University of Colorado. All rights
! reserved.
! Copyright 2016 Regents of the Colorado State University. All rights
! reserved.

subroutine pre_trisol(nvars,ims,ime,ips,ipe,wr,tp,trp,ca4k,ca4p,z,zm,tebb,bedgvar)
!**************************************************************************
!       Loops prior to call to trisol()
!       Jin Lee                  December, 2008
!**************************************************************************
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent(IN   ) :: nvars,ims,ime,ips,ipe
real   ,intent(IN   ) :: wr     (0:NZ,    ims:ime      )
real   ,intent(IN   ) :: tp     (  NZ,    ims:ime      )
real   ,intent(IN   ) :: trp    (  NZ,    ims:ime      )
real   ,intent(IN   ) :: ca4k   (  NZ                  )
real   ,intent(IN   ) :: ca4p   (  NZ                  )
real   ,intent(IN   ) :: z      (0:NZ,    ims:ime      )
real   ,intent(IN   ) :: zm     (  NZ,    ims:ime      )
real   ,intent(IN   ) :: tebb   (0:NZ,    ims:ime      )
real   ,intent(INOUT) :: bedgvar(0:NZ,    ims:ime,nvars)
integer ipn,k,kp1

!$acc parallel vector_length(96)

!$OMP PARALLEL DO PRIVATE(k,kp1) SCHEDULE(runtime)
!$acc loop gang
do ipn=ips,ipe
!JR Does this vectorize? It should.
!$acc loop vector
  do k=1,NZ
    kp1=min(NZ,k+1)
    bedgvar(k,ipn,3)=.5_rt*(wr(k,ipn)+wr(kp1,ipn))
    bedgvar(k,ipn,5)=(ca4k(k)*tp(k,ipn)+ca4p(k)*tp(kp1,ipn))
    bedgvar(k,ipn,6)=ca4k(k)* trp(k,ipn)+ca4p(k)* trp(kp1,ipn)
  end do
  bedgvar(0,ipn,3) = .5_rt*(wr(0,ipn)+wr(1  ,ipn) )
  bedgvar(0,ipn,5) = (zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*bedgvar(1,ipn,5) &
                    +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*bedgvar(2,ipn,5)
  bedgvar(0,ipn,6) = (zm(1,ipn)-zm(3,ipn))/(zm(2,ipn)-zm(3,ipn))*bedgvar(1,ipn,6) &
                    +(zm(1,ipn)-zm(2,ipn))/(zm(3,ipn)-zm(2,ipn))*bedgvar(2,ipn,6)
  bedgvar(0,ipn,5)=tebb(0,ipn)+bedgvar(0,ipn,5)
!$acc loop vector
  do k=1,NZ
    bedgvar(k,ipn,5)=tebb(k,ipn)+bedgvar(k,ipn,5)
  end do
enddo !ipn-loop
!$acc end parallel
!$OMP END PARALLEL DO

return
end subroutine pre_trisol
!**************************************************************************
!
!       Nonhydrostatic Icosahedral Model (NIM)
!
!       Design:  Jin-luen Lee and Alexander E. MacDonald (2007-2010)
!       Development Leader:  Jin-luen Lee (2008-2011)
!       Computational Development:  Jacques Middlecoff (2009-2011)
!       GPU Design:  Mark Govett (2009-2011)
!       GPU Development: Mark Govett and Tom Henderson (2009-2011)
!       Documentation:  Alexander E. MacDonald (2010)!
!
!**************************************************************************
