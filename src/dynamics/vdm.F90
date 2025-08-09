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

subroutine vdm
!**************************************************************************
!
!       Obtain values at cellular edges from those at centers
!       Jin Lee                  December, 2008
!
!**************************************************************************

  use module_control
  use module_variables
  use module_constants
  implicit none
#include <gptl.inc>
  integer ret

  ret = gptlstart('Vdtotal')
  ret = gptlstart('Vdm')
  if (vdmint_combine==5) then
! Combine u,v,rp,tp,trp into one call
! note:  Vdmints and Vdmintv times are combined for this option
    ret = gptlstart('Vdmints+Vdmintv')
    ret = gptlstart('Vdmintsv5')
    call vdmintsv5(u,v,rp,tp,trp,teb,sedgvar,bedgvar,nvars,4,amtx1,nbf,npp, &
         ims,ime,nob,nprox,prox,cs,sn,xyc,zc,zm,z,ca4k,ca4p,ips,ipe)
    ret = gptlstop ('Vdmintsv5')
    ret = gptlstart('Vdmints0')
    call vdmints0(w,sedgvar,bedgvar,nvars,3,amtx2,nbf,npp,ims,ime,nob, &
         nprox,prox,xyc,zc,z ,ca4k,ca4p,ips,ipe)
    ret = gptlstop('Vdmints0')
    ret = gptlstop ('Vdmints+Vdmintv')
  else
! Combine u,v into one call
! Combine rp,tp,trp into a second call
    ret = gptlstart('Vdmints+Vdmintv')
    ret = gptlstart('Vdmintv')
    call vdmintv(u,v,sedgvar,bedgvar,nvars,amtx1,nbf,npp,ims,ime,nob, &
         nprox,prox,cs,sn,xyc,zc,zm,ca4k,ca4p,z,ips,ipe)
    ret = gptlstop ('Vdmintv')
    ret = gptlstart('Vdmints0')
    call vdmints0(w,sedgvar,bedgvar,nvars,3,amtx2,nbf,npp,ims,ime, &
         nob,nprox,prox,xyc,zc,z ,ca4k,ca4p,ips,ipe)
    ret = gptlstop('Vdmints0')
    ret = gptlstart('Vdmints3')
    call vdmints3(rp,tp,trp,teb,sedgvar,bedgvar,nvars,4,amtx1,nbf,npp, &
         ims,ime,nob,nprox,prox,xyc,zc,zm,ca4k,ca4p,ips,ipe)
    ret = gptlstop ('Vdmints3')
    ret = gptlstop ('Vdmints+Vdmintv')
  endif
!sms$compare_var(bedgvar,"vdm.F90 - bedgvar5")
!sms$compare_var(sedgvar,"vdm.F90 - sedgvar5")
  doBarrier = preExchangeBarrier.and.noOverlap
  call preExchange ('VdmExchbegPreBar' ,'VdmExchangebeg')
!sms$exchange_begin (sedgvar)
  doBarrier = postExchangeBarrier.and.noOverlap
  call postExchange('VdmExchbegPostBar','VdmExchangebeg')
  ret = gptlstart('vdmfinish')
  call vdmfinish(npp,nvars,ims,ime,ips,ipe,nprox,rebb,tebb,reb, &
                ca4k,ca4p,bedgvar,sedgvar,uw8s,uw8b)
  ret = gptlstop('vdmfinish')
  ret = gptlstop ('Vdm')
  ret = gptlstop ('Vdtotal')
!sms$compare_var(bedgvar,"vdm.F90 - bedgvar6")
!sms$compare_var(sedgvar,"vdm.F90 - sedgvar6")
!sms$compare_var(uw8b,"vdm.F90 - uw8b6")
!sms$compare_var(uw8s,"vdm.F90 - uw8s6")

  return
end subroutine vdm
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
