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

subroutine dyntnc

!****************************************************************************
!
!       Modify nim.F90 to calculate dynamics tendency for Runge-Kutta schemes
!       Jin-Luen Lee Feb, 2010
!       Jacques Middlecoff July 2009 Put variables in calling sequences.
!       Jacques Middlecoff May  2010 Split caldiv into caldiv and ApplyDiv
!       Jacques Middlecoff May  2010 Combined exchanges of div and vdns,sa,sedgvar
!       Jacques Middlecoff Sep  2010 Removed exchange of sa.
!
!****************************************************************************

use module_control
use module_variables
use module_constants
use globals, only: myrank
implicit none
#include <gptl.inc>
integer :: ret

             ret = gptlstart('Dyntnc')
             doBarrier = preExchangeBarrier.and.noOverlap
             call preExchange ('DynExch1begPreBar' ,'DynExch1beg')
!sms$exchange_begin (ur,vr,wr,u,v,w,tp,trp,rp)
             doBarrier = postExchangeBarrier.and.noOverlap
             call postExchange('DynExch1begPostBar','DynExch1beg')
             ret = gptlstart('RHStendencies')
call RHStendencies(ims,ime,ips,ipe,spncef,su,sv,tfr,tfur,tfvr,tfwr,tftr)
             ret = gptlstop ('RHStendencies')
             doBarrier = preExchangeBarrier.and.(.not.noOverlap)
             call preExchange ('DynExch1endPreBar' ,'DynExch1end')
!sms$exchange_end(ur,vr,wr,u,v,w,tp,trp,rp)
             doBarrier = postExchangeBarrier.and.(.not.noOverlap)
             call postExchange('DynExch1endPostBar','DynExch1end')
call vdm
             doBarrier = preExchangeBarrier.and.noOverlap
             call preExchange ('DynExch2PreBar' ,'DynExchange2beg')
!sms$exchange_begin (uw8s)
             doBarrier = postExchangeBarrier.and.noOverlap
             call postExchange('Dyn1Exch2PostBar','DynExchange2beg')
             ret = gptlstart('PreForce')
call preForce(ims,ime,ips,ipe,npp,nvars,gamma,rd,cfh,nprox,prox,cs,sn,nvecs,sedgvar,sa,e,saoprx,ur,vr,wr,fu0,fv0,fw0)
             ret = gptlstop ('PreForce')
             doBarrier = preExchangeBarrier.and.(.not.noOverlap)
             call preExchange ('DynExch2PreBar' ,'DynExchange2end')
!sms$exchange_end (uw8s)
             doBarrier = postExchangeBarrier.and.(.not.noOverlap)
             call postExchange('Dyn2Exch2PostBar','DynExchange2end')
             ret = gptlstart('Vdtotal')
             ret = gptlstart('Vdn')
call vdn(npp,ims,ime,ips,ipe,nprox,prox,proxs,nvecs,nvecb,uw8b,uw8s,vdnb,vdns)
             ret = gptlstop ('Vdn')
             ret = gptlstop ('Vdtotal')
             doBarrier = preExchangeBarrier
             call preExchange ('DynExchan3PreBar' ,'DynExch3')
!sms$exchange (vdns)
             doBarrier = postExchangeBarrier
             call postExchange('DynExch3PostBar','DynExch3')
!sms$compare_var(vdnb,"dyntnc.F90 - vdnb1")
!sms$compare_var(vdns,"dyntnc.F90 - vdns1")
!sms$compare_var(sedgvar,"dyntnc.F90 - sedgvar1")
             doBarrier = preExchangeBarrier.and.(.not.noOverlap)
             call preExchange ('VdmExchEndPreBar' ,'VdmExchEnd')
!sms$exchange_end(sedgvar)
             doBarrier = postExchangeBarrier.and.(.not.noOverlap)
             call postExchange('VdmExchEndPostBar','VdmExchEnd')
             ret = gptlstart('Flux')
call flux    (npp,ims,ime,ips,ipe,nvars,nprox,prox,proxs,vdns,vdnb,sa,area,vol,ca4k,ca4p,sedgvar,bedgvar,tfur,tfvr,tfwr,tftr,tefr)
             ret = gptlstop ('Flux')
!sms$compare_var(tfur,"dyntnc.F90 - tfur1")
!sms$compare_var(tfvr,"dyntnc.F90 - tfvr1")
!sms$compare_var(tfwr,"dyntnc.F90 - tfwr1")
!sms$compare_var(tftr,"dyntnc.F90 - tftr1")
!sms$compare_var(tefr,"dyntnc.F90 - tefr1")
             ret = gptlstart('Force')
call force   (ims,ime,ips,ipe,vol,f,ur,vr,fu0,fv0,fw0,tfur,tfvr,tfwr)
             ret = gptlstop ('Force')
!sms$compare_var(tfur,"dyntnc.F90 - tfur2")
!sms$compare_var(tfvr,"dyntnc.F90 - tfvr2")
!sms$compare_var(tfwr,"dyntnc.F90 - tfwr2")
             ret = gptlstop ('Dyntnc')
return
end subroutine dyntnc

!****************************************************************************
!       Nonhydrostatic Icosahedral Model (NIM)
!
!       Design:  Jin-luen Lee and Alexander E. MacDonald (2007-2010)
!       Development Leader:  Jin-luen Lee (2008-2011)
!       Computational Development:  Jacques Middlecoff (2009-2011)
!       GPU Design:  Mark Govett (2009-2011)
!       GPU Development: Mark Govett and Tom Henderson (2009-2011)
!       Documentation:  Alexander E. MacDonald (2010)!
!
!****************************************************************************
