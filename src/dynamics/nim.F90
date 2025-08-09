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

!******************************************************************************
!       Nonhydrostatic Icosahedral Model (NIM)
!
!       Design:  Jin-luen Lee and Alexander E. MacDonald (2007-2010)
!       Development Leader:  Jin-luen Lee (2008-2011)
!       Computational Development:  Jacques Middlecoff (2009-2011)
!       GPU Design:  Mark Govett (2009-2011)
!       GPU Development: Mark Govett and Tom Henderson (2009-2011)
!       Documentation:  Alexander E. MacDonald (2010)!
!
!       NIM is an ultra-high resolution non-hydrosatic global model designed for
!       cloud resolving simulations.
!       Jin-Luen Lee December, 2008
!       Jacques Middlecoff July 2009 Put variables in calling sequences.
!       Jacques Middlecoff May  2010 Moved calculations to called routines
!       Jacques Middlecoff May  2010 Split caldiv into CalDiv and ApplyDiv
!       Mark Govett      August 2010 Added initialization for SMS GPU parallelization
!       Jacques Middlecoff Oct. 2010 Removed SMS from GPU routines
!       Jacques Middlecoff Feb. 2011 Merged mediation branch with Jins trunk
!       Jacques Middlecoff Sep. 2011 Merged mediation branch with Jins trunk
!       Jacques Middlecoff Jan. 2012 Ported to the GPU
!******************************************************************************

program nim
use module_control
use module_variables
use module_constants
!
!
use icosio    , only: icosio_stop
use core_setup, only: iam_nim_task, cpu_cores_per_node, max_compute_tasks_per_node
use globals   , only: myrank, npes, GPUrun
!use ReadNamelist, only: nz
implicit none
#include <gptl.inc>
integer :: ret
integer :: num_mic_tasks = 0
character(len=40) :: str
#ifdef ENABLE_AFFINITY
integer, external :: set_affinity
#endif
!sms$start

ret = gptlsetoption (gptlverbose, 1)
!ret = gptlsetoption (gptlmaxthreads_gpu, 262208)
!ret = gptlsetoption (gptlmaxwarps_gpu, 448)   ! for flatten
! This one is the highest number that works
!ret = gptlsetoption (gptlmaxwarps_gpu, 983187)
ret = gptlsetoption (gptlmaxwarps_gpu, 30724)
ret = gptlsetoption (gptlmaxtimers_gpu, 20)
! This one seems to hang
!ret = gptlsetoption (gptlmaxthreads_gpu, 1573248)
!ret = gptlsetoption (gptlmaxthreads_gpu, 65552)

ret = gptlinitialize ()
ret = gptlstart ('Total')

call start

! If MICs are present, powDouble must be false because then calling diagexact does nothing 
! but slow down the code
#ifdef XEONPHI
num_mic_tasks = 1
#endif
!SMS$REDUCE(num_mic_tasks,sum)
write(6,*)num_mic_tasks,' MIC tasks are participating'
! MWG: removed, but Jim needs to test
!if (num_mic_tasks > 0 .and. powDouble) then
!  stop 'NIM error: powDouble must be false when MICs are involved'
!end if

#ifdef ENABLE_AFFINITY
ret = set_affinity (myrank,cpu_cores_per_node,max_compute_tasks_per_node,pin_to_single_core,root_on_socket1)
if (ret /= 0) then
  write(6,*)'NIM error: set_affinity gave non-zero return value=', ret
  stop
end if
call print_affinity (myrank)
#endif

if (myrank == 0 .or. myrank == npes-1) then
  str = ' '
  write(str,'("Before input myrank=",i3)') myrank
  ret = gptlprint_memusage (str)
end if

call input

if (myrank == 0 .or. myrank == npes-1) then
  str = ' '
  write(str,'("Before init myrank=",i3)') myrank
  ret = gptlprint_memusage (str)
end if

call init

!sms$compare_var(trc1,"nim.F90 - trc1_1")

ret = gptlstart('copytoGPUnim')
!$acc update device(ct4w,caw4tk,caw4tp,f,amtx1,amtx2,uw8s,uw8b,zc,zcs,zb3d,reb,rebb,teb,tebb,xys,xyb,xyc,prox,saoprx,div,sedgvar,st,su,sv,sw)
ret = gptlstop ('copytoGPUnim')

ret = gptlstart('InitExchange')
!sms$exchange (u,v,w,ur,vr,wr,tp,trp,rp,trb,trc1)
ret = gptlstop ('InitExchange')

!sms$compare_var(trc1,"nim.F90 - trc1_2")
if (itsbeg == 1 .and. writeOutput) then
  if (myrank == 0 .or. myrank == npes-1) then
    str = ' '
    write(str,'("Before 1st output myrank=",i3)') myrank
    ret = gptlprint_memusage (str)
  end if
  call output (nip,itsbeg-1,minmaxPrintInterval,itsend,numvars,var_list,r,wr,trp,u,v,w,t,tp,p,rp,ur,vr)
endif

!sms$compare_var(u,"nim.F90 - u1")
!sms$compare_var(v,"nim.F90 - v1")
!sms$compare_var(w,"nim.F90 - w1")
!sms$compare_var(ur,"nim.F90 - ur1")
!sms$compare_var(vr,"nim.F90 - vr1")
!sms$compare_var(wr,"nim.F90 - wr1")
!sms$compare_var(trp,"nim.F90 - trp1")
!sms$compare_var(rp,"nim.F90 - rp1")
!sms$compare_var(trb,"nim.F90 - trb1")
!sms$compare_var(trc1,"nim.F90 - trc1_3")

if (iam_nim_task) then
  if(TimeInitExchanges) then
!sms$unstructured_print_timers
  endif
  if(zeroSMStimers) then
!SMS$zerotimers
  endif
endif

ret = gptlstart ('MainLoop')
ret = gptlstart ('MainLoop_nooutput')
ret = gptlstart ('MainLoopCompute')

do its=itsbeg,itsend
  print*,'its=',its
  ret = gptlstart('ZeroTendencies')
  call ZeroTendencies(ims,ime,ips,ipe,fr,fur,fvr,fwr,ftr)
  ret = gptlstop ('ZeroTendencies')
! Normally would call physics here but this is a dynamics-only code

!sms$compare_var(trc1,"nim.F90 - trc1_4")
!sms$compare_var(ur,"nim.F90 - ur2")
!sms$compare_var(vr,"nim.F90 - vr2")
!sms$compare_var(wr,"nim.F90 - wr2")
!sms$compare_var(trp,"nim.F90 - trp2")
!sms$compare_var(r,"nim.F90 - r2")
  ret = gptlstart('SaveFlux')
  call SaveFlux  (ims,ime,ips,ipe,ur,vr,wr,trp,r,urs,vrs,wrs,trs,rs)
  ret = gptlstop ('SaveFlux')
!sms$compare_var(urs,"nim.F90 - urs3")
!sms$compare_var(vrs,"nim.F90 - vrs3")
!sms$compare_var(wrs,"nim.F90 - wrs3")
!sms$compare_var(trs,"nim.F90 - trs3")
!sms$compare_var(rs,"nim.F90 - rs3")
!sms$compare_var(u,"nim.F90 - u3")
!sms$compare_var(ur,"nim.F90 - ur3")
!sms$compare_var(tfur,"nim.F90 - tfur3")
  call dyntnc
!sms$compare_var(trc1,"nim.F90 - trc1_5")
  ret = gptlstart('RKdiff')
  call RKdiff    (npp,ims,ime,ips,ipe,FullStep,frk(1),dt,rhomin,nprox,sa,vol,tefr,efr,tfr,tfur,tfvr,tfwr,tftr,fr,fur,fvr,fwr,ftr)
  ret = gptlstop ('RKdiff')
!sms$compare_var(tfr,"nim.F90 - tfr4")
!sms$compare_var(fr,"nim.F90 - fr4")
!sms$compare_var(tfur,"nim.F90 - tfur4")
!sms$compare_var(tfvr,"nim.F90 - tfvr4")
!sms$compare_var(tfwr,"nim.F90 - tfwr4")
!sms$compare_var(tftr,"nim.F90 - tftr4")
!sms$compare_var(efr,"nim.F90 - efr4")
!sms$compare_var(amtx1,"nim.F90 - amtx1_4")
!sms$compare_var(trc1,"nim.F90 - trc1_6")
  do irk=1,nrkl-1
    ret = gptlstart('TimeDiff')
    call TimeDiff(ims,ime,ips,ipe,PartialStep,fdt(irk), fr,fur,fvr,fwr,ftr, tfr,tfur,tfvr,tfwr,tftr, rs,urs,vrs,wrs,trs, r,ur,vr,wr,trp)
    ret = gptlstop ('TimeDiff')
!sms$compare_var(r,"nim.F90 - r7")
!sms$compare_var(ur,"nim.F90 - ur7")
!sms$compare_var(vr,"nim.F90 - vr7")
!sms$compare_var(wr,"nim.F90 - wr7")
!sms$compare_var(trp,"nim.F90 - trp7")
    ret = gptlstart('Sponge')
    call sponge    (ims,ime,ips,ipe,dt,ca4k,ca4p,spncef,r,ur,vr,wr,u,v,w)
    ret = gptlstop ('Sponge')
    ret = gptlstart('Diag')
    call diag    (ntr1,ims,ime,ips,ipe,powDouble,p1000,rd,gamma, &
                  kapa,g,qvmin,qwmin,qrmin,area,vol,pm1d,spncef,ca4k, &
                  ca4p,r,rb,trp,trb,tb,trc1,ur,vr,wr, &
                  rp,tp,tr,t,p,qv,qw,qr,u,v, &
                  e,w,tkv)
    ret = gptlstop ('Diag')
!sms$compare_var(rp,"nim.F90 - rp8")
!sms$compare_var(tp,"nim.F90 - tp8")
!sms$compare_var(tr,"nim.F90 - tr8")
!sms$compare_var(t,"nim.F90 - t8")
!sms$compare_var(p,"nim.F90 - p8")
!sms$compare_var(qv,"nim.F90 - qv8")
!sms$compare_var(qw,"nim.F90 - qw8")
!sms$compare_var(u,"nim.F90 - u8")
!sms$compare_var(v,"nim.F90 - v8")
!sms$compare_var(e,"nim.F90 - e8")
!sms$compare_var(w,"nim.F90 - w8")
    call dyntnc
    ret = gptlstart('RKdiff')
    call RKdiff  (npp,ims,ime,ips,ipe,PartialStep,frk(irk+1),dt,rhomin,nprox,sa,vol,tefr,efr,tfr,tfur,tfvr,tfwr,tftr,fr,fur,fvr,fwr,ftr)
    ret = gptlstop ('RKdiff')
!sms$compare_var(tfr,"nim.F90 - tfr9")
!sms$compare_var(fr,"nim.F90 - fr9")
!sms$compare_var(tfur,"nim.F90 - tfur9")
!sms$compare_var(tfvr,"nim.F90 - tfvr9")
!sms$compare_var(tfwr,"nim.F90 - tfwr9")
!sms$compare_var(tftr,"nim.F90 - tftr9")
!sms$compare_var(efr,"nim.F90 - efr9")
  enddo !irk loop
  ret = gptlstart('TimeDiff')
  call TimeDiff  (ims,ime,ips,ipe,FullStep,fdt(1), fr,fur,fvr,fwr,ftr, tfr,tfur,tfvr,tfwr,tftr, rs,urs,vrs,wrs,trs, r,ur,vr,wr,trp)
  ret = gptlstop ('TimeDiff')
!sms$compare_var(r,"nim.F90 - r10")
!sms$compare_var(ur,"nim.F90 - ur10")
!sms$compare_var(vr,"nim.F90 - vr10")
!sms$compare_var(wr,"nim.F90 - wr10")
!sms$compare_var(trp,"nim.F90 - trp10")
  ret = gptlstart('Sponge')
  call sponge    (ims,ime,ips,ipe,dt,ca4k,ca4p,spncef,r,ur,vr,wr,u,v,w)
  ret = gptlstop ('Sponge')
  ret = gptlstart('Diag')
  call diag   (ntr1,ims,ime,ips,ipe,powDouble,p1000,rd,gamma, &
                  kapa,g,qvmin,qwmin,qrmin,area,vol,pm1d,spncef,ca4k, &
                  ca4p,r,rb,trp,trb,tb,trc1,ur,vr,wr, &
                  rp,tp,tr,t,p,qv,qw,qr,u,v, &
                  e,w,tkv)
  ret = gptlstop ('Diag')

  doBarrier = preExchangeBarrier.and.noOverlap
  call preExchange ('MLexch1begPreBarrier' ,'MLexch1beg')
!sms$exchange_begin (trc1)
  doBarrier = postExchangeBarrier.and.noOverlap
  call postExchange('MLexch1begPostBar','MLexch1beg')

!sms$compare_var(rp,"nim.F90 - rp11")
!sms$compare_var(tp,"nim.F90 - tp11")
!sms$compare_var(tr,"nim.F90 - tr11")
!sms$compare_var(t,"nim.F90 - t11")
!sms$compare_var(p,"nim.F90 - p11")
!sms$compare_var(qv,"nim.F90 - qv11")
!sms$compare_var(qw,"nim.F90 - qw11")
!sms$compare_var(u,"nim.F90 - u11")
!sms$compare_var(v,"nim.F90 - v11")
!sms$compare_var(e,"nim.F90 - e11")
!sms$compare_var(w,"nim.F90 - w11")
  ret = gptlstart('pre_trisol')
  call pre_trisol(nvars,ims,ime,ips,ipe,wr,tp,trp,ca4k,ca4p,z,zm,tebb,bedgvar)
  ret = gptlstop('pre_trisol')
!sms$compare_var(trp,"nim.F90 - trp11.2")
!sms$compare_var(wr,"nim.F90 - wr11.2")
!sms$compare_var(r,"nim.F90 - r11.2")
!sms$compare_var(bedgvar,"nim.F90 - bedgvar11.2")
                 ret = gptlstart('Trisol')
  call trisol(ims,ime,ips,ipe,nvars,cn0,cn1,dt,rd,cv,g,ct4w,caw4tk,caw4tp,ca4p,ca4k,spncef,e,nvecb,tb,rb,st,rp,trp,wr,r,uw8b,vdnb,bedgvar)
                 ret = gptlstop ('Trisol')
!sms$compare_var(trp,"nim.F90 - trp11.3")
!sms$compare_var(wr,"nim.F90 - wr11.3")
!sms$compare_var(r,"nim.F90 - r11.3")
!sms$compare_var(uw8b,"nim.F90 - uw8b11.3")
!sms$compare_var(vdnb,"nim.F90 - vdnb11.3")
!sms$compare_var(bedgvar,"nim.F90 - bedgvar11.3")
!sms$compare_var(trc1,"nim.F90 - trc111.4")
!sms$compare_var(sedgvar,"nim.F90 - sedgvar11.4")
!sms$compare_var(bedgvar,"nim.F90 - bedgvar11.4")
!sms$compare_var(uw8b,"nim.F90 - uw8b11.4")

  doBarrier = preExchangeBarrier.and.(.not.noOverlap)
  call preExchange ('MLexch1endPreBarrier' ,'MLexch1end')
!sms$exchange_end(trc1)
  doBarrier = postExchangeBarrier.and.(.not.noOverlap)
  call postExchange('MLexch1endPostBar','MLexch1end')

!TODO:  Optimize these vdmints() calls using approach of vdmintsG5 to improve re-use 
!TODO:  of amtx1.  
  ret = gptlstart('Vdmints')
  call vdmints(trc1,sedgvar,bedgvar,nvars,amtx1,nbf,npp,ims,ime,nob,ntr1,nprox,prox,xyc,zc,zm,z,ca4k,ca4p,ips,ipe)
  ret = gptlstop ('Vdmints')
!sms$compare_var(vdns,"nim.F90 - vdns11.5")
!sms$compare_var(vdnb,"nim.F90 - vdnb11.5")
!sms$compare_var(sedgvar,"nim.F90 - sedgvar11.5")
!sms$compare_var(bedgvar,"nim.F90 - bedgvar11.5")
!sms$compare_var(sa,"nim.F90 - sa11.5")

  doBarrier = preExchangeBarrier.and.noOverlap
  call preExchange ('MLexchange2begPreBarrier' ,'MLexchange2beg')
!sms$exchange_begin (sedgvar)
  doBarrier = postExchangeBarrier.and.noOverlap
  call postExchange('MLexchange2begPostBar','MLexchange2beg')

  ret = gptlstart('post_trisol')
  call post_trisol(nvars,ims,ime,ips,ipe,wr,wrs,ca4k,ca4p,z,nvecb,vdnb,uw8b,bedgvar)
  ret = gptlstop('post_trisol')

  doBarrier = preExchangeBarrier.and.(.not.noOverlap)
  call preExchange ('MLexchange2endPreBar' ,'MLexchange2end')
!sms$exchange_end(sedgvar)
  doBarrier = postExchangeBarrier.and.(.not.noOverlap)
  call postExchange('MLexchange2endPostBar','MLexchange2end')

  ret = gptlstart('Pstadv')
  call pstadv(ntr1,ims,ime,ips,ipe,npp,nvars,dt,nprox,prox,proxs,vdns,vdnb,sa,area,vol,ca4k,ca4p,sedgvar,trc1)
  ret = gptlstop ('Pstadv')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!sms$compare_var(rp,"nim.F90 - rp12")
!sms$compare_var(trp,"nim.F90 - trp12")
!sms$compare_var(wr,"nim.F90 - wr12")
!sms$compare_var(r,"nim.F90 - r12")
!sms$compare_var(uw8b,"nim.F90 - uw8b12")
!sms$compare_var(vdnb,"nim.F90 - vdnb12")
!sms$compare_var(bedgvar,"nim.F90 - bedgvar12")
  ret = gptlstart('Diag')
  call diag    (ntr1,ims,ime,ips,ipe,powDouble,p1000,rd,gamma, &
                  kapa,g,qvmin,qwmin,qrmin,area,vol,pm1d,spncef,ca4k, &
                  ca4p,r,rb,trp,trb,tb,trc1,ur,vr,wr, &
                  rp,tp,tr,t,p,qv,qw,qr,u,v, &
                  e,w,tkv)
  ret = gptlstop ('Diag')
!sms$compare_var(rp,"nim.F90 - rp13")
!sms$compare_var(tp,"nim.F90 - tp13")
!sms$compare_var(tr,"nim.F90 - tr13")
!sms$compare_var(t,"nim.F90 - t13")
!sms$compare_var(p,"nim.F90 - p13")
!sms$compare_var(qv,"nim.F90 - qv13")
!sms$compare_var(qw,"nim.F90 - qw13")
!sms$compare_var(u,"nim.F90 - u13")
!sms$compare_var(v,"nim.F90 - v13")
!sms$compare_var(e,"nim.F90 - e13")
!sms$compare_var(w,"nim.F90 - w13")
  if(writeOutput) then
    ret = gptlstop ('MainLoop_nooutput')
    ret = gptlstop ('MainLoopCompute')
    call output(nip,its,minmaxPrintInterval,itsend,numvars,var_list,r,wr,trp,u,v,w,t,tp,p,rp,ur,vr)
    ret = gptlstart ('MainLoopCompute')
    ret = gptlstart ('MainLoop_nooutput')
  endif
enddo ! time steps

if (myrank == 0 .or. myrank == npes-1) then
  str = ' '
  write(str,'("End time loop myrank=",i3)') myrank
  ret = gptlprint_memusage (str)
end if

ret = gptlstop  ('MainLoopCompute')
ret = gptlstop  ('MainLoop_nooutput')
ret = gptlstop  ('MainLoop')
call icosio_stop(its-1)
call stop
write (*,'(a)') 'NIM completed.'
call flush(6)

!sms$stop

end program nim
