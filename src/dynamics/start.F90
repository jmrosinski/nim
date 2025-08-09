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

subroutine start
!**************************************************************************
!
!       Starts MPI and SMS and reads the initial variables and constants.
!       Also calls core_setup_nim() to split the MPI communicator (if MPI 
!       is being used) into "compute", "write", and "do-nothing"
!       sub-communicators.  Any "write" tasks and "do-nothing" tasks will 
!       not return from this subroutine.  All "compute" tasks will return 
!       from this subroutine.  
!       Alexander E. MacDonald  11/27/04
!       J. Lee                  September, 2005
!       Jacques Middlecoff      7/09 Separated start into start, input and init.
!
!**************************************************************************
!use namelistdata,   only: glvl
use module_variables
use module_constants
use module_control
use core_setup,        only: core_setup_nim, iam_nim_task, iam_write_task, &
                             my_comm, write_intercomm, use_write_tasks,    &
                             icosio_debugmsg_on,max_accelerators_per_node
use module_alloc,      only: alloc
use globals,           only: myrank, npes, GPUrun, parallelBuild
use print_taskinfo,    only: printmem
use icosio,            only: icosio_setup
use units,             only: getunit,initunit
use kinds,             only: rt,sp
use ReadNamelist, only: nz
#ifndef SERIAL
use mpi
#endif

IMPLICIT NONE
#include <gptl.inc>

integer,allocatable :: icos_prox (:,:) !icos global index of proximity points
integer,allocatable :: icos_nprox(  :) !# pts in dim 1 of icos_prox; 5 or 6 in NIM
integer :: ret
! non-decomposed versions of various grid-related index arrays used during start-up only
integer :: ipn,ipns,isn,j,iprox
integer :: unitno              !Unit number for I/O
integer :: icosio_comm
integer :: ij,k,patchsize,tilesize,numlargertiles
integer :: ips_set,ipe_set,ierr
real(sp),allocatable :: tmpr4 (:)
CHARACTER(80) :: FileName
integer,allocatable :: regionSize(:)!Region size for each PE

ret = gptlstart ('Start')

!TODO:  pass in list of reserved unit numbers if needed
call initunit ()

! When MPI is used, set up communicators for compute tasks and 
! optional write tasks.  Mirrors NEMS approach.  
! NOTE:  Executable SMS directives must not be placed before 
! NOTE:  this call.  
! NOTE:  This includes writes or prints without !SMS$ignore because they 
! NOTE:  cause SMS to generate code.

call core_setup_nim ()

!$$$DEBUG BEGIN
if (.not.(iam_nim_task .or. iam_write_task)) then
!SMS$IGNORE BEGIN
  print *,'ERROR:  do-nothing task still in start() after core_setup_nim()'
  call flush(6)
!SMS$IGNORE END
  stop
endif
!$$$DEBUG END

if (iam_nim_task) then
  !JR Cannot use printmem here because create_decomp has not yet been done
  if (myrank == 0) then
    ret = gptlprint_memusage ('Start of run')
  end if
endif

call control !Read the namelist

if (iam_nim_task) then

  allocate(icos_prox(npp,nip),icos_nprox(nip))

!SMS$SERIAL(default=ignore) BEGIN
  unitno = getunit (28)
  if (unitno < 0) then
    print*,'dyn_init: getunit failed for icosio_setup. Stopping in start'
    call flush(6)
    stop
  end if
  ! Get prox
  !TODO:  Remove duplication with input.F90 !!  
  FileName="g3d.dat"
  open(unit=unitno, file=TRIM(FileName), form="unformatted", status='old', action='read', iostat=ierr)
  if (ierr.ne.0) then
    write (*,'(a,a)') 'start: error opening ',trim(filename)
    stop
  endif
  print *,'  Reading icos_prox & icos_nprox from ',TRIM(FileName),' ...'

  allocate(tmpr4(nip))
  call read1Dr4  ('lat',unitno,nip,tmpr4)
  call read1Dr4  ('lon',unitno,nip,tmpr4)
  do k=0,nz
    call read1Dr4('z'  ,unitno,nip,tmpr4)
  enddo
  do k=1,nz
    call read1Dr4('zm' ,unitno,nip,tmpr4)
  enddo
  deallocate(tmpr4)
  call read2Di4('icos_prox' ,unitno,npp,nip,icos_prox )
  call read1Di4('icos_nprox',unitno,    nip,icos_nprox)
  close(unitno)
  do ipn = 1, nip
    if (icos_nprox(ipn)==5) then
      icos_prox(6,ipn)=-1
    endif
  enddo

!SMS$SERIAL END

  allocate(regionSize(nPEs))
!SMS$barrier
ret = gptlstart ('order_grid')
!SMS$order_grid(npp,nip,icos_nprox,icos_prox,regionSize,haloSize)
ret = gptlstop  ('order_grid')
ret = gptlstart ('CREATE_DECOMP')
!SMS$CREATE_DECOMP(dh,<nip>,<haloSize>:regionsize=regionSize)
ret = gptlstop  ('CREATE_DECOMP')
  deallocate(RegionSize)

! Distribute proximity grid index arrays.  Note that SERIAL can correctly 
! handle distributed arrays only after CREATE_DECOMP (above).  
! Allocate proximity grid index arrays.  
allocate(prox           (npp, nip)) ! holds index of proximity points
allocate(proxs          (npp, nip)) ! holds index of proximity sides
allocate(nprox          (     nip)) ! holds number of proximity points
allocate(perm           (     nip)) ! permutation of the grid
if (myrank == 0) then
  allocate(globalperm     (nip)) ! global permutation of the grid
  allocate(global_inv_perm(nip)) ! global inverse permutation of the grid
else
  allocate(globalperm     (1)) ! global permutation of the grid
  allocate(global_inv_perm(1)) ! global inverse permutation of the grid
endif

if(.not.parallelBuild .or. nPEs==1) then! serial or 1 processor
  do ipn=1,nip
    perm           (ipn) = ipn
    globalperm     (ipn) = ipn
    global_inv_perm(ipn) = ipn
    nprox          (ipn) = icos_nprox(ipn)
    do isn=1,nprox(ipn)
      prox(isn,ipn) = icos_prox(isn,ipn)
    enddo
    if(nprox(ipn)==5)then
      prox(6,ipn) = prox(5,ipn)
    endif
  enddo
  proxs=-99
  do ipn=1,nip
    ISNLOOP: do isn=1,nprox(ipn)
      iprox=prox(isn,ipn)
      do j=1,nprox(iprox)
        if(prox(j,iprox).eq.ipn) then
          proxs(isn,ipn)=j
          cycle ISNLOOP
        end if
      enddo
      print*,'Error in start: Could not find prox.'
      print"(11i7)",ipn,isn,iprox,perm(ipn),perm(iprox),prox(1:nprox(iprox),iprox)
      print*,'Stopping in start.'
      stop
    enddo ISNLOOP
  enddo !ipn
endif
deallocate(icos_prox,icos_nprox)

ret = gptlstart ('set_Prox')
!SMS$set_Prox_And_Proxs(dh,perm,globalPerm,global_inv_perm,nprox,prox,proxs)
ret = gptlstop  ('set_Prox')

  ips_set=1
  ipe_set=nip
!SMS$TO_LOCAL(dh:<1,ips_set:lbound>,<1,ipe_set:ubound>) BEGIN
  ips = ips_set
  ipe = ipe_set
!SMS$TO_LOCAL END
!sms$ignore begin
  ims = lbound(nprox,1)
  ime = ubound(nprox,1)
!sms$ignore end
  ihs = ipe+1
  ihe = ipe
!sms$parallel (dh,ipn) begin
!sms$halo_comp(<1,1>) begin
  do ipn=1,nip
    ipns = ipn
  enddo
!sms$halo_comp end
!sms$parallel end
  ihe=ipns-1
! ihe = ipe+sum(recvCount(1:numRecvs(myrank+1),myrank+1))
endif !(iam_nim_task)

!JR Tell icosio to use unit 50. Since getunit() with no arguments may return 
!JR different numbers for different MPI tasks, want to be safe.
unitno = getunit (50)
if (unitno < 0) then
  print*,'dyn_init: getunit failed for icosio_setup. Stopping'
  call flush(6)
  stop
end if

icosio_comm = my_comm
if (use_write_tasks) icosio_comm = write_intercomm

! Send icosio required values. Since client/server io is enabled, write tasks 
! will not return from the icosio_setup() call.
if (iam_write_task) then
  call icosio_setup(                     &
    client_server_io_in=.true.,          &
    debugmsg_on_in=icosio_debugmsg_on,   &
    comm_in=icosio_comm,                 &
    i_am_write_task_in=iam_write_task    &
    )
else
  call icosio_setup(                     &
    binout_in=.true.,                    &
    client_server_io_in=.true.,          &
    comm_in=icosio_comm,                 &
    debugmsg_on_in=icosio_debugmsg_on,   &
    gribout_in=.false.,                  &
    i_am_write_task_in=iam_write_task,   &
    inv_perm_global_in=global_inv_perm,  &
    ipe_in=ipe,                          &
    ips_in=ips,                          &
    lunout_in=unitno,                    &
    nip_in=nip,                          &
    permute_in=.true.,                   &
    print_diags_in=.false.,              &
    using_write_tasks_in=use_write_tasks &
    ) 
endif

!$$$DEBUG BEGIN
if (iam_write_task) then
!SMS$IGNORE BEGIN
  print *,'ERROR:  write task still in start() after icosio_setup!'
  call flush(6)
!SMS$IGNORE END
  stop
endif
!$$$DEBUG END

! Phony loop to verify threading configuration--needed since default timing
! is zero calls inside of threaded loops
!$OMP PARALLEL DO PRIVATE (ret)
do ipn=ips,ipe
  ret = gptlstart ('threading')
  ret = gptlstop ('threading')
end do
!$OMP END PARALLEL DO

!JR TODO: Fix this:
!JR Where is num_tiles coming from? It lives in ReadNamelist so that must be globally "used" in one
!JR of the global "use" stmts above.
allocate(i_start(num_tiles))
allocate(i_end(num_tiles))
patchsize = ipe - ips + 1
if (num_tiles > patchsize) then
  print*,'Error start.F90 num_tiles is too large, reduce in NIMnamelist'
  stop
endif
tilesize = patchsize/num_tiles
numlargertiles = patchsize - (tilesize*num_tiles)
i_start(1) = ips
do ij=1,num_tiles-1
  i_end(ij) = i_start(ij) + tilesize - 1
  if (ij <= numlargertiles) then
    i_end(ij) = i_end(ij) + 1
  endif
  i_start(ij+1) = i_end(ij) + 1
enddo
i_end(num_tiles) = ipe
call GPUinit(myRank,max_accelerators_per_node,GPUrun,ret)
if(ret /= 0) then
!SMS$ignore begin
  print*,'Error start.F90 calling GPUinit.cu',ret,GPUrun
  call flush(6)
!SMS$ignore end
  stop
endif
call printmem ('start before alloc')
call alloc ()   ! Allocate the dynamics variables
call printmem ('start after alloc')

! allocate physics variables
call printmem ('start before alloc physics vars')
if(physics == 'NONE') then
  !Physics variables initialized here so NIM can be run without physics.
  st = 0.0_rt
  su = 0.0_rt
  sv = 0.0_rt
else
  print*,'Only physics = "NONE" is supported'
  stop
endif
call printmem ('start after alloc physics vars')

!SMS$SERIAL BEGIN
open(unit=37,file="nim_diat.dat",form="unformatted", action='write', iostat=ierr)
if (ierr.ne.0) then
  write (*,'(a)') 'start: error opening nim_diat.dat'
  stop
endif
open(unit=39,file="nim_diag.dat",form="unformatted", action='write', iostat=ierr)
if (ierr.ne.0) then
  write (*,'(a)') 'start: error opening nim_diag.dat'
  stop
endif
open(unit=89,file="physics_diag.dat",form="unformatted", action='write', iostat=ierr)
if (ierr.ne.0) then
  write (*,'(a)') 'start: error opening physics_diag.dat'
  stop
endif
!SMS$SERIAL END

ret = gptlstop  ('Start')
return

end subroutine start 

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
