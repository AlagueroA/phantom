!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module centreofmass
!
! Utilities for computing the centre of mass on the particles
!   and correcting the bulk motion (used for turbulent driving)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, io, mpiutils, part, vectorutils
!
 implicit none
 public :: reset_centreofmass,get_centreofmass,correct_bulk_motion,get_total_angular_momentum
 public :: get_centreofmass_accel

 private

contains

!----------------------------------------------------------------
!+
!  routine to reset the centre of mass in the initial conditions
!  (assuming equal mass particles)
!+
!----------------------------------------------------------------
subroutine reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,only_sinks,sinks_index)
 use io, only:iprint
 integer, intent(in)    :: npart
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(in),    optional :: nptmass
 real,    intent(inout), optional :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, optional :: sinks_index(:)
 logical, optional :: only_sinks
 logical :: only_sinks_
 integer, allocatable :: sinks_index_(:)
 real :: xcom(3),vcom(3),xcomold(3),vcomold(3)
 integer :: i

! Check if we consider all the particles or only the sinks
 if (present(only_sinks)) then
    only_sinks_ = only_sinks
 else
    only_sinks_ = .false.
 endif
! Construction of the sinks list if not present
 if (present(sinks_index)) then
    allocate(sinks_index_(SIZE(sinks_index)))
    sinks_index_ = sinks_index
 else
    allocate(sinks_index_(nptmass))
    do i = 1, nptmass
        sinks_index_(i) = i
    enddo
 endif

! Get CoM
 if (present(xyzmh_ptmass) .and. present(vxyz_ptmass) .and. present(nptmass)) then
    call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,&
                          only_sinks=only_sinks_,sinks_index=sinks_index_)
 else
    call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu)
 endif

! Move particles & sinks to put the CoM at the origin
 xcomold = xcom
 vcomold = vcom
 do i=1,npart
    xyzh(1:3,i) = xyzh(1:3,i) - xcom(1:3)
    vxyzu(1:3,i) = vxyzu(1:3,i) - vcom(1:3)
 enddo
 if (present(xyzmh_ptmass) .and. present(vxyz_ptmass) .and. present(nptmass)) then
    do i=1,nptmass
       xyzmh_ptmass(1:3,i) = xyzmh_ptmass(1:3,i) - xcom(1:3)
       vxyz_ptmass(1:3,i)  = vxyz_ptmass(1:3,i) - vcom(1:3)
    enddo
    call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,&
                          only_sinks=only_sinks_,sinks_index=sinks_index_)
 else
    call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu)
 endif
 write(iprint,"(' reset CofM: (',3(es9.2,1x),') -> (',3(es9.2,1x),')')") xcomold,xcom

 return
end subroutine reset_centreofmass

!----------------------------------------------------------------
!+
! Routine returns the centre of mass and centre of mass velocity
! Accounting for sink particles is optional
!+
!----------------------------------------------------------------
subroutine get_centreofmass(xcom,vcom,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,mass,only_sinks,sinks_index)
 use io,       only:id,master
 use dim,      only:maxphase,maxp
 use part,     only:massoftype,iamtype,iphase,igas,isdead_or_accreted
 use mpiutils, only:reduceall_mpi
 real,         intent(out) :: xcom(3),vcom(3)
 integer,      intent(in)  :: npart
 real,         intent(in)  :: xyzh(:,:),vxyzu(:,:)
 integer,      intent(in),  optional :: nptmass
 real,         intent(in),  optional :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, optional, intent(in) :: sinks_index(:)
 logical, optional :: only_sinks
 integer, allocatable :: sinks_index_(:)
 logical :: only_sinks_
 real,         intent(out), optional :: mass
 integer :: i,itype,n_sinks,i_sink
 real :: xi,yi,zi,hi
 real(kind=8) :: xpos,ypos,zpos,vxpos,vypos,vzpos
 real(kind=8) :: dm,pmassi,totmass

 xpos  = 0.d0
 ypos  = 0.d0
 zpos  = 0.d0
 vxpos = 0.d0
 vypos = 0.d0
 vzpos = 0.d0
 totmass = 0.d0
 pmassi = massoftype(igas)

! Check if we consider all the particles or only the sinks
 if (present(only_sinks)) then
    only_sinks_ = only_sinks
 else
    only_sinks_ = .false.
 endif
! Construction of the sinks list if not present
 if (present(sinks_index)) then
    allocate(sinks_index_(SIZE(sinks_index)))
    sinks_index_ = sinks_index
 else
    allocate(sinks_index_(nptmass))
    do i = 1, nptmass
        sinks_index_(i) = i
    enddo
 endif
 n_sinks = SIZE(sinks_index_)

!$omp parallel default(none) &
!$omp shared(only_sinks_) &
!$omp shared(maxphase,maxp) &
!$omp shared(npart,xyzh,vxyzu,iphase,massoftype) &
!$omp private(i,itype,xi,yi,zi,hi) &
!$omp firstprivate(pmassi) &
!$omp reduction(+:xpos,ypos,zpos,vxpos,vypos,vzpos,totmass)
 if (.not.only_sinks_) then
    !$omp do
     do i=1,npart
        xi = xyzh(1,i)
        yi = xyzh(2,i)
        zi = xyzh(3,i)
        hi = xyzh(4,i)
        if (.not.isdead_or_accreted(hi)) then
           if (maxphase==maxp) then
              itype = iamtype(iphase(i))
              if (itype > 0) then ! avoid problems if called from ICs
                 pmassi = massoftype(itype)
              else
                 pmassi = massoftype(igas)
              endif
           endif
           totmass = totmass + pmassi
           xpos    = xpos  + pmassi*xi
           ypos    = ypos  + pmassi*yi
           zpos    = zpos  + pmassi*zi
           vxpos   = vxpos + pmassi*vxyzu(1,i)
           vypos   = vypos + pmassi*vxyzu(2,i)
           vzpos   = vzpos + pmassi*vxyzu(3,i)
        endif
     enddo
    !$omp enddo
 endif
!$omp end parallel

 if (id==master .and. present(xyzmh_ptmass) .and. present(vxyz_ptmass) .and. present(nptmass)) then
    do i=1,n_sinks
       i_sink = sinks_index_(i)
       pmassi = xyzmh_ptmass(4,i_sink)
       if (pmassi > 0.) then
          totmass = totmass + pmassi
          xpos  = xpos  + pmassi*xyzmh_ptmass(1,i_sink)
          ypos  = ypos  + pmassi*xyzmh_ptmass(2,i_sink)
          zpos  = zpos  + pmassi*xyzmh_ptmass(3,i_sink)
          vxpos = vxpos + pmassi*vxyz_ptmass(1,i_sink)
          vypos = vypos + pmassi*vxyz_ptmass(2,i_sink)
          vzpos = vzpos + pmassi*vxyz_ptmass(3,i_sink)
       endif
    enddo
 endif
 xcom = (/xpos,ypos,zpos/)
 vcom = (/vxpos,vypos,vzpos/)
 xcom = reduceall_mpi('+',xcom)
 vcom = reduceall_mpi('+',vcom)
 totmass = reduceall_mpi('+',totmass)

 if (totmass > tiny(totmass)) then
    dm = 1.d0/totmass
 else
    dm = 0.d0
 endif
 xcom = xcom*dm
 vcom = vcom*dm

 if (present(mass)) mass = totmass

 return
end subroutine get_centreofmass

!----------------------------------------------------------------
!+
!  Subroutine to compute com acceleration
!+
!----------------------------------------------------------------
subroutine get_centreofmass_accel(acom,npart,xyzh,fxyzu,fext,nptmass,xyzmh_ptmass,fxyz_ptmass,only_sinks,sinks_index)
 use io,       only:id,master
 use dim,      only:maxphase,maxp
 use part,     only:massoftype,iamtype,iphase,igas,isdead_or_accreted
 use mpiutils, only:reduceall_mpi
 real,         intent(out) :: acom(3)
 integer,      intent(in)  :: npart
 real,         intent(in)  :: xyzh(:,:),fxyzu(:,:),fext(:,:)
 integer,      intent(in),  optional :: nptmass
 real,         intent(in),  optional :: xyzmh_ptmass(:,:), fxyz_ptmass(:,:)
 integer, optional, intent(in) :: sinks_index(:)
 logical, optional, intent(in) :: only_sinks
 integer :: i,n_sinks,i_sink
 real :: hi
 real(kind=8) :: dm,pmassi,totmass


 acom(:) = 0.
 totmass = 0.


!$omp parallel default(none) &
!$omp shared(only_sinks,sinks_index,n_sinks) &
!$omp shared(maxphase,maxp,id) &
!$omp shared(xyzh,fxyzu,fext,npart) &
!$omp shared(massoftype,iphase,nptmass) &
!$omp shared(xyzmh_ptmass,fxyz_ptmass) &
!$omp private(i,pmassi,hi,i_sink) &
!$omp reduction(+:acom) &
!$omp reduction(+:totmass)
 if (.not.only_sinks) then
    !$omp do
     do i=1,npart
        hi = xyzh(4,i)
        if (.not.isdead_or_accreted(hi)) then
           if (maxphase==maxp) then
              pmassi = massoftype(iamtype(iphase(i)))
           else
              pmassi = massoftype(igas)
           endif
           totmass = totmass + pmassi
           acom(1) = acom(1) + pmassi*(fxyzu(1,i) + fext(1,i))
           acom(2) = acom(2) + pmassi*(fxyzu(2,i) + fext(2,i))
           acom(3) = acom(3) + pmassi*(fxyzu(3,i) + fext(3,i))
        endif
     enddo
    !$omp enddo
 endif

!
! add acceleration from sink particles
!
! Check if some sinks are specified or if all sinks are considered
 if (present(sinks_index)) then
    n_sinks = SIZE(sinks_index)
 else
    n_sinks = nptmass
 endif

 if (id==master) then
    !$omp do
    do i=1,n_sinks
       if (present(sinks_index)) then
            i_sink = sinks_index(i)
       else
            i_sink = i
       endif
       pmassi = xyzmh_ptmass(4,i_sink)
       totmass = totmass + pmassi
       acom(1) = acom(1) + pmassi*fxyz_ptmass(1,i_sink)
       acom(2) = acom(2) + pmassi*fxyz_ptmass(2,i_sink)
       acom(3) = acom(3) + pmassi*fxyz_ptmass(3,i_sink)
    enddo
    !$omp enddo
 endif
!$omp end parallel

 acom = reduceall_mpi('+',acom)
 totmass = reduceall_mpi('+',totmass)
 if (totmass > tiny(totmass)) then
    dm = 1.d0/totmass
 else
    dm = 0.d0
 endif
 acom = acom*dm

end subroutine get_centreofmass_accel

!----------------------------------------------------------------
!+
!  Subroutine to compute and correct the bulk motion
!
!  This is really for use with turbulent driving routines
!  which give a net motion to the flow.
!+
!----------------------------------------------------------------
subroutine correct_bulk_motion()
 use dim,      only:maxp,maxphase
 use part,     only:npart,xyzh,vxyzu,fxyzu,iamtype,igas,iphase,&
                    nptmass,xyzmh_ptmass,vxyz_ptmass,isdead_or_accreted,&
                    massoftype
 use mpiutils, only:reduceall_mpi
 use io,       only:iprint,iverbose,id,master
 real    :: totmass,pmassi,hi,xmom,ymom,zmom
 real    :: fmeanx,fmeany,fmeanz,totmass1
 integer :: i

 xmom   = 0.
 ymom   = 0.
 zmom   = 0.
 fmeanx = 0.
 fmeany = 0.
 fmeanz = 0.
 totmass = 0.
 pmassi = massoftype(igas)
!$omp parallel default(none) &
!$omp shared(maxphase,maxp) &
!$omp shared(xyzh,vxyzu,fxyzu,npart) &
!$omp shared(massoftype,iphase) &
!$omp shared(xyzmh_ptmass,vxyz_ptmass,nptmass) &
!$omp private(i,hi) &
!$omp firstprivate(pmassi) &
!$omp reduction(+:fmeanx,fmeany,fmeanz) &
!$omp reduction(+:xmom,ymom,zmom,totmass)
!$omp do
 do i=1,npart
    hi = xyzh(4,i)
    if (.not.isdead_or_accreted(hi)) then
       if (maxphase==maxp) then
          pmassi = massoftype(iamtype(iphase(i)))
       else
          pmassi = massoftype(igas)
       endif
       totmass = totmass + pmassi

       xmom = xmom + pmassi*vxyzu(1,i)
       ymom = ymom + pmassi*vxyzu(2,i)
       zmom = zmom + pmassi*vxyzu(3,i)

       fmeanx = fmeanx + pmassi*fxyzu(1,i)
       fmeany = fmeany + pmassi*fxyzu(2,i)
       fmeanz = fmeanz + pmassi*fxyzu(3,i)
    endif
 enddo
!$omp enddo
!
!--add contribution from sink particles
!
!$omp do
 do i=1,nptmass
    pmassi  = xyzmh_ptmass(4,i)
    if (pmassi > 0.) then
       totmass = totmass + pmassi
       xmom    = xmom + pmassi*vxyz_ptmass(1,i)
       ymom    = ymom + pmassi*vxyz_ptmass(2,i)
       zmom    = zmom + pmassi*vxyz_ptmass(3,i)
    endif
 enddo
!$omp enddo
!$omp end parallel

 xmom = reduceall_mpi('+',xmom)
 ymom = reduceall_mpi('+',ymom)
 zmom = reduceall_mpi('+',zmom)
 fmeanx = reduceall_mpi('+',fmeanx)
 fmeany = reduceall_mpi('+',fmeany)
 fmeanz = reduceall_mpi('+',fmeanz)
 totmass = reduceall_mpi('+',totmass)

 totmass1 = 1./totmass
 xmom     = xmom*totmass1
 ymom     = ymom*totmass1
 zmom     = zmom*totmass1
 fmeanx   = fmeanx*totmass1
 fmeany   = fmeany*totmass1
 fmeanz   = fmeanz*totmass1

 if (id==master .and. iverbose >= 2) then
    write(iprint,*) ' correcting bulk motion ',xmom,ymom,zmom
 endif
 !$omp parallel default(none) &
 !$omp shared(npart,xyzh,vxyzu,fxyzu,xmom,ymom,zmom,fmeanx,fmeany,fmeanz) &
 !$omp shared(nptmass,vxyz_ptmass) &
 !$omp private(i)
 !$omp do schedule(static)
 do i=1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       vxyzu(1,i) = vxyzu(1,i) - xmom
       vxyzu(2,i) = vxyzu(2,i) - ymom
       vxyzu(3,i) = vxyzu(3,i) - zmom
       fxyzu(1,i) = fxyzu(1,i) - fmeanx
       fxyzu(2,i) = fxyzu(2,i) - fmeany
       fxyzu(3,i) = fxyzu(3,i) - fmeanz
    endif
 enddo
 !$omp enddo
 !$omp do
 do i=1,nptmass
    vxyz_ptmass(1,i) = vxyz_ptmass(1,i) - xmom
    vxyz_ptmass(2,i) = vxyz_ptmass(2,i) - ymom
    vxyz_ptmass(3,i) = vxyz_ptmass(3,i) - zmom
 enddo
 !$omp enddo
 !$omp end parallel

end subroutine correct_bulk_motion

!------------------------------------------------------------------------
!
! Small routine to calculate the total angular momentum vector of
! the system.
! Sinks index can be specify and not taking particles into account is optional.
!
!------------------------------------------------------------------------
subroutine get_total_angular_momentum(xyzh,vxyz,npart,L_tot,xyzmh_ptmass,vxyz_ptmass,npart_ptmass,only_sinks,sinks_index)
 use vectorutils, only:cross_product3D
 use part,        only:iphase,iamtype,massoftype,isdead_or_accreted
 use mpiutils,    only:reduceall_mpi
 real, intent(in)  :: xyzh(:,:),vxyz(:,:)
 real, optional, intent(in):: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in) :: npart
 integer, optional, intent(in) :: npart_ptmass
 integer, optional, intent(in) :: sinks_index(:)
 logical, optional :: only_sinks
 integer, allocatable :: sinks_index_(:)
 logical :: only_sinks_
 real, intent(out) :: L_tot(3)
 integer           :: ii,itype,n_sinks,i_sink
 real              :: temp(3),pmassi

 L_tot(:) = 0.
! Check if we consider all the particles or only the sinks
 if (present(only_sinks)) then
    only_sinks_ = only_sinks
 else
    only_sinks_ = .false.
 endif
! Construction of the sinks index list if not present
 if (present(sinks_index)) then
    allocate(sinks_index_(SIZE(sinks_index)))
    sinks_index_ = sinks_index
 else
    allocate(sinks_index_(npart_ptmass))
    do ii = 1, npart_ptmass
        sinks_index_(ii) = ii
    enddo
 endif
 n_sinks = SIZE(sinks_index_)

 ! Calculate the angular momentum from all the particles
 ! Check if particles are dead or have been accreted first
!$omp parallel default(none) &
!$omp shared(only_sinks_,sinks_index_,n_sinks) &
!$omp shared(xyzh,vxyz,npart) &
!$omp shared(massoftype,iphase) &
!$omp shared(xyzmh_ptmass,vxyz_ptmass,npart_ptmass) &
!$omp private(ii,itype,pmassi,temp,i_sink) &
!$omp reduction(+:L_tot)
 if (.not.only_sinks_) then
    !$omp do
     do ii = 1,npart
        if (.not.isdead_or_accreted(xyzh(4,ii))) then
           itype = iamtype(iphase(ii))
           pmassi = massoftype(itype)
           call cross_product3D(xyzh(1:3,ii),vxyz(1:3,ii),temp)
           L_tot = L_tot + temp*pmassi
        endif
     enddo
    !$omp enddo
 endif

 ! Calculate from the sinks
 if (present(npart_ptmass)) then
    !$omp do
     do ii = 1,n_sinks
        i_sink = sinks_index_(ii)
        call cross_product3D(xyzmh_ptmass(1:3,i_sink),vxyz_ptmass(1:3,i_sink),temp)
        L_tot = L_tot + temp*xyzmh_ptmass(4,i_sink)
     enddo
    !$omp enddo
 endif
!$omp end parallel

 L_tot = reduceall_mpi('+',L_tot)

end subroutine get_total_angular_momentum

end module centreofmass
