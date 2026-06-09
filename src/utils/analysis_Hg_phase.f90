!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Get gas scale height and binary orbital phase
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, discanalysisutils, infile_utils, io, part, physcon
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'disc'
 public :: do_analysis

 integer, parameter :: nr = 300
 integer, parameter :: nphi = 36
 real,dimension(nr) :: twist,twistprev

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 use io,      only:fatal
 use physcon, only:pi
 use part,    only:xyzmh_ptmass,vxyz_ptmass,nptmass,dustfrac
 use dim,     only:gr
 use infile_utils, only:open_db_from_file,read_inopt,close_db,inopts
 character(len=*), intent(in) :: dumpfile
 real,             intent(inout) :: xyzh(:,:),vxyz(:,:)
 real,             intent(inout) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile
 character(len=9)   :: output
 character(len=20)  :: filename
 character(len=20)  :: discprefix
 character(len=13)  :: phasefile
 character(len=7)   :: Hgfile
 character(len=7)   :: Hdfile
 integer :: i,ierr,iline
 real :: rmin,rmax
 type(inopts), allocatable :: db(:)

 integer :: ibin
 integer :: counts(nphi),countsd(nphi)
 real :: r,phi
 real :: dphi
 real :: theta
 real :: zsum(nphi),zsumd(nphi)
 real :: z2sum(nphi),z2sumd(nphi)
 real :: zmean(nphi),zmeand(nphi)
 real :: sigma(nphi),sigmad(nphi)
 real :: phase,ddg


! User
 rmin = 10 !5     ! radial range (code units)
 rmax = 25 !20

! Name files
 write(phasefile,"(a8,i5.5)") 'binphase',numfile
 write(Hgfile,"(a2,i5.5)") 'Hg',numfile
 write(Hdfile,"(a2,i5.5)") 'Hd',numfile

! ==========================
! BINARY ORBITAL PHASE
! ==========================
 phase = atan2(xyzmh_ptmass(2,2)-xyzmh_ptmass(2,1), xyzmh_ptmass(1,2)-xyzmh_ptmass(1,1))

 if (phase < 0.0) then
    phase = phase + 2.0*pi
 end if

! ==========================
! BIN INITIALIZATION
! ==========================
 dphi = 2.0*pi / real(nphi)

 counts = 0  ! gas
 zsum   = 0.0
 z2sum  = 0.0

 countsd = 0 ! dust
 zsumd   = 0.0
 z2sumd  = 0.0

! ==========================
! LOOP OVER PARTICLES
! ==========================
 do i = 1,npart

    r = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2)

    ! Radial selection
    if (r >= rmin .and. r < rmax) then

        theta = atan2(xyzh(2,i), xyzh(1,i))

        if (theta < 0.0) then
            theta = theta + 2.0*pi
        end if

        ibin = int(theta / dphi) + 1

        if (ibin > nphi) ibin = nphi

        counts(ibin) = counts(ibin) + 1

        zsum(ibin)  = zsum(ibin)  + xyzh(3,i)
        z2sum(ibin) = z2sum(ibin) + xyzh(3,i)**2
        
        ddg = dustfrac(1,i) / (1-dustfrac(1,i)) !dust-to-gas ratio
        if (ddg>0.01) then
            countsd(ibin) = countsd(ibin) + 1
            zsumd(ibin)   = zsumd(ibin)  + xyzh(3,i)
            z2sumd(ibin)  = z2sumd(ibin) + xyzh(3,i)**2
        endif

    endif

 end do

! ==========================
! COMPUTE GAUSSIAN PARAMETERS
! ==========================
 do ibin = 1,nphi

    if (counts(ibin) > 1) then

        zmean(ibin) = zsum(ibin) / real(counts(ibin))

        sigma(ibin) = sqrt( &
            z2sum(ibin)/real(counts(ibin)) &
            - zmean(ibin)**2 )

    else
        zmean(ibin)  = 0.0
        sigma(ibin)  = 0.0
    endif

    if (countsd(ibin) > 1) then

        zmeand(ibin) = zsumd(ibin) / real(countsd(ibin))

        sigmad(ibin) = sqrt( &
            z2sumd(ibin)/real(counts(ibin)) &
            - zmeand(ibin)**2 )
    else
        zmeand(ibin) = 0.0
        sigmad(ibin) = 0.0

    endif

 enddo

! ==========================
! OUTPUT BINARY PHASE
! ==========================
 open(unit=10, file=phasefile, status="replace")

 write(10,'(F20.10)') phase

 close(10)

! ==========================
! OUTPUT GAUSSIAN RESULTS
! ==========================
! gas
 open(unit=20, file=Hgfile, status="replace")

 !write(20,'(A)') &
 !   '# phi_center(rad)   z0(mean)   sigma   Nparticles'

 do ibin = 1, nphi

    phi = (real(ibin)-0.5)*dphi

    write(20,'(4ES20.10)') &
        phi, zmean(ibin), sigma(ibin), &
        real(counts(ibin))

 enddo

 close(20)

! dust
 open(unit=20, file=Hdfile, status="replace")

!write(20,'(A)') &
!   '# phi_center(rad)   z0(mean)   sigma   Nparticles'

 do ibin = 1, nphi

   phi = (real(ibin)-0.5)*dphi

   write(20,'(4ES20.10)') &
       phi, zmeand(ibin), sigmad(ibin), &
       real(countsd(ibin))

 enddo

 close(20)

 !print *, 'Binary phase written to ', phasefile
 !print *, 'Gaussian results written to ', gaussfile




end subroutine do_analysis

end module analysis
