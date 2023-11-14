!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for discs
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, bindiscanalysisutils, infile_utils, io, part, physcon
!
 use discanalysisutils, only:bindisc_analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'disc'
 public :: do_analysis

 integer, parameter :: nr = 40
 real,dimension(nr) :: twist,twistprev

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyz,pmass,npart,time,iunit)
 use io,      only:fatal
 use physcon, only:pi
 use part,    only:xyzmh_ptmass,vxyz_ptmass,nptmass
 use dim,     only:gr
 use infile_utils, only:open_db_from_file,read_inopt,close_db,inopts
 character(len=*), intent(in) :: dumpfile
 real,             intent(inout) :: xyzh(:,:),vxyz(:,:)
 real,             intent(inout) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile
 character(len=9) :: output
 character(len=20) :: filename
 character(len=20) :: discprefix
 integer :: i,j,ierr,iline
 real :: R_in,R_out,H_R,p_index,q_index,M_star
 real :: G,rmin,rmax
 real :: tilt(nr),Lx(nr),Ly(nr),Lz(nr)
 real :: Lx_mean(nr),Ly_mean(nr),Lz_mean(nr)
 real :: rad(nr),h_smooth(nr),sigma(nr),H(nr)
 real :: unitlx(nr),unitly(nr),unitlz(nr),ecc(nr)
 real :: psi(nr),tilt_acc(nr)
 real :: incl(nr),Omega(nr),omega_(nr)
 integer :: ninbin(nr),iexternalforce_read
 logical :: assume_Ltot_is_same_as_zaxis,iexist,ref_from_only_sinks
 type(inopts), allocatable :: db(:)

 integer, parameter :: iparams = 10
 integer, parameter :: iprec   = 24
 logical :: do_precession,ifile

 !real, dimension(SIZE(xyzmh_ptmass, DIM=1),nptmass-1) :: xyzmh_bin_ptmass  !mod
 !real, dimension(SIZE(vxyz_ptmass, DIM=1),nptmass-1) :: vxyz_bin_ptmass
 integer, allocatable :: sinks_CoM(:), sinks_L(:)
 allocate(sinks_CoM(3))
 allocate(sinks_L(3))

! This variable should be set to false for any discs that use sink particles to set
! the potential, any discs that have a warp or any time precession is measured
! For any setup that uses iexternalforce and assumes that the vast majority of the angular
! momentum is held by the central potential, this should be set to true
 assume_Ltot_is_same_as_zaxis = .false.

! Option for if you want precession files printed
 do_precession = .false.

! Sinks to be considered
 sinks_CoM = (/1,3/)
 sinks_L = (/1,3/)

! Need the CoM and the L reference plane be calculated from the sinks only ?
 ref_from_only_sinks = .true.


! Construction of the sinks list if input is empty
! if (SIZE(sinks_index_)==0) then
!    allocate(sinks_index(nptmass))
!    do i = 1, nptmass
!        sinks_index(i) = i
!    enddo
! else
!    sinks_index=sinks_index_
! endif

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 write(output,"(a4,i5.5)") 'angm',numfile
 write(*,'("Output file name is ",A)') output

! Assuming G=1
 write(*,*)
 write(*,'("ASSUMING G==1")')
 G = 1.0

 iline = index(dumpfile,'_')
 discprefix = dumpfile(1:iline-1)
 inquire(file=trim(discprefix)//'.discparams', exist=ifile)
 if (ifile) then
    call read_discparams(trim(discprefix)//'.discparams',R_in,R_out,H_R,p_index,q_index,M_star,iparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read .discparams file')
 else
    call read_discparams('discparams.list',R_in,R_out,H_R,p_index,q_index,M_star,iparams,ierr)
    if (ierr /= 0) call fatal('analysis','could not open/read discparams.list')
 endif

! Print out the parameters
 write(*,*)
 write(*,'("Parameters are:")')
 write(*,*) 'R_in    = ',R_in
 write(*,*) 'R_out   = ',R_out
 write(*,*) 'H/R_ref = ',H_R
 write(*,*) 'p_index = ',p_index
 write(*,*) 'q_index = ',q_index
 write(*,*) 'M_star  = ',M_star
 write(*,*)
 write(*,*)

! Setup rmin and rmax for the analysis
 rmin = 0.1*R_in
 rmax = R_out

! If do_precession is true and then this variable should be false, so do a check
 if (do_precession) assume_Ltot_is_same_as_zaxis = .false.

 ! Check, if iexternalforce > 0 from the *.in file
 ! this value should be set to true (or if GR is used)
 ! Open *.in file and read the iexternalforce variable
 iexternalforce_read = 0
 inquire(file=trim(discprefix)//'.in', exist=ifile)
 if (ifile) then
    call open_db_from_file(db,trim(discprefix)//'.in',iparams,ierr)
    call read_inopt(iexternalforce_read,'iexternalforce',db,ierr)
    call close_db(db)
 endif
 if (iexternalforce_read > 0 .or. gr) then
    assume_Ltot_is_same_as_zaxis = .true.
    print*,'Resetting assume_Ltot_is_same_as_zaxis=.true. in analysis'
 endif

!consider only the inner binary for the analysis : put the inner binary as the center of the system and ignore the tertiary body !mod
! write(*,*) SIZE(xyzmh_ptmass, DIM=1), SIZE(vxyz_ptmass, DIM=1), nptmass
!
! do i=1,SIZE(xyzmh_ptmass, DIM=1)
!     do j=1,nptmass
!         if (j<2) then
!            xyzmh_bin_ptmass(i,j) = xyzmh_ptmass(i,j)
!         endif
!         if (j>2) then
!            xyzmh_bin_ptmass(i,j-1) = xyzmh_ptmass(i,j)
!         endif
!     enddo
! enddo
!
! do i=1,SIZE(vxyz_ptmass, DIM=1)
!     do j=1,nptmass
!         if (j<2) then
!             vxyz_bin_ptmass(i,j) = vxyz_ptmass(i,j)
!         endif
!         if (j>2) then
!             vxyz_bin_ptmass(i,j-1) = vxyz_ptmass(i,j)
!         endif
!     enddo
! enddo
!
!
! nptmass = 2

 call bindisc_analysis(xyzh,vxyz,npart,pmass,time,nr,rmin,rmax,G,M_star,&             !mod
                     tilt,tilt_acc,twist,twistprev,psi,H,rad,h_smooth,sigma,unitlx,unitly,unitlz,&
                     Lx,Ly,Lz,Lx_mean,Ly_mean,Lz_mean,&
                     ecc,ninbin,assume_Ltot_is_same_as_zaxis,xyzmh_ptmass,vxyz_ptmass,nptmass,&
                     incl=incl,Omega=Omega,omega_=omega_,&
                     ref_is_only_sinks=ref_from_only_sinks,sinks_CoM=sinks_CoM,sinks_L=sinks_L)

 open(iunit,file=output)
 write(iunit,'("# Analysis data at t = ",es20.12)') time
 write(iunit,"('#',18(1x,'[',i2.2,1x,a11,']',2x))") &
       1,'radius', &
       2,'sigma', &
       3,'<h>/H', &
       4,'lx', &
       5,'ly', &
       6,'lz', &
       7,'tilt', &
       8,'twist', &
       9,'psi', &
       10,'H/R', &
       11,'|e|', &
       12,'i', &
       13,'Omega', &
       14,'omega', &
       15,'Lx', &
       16,'Ly', &
       17,'Lz', &
       18,'N'

 do i = 1,nr
    if (ninbin(i) > 0) then
       write(iunit,'(18(es18.10,1X))') rad(i),sigma(i),h_smooth(i),unitlx(i),unitly(i),unitlz(i),&
                                         tilt(i),twistprev(i),psi(i),H(i)/rad(i),ecc(i),&
                                         incl(i),Omega(i),omega_(i),&
                                         Lx_mean(i),Ly_mean(i),Lz_mean(i),&
                                         real(ninbin(i))
    endif



! Printing time and twist for each radius bin
    if (do_precession) then
       write(filename,"(a,i3.3)")"precess",i
       inquire(file=filename,exist=iexist)
       if (.not.iexist .or. numfile==0) then
          open(unit=iprec,file=filename,status="replace")
          write(iprec,'("# tilt and twist with time for r = ",es18.10)') rad(i)
          write(iprec,"('#',6(1x,'[',i2.2,1x,a11,']',2x))") &
               1,'rad', &
               2,'time', &
               3,'tilt', &
               4,'twist', &
               5,'tot twist', &
               6,'|e|'
       else
          open(unit=iprec,file=filename,status="old",position="append")
       endif
       write(iprec,'(6(es18.10,1X))') rad(i),time,tilt(i),twist(i),twistprev(i),ecc(i)
       close(unit=iprec)
    endif

 enddo

 close(iunit)

end subroutine do_analysis

!----------------------------------------------------------------
!+
!  Read disc information from discparams.list file
!+
!----------------------------------------------------------------
subroutine read_discparams(filename,R_in,R_out,H_R,p_index,q_index,M_star,iunit,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 real,             intent(out) :: R_in,R_out,H_R,p_index,q_index,M_star
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 type(inopts), allocatable :: db(:)

! Read in parameters from the file discparams.list
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(R_in,'R_in',db,ierr)
 if (ierr /= 0) return
 call read_inopt(R_out,'R_out',db,ierr)
 if (ierr /= 0) return
 call read_inopt(H_R,'H/R_ref',db,ierr)
 if (ierr /= 0) return
 call read_inopt(p_index,'p_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(q_index,'q_index',db,ierr)
 if (ierr /= 0) return
 call read_inopt(M_star,'M_star',db,ierr)
 if (ierr /= 0) return
 call close_db(db)

end subroutine read_discparams

end module analysis
