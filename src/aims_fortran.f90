!------------------------------------------------------------------------------
! This method takes two sets of modes and produces a third set by combining
! the two first sets using the interpolation weights coef1 and coef2. 
! Only modes which are common to both sets will be present in the output set.
!------------------------------------------------------------------------------
! INPUT:
! coef1 = first interpolation coefficient
! narr1 = n values of the first set of modes
! larr1 = l values of the first set of modes
! farr1 = frequencies of the first set of modes
! iarr1 = inertias of the first set of modes
! n1    = number of modes in the first set of modes
! coef2 = second interpolation coefficient
! narr2 = n values of the second set of modes
! larr2 = l values of the second set of modes
! farr2 = frequencies of the second set of modes
! iarr2 = inertias of the second set of modes
! n2    = number of modes in the second set of modes
!
! OUTPUT: 
! narr3 = n values of the resultant set of modes
! larr3 = l values of the resultant set of modes
! farr3 = frequencies of the resultant set of modes
! iarr3 = inertias of the resultant set of modes
! n3    = number of modes in the resultant set of modes
!------------------------------------------------------------------------------
      subroutine combine_modes(coef1,narr1,larr1,farr1,iarr1, n1, &
                               coef2,narr2,larr2,farr2,iarr2, n2, &
                               narr3,larr3,farr3,iarr3, n3)

      implicit none

      ! input arguments
      integer, intent(in) :: n1, n2
      integer(kind=2), intent(in) :: narr1(0:n1-1), narr2(0:n2-1)  ! np.int16
      integer(kind=1), intent(in) :: larr1(0:n1-1), larr2(0:n2-1)  ! np.int8
      real(kind=8), intent(in) :: coef1, coef2
      real(kind=8), intent(in) :: farr1(0:n1-1),iarr1(0:n1-1)
      real(kind=8), intent(in) :: farr2(0:n2-1),iarr2(0:n2-1)

      ! output arguments
      integer, intent(inout) :: n3
      !f2py intent(in,out) n3
      integer(kind=2), intent(inout) :: narr3(0:n3-1)
      !f2py intent(in,out) narr3
      integer(kind=1), intent(inout) :: larr3(0:n3-1)
      !f2py intent(in,out) larr3
      real(kind=8), intent(inout) :: farr3(0:n3-1),iarr3(0:n3-1)
      !f2py intent(in,out) farr3, iarr3

      integer i1, i2, i3

  
      i1 = 0
      i2 = 0
      i3 =-1 
      do while ((i1.lt.n1).and.(i2.lt.n2))
        if (larr1(i1).lt.larr2(i2)) then
          i1 = i1+1
          cycle
        endif
        if (larr1(i1).gt.larr2(i2)) then
          i2 = i2+1
          cycle
        endif
        if (narr1(i1).lt.narr2(i2)) then
          i1 = i1+1
          cycle
        endif
        if (narr1(i1).gt.narr2(i2)) then
          i2 = i2+1
          cycle
        endif
        ! now the two modes have the same n and l values:

        ! first increment i3
        i3 = i3 + 1
        !! sanity check (to be removed eventually)
        !if (i3.gt.n3) then
        !  print*,i3,n3
        !  stop "array3 too small in aims_fortran.combine"
        !endif

        narr3(i3) = narr1(i1)
        larr3(i3) = larr1(i1)
        farr3(i3) = coef1*farr1(i1) + coef2*farr2(i2) 
        iarr3(i3) = coef1*iarr1(i1) + coef2*iarr2(i2) 

        i1 = i1+1
        i2 = i2+1
      enddo
      n3 = i3+1

      end subroutine combine_modes

!------------------------------------------------------------------------------
! This method finds a mapping between an observed set of modes and a
! theoretical set of modes from a model.  The mapping is obtained from
! matching (l,n) values.
!------------------------------------------------------------------------------
! INPUT:
! nobs     = n values of the observed modes
! lobs     = l values of the observed modes
! size_obs = number of observed modes
! nmod     = n values of the theoretical modes
! lmod     = l values of the theoretical modes
! size_mod = number of theoretical modes
!
! OUTPUT: 
! mode_map = mapping between observed and theoretical modes
! nmissing = number of unmatched observed modes
!------------------------------------------------------------------------------
      subroutine find_map_n(nobs, lobs, size_obs, nmod, lmod, &
                            size_mod, mode_map, nmissing)

      implicit none

      ! input arguments
      integer, intent(in) :: size_obs, size_mod
      integer(kind=2), intent(in) :: nmod(0:size_mod-1), nobs(0:size_obs-1)
      integer(kind=1), intent(in) :: lmod(0:size_mod-1), lobs(0:size_obs-1)

      ! output arguments
      integer, intent(out) :: nmissing
      integer, intent(inout) :: mode_map(0:size_obs-1)
      !f2py intent(in,out) mode_map

      integer iobs, imod

      !initialisation
      iobs = 0
      imod = 0
      nmissing = 0
      mode_map = -1

      ! NOTE: this assumes the observed modes are sorted according to (l,n)
      do while((iobs.lt.size_obs).and.(imod.lt.size_mod))

        if (lmod(imod).lt.lobs(iobs)) then
           imod = imod + 1
           cycle
        endif

        if (lmod(imod).gt.lobs(iobs)) then
           iobs = iobs + 1
           nmissing = nmissing + 1
           cycle
        endif

        if (nmod(imod).lt.nobs(iobs)) then
           imod = imod + 1
           cycle
        endif

        if (nmod(imod).gt.nobs(iobs)) then
           iobs = iobs + 1
           nmissing = nmissing + 1
           cycle
        endif

        mode_map(iobs) = imod

        imod = imod + 1
        iobs = iobs + 1
      enddo

      ! keep track of the number of missing modes:
      nmissing = nmissing + (size_obs-iobs)
      end subroutine find_map_n
        
!------------------------------------------------------------------------------
! This method finds a mapping between an observed set of modes and a
! theoretical set of modes from a model.  The mapping is obtained from
! l values and frequency proximity.
!------------------------------------------------------------------------------
! INPUT:
! fobs     = frequencies of the observed modes (units = muHz)
! lobs     = l values of the observed modes
! size_obs = number of observed modes
! fmod     = frequencies of the theoretical modes (units = muHz)
! lmod     = l values of the theoretical modes
! ind      = index array for which the theoretical modes are sorted according
!            to (l,freq)
! size_mod = number of theoretical modes
!
! OUTPUT: 
! mode_map = mapping between observed and theoretical modes
! nmissing = number of unmatched observed modes
!------------------------------------------------------------------------------
      subroutine find_map_freq(fobs, lobs, size_obs, fmod, lmod, ind, &
                            size_mod, mode_map, nmissing) 

      implicit none

      ! input arguments
      integer, intent(in) :: size_obs, size_mod
      integer, intent(in) :: ind(0:size_mod-1)
      integer(kind=1), intent(in) :: lmod(0:size_mod-1), lobs(0:size_obs-1)
      real(kind=8), intent(in) ::    fmod(0:size_mod-1), fobs(0:size_obs-1)

      ! output arguments
      integer, intent(out) :: nmissing
      integer, intent(inout) :: mode_map(0:size_obs-1)
      !f2py intent(in,out) mode_map

      real(kind=8) :: diff0, diff1, diffmin
      real(kind=8), allocatable :: diffs(:)
      integer :: iobs, imod, imod0, imod1, val, j, jmin

      !initialisation
      iobs = 0
      imod = 0
      nmissing = 0
      mode_map = -1
      allocate(diffs(0:size_obs-1))

      ! This is a two step process:
      !   1. The first step finds the nearest theoretical mode to
      !      each observed mode.
      !   2. When a theoretical mode has several observed modes that
      !      correspond to it, only the nearest observed mode is kept.
      !      The other observed modes are no longer mapped to the
      !      theoretical mode and become "missing".

      ! Step 1: find nearest theoretical modes:
      do while((iobs.lt.size_obs).and.(imod.lt.size_mod))
        imod1 = ind(imod)

        if (lmod(imod1).lt.lobs(iobs)) then
           imod = imod + 1
           cycle
        endif

        if (lmod(imod1).gt.lobs(iobs)) then
           iobs = iobs + 1
           nmissing = nmissing + 1
           cycle
        endif
        
        if (fmod(imod1).lt.fobs(iobs)) then
           imod = imod + 1
           cycle
        endif

        diff1 = fmod(imod1) - fobs(iobs)
        if (imod.gt.0) then
          imod0 = ind(imod-1)
          if (lmod(imod0).eq.lobs(iobs)) then
            diff0 = fobs(iobs) - fmod(imod0)
            if (diff0 < diff1) then
              mode_map(iobs) = imod0
              diffs(iobs) = diff0
            else
              mode_map(iobs) = imod1
              diffs(iobs) = diff1
            endif
          else
            mode_map(iobs) = imod1
            diffs(iobs) = diff1
          endif
        else
          mode_map(iobs) = imod1
          diffs(iobs) = diff1
        endif

        iobs = iobs + 1
      enddo

      ! keep track of the number of missing modes:
      nmissing = nmissing + (size_obs-iobs)

      ! Step 2: filter out cases where more than one observed mode 
      !         corresponds to the same theoretical mode.
      do iobs = 0, size_obs-1
        if (mode_map(iobs).eq.-1) cycle

        val = mode_map(iobs)
        jmin = iobs
        diffmin = diffs(iobs)
        do j=iobs+1,size_obs-1
          if (mode_map(j).ne.val) exit
          if (diffs(j).lt.diffmin) then
            jmin = j
            diffmin = diffs(j)
          endif
        enddo
        mode_map(iobs:j-1) = -1
        mode_map(jmin) = val
        nmissing = nmissing + j-iobs-1
      enddo

      deallocate(diffs)

      end subroutine find_map_freq

!------------------------------------------------------------------------------
! This method calculates differences between observed and theoretical frequency
! combinations.
!------------------------------------------------------------------------------
! INPUT:
! freq           = theoretical frequencies
! mode_map       = mapping between observed and theoretical frequencies
! ncoeff         = number of terms in the frequency combinations
! coeff          = coefficients in the frequency combinations
! indices        = mode indices in the frequency combinations
! nmodes         = number of theoretical modes
! nobs           = number of observed modes
! ncomb          = number of frequency combinations
! nmax           = maximum number of terms in a frequency combination
!
! OUTPUT:
! values         = differences between observed and theoretical frequency
!                  combinations
!------------------------------------------------------------------------------
      subroutine compare_frequency_combinations(freq,mode_map,&
                    ncoeff,coeff,indices,nmodes,nobs,ncomb,nmax,values)

      implicit none

      integer, intent(in) :: nmodes, nobs, ncomb, nmax
      real(kind=8), intent(in) :: freq(0:nmodes-1),&
                                  coeff(0:nmax-1,0:ncomb-1)
      integer, intent(in) :: mode_map(0:nobs-1), ncoeff(0:ncomb-1),&
                             indices(0:nmax-1,0:ncomb-1)
      real(kind=8), intent(inout) :: values(0:ncomb-1)
      !f2py intent(in,out) values

      integer i, j

      ! NOTE: values has already been initialised to 0 in AIMS.py
      do i=0, ncomb-1
        do j=0,ncoeff(i)-1
          values(i) = values(i) + coeff(j,i)*freq(mode_map(indices(j,i)))
        enddo
      enddo

      end subroutine compare_frequency_combinations

!------------------------------------------------------------------------------
! This method reads a set of modes from a file in "agsm" format, a fortran
! binary format used by ADIPLS.
!------------------------------------------------------------------------------
! INPUT:
! filename       = name of the agsm file
! npositive      = specifies whether to retain only modes with n >= 0
! below_cutoff   = specifies whether to only retain modes below the cutoff frequency
! freqlim        = upper frequency bound. Mode with frequencies above this bound
!                  are discarded
!
! OUTPUT:
! narr           = n values of the modes which are read
! larr           = l values of the modes which are read
! farr           = frequencies of the modes which are read
! iarr           = inertias of the modes which are read
! nn             = number of modes which are read
! exceed_freqlim = specifies whether a mode exceeded the frequency limit
!------------------------------------------------------------------------------
      subroutine read_file_agsm(filename,npositive,below_cutoff, &
                                freqlim,narr,larr,farr,iarr,nn,  &
                                exceed_freqlim)

      implicit none
      integer, parameter :: nmax = 1000
      character(len=*), intent(in) :: filename
      logical, intent(in) :: npositive, below_cutoff
      real(kind=8), intent(in) :: freqlim
      real(kind=8), intent(out) :: farr(nmax), iarr(nmax)
      integer(kind=2), intent(out) :: narr(nmax)
      integer(kind=1), intent(out) :: larr(nmax)
      logical, intent(out) :: exceed_freqlim
      integer, intent(out) :: nn
      real(kind=8) :: cs(50)
      integer :: ics(8), ntot, i, j

      ! transform the last part of the cs array into integers:
      equivalence(ics(1), cs(39))

      ntot = 0
      open(unit=31,file=filename,status='old',form='unformatted')
      ! note: the convert keyword can be used to change the endian
      do
        read(31,end=10) (cs(i),i=1,50)
        ntot = ntot + 1
      enddo
10    close(31)
      
      ! sanity check:
      if (ntot.gt.nmax) then
        stop "Please increase nmax in aims_fortran.read_file_agsm"
      endif

      exceed_freqlim = .False.
      open(unit=31,file=filename,status='old',form='unformatted')
      nn = 1 
      do i=1,ntot
        read(31) (cs(j),j=1,50)
        if (below_cutoff.and.(ics(5).ne.10010)) cycle
        if (npositive.and.(cs(19).lt.-0.1d0)) cycle
        if (cs(37)*1d3.gt.freqlim) then
          exceed_freqlim = .True.
          cycle
        endif
        larr(nn) = nint(cs(18))
        narr(nn) = nint(cs(19))
        farr(nn) = cs(37)*1d3 ! Richardson frequency in muHz
        !farr(nn) = cs(27)*1d3 ! variational frequency in muHz
        iarr(nn) = cs(24)
        nn = nn + 1
      enddo
      close(31)
      nn = nn - 1

      end subroutine read_file_agsm

!------------------------------------------------------------------------------
! This method finds the age parameter for an input physical age, using a
! a weighted combination of evolutionary tracks (which contain a piecewise
! affine relation between the age parameter and the physical age).
!------------------------------------------------------------------------------
! IMPORTANT: this method assumes that each track has its own set of age
!            parameters, i.e. they need not be the same.
!------------------------------------------------------------------------------
! INPUT:
! tau_array = array with age parameters for each track
! age_array = array with physical ages for each track
! coef      = array with interpolation coefficients (or weights) for the tracks
! sze       = array with sizes of each evolutionary track
! n         = maximum track size
! ntracks   = number of tracks
! trget     = target value on age
!
! OUTPUT:
! tau       = output age parameter
! indices   = indices right below relevant age and additional work space
! weights   = weights assigned to lower indices along tracks
!------------------------------------------------------------------------------
       subroutine find_tau(tau_array, age_array, coef, sze, n, ntracks, trget, &
                           tau, indices, weights)

       implicit none

       ! input arguments
       integer, intent(in) :: n, ntracks
       integer, intent(in) :: sze(ntracks)
       real(kind=8), intent(in) :: tau_array(n,ntracks), age_array(n,ntracks)
       real(kind=8), intent(in) :: coef(ntracks), trget

       ! output arguments
       real(kind=8), intent(inout) :: tau
       !f2py intent(in,out) tau
       integer, intent(inout) :: indices(ntracks,3)
       !f2py intent(in,out) indices
       real(kind=8), intent(inout) :: weights(ntracks)
       !f2py intent(in,out) weights

       real(kind=8) :: age_min, age_max, age_mid, weight
       real(kind=8) :: tau_min, tau_max, tau_mid, tau_new, eta
       real(kind=8), external :: find_age_single
       logical, external :: test_indices
       integer :: j, ndx

       ! initialisation
       indices = -1
       weights = 0d0

       ! test lower bound
       tau_min = tau_array(1,1)
       do j=2,ntracks
         if (tau_array(1,j).gt.tau_min) tau_min = tau_array(1,j)
       enddo
       age_min = 0d0
       do j=1,ntracks
         age_min = age_min + coef(j)*find_age_single(tau_array(:,j), &
                   age_array(:,j),n,1,sze(j),tau_min,ndx,weight)
         indices(j,1) = ndx
       enddo
       if (trget.lt.age_min) return

       ! test upper bound
       tau_max = tau_array(sze(1),1)
       do j=2,ntracks
         if (tau_array(sze(j),j).lt.tau_max) tau_max = tau_array(sze(j),j)
       enddo
       age_max = 0d0
       do j=1,ntracks
         age_max = age_max + coef(j)*find_age_single(tau_array(:,j), &
                   age_array(:,j),n,1,sze(j),tau_max,ndx,weight)
         indices(j,2) = ndx+1
       enddo
       if (trget.gt.age_max) return

       do while(.not.test_indices(indices,ntracks))
         tau_new = (tau_min+tau_max)/2.0

         ! test degenerate case where solution lies on one of the
         ! values of tau for at least one of the tracks:
         if (tau_new.eq.tau_mid) then
           !print*,"Dicho degenerate"
           ! interpolating between tau_min and tau_max may be risky, hence
           ! we simply use the mid point as it probably is more robust:
           tau = tau_new
           age_mid  = 0d0
           do j=1,ntracks
             age_mid = age_mid + coef(j)*find_age_single(tau_array(:,j), &
               age_array(:,j),n,indices(j,1),indices(j,2),tau,ndx,weight)
             indices(j,1) = ndx
             indices(j,2) = ndx+1
             weights(j) = weight 
           enddo
           return
         endif

         tau_mid = tau_new
         age_mid = 0d0
         do j=1,ntracks
           age_mid = age_mid + coef(j)*find_age_single(tau_array(:,j), &
             age_array(:,j),n,indices(j,1),indices(j,2),tau_mid,ndx,weight)
           indices(j,3) = ndx
         enddo
         if (trget.gt.age_mid) then
           age_min = age_mid
           tau_min = tau_mid
           indices(:,1) = indices(:,3)
         else
           age_max = age_mid
           tau_max = tau_mid
           indices(:,2) = indices(:,3)+1
         endif
       enddo

       ! save results
       eta = (age_max-trget)/(age_max-age_min)
       tau = eta*tau_min + (1d0-eta)*tau_max
       age_mid  = 0d0
       do j=1,ntracks
         age_mid = age_mid + coef(j)*find_age_single(tau_array(:,j), &
           age_array(:,j),n,indices(j,1),indices(j,2),tau,ndx,weight)
         indices(j,1) = ndx
         indices(j,2) = ndx+1
         weights(j) = weight 
       enddo

       end subroutine find_tau

!------------------------------------------------------------------------------
! Function which tests if series of lower and upper indices differ by one.
! This provides the stop condition on the dichotomy loop in find_tau.
!------------------------------------------------------------------------------
! INPUT:
! indices = array with lower and upper indices
! ntracks = number of tracks
!
! NOTE: the second index on indices ranges from 1 to 3:
!       1 corresponds to the lower indices
!       2 corresponds to the upper indices
!       3 is additional work space needed by find_tau (and not used here)
!------------------------------------------------------------------------------
       logical function test_indices(indices,ntracks)

       implicit none
       integer, intent(in) :: ntracks
       integer, intent(in) :: indices(ntracks,3)
       integer :: i

       test_indices = .false.
       do i=1,ntracks
         if (indices(i,2)-indices(i,1).gt.1) return
       enddo
       test_indices = .true.
       return
       end function test_indices

!------------------------------------------------------------------------------
! Find physical age for a given input age parameter, for a single track.
!------------------------------------------------------------------------------
! INPUT:
! tau    = array with age parameter
! age    = array with physical age
! n      = size of tau and age arrays
! nstart = starting index at which to look for relevant age
! nstop  = stoping index at which to look for relevant age
! trget  = target value on tau
!
! OUTPUT (in addition to physical age):
! ndx    = index right below relevant age
! weight = weight assigned to tau(ndx)
!------------------------------------------------------------------------------
! IMPORTANT: the subroutine assumes tau and age are sorted
!------------------------------------------------------------------------------
       real(kind=8) function find_age_single(tau, age, n, nstart, &
                                       nstop, trget, ndx, weight)

       implicit none
       integer, intent(in) :: n, nstart, nstop
       real(kind=8), intent(in) :: tau(n), age(n), trget
       integer, intent(out) :: ndx
       real(kind=8), intent(out) :: weight
       integer :: istart, istop, imid

       ! initialisation
       ndx = -1
       weight = 0d0
       find_age_single = 0d0 

       ! easy exit
       if (trget.lt.tau(nstart)) return
       if (trget.gt.tau(nstop))  return

       istart = nstart
       istop = nstop
       do while((istop-istart).gt.1)
         imid = (istart+istop)/2
         if (trget.lt.tau(imid)) then
           istop = imid
         else
           istart = imid
         endif
       enddo

       ndx = istart
       weight = (tau(istop)-trget)/(tau(istop)-tau(istart))
       find_age_single = weight*age(istart) + (1d0-weight)*(age(istop))
       return
       end function find_age_single
