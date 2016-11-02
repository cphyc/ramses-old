! Hooks for RAMSES
module hooks
  implicit none

contains

  !---------------------------------------------------------------------------------
  ! Hook after flux computation
  !---------------------------------------------------------------------------------
  subroutine post_flux_hook(ind_grid, flux, ilevel)
    use amr_commons
    use hydro_commons

    integer, dimension(1:nvector), intent(in)                                           :: ind_grid
    real(dp),dimension(1:nvector, if1:if2, jf1:jf2, kf1:kf2, 1:nvar, 1:ndim),intent(in) :: flux
    integer, intent(in)                                                                 :: ilevel

    !######################
    ! Customize here
    !######################
    if (MC_tracer) then
       call move_tracers_oct(ind_grid, flux, ilevel)
    end if
  end subroutine post_flux_hook

  !---------------------------------------------------------------------------------
  ! Hooks before and after grid creation
  !---------------------------------------------------------------------------------
  subroutine pre_kill_grid_hook(ind_cell, ilevel, nn, ibound, boundary_region)
    use amr_commons
    use pm_commons

    integer, intent(in) :: nn, ilevel, ibound
    logical, intent(in) :: boundary_region
    integer, intent(in), dimension(1:nvector) :: ind_cell

    !######################
    ! Customize here
    !######################

    integer :: ipart, igrid, i, j, dim
    logical, save :: firstCall = .true.
    real(dp), save :: skip_loc(3), nx_loc
    real(dp) :: dx, scale

    if (firstCall) then
       skip_loc=(/0.0d0, 0.0d0, 0.0d0/)
       if(ndim>0) skip_loc(1) = dble(icoarse_min)
       if(ndim>1) skip_loc(2) = dble(jcoarse_min)
       if(ndim>2) skip_loc(3) = dble(kcoarse_min)

       nx_loc=(icoarse_max-icoarse_min+1)

       firstCall = .false.
    end if

    scale = boxlen/dble(nx_loc)
    dx = 0.5D0**ilevel

    if (MC_tracer) then
       ! For all particles, recenter them in the center of the grid that's being deleted
       ! (becoming a cell)
       do j = 1, nn
          igrid = son(ind_cell(j))
          ipart = headp(igrid)

          do i = 1, numbp(igrid)
             if (mp(ipart) == 0d0) then
                ! print*, 'pre_kill_grid_hook', ipart, igrid, ilevel
                ! print*, 'pre_kill_grid_hook', xp(ipart, :)
                if(ipart == -1) print*, 'pre_kill_grid_hook'
                do dim = 1, ndim
                   xp(ipart, dim) = (xg(igrid, dim) - skip_loc(dim)) * scale
                end do
                ! print*, 'pre_kill_grid_hook', xp(ipart, :)
             end if
             ipart = nextp(ipart)
          end do
       end do
    end if
  end subroutine pre_kill_grid_hook

  subroutine post_kill_grid_hook(ind_cell, ilevel, nn, ibound, boundary_region)
    use amr_commons

    integer, intent(in) :: nn, ilevel, ibound
    logical, intent(in) :: boundary_region
    integer, intent(in), dimension(:) :: ind_cell

    !######################
    ! Customize here
    !######################

  end subroutine post_kill_grid_hook

  !---------------------------------------------------------------------------------
  ! Hooks before and after grid creation
  !---------------------------------------------------------------------------------
  subroutine pre_make_grid_fine_hook(ind_grid, ind_cell, ind, &
       ilevel, nn, ibound, boundary_region)
    use amr_commons
    integer, intent(in) :: nn, ind, ilevel, ibound
    logical, intent(in) :: boundary_region
    integer, dimension(1:nvector), intent(in) :: ind_grid, ind_cell

    !######################
    ! Customize here
    !######################
  end subroutine pre_make_grid_fine_hook

  subroutine post_make_grid_fine_hook(ind_grid, ind_cell, ind, &
       ilevel, nn, ibound, boundary_region)
    use amr_commons
    use hydro_commons
    use pm_commons
    use random
    integer, intent(in) :: nn, ind, ilevel, ibound
    logical, intent(in) :: boundary_region
    integer, dimension(1:nvector), intent(in) :: ind_grid, ind_cell

    !######################
    ! Customize here
    !######################
    logical, save :: firstCall = .true.
    real(dp), save :: skip_loc(3), nx_loc
    real(dp) :: dx, dxcoarse, scale

    real(dp) :: x, rand, mass(0:twotondim), xc(1:3)
    integer :: j, i, dim, fgrid, igrid, ipart, part_dir, ison, iskip, icell
    integer :: ix, iy, iz, loc(1:3)
    logical :: ok


    ! Save some time here...
    if (firstCall) then
       skip_loc=(/0.0d0, 0.0d0, 0.0d0/)
       if(ndim>0) skip_loc(1) = dble(icoarse_min)
       if(ndim>1) skip_loc(2) = dble(jcoarse_min)
       if(ndim>2) skip_loc(3) = dble(kcoarse_min)

       nx_loc=(icoarse_max-icoarse_min+1)

       firstCall = .false.
    end if

    scale = boxlen/dble(nx_loc)
    dx = 0.5D0**ilevel           ! dx of the new level
    dxcoarse = 0.5D0**(ilevel-1) ! dx of the previous level

    if (MC_tracer) then
       ! Compute the expected location of particles relative to xg in dx units
       loc(3) = (ind-1) / 4
       loc(2) = (ind-1-loc(3)*4) / 2
       loc(1) = (ind-1-loc(3)*4-loc(2)*2)

       do j = 1, nn
          fgrid = ind_grid(j)
          igrid = son(ind_cell(j))
          ipart = headp(fgrid)

          ! Load masses
          mass(0) = uold(ind_cell(j), 1)
          do ison = 1, twotondim
             iskip = ncoarse + (ison-1)*ngridmax
             icell = iskip + igrid
             mass(ison) = uold(icell, 1)
          end do

          do i = 1, numbp(fgrid)
             if (mp(ipart) == 0d0) then

                ! print*, 'post_make_grid_fine_hook', ipart, igrid, ilevel

                ok = .true.
                ! Check whether the particle was in the refined cell
                do dim = 1, ndim
                   x = (xp(ipart, dim) / scale - xg(fgrid, dim) + skip_loc(dim)) / dxcoarse &
                        + 0.5d0
                   ok = ok .and. (int(x) == loc(dim))

                   ! It is possible that x is not 0 or 1 if the particle has already been
                   ! moved by the routine. In this case, we don't need to move it again
                   ! print*, x, loc(dim)
                end do

                ! If the particle is in refined cell, spread it accordingly
                if (ok) then
                   ! print*, 'post_make_grid_fine_hook', ipart, igrid, ilevel
                   ! print*, 'post_make_grid_fine_hook', xp(ipart, :)

                   ! Pick a random direction
                   call ranf(localseed, rand)

                   do ison = 1, twotondim
                      if (rand < mass(ison) / mass(0)) then
                         ! Move particle to center of new cells
                         iz=(ison-1)/4
                         iy=(ison-1-4*iz)/2
                         ix=(ison-1-2*iy-4*iz)
                         xc(1) = (dble(ix)-0.5D0)*dx
                         xc(2) = (dble(iy)-0.5D0)*dx
                         xc(3) = (dble(iz)-0.5D0)*dx
                         if(ipart == -1) print*, 'post_make_grid_fine_hook'
                         do dim = 1, ndim
                            xp(ipart, dim) = (xg(igrid, dim) + xc(dim) - skip_loc(dim)) * scale
                         end do
                         rand = 1
                      else
                         rand = rand - mass(ison) / mass(0)
                      end if
                   end do
                   ! print*, 'post_make_grid_fine_hook', xp(ipart, :)

                end if
             end if
             ipart = nextp(ipart)
          end do
       end do
    end if
  end subroutine post_make_grid_fine_hook

end module hooks

!################################################################
! Write your hooks here
!################################################################
subroutine move_tracers_oct(ind_grid, fluxes, ilevel)
  use amr_parameters      ! nx, ny, nz, dp, ngridmax, nvector, …
  use amr_commons, only   : ncoarse, father, xg, son, myid
  use pm_commons, only    : xp, headp, mp, localseed, numbp, nextp, move_flag
  use random, only        : ranf
  use hydro_commons, only : uold, if1, if2, jf1, jf2, kf1, kf2, nvar

  implicit none

  integer,dimension(1:nvector), intent(in)   :: ind_grid
  real(dp),dimension(1:nvector, if1:if2, jf1:jf2, kf1:kf2, 1:nvar, 1:ndim),intent(in)::fluxes
  integer, intent(in) :: ilevel

  !-------------------------------------------------------------------
  ! Implement a Monte Carlo sampling of the gas properties following
  ! the fluxes of mass from oct to oct.
  !-------------------------------------------------------------------
  integer,dimension(1:nvector,0:twondim) :: ind_ngrid
  real(dp) :: flux, Fout, mass, rand, masses(0:twotondim/2)
  integer,dimension(1:nvector, 1:twotondim) :: ind_cell
  integer,dimension(1:nvector, 1:twotondim, 0:twondim) :: ind_ncell
  integer,dimension(1:nvector, 0:twondim) :: tmp_ncell
  integer, dimension(1:nvector) :: tmp_ind_cell

  integer :: i, j, k, dir, dim, ison, iskip, icell, igrid_cell, igrid_ncell, ficell, ncell, ipart, new_icell, ix, iy, iz, ii, jj, kk, indn
  integer :: ixn, iyn, izn, ixm, iym, izm
  integer :: dim0, nxny, dirm2, diro2, nx_loc
  real(kind=dp) :: scale, dx, dxcoarse
  real(kind=dp), dimension(1:3) :: x, xbound, skip_loc, xc
  integer, dimension(1:ndim) :: loc
  integer :: npos, ind_ngrid_ncell, twotondimo2
  logical :: ok

  integer :: relLvl
  ! Maximum indexes when looping on cells
#IF NDIM == 1
  ixm = 2; iym = 1; izm = 1
#ENDIF
#IF NDIM == 2
  ixm = 2; iym = 2; izm = 1
#ENDIF
#IF NDIM == 3
  ixm = 2; iym = 2; izm = 2
#ENDIF
  ! integer, save :: ncall = 0

  twotondimo2 = twotondim/2
  ! ncall = ncall + 1
  ! write(*, '(a,i3,a)') '===', ncall, "'th call of move_tracers_oct"
  nxny = nx*ny

  nx_loc=(icoarse_max-icoarse_min+1)
  dx = 0.5D0**ilevel
  dxcoarse = 2d0 * dx
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  skip_loc=(/0.0d0, 0.0d0, 0.0d0/)
  if(ndim>0) skip_loc(1) = dble(icoarse_min)
  if(ndim>1) skip_loc(2) = dble(jcoarse_min)
  if(ndim>2) skip_loc(3) = dble(kcoarse_min)
  scale = boxlen/dble(nx_loc)

  ! Get neighbour grids
  call getnborgrids(ind_grid, ind_ngrid, nvector)

  ! Precompute the indexes of the cells of each grid
  do j = 1, nvector
     do i = 1, twotondim
        iskip = ncoarse + (i - 1)*ngridmax
        ind_cell(j, i) = ind_grid(j) + iskip
     end do
  end do

  ! Get each cell's neighbour in all 6 directions
  do i = 1, twotondim
     do j = 1, nvector
        tmp_ind_cell(j) = ind_cell(j, i)
     end do
     ! print*, 'direction   ', i
     ! print*, 'grid        ', ind_grid
     ! print*, 'tmp_ind_cell', tmp_ind_cell
     call getnborfather(tmp_ind_cell, tmp_ncell, nvector, ilevel)
     do dir = 0, twondim
        do j = 1, nvector
           ind_ncell(j, i, dir) = tmp_ncell(j, dir)
        end do
     end do
  end do

  !--------------------------------------------------------------------------
  ! Move particles from neighbor cell into cell following fluxes
  !--------------------------------------------------------------------------
  ! Here, the neighbor grids are iterated through. Each unrefined one (the cell contained
  ! are at a coarser level) are treated and their particles moved into the current grid.
  do j = 1, nvector
     do dir = 1, twondim
        if (ind_ngrid(j, dir) == 0) then
           ! Take cell in grid at boundary
           if (mod(dir, 2) == 1) then
              ! Upper left (front) cell in 2D (3D)
              ison = 1
           else
              ! Lower right (back) cell n 2D (3D)
              ison = twotondim
           end if
           ! Get its neighbour in direction dir
           ncell = ind_ncell(j, ison, dir)

           ! Get the neighbour grid
           npos = (ncell-ncoarse-1)/ngridmax+1
           ind_ngrid_ncell = ncell-ncoarse-(npos-1)*ngridmax

           ! Compute the flux from ngrid to grid
           flux = 0

           ! Iterate over cells on the face in direction
           ! compute the flux + store masses
           tmp_ind_cell(1:twotondimo2) = getCells(dir)
           mass = uold(ncell, 1) ! ncell mass
           masses(0) = 0d0       ! total mass
           do i = 1, twotondimo2
              ison = tmp_ind_cell(i)
              iz = (ison-1)/4
              iy = (ison-1-4*iz)/2
              ix = (ison-1-4*iz-2*iy)
              call getFlux(dir, ix+1, iy+1, iz+1, j, fluxes, flux)
              masses(i) = uold(ncoarse + (ison-1)*ngridmax + ind_grid(j), 1)
              masses(0) = masses(0) + masses(i) ! total mass
           end do

           ! Skip to next direction if flux<0 (outflow from cell to ncell, not interested)
           if (flux < 0) then
              exit
           end if
           ! print*, '================================================'
           ! print*, ind_grid(j), dir, numbp(ind_ngrid_ncell)


           ! Iterate over particles
           ipart = headp(ind_ngrid_ncell)

           ! Compute neighbor cell center
           x(1:ndim) = cellCenter(npos, ind_ngrid_ncell, dxcoarse)

           ! For the tracer particles not already moved in this cell
           do i = 1, numbp(ind_ngrid_ncell)
              ! print*, 'ipart move_flag', ilevel
              ! print*, ipart, move_flag(ipart)
              ! print*, x(1:ndim)/dx
              ! print*, xp(ipart, 1:ndim)/dx
              if (mp(ipart) == 0d0 .and. move_flag(ipart) .and. &
                   all(x(1:ndim) == xp(ipart, 1:ndim))) then
                 call ranf(localseed, rand)
                 ! print*, 'ipart flux/mass rand'
                 ! print*, ipart, flux/mass, rand
                 if (rand < flux/mass) then
                    ! Distribute randomly the particle

                    ! print*, 'moving', ipart
                    call ranf(localseed, rand)
                    do k = 1, twotondimo2
                       if (rand < masses(k) / masses(0)) then
                          ! Put particle at center of cell
                          ! print*,'moved from to', ncell, ind_ngrid_ncell+ncoarse+(tmp_ind_cell(k)-1)*ngridmax+1, ind_grid(j), dir

                          ! print*, xp(ipart, 1:ndim)
                          xp(ipart, 1:ndim) = cellCenter(tmp_ind_cell(k), ind_grid(j), dx)
                          tp(ipart) = masses(k)
                          ! print*, xp(ipart, 1:ndim)
                          rand = 1
                       else
                          rand = rand - masses(k) / masses(0)
                       end if
                    end do

                    ! Mark particle as treated
                    move_flag(ipart) = .false.
                 end if
              end if

              ipart = nextp(ipart)
           end do
        end if
     end do
  end do

  !--------------------------------------------------------------------------
  ! Move the particles in cell following the flux
  !--------------------------------------------------------------------------
  do j = 1, nvector
     ! Loop over all particles
     ipart = headp(ind_grid(j))
     if(debug) print*, 'particle', ipart, 'level', ilevel
     do i = 1, numbp(ind_grid(j))
        if(ipart == -1) print*, 'enter i loop', move_flag(ipart), myid


        if (mp(ipart) == 0d0 .and. move_flag(ipart)) then

           ! Find cell in which the particle is
           ! see amr/refine_utils.f90:202
           ix = 1; iy = 1; iz = 1
           do dim = 1, ndim
              x(dim) = (xp(ipart, dim) / scale + skip_loc(dim) - xg(ind_grid(j), dim))/dx + 1.5d0
           end do

           ok = .true.
           ix = x(1)
           if (real(ix, dp) /= x(1)) ok = .false.
#IF NDIM > 1
           iy = x(2)
           if (real(iy, dp) /= x(2)) ok = .false.
#ENDIF
#IF NDIM > 2
           iz = x(3)
           if (real(iz, dp) /= x(3)) ok = .false.
#ENDIF

           !--------------------------------------------------------------------------
           ! Check that the particles is in the center of a cell
           !--------------------------------------------------------------------------
           ! Three cases for the location of the particle:
           ! 1. in the center of a cell
           !    -> do nothing
           ! 2. in the center of the grid (the grid just got refined)
           !    -> distribute the particle randomly into the grid (following the mass)
           ! 3. in the center of a deleted cell
           !    -> recenter the particle in the enclosing cell

           if (.not. ok) then ! cases 2 and 3
              if (all(x(1:ndim) == 1.5d0)) then ! case 2
                 ! print*, 'case 2'
                 ! compute grid mass
                 mass = 0d0
                 do ison = 1, twotondim
                    iskip = ncoarse + (ison-1)*ngridmax
                    icell = iskip + ind_grid(j)
                    mass = mass + uold(icell, 1)
                 end do

                 ! randomly pick a destination for the particle
                 call ranf(localseed, rand)
                 ! print*, 'before', xp(ipart, :) / dx

                 do ison = 1, twotondim
                    iskip = ncoarse + (ison-1)*ngridmax
                    icell = iskip + ind_grid(j)
                    if (rand < uold(icell, 1) / mass) then
                       if(ipart == -1) print*, 'projecting'
                       xp(ipart, 1:ndim) = cellCenter(ison, ind_grid(j), dx)
                       tp(ipart) = uold(icell, 1)
                       kk = (ison-1) / 4
                       jj = (ison-1-4*kk) / 2
                       ii = (ison-1-4*kk-2*jj)

                       ix = ii + 1
                       iy = jj + 1
                       iz = kk + 1
                       exit
                    else
                       rand = rand - uold(icell, 1) / mass
                    end if
                 end do
                 ! print*, 'after ', xp(ipart, :) / dx
                 ison = ison
                 icell = icell
              else                                                ! case 3
                 ! print*, 'case 3'

                 ! print*, 'before', xp(ipart, :) / dx
                 ison = 1
                 do dim = 1, ndim
                    if(ipart == -1) print*, 'recentering'
                    if (x(dim) < 1.5) then
                       xp(ipart, dim) = (xg(ind_grid(j), dim) - skip_loc(dim) - dx/2)*scale
                    else
                       xp(ipart, dim) = (xg(ind_grid(j), dim) - skip_loc(dim) + dx/2)*scale
                       ison = ison + 2**(dim-1)
                    end if
                 end do
                 ! print*, 'after ', xp(ipart, :) / dx
                 iz = (ison-1) / 4 + 1
                 iy = (ison-1-4*(iz - 1)) / 2 + 1
                 ix = (ison-1-4*(iz - 1)-2*(iy - 1)) + 1

                 ison = ison
                 iskip = ncoarse + (ison-1)*ngridmax
                 icell = iskip + ind_grid(j)

                 tp(ipart) = uold(icell, 1)
                 ! print*, 'recentering', ipart, move_flag(ipart), icell, ind_grid(j)
              end if
           else
              ison = 1 + (ix-1) + 2*(iy-1) + 4*(iz-1)

              iskip = ncoarse + (ison-1)*ngridmax
              icell = iskip + ind_grid(j)
           end if
           ! At this stage, we have icell and ison for the particle

           ! Compute the outflux (<0)
           Fout = 0
           do dir = 1, twondim
              ! Check that the neighbor cell in dir is not finer (in this case,
              ! ignore flux (already treated above))
              if (neighborRelativeLevel(icell, ind_ngrid(j,&
                   & 0:twotondim), ind_ncell(j, ison, 0:twotondim),&
                   & dir) < 1) then
                 call getFlux(dir, ix, iy, iz, j, fluxes, flux)
                 if (flux < 0) Fout = Fout + flux
              end if

           end do

           mass = uold(icell, 1) ! mass of cell

           call ranf(localseed, rand)

           ! Move the particle
           if (rand < -Fout/mass) then
              ! Pick a direction
              call ranf(localseed, rand)

              do dir = 1, twondim
                 if(ipart == -1) print*, 'enter dir loop', dir
                 call getFlux(dir, ix, iy, iz, j, fluxes, flux)
                 if (rand < flux/Fout) then ! Move particle
                    ! 3 cases:
                    ! --------
                    !
                    ! 1. the neighbour cell is refined (a.k.a. has a
                    ! son)
                    ! 2. the neighbour cell is unrefined, of same
                    ! level as current cell
                    ! (their fathers are neighbors, or are the same)
                    ! 3. the neighbour cell is unrefined, of coarser
                    ! level as current cell
                    ! (neighbour grid == 0)

                    ! Case 1 has already been treated at a finer
                    ! level, do nothing
                    ! Case 2 is easy, shift the particle in the right
                    ! direction
                    ! Case 3 we juste move the particle in the center
                    ! of the neighbour cell
                    relLvl = neighborRelativeLevel(icell, ind_ngrid(j&
                         &, 0:twondim), ind_ncell(j, ison, 0:twondim)&
                         &, dir)
                    if (relLvl == 1) then                            &
                         &                     ! case 1
                       ! Nothing to do, it should have already been
                         ! done at a finer level
                       print*, 'Too lazy for da job'
                       stop
                    else if (relLvl == 0) then                       &
                         &                     ! case 2
                       dirm2 = (dir - 1) / 2 + 1 ! 1,2->1 | 3,4->2 |
                       ! 5,6->3

                       if(ipart == -1) print*, 'cell to cell (same&
                            & level)'
                       ! print*, 'case  2'
                       if (mod(dir, 2) == 1) then ! going left
                          xp(ipart, dirm2) = xp(ipart, dirm2) - dx
                       else                       ! going right
                          xp(ipart, dirm2) = xp(ipart, dirm2) + dx
                       end if
                       tp(ipart) = uold(ind_ncell(j, ison, dir), 1)

                    else if (relLvl == -1) then                      &
                         &                     ! case 3
                       ! get the center of the neighboring cell
                       ! print*, 'case   3', dir

                       ! Get the index of neigh. in its grid
                       indn = (ind_ncell(j, ison, dir) - ncoarse - 1)&
                       &/ngridmax + 1

                       ! Get the grid number
                       ind_ngrid_ncell = ind_ncell(j, ison, dir)&
                            &-ncoarse-(indn-1)*ngridmax

                       if(ipart == -1) then
                          print*, 'cell to cell (fine to&
                            & coarse)'
                          print*, xp(ipart, 1:ndim)
                          print*, xg(ind_grid(j), :)
                          print*, xg(ind_ngrid_ncell, :)
                          write(*, '(5(a12, 2x))') 'indn', 'dir', 'ind_ngrid_ncell'
                          write(*, '(5(i12, 2x))') indn, dir, ind_ngrid_ncell
                       end if

                       ! Set particle in cell (at a coarser level)
                       x(1:ndim) = cellCenter(indn,&
                            & ind_ngrid_ncell, dxcoarse)

                       do dim = 1, ndim
                          ! Detect when moving particles of more than half the box size
                          if (x(dim) - xp(ipart, dim) > 0.5d0) then
                             print*, 'fucking patch -1 ', ipart
                             xp(ipart, dim) = x(dim) - 1d0
                          else if (x(dim) - xp(ipart, dim) < -0.5d0) then
                             print*, 'fucking patch +1 ', ipart
                             xp(ipart, dim) = x(dim) + 1d0
                          else
                             xp(ipart, dim) = x(dim)
                          end if
                       end do

                       tp(ipart) = uold(ind_ncell(j, ison, dir), 1)

                       if (ipart == -1) then
                          print*, xp(ipart, 1:ndim)
                          print*, xp(ipart, 1:ndim) / dx
                       end if
                    else
                       print*, 'should not happen'
                       stop
                    end if

                    if(ipart == -1) print*, 'exit dir loop'
                    exit
                 else ! Increase proba for moving in next direction
                    if (flux < 0) rand = rand - flux / Fout
                 end if
              end do
           end if
        end if

        if(ipart == -1) print*, 'set move_flag=.false.'
        ! Next one
        move_flag(ipart) = .false.
        ipart = nextp(ipart)
     end do
  end do

contains
  subroutine getFlux(dir, ii0, jj0, kk0, j, fluxes, flux)
    integer, intent(in)  :: dir, j
    integer, intent(in) :: ii0, jj0, kk0
    real(dp), intent(in) :: fluxes(:,:,:,:,:,:)

    real(dp), intent(out) :: flux

    integer :: sign ! this is the sign to get fluxes from cell out of it
    integer :: dim, nxny, ii, jj, kk

    ii = ii0
    jj = jj0
    kk = kk0
    ! the sign is +1 when the direction is 1, 3, 5
    ! the sign is -1 when the direction is 2, 4, 6
    sign = 1
    if (mod(dir, 2) == 0) sign = -1
    dim = dir / 3 + 1 ! so that 1,2 -> 1 & 3,4 -> 2 & 4,5 -> 3

    ! Select the location given ii, jj, kk and the direction
    if (dim == 1)      then
       if (sign == -1) ii = ii + 1
    else if (dim == 2) then
       if (sign == -1) jj = jj + 1
       kk = 1
    else
       if (sign == -1) kk = kk + 1
    end if

    flux = sign*fluxes(j, ii, jj, kk, 1, dim)
  end subroutine getFlux

  ! Compute the relative level of the neighboring cell in direction, returning
  ! +1 if the cell in dir is refined, 0 if the cell is unrefined of same level
  ! as the input cell and -1 if the cell is unrefined at level-1
  function neighborRelativeLevel(ind_cell, ind_ngrid, ind_ncell, direction)
    integer, intent(in) :: ind_cell, direction
    integer, intent(in), dimension(0:twondim) :: ind_ncell ! the neighboring cells
    integer, intent(in), dimension(0:twondim) :: ind_ngrid ! the neighboring parent grids
    integer :: neighborRelativeLevel

    integer :: pos, ind_ngrid_ncell

    ! If the neighbor cell is a grid
    if (son(ind_ncell(direction)) > 0 ) then
       neighborRelativeLevel = +1
    else
       ! If the grid of the neighbour cell is a neighbor of the cell's grid

       ! get the location of neighbor in grid
       pos = (ind_ncell(direction)-ncoarse-1)/ngridmax+1
       ind_ngrid_ncell = ind_ncell(direction)-ncoarse-(pos-1)*ngridmax

       ! If the neighbor cell is in the same/a neighbor grid
       if (ind_ngrid_ncell == ind_ngrid(direction) .or. ind_ngrid_ncell == ind_ngrid(0)) then
          neighborRelativeLevel = 0
       else
          neighborRelativeLevel = -1
       end if

    end if

  end function neighborRelativeLevel

  function getCells(direction) result (locs)
    integer, intent(in) :: direction
    integer, dimension(1:twotondim/2) :: locs

    integer, save, dimension(1:6, 1:4) :: mapping
    logical, save :: firstCall = .true.

    if (firstCall) then
       mapping(1, 1:4) = (/1, 3, 5, 7/) ! left cells
       mapping(2, 1:4) = (/2, 4, 6, 8/) ! right cells
       mapping(3, 1:4) = (/1, 2, 5, 6/) ! top cells
       mapping(4, 1:4) = (/3, 4, 7, 8/) ! bottom cells
       mapping(5, 1:4) = (/1, 2, 3, 4/) ! front cells
       mapping(6, 1:4) = (/5, 6, 7, 8/) ! back cells

       firstCall = .false.
    end if

    locs(1:twotondimo2) = mapping(direction, 1:twotondimo2)

  end function getCells

  function cellCenter(ison, ind_grid, dx) result (x)
    integer, intent(in)  :: ison, ind_grid
    real(dp), intent(in) :: dx
    real(dp), dimension(1:ndim) :: x

    real(dp), dimension(1:3) :: xc

    integer :: ixn, iyn, izn, dim
    ! Get the location of neighbor cell in its grid
    izn = (ison-1)/4
    iyn = (ison-1-4*izn)/2
    ixn = (ison-1-4*izn-2*iyn)

    ! Compute the expected location of the particles
    xc(1) = (dble(ixn)-0.5d0)*dx
    xc(2) = (dble(iyn)-0.5d0)*dx
    xc(3) = (dble(izn)-0.5d0)*dx

    do dim = 1, ndim
       x(dim) = (xg(ind_grid, dim) + xc(dim) - skip_loc(dim)) * scale
    end do

    ! print*, 'cellCenter', ison, ind_grid, dx, x(1:ndim)

  end function cellCenter

end subroutine move_tracers_oct
