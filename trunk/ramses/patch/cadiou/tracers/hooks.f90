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
    integer, intent(in), dimension(:) :: ind_cell

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
             print*, 'pre_kill_grid_hook', ipart, igrid

             if (mp(ipart) == 0d0) then
                do dim = 1, ndim
                   xp(ipart, dim) = (xg(igrid, dim) - skip_loc(dim)) * scale
                end do
             end if
             ipart = nextp(igrid)
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
                ok = .true.
                ! Check whether the particle was in the refined cell
                do dim = 1, ndim
                   x = (xp(ipart, dim) / scale - xg(fgrid, dim) + skip_loc(dim)) / dxcoarse &
                        + 0.5d0
                   ok = ok .and. (int(x) == loc(dim))

                   ! It is possible that x is not 0 or 1 if the particle has already been
                   ! moved by the routine. In this case, we don't need to move it again
                   print*, x, loc(dim)
                end do

                ! If the particle is in refined cell, spread it accordingly
                if (ok) then
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
                         do dim = 1, ndim
                            xp(ipart, dim) = (xg(igrid, dim) + xc(dim) - skip_loc(dim)) * scale
                         end do
                         rand = 1
                      else
                         rand = rand - mass(ison) / mass(0)
                      end if
                   end do
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
  use amr_commons, only   : ncoarse, father, xg, son
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
  real(dp) :: flux, Fout, mass, rand
  integer,dimension(1:nvector, 1:twotondim) :: ind_cell
  integer,dimension(1:nvector, 1:twotondim, 0:twondim) :: ind_ncell
  integer,dimension(1:nvector, 0:twondim) :: tmp_ncell
  integer, dimension(1:nvector) :: tmp_ind_cell

  integer :: i, j, k, dir, dim, ison, iskip, icell, igrid_cell, ficell, ncell, ipart, new_icell, ix, iy, iz, ii, jj, kk
  integer :: ixnc, iync, iznc
  integer :: dim0, nxny, dirm2, nx_loc
  real(kind=dp) :: scale, dx
  real(kind=dp), dimension(1:3) :: x, xbound, skip_loc, disp
  integer :: npos, ind_fathergrid_ncell
  logical :: ok


  ! Maximum indexes when looping on cells
#IF NDIM == 1
  ixnc = 2; iync = 1; iznc = 1
#ENDIF
#IF NDIM == 2
  ixnc = 2; iync = 2; iznc = 1
#ENDIF
#IF NDIM == 3
  ixnc = 2; iync = 2; iznc = 2
#ENDIF
  ! integer, save :: ncall = 0

  ! ncall = ncall + 1
  ! write(*, '(a,i3,a)') '===', ncall, "'th call of move_tracers_oct"
  nxny = nx*ny

  nx_loc=(icoarse_max-icoarse_min+1)
  dx = 0.5D0**ilevel
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
     do k = 0, twondim
        do j = 1, nvector
           ind_ncell(j, i, k) = tmp_ncell(j, k)
        end do
     end do
  end do

  ! do j = 1, nvector
  !    do i = 1, twotondim
  !       print*, 'Neighbour of ', ind_cell(j, i), ind_grid(j)
  !       print*, ind_ncell(j, i, :)
  !    end do
  ! end do

  do j = 1, nvector
     ! Loop over all particles
     ipart = headp(ind_grid(j))
     if(debug) print*, 'particle', ipart, 'level', ilevel

     do i = 1, numbp(ind_grid(j))

        if (mp(ipart) == 0d0 .and. move_flag(ipart)) then

           ! Find cell in which the particle is
           ! see amr/refine_utils.f90:202
           ix = 1; iy = 1; iz = 1
           do dim = 1, ndim
              x(dim) = (xp(ipart, dim) / scale + skip_loc(dim) - xg(ind_grid(j), dim))/dx + 1.5d0
           end do

           ! print*, 'xp', (xp(ipart, :) / scale + skip_loc(:)) / dx + 1.5d0
           ! print*, 'xg', (xg(ind_grid(j), :)) / dx
           ! print*, 'x ', x

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
           ! Three cases for the location of the particle:
           ! 1. in the center of a cell
           !    -> do nothing
           ! 2. in the center of the grid (the grid just got refined)
           !    -> distribute the particle randomly into the grid (following the mass)
           ! 3. in the center of a deleted cell
           !    -> recenter the particle in the enclosing cell

           if (.not. ok) then ! cases 2 and 3
              !TODO : project to nearest randomly
              ! print*, x
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

                 do ii = 1, ixnc
                    do jj = 1, iync
                       do kk = 1, iznc
                          ison = 1 + (ii-1) + 2*(jj-1) + 4*(kk-1)
                          iskip = ncoarse + (ison-1)*ngridmax
                          icell = iskip + ind_grid(j)
                          if (rand < uold(icell, 1) / mass) then
                             ! print*, 'projecting', ipart, ilevel, 1+ii+2*jj+4*kk

                             disp = ((/ii, jj, kk/)*2 - 3) * dx / 2d0
                             do dim = 1, ndim
                                xp(ipart, dim) = xp(ipart, dim) + disp(dim)
                             end do
                             ix = ii
                             iy = jj
                             iz = kk
                             rand = 1
                          else
                             rand = rand - uold(icell, 1) / mass
                          end if
                       end do
                    end do
                 end do
                 ! print*, 'after ', xp(ipart, :) / dx
                 ison = 1 + (ix-1) + 2*(iy-1) + 4*(iz-1)

                 iskip = ncoarse + (ison-1)*ngridmax
                 icell = iskip + ind_grid(j)
              else                                                ! case 3
                 ! print*, 'case 3'
                 ! print*, 'recentering', ipart, ilevel
                 ! print*, 'before', xp(ipart, :) / dx
                 ison = 1
                 do dim = 1, ndim
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

                 iskip = ncoarse + (ison-1)*ngridmax
                 icell = iskip + ind_grid(j)
              end if
           else
              ison = 1 + (ix-1) + 2*(iy-1) + 4*(iz-1)

              iskip = ncoarse + (ison-1)*ngridmax
              icell = iskip + ind_grid(j)
           end if
           ! print*, ison, icell, ix, iy, iz

           ! print*, 'computing flux', ipart

           ! Compute the outflux (<0)
           Fout = 0
           do dir = 1, twondim
              ! print*, ix, iy, iz, j, ilevel
              call getFlux(dir, ix, iy, iz, j, fluxes, flux)
              if (flux < 0) Fout = Fout + flux
           end do

           mass = uold(icell, 1) ! mass of cell

           call ranf(localseed, rand)

           ! Move the particle
           if (rand < -Fout/mass) then
              ! Pick a direction
              call ranf(localseed, rand)

              ! Two cases:
              ! 1. stay in the same grid
              ! 2. move to a neighbor grid
              do dir = 1, twondim
                 call getFlux(dir, ix, iy, iz, j, fluxes, flux)
                 if (rand < flux/Fout) then ! Move particle

                    ! 4 cases:
                    ! 1. the neighbour cell is refined (a.k.a. has a son)
                    ! 2. the neighbour cell is unrefined, of same level as current cell (their fathers are neighbors, or are the same)
                    ! 3. the neighbour cell is unrefined, of coarser level as current cell (neighbour grid == 0)
                    if (son(ind_ncell(j, ison, dir)) > 0) then                      ! case 1
                       ! The algorithm has already treated it, hasn't it?
                       ! print*, 'case 1'
                    else
                       if (ind_ngrid(j, dir) == 0) then                             ! case 3
                          ! get the center of the cell
                          ! print*, 'case 3'
                       else                                                         ! case 2
                          dirm2 = (dir - 1) / 2 + 1 ! 1,2->1 | 3,4->2 | 5,6->3

                          ! print*, 'case 2'
                          if (mod(dir, 2) == 1) then ! going left
                             xp(ipart, dirm2) = xp(ipart, dirm2) - dx
                          else                       ! going right
                             xp(ipart, dirm2) = xp(ipart, dirm2) + dx
                          end if
                       end if
                    end if

                    ! do dim = 1, ndim
                    !    if (dirm2 == dim) then
                    !       if (mod(dir, 2) == 1) then ! going left
                    !          xp(ipart, dim) = xp(ipart, dim) - dx
                    !       else                       ! going right
                    !          xp(ipart, dim) = xp(ipart, dim) + dx
                    !       end if
                    !       ! print*, 'moved particle', ipart, dir, dx
                    !    end if
                    ! end do
                    exit
                 else ! Increase proba for moving in next direction
                    if (flux < 0) rand = rand - flux / Fout
                 end if
              end do
           end if
        end if

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
    integer, intent(in), dimension(1:twondim) :: ind_ncell ! the neighboring cells
    integer, intent(in), dimension(1:twondim) :: ind_ngrid ! the neighboring parent grids
    integer :: neighborRelativeLevel

    integer :: pos, ind_fathergrid_ncell

    if (son(ind_ncell(direction)) > 0 ) then
       ! If the neighbor cell is a grid, then because of how ind_ncell is built
       neighborRelativeLevel = +1
    else
       ! If the grid of the neighbour cell is a neighbor of the cell's grid

       ! get the location of neighbor in grid
       pos = (ind_ncell(direction)-ncoarse-1)/ngridmax+1
       ind_fathergrid_ncell = ind_ncell(direction)-ncoarse-(pos-1)*ngridmax

       if (ind_fathergrid_ncell == ind_ngrid(direction)) then
          neighborRelativeLevel = 0
       else
          neighborRelativeLevel = -1
       end if

    end if

  end function neighborRelativeLevel


end subroutine move_tracers_oct
