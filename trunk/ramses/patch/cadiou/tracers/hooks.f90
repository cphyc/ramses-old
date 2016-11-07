! Hooks for RAMSES
module hooks
  use amr_parameters      ! nx, ny, nz, dp, ngridmax, nvector, …
  use amr_commons, only   : ncoarse, father, xg, son, myid
  use pm_commons, only    : xp, tp, idp, levelp, headp, mp, localseed, numbp, nextp, move_flag
  use random, only        : ranf
  use hydro_commons, only : uold, if1, if2, jf1, jf2, kf1, kf2, nvar

  implicit none

  real(kind=dp), dimension(1:3) :: skip_loc
  real(kind=dp) :: scale
  integer :: nx_loc

  logical :: ddebug = .false.
  integer :: dbpart = -1

contains

  subroutine initialize_skip_loc
    logical, save :: firstCall = .true.

    if (firstCall) then
       skip_loc=(/0.0d0, 0.0d0, 0.0d0/)
       if(ndim>0) skip_loc(1) = dble(icoarse_min)
       if(ndim>1) skip_loc(2) = dble(jcoarse_min)
       if(ndim>2) skip_loc(3) = dble(kcoarse_min)

       nx_loc=(icoarse_max-icoarse_min+1)
       scale = boxlen/dble(nx_loc)

       firstCall = .false.
    end if

  end subroutine initialize_skip_loc


  !---------------------------------------------------------------------------------
  ! Hook after flux computation
  !---------------------------------------------------------------------------------
  subroutine post_flux_hook(ind_grid, flux, ilevel)
    use amr_commons
    use hydro_commons

    integer, dimension(1:nvector), intent(in)                                           :: ind_grid
    real(dp),dimension(1:nvector, if1:if2, jf1:jf2, kf1:kf2, 1:nvar, 1:ndim),intent(in) :: flux
    integer, intent(in)                                                                 :: ilevel

    call initialize_skip_loc

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
    real(dp) :: dx, prevxp(1:ndim)

    call initialize_skip_loc

    dx = 0.5D0**ilevel

    if (MC_tracer) then
       ! For all particles, recenter them in the center of the grid that's being deleted
       ! (becoming a cell)
       do j = 1, nn
          igrid = son(ind_cell(j))
          ipart = headp(igrid)

          do i = 1, numbp(igrid)
             if (mp(ipart) == 0d0) then
                if (ddebug) print*,'kill_grid', ipart, igrid

                if (idp(ipart) == dbpart) print*, 'DEBUG: ', 'kill_grid'
                if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)

                prevxp(:) = xp(ipart, :)
                do dim = 1, ndim
                   xp(ipart, dim) = (xg(igrid, dim) - skip_loc(dim)) * scale
                end do
                call checkBadMove(xp(ipart, :), prevxp)
                if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)
                ! print*, xp(ipart, :) * 2**8
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

    call initialize_skip_loc

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
    call initialize_skip_loc

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
    real(dp) :: dx, dxcoarse

    real(dp) :: x(1:ndim), rand, mass(0:twotondim), prevxp(1:ndim)
    integer :: j, i, dim, fgrid, igrid, ipart, part_dir, ison, iskip, icell
    integer :: loc(1:3)
    logical :: ok

    call initialize_skip_loc

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
          do ison = 1, twotondim
             iskip = ncoarse + (ison-1)*ngridmax
             icell = iskip + igrid
             mass(ison) = uold(icell, 1)
          end do

          mass(0) = sum(mass(1:twotondim))

          do i = 1, numbp(fgrid)
             if (mp(ipart) == 0d0) then

                ! Check whether the particle was in the refined cell
                x(1:ndim) = cellCenter(ind, fgrid, dxcoarse)

                ok = all(xp(ipart, 1:ndim) == x(1:ndim))
                ! If the particle is in refined cell, spread it accordingly
                if (ok) then

                   ! Pick a random direction
                   call ranf(localseed, rand)

                   do ison = 1, twotondim
                      if (rand < mass(ison) / mass(0)) then
                         if (ddebug) print*,'make_grid', idp(ipart), igrid
                         ! Move particle to center of new cells
                         if (idp(ipart) == dbpart) print*, 'DEBUG: ', 'make_grid'
                         if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)
                         prevxp(:) = xp(ipart, :)
                         xp(ipart, :) = cellCenter(ison, igrid, dx)
                         call checkBadMove(xp(ipart, :), prevxp(:))
                         if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)
                         exit
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

  function cellCenter(ison, ind_grid, dx) result (x)
    implicit none

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

  subroutine move_tracers_oct(ind_grid, fluxes, ilevel)

    integer,dimension(1:nvector), intent(in)   :: ind_grid
    real(dp),dimension(1:nvector, if1:if2, jf1:jf2, kf1:kf2, 1:nvar, 1:ndim),intent(in)::fluxes
    integer, intent(in) :: ilevel

    !-------------------------------------------------------------------
    ! Implement a Monte Carlo sampling of the gas properties following
    ! the fluxes of mass from oct to oct.
    !-------------------------------------------------------------------
    integer,dimension(1:nvector,0:twondim) :: ind_ngrid
    real(dp) :: flux, Fout, mass, rand, masses(0:twotondim/2), x_tmp
    integer,dimension(1:nvector, 1:twotondim) :: ind_cell
    integer,dimension(1:nvector, 1:twotondim, 0:twondim) :: ind_ncell
    integer,dimension(1:nvector, 0:twondim) :: tmp_ncell
    integer, dimension(1:nvector) :: tmp_ind_cell

    integer :: i, j, k, dir, dim, ison, iskip, icell, igrid_cell, igrid_ncell, ficell, ncell, ipart, new_icell, ix, iy, iz, ii, jj, kk, indn
    integer :: ixn, iyn, izn, ixm, iym, izm
    integer :: dim0, nxny, dirm2, diro2, nx_loc
    real(kind=dp) :: dx, dxcoarse
    real(kind=dp), dimension(1:3) :: x, xbound, xc
    real(dp), dimension(1:ndim) :: prevxp
    integer, dimension(1:ndim) :: loc
    integer :: npos, ind_ngrid_ncell, twotondimo2
    logical :: ok

    integer, dimension(1:twondim) :: relLvl

    call initialize_skip_loc

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

    dx = 0.5D0**ilevel
    dxcoarse = 2d0 * dx
    xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)

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
          if (ind_ngrid(j, dir) == 0) then ! No neighbor grid -> neighbor is an unrefined cell
             ! Take cell in grid at boundary
             if (mod(dir, 2) == 1) then
                ! Upper left (front) cell in 2D (3D)
                ison = 1
             else
                ! Lower right (back) cell in 2D (3D)
                ison = twotondim
             end if
             ! Get its neighbour in direction dir
             ncell = ind_ncell(j, ison, dir)

             ! Get the neighbour grid
             npos = (ncell-ncoarse-1)/ngridmax+1
             ind_ngrid_ncell = ncell-ncoarse-(npos-1)*ngridmax

             ! Compute the flux from ngrid to grid
             Fout = 0

             ! Iterate over cells on the face in direction
             ! compute the flux + store masses
             tmp_ind_cell(1:twotondimo2) = getCells(dir)
             mass = uold(ncell, 1) * twotondim ! ncell mass (in units of current level)
             masses(0) = 0d0                   ! total mass
             do k = 1, twotondimo2
                ison = tmp_ind_cell(k)
                iz = (ison-1)/4
                iy = (ison-1-4*iz)/2
                ix = (ison-1-4*iz-2*iy)
                ! Compute the inflow in the cell
                call getFlux(dir, ix+1, iy+1, iz+1, j, fluxes, flux)
                Fout = Fout + flux
                if (flux > 0) then
                   masses(k) = uold(ncoarse + (ison-1)*ngridmax + ind_grid(j), 1)
                else
                   ! Store a 0 value so that we don't move particles into this one (flux<0)
                   masses(k) = 0
                end if

                masses(0) = masses(0) + masses(k)
             end do

             ! Skip to next direction if total flux < 0
             if (Fout < 0) then
                exit
             end if

             ! Iterate over particles
             ipart = headp(ind_ngrid_ncell)

             ! Compute neighbor cell center
             x(1:ndim) = cellCenter(npos, ind_ngrid_ncell, dxcoarse)

             ! For the tracer particles in the neighbor cell
             do i = 1, numbp(ind_ngrid_ncell)
                if (mp(ipart) == 0d0 .and. move_flag(ipart) .and. &
                     all(x(1:ndim) == xp(ipart, 1:ndim))) then

                   ! Decide whether or not to move it depending on the total flux
                   call ranf(localseed, rand)

                   if (idp(ipart) == dbpart) print*, 'DEBUG: ', 'move from coarser?', dir, (dir-1)/2+1
                   if (idp(ipart) == dbpart) print*, 'DEBUG: ', '                  ', ix, iy
                   if (idp(ipart) == dbpart) then
                      block
                        integer :: ii, jj
                        do ii = 1, 3
                           print*, 'DEBUG: F ', fluxes(j, ii, 1:3, 1, 1, dir/3+1)
                        end do
                        print*, 'DEBUG: mass ', mass
                        print*, 'DEBUG: masses ', masses(1:)
                      end block
                   end if

                   if (rand < Fout/mass) then
                      ! Distribute randomly the particle in cells
                      call ranf(localseed, rand)
                      do k = 1, twotondimo2
                         if (rand < masses(k) / masses(0)) then
                            ! Put particle at center of cell
                            if (idp(ipart) == dbpart) print*, 'DEBUG: ', 'move from coarser!'
                            if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)

                            prevxp(1:ndim) = xp(ipart, 1:ndim)
                            xp(ipart, 1:ndim) = cellCenter(tmp_ind_cell(k), ind_grid(j), dx)
                            call checkBadMove(xp(ipart, :), prevxp(:))

                            tp(ipart) = masses(k)
                            exit
                         else
                            rand = rand - masses(k) / masses(0)
                         end if
                      end do

                      ! Take note that the particle has been moved
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
       do i = 1, numbp(ind_grid(j))

          ! Store old position
          prevxp(1:ndim) = xp(ipart, 1:ndim)

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

             if (ok) then
                !~~~~~~~~~~~~~~~~~ Don't move particles ~~~~~~~~~~~~~~~~~~~~!
                ison = 1 + (ix-1) + 2*(iy-1) + 4*(iz-1)

                iskip = ncoarse + (ison-1)*ngridmax
                icell = iskip + ind_grid(j)

                if (ison < 1)  then
                   print*, 'ouch, icell < 0', ison, idp(ipart), ix, iy, iz
                   print*, xp(ipart, :)
                   print*, x(1:ndim)
                end if
             else
                if (all(x(1:ndim) == 1.5d0)) then
                   !~~~~~~~~~~~~~~~~~ Project particles ~~~~~~~~~~~~~~~~~~~~!

                   ! compute grid mass
                   mass = 0d0
                   do ison = 1, twotondim
                      iskip = ncoarse + (ison-1)*ngridmax
                      icell = iskip + ind_grid(j)
                      mass = mass + uold(icell, 1)
                   end do

                   ! randomly pick a destination for the particle
                   call ranf(localseed, rand)

                   do ison = 1, twotondim
                      iskip = ncoarse + (ison-1)*ngridmax
                      icell = iskip + ind_grid(j)
                      if (rand < uold(icell, 1) / mass) then
                         if (ddebug) write(*, '(a12, 1x, i8, 1x, i5, 1x, 3(f10.8, 1x))') &
                              'projecting', idp(ipart), ind_grid(j), xp(ipart, 1:ndim)

                         if (idp(ipart) == dbpart) print*, 'DEBUG: ', 'projecting', ind_grid(j)
                         if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)

                         xp(ipart, 1:ndim) = cellCenter(ison, ind_grid(j), dx)
                         tp(ipart) = uold(icell, 1)

                         kk = (ison-1) / 4
                         jj = (ison-1-4*kk) / 2
                         ii = (ison-1-4*kk-2*jj)

                         ix = ii + 1
                         iy = jj + 1
                         iz = kk + 1
                         if (ddebug) write(*, '(a12, 1x, i8, 1x, i5, 1x, 3(f10.8, 1x))') &
                              'projected', idp(ipart), ind_grid(j), xp(ipart, 1:ndim)
                         exit
                      else
                         rand = rand - uold(icell, 1) / mass
                      end if
                   end do
                   ison = ison
                   icell = icell
                else
                   !~~~~~~~~~~~~~~~~~ Center particles ~~~~~~~~~~~~~~~~~~~~!
                   ison = 1

                   if (ddebug) write(*, '(a12, 1x, i8, 1x, i5, 1x, 2(f10.8, 1x), i4)') &
                        'recentering', idp(ipart), ind_grid(j), xp(ipart, 1:ndim), ilevel
                   if (idp(ipart) == dbpart) print*, 'DEBUG: ', 'recentering', ind_grid(j)
                   if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)

                   ! DEBGU
                   if (ddebug) then
                      block
                        integer :: fgrid, fpos, fcell

                        fcell = father(ind_grid(j))
                        fpos = (fcell-ncoarse-1)/ngridmax+1
                        fgrid = fcell-ncoarse-(fpos-1)*ngridmax
                        print*, '|        ', xp(ipart, 1:ndim) / dx
                        print*, '|--      ', (xg(fgrid, 1:ndim)-skip_loc(1:ndim))*scale / dx
                        print*, '\--      ', (xg(ind_grid(j), 1:ndim)-skip_loc(1:ndim))*scale / dx
                        if (son(father(ind_grid(j))) /= ind_grid(j)) then
                           print*, 'sadfj;aklsdfjas;kljdfa'
                           stop
                        end if

                      end block
                   end if

                   do dim = 1, ndim
                      if (x(dim) < 1.5) then
                         xp(ipart, dim) = (xg(ind_grid(j), dim) - skip_loc(dim) - dx/2)*scale
                      else
                         xp(ipart, dim) = (xg(ind_grid(j), dim) - skip_loc(dim) + dx/2)*scale
                         ison = ison + 2**(dim-1)
                      end if
                   end do
                   if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)
                   iz = (ison-1) / 4 + 1
                   iy = (ison-1-4*(iz - 1)) / 2 + 1
                   ix = (ison-1-4*(iz - 1)-2*(iy - 1)) + 1

                   ison = ison
                   iskip = ncoarse + (ison-1)*ngridmax
                   icell = iskip + ind_grid(j)

                   tp(ipart) = uold(icell, 1)

                end if
             end if
             ! At this stage, we have icell and ison for the particle

             !~~~~~~~~~~~~~~~~~ Compute flux ~~~~~~~~~~~~~~~~~~~~!
             ! Compute the outflux (<0)
             Fout = 0
             do dir = 1, twondim
                ! Only use flux to same or coarser cells
                relLvl(dir) = neighborRelativeLevel(icell, ind_ngrid(j, 0:twotondim), &
                     ind_ncell(j, ison, 0:twotondim), &
                     & dir)
                if (relLvl(dir) < 1) then
                   call getFlux(dir, ix, iy, iz, j, fluxes, flux)
                   if (flux < 0) Fout = Fout + flux
                end if

             end do

             mass = uold(icell, 1) ! mass of cell

             call ranf(localseed, rand)

             !~~~~~~~~~~~~~~~~~ Move along flux ~~~~~~~~~~~~~~~~~~~~!
             if (rand < -Fout/mass) then
                ! Pick a direction
                call ranf(localseed, rand)

                do dir = 1, twondim
                   if (relLvl(dir) == 1) cycle
                   call getFlux(dir, ix, iy, iz, j, fluxes, flux)
                   if (rand < flux/Fout) then ! Move particle
                      ! 3 cases:
                      ! --------
                      ! 1. the neighbour cell is refined
                      ! 2. the neighbour cell is unrefined, as coarse as cell
                      ! 3. the neighbour cell is unrefined, coarser than cell

                      if (relLvl(dir) == 1) then                            &
                           &                     ! case 1
                           ! Nothing to do, it should have already been
                           ! done at a finer level
                           print*, 'should not happen'
                         stop
                      else if (relLvl(dir) == 0) then                       &
                           &                     ! case 2
                           dirm2 = (dir - 1) / 2 + 1 ! 1,2->1 | 3,4->2 |
                         ! 5,6->3

                         if (idp(ipart) == dbpart) print*, 'DEBUG: ', 'move along flux 2', ind_grid(j)
                         if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)
                         if (mod(dir, 2) == 1) then ! going left
                            xp(ipart, dirm2) = xp(ipart, dirm2) - dx
                         else                       ! going right
                            xp(ipart, dirm2) = xp(ipart, dirm2) + dx
                         end if
                         if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)
                         tp(ipart) = uold(ind_ncell(j, ison, dir), 1)

                         exit
                      else if (relLvl(dir) == -1) then ! case 3
                         ! get the center of the neighboring cell

                         ! Get the index of neigh. in its grid
                         indn = (ind_ncell(j, ison, dir) - ncoarse - 1)&
                              &/ngridmax + 1

                         ! Get the grid number
                         ind_ngrid_ncell = ind_ncell(j, ison, dir)&
                              &-ncoarse-(indn-1)*ngridmax

                         ! Set particle in cell (at a coarser level)
                         x(1:ndim) = cellCenter(indn,&
                              & ind_ngrid_ncell, dxcoarse)

                         if (idp(ipart) == dbpart) print*, 'DEBUG: ', 'move along flux 3'
                         if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)

                         do dim = 1, ndim
                            ! Detect when moving particles of more than half the box size
                            if (x(dim) - xp(ipart, dim) > 0.5d0) then
                               xp(ipart, dim) = x(dim) - 1d0
                            else if (x(dim) - xp(ipart, dim) < -0.5d0) then
                               xp(ipart, dim) = x(dim) + 1d0
                            else
                               xp(ipart, dim) = x(dim)
                            end if
                         end do
                         if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)

                         tp(ipart) = uold(ind_ncell(j, ison, dir), 1)
                         exit
                      end if

                      ! Important: do not exit here, else particles

                   else ! Increase proba for moving in next direction
                      if (flux < 0) rand = rand - flux / Fout
                   end if
                end do
             end if
          end if

          ! Check bad moves
          call checkBadMove(xp(ipart, :), prevxp(:))
          if (idp(ipart) == dbpart) print*, 'DEBUG: ', xp(ipart, :)

          ! Mark particle as moved
          move_flag(ipart) = .false.

          ipart = nextp(ipart)
       end do
    end do

  contains
    ! This function prevents the particles from jumping to the other side of the box
    ! by limiting the maximum distance traveled to less than 0.5

    ! Compute the influx in direction at a given location
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
      ! The sign is given so that a flux > 0 is an inflow into the cell
      ! This is how the flux is stored:
      !          3
      !
      !      +---|---+
      !      |   v   |
      !  1  ->       ->   2
      !      |       |
      !      +---|---+
      !          v
      !
      !          4
      !
      ! the sign is +1 when the direction is 1, 3, 5
      ! the sign is -1 when the direction is 2, 4, 6
      sign = 1
      if (mod(dir, 2) == 0) sign = -1
      dim = (dir-1) / 2 + 1 ! so that 1,2 -> 1 & 3,4 -> 2 & 4,5 -> 3

      ! Select the location given ii, jj, kk and the direction
      if (dim == 1)      then
         if (sign == -1) ii = ii + 1
         jj = 1; kk = 1
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

    ! Get the cells in direction (for example, direction=1 for cells on the left of grid)
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

  end subroutine move_tracers_oct

  subroutine checkBadMove(xp, prevxp)
    real(dp), dimension(1:ndim), intent(inout) :: xp
    real(dp), dimension(1:ndim), intent(in) :: prevxp

    integer :: dim

    do dim = 1, ndim
       if (xp(dim) - prevxp(dim) > 0.5d0) then
          ! print*,'correcting -', xp(dim)
          xp(dim) = xp(dim) - 1d0
       else if (xp(dim) - prevxp(dim) < -0.5d0) then
          ! print*,'correcting +', xp(dim)
          xp(dim) = xp(dim) + 1d0
       end if
    end do

  end subroutine checkBadMove
end module hooks
