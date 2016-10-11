subroutine move_fine(ilevel, flux)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in)  :: ilevel
  real(dp), intent(in) :: flux(:, :, :, :, :, :)
  !----------------------------------------------------------------------
  ! Update particle position and time-centred velocity at level ilevel.
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel - 1) force.
  !----------------------------------------------------------------------
  integer :: igrid, jgrid, ipart, jpart, next_part, ig, ip, npart1, info
  integer, dimension(1:nvector), save :: ind_grid, ind_part, ind_grid_part

  if(numbtot(1, ilevel) == 0) return
  if(verbose) write(*, 111)ilevel

  ! Update particles position and velocity
  ig = 0
  ip = 0
  ! Loop over grids
  igrid = headl(myid, ilevel)
  do jgrid = 1, numbl(myid, ilevel)
     npart1 = numbp(igrid)  ! Number of particles in the grid
     if (npart1 > 0) then
        ig = ig + 1
        ind_grid(ig) = igrid
        ipart = headp(igrid)
        ! Loop over particles
        do jpart = 1, npart1
           ! Save next particle  < ---- Very important !!!
           next_part = nextp(ipart)
           if (ig == 0) then
              ig = 1
              ind_grid(ig) = igrid
           end if
           ip = ip + 1
           ind_part(ip) = ipart
           ind_grid_part(ip) = ig
           if (ip == nvector) then
              call move1(ind_grid, ind_part, ind_grid_part, ig, ip, ilevel)
              ip = 0
              ig = 0
           end if
           ipart = next_part  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid = next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if (ip>0) call move1(ind_grid, ind_part, ind_grid_part, ig, ip, ilevel)

111 format('   Entering move_fine for level ', I2)

end subroutine move_fine
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine move_fine_static(ilevel, flux)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in)  :: ilevel
  real(dp), intent(in) :: flux(:, :, :, :, :, :)
  !----------------------------------------------------------------------
  ! Update particle position and time - centred velocity at level ilevel.
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel - 1) force.
  !----------------------------------------------------------------------
  integer :: igrid, jgrid, ipart, jpart, next_part, ig, ip, npart1, npart2, info
  integer, dimension(1:nvector), save :: ind_grid, ind_part, ind_grid_part

  if (numbtot(1, ilevel) == 0) return
  if (verbose) write(*, 111) ilevel

  ! Update particles position and velocity
  ig = 0
  ip = 0
  ! Loop over grids
  igrid = headl(myid, ilevel)
  do jgrid = 1, numbl(myid, ilevel)
     npart1 = numbp(igrid)  ! Number of particles in the grid
     npart2 = 0

     ! Count particles
     if (npart1>0) then
        ipart = headp(igrid)
        ! Loop over particles
        do jpart = 1, npart1
           ! Save next particle   <--- Very important !!!
           next_part = nextp(ipart)
           if(star) then
              if((.not.static_dm.and.tp(ipart).eq.0).or.(.not.static_stars.and.tp(ipart).ne.0)) then
                 npart2 = npart2 + 1
              endif
           else
              if(.not.static_dm) then
                 npart2 = npart2 + 1
              endif
           endif
           ipart = next_part  ! Go to next particle
        end do
     endif

     ! Gather star particles
     if (npart2 > 0) then
        ig = ig + 1
        ind_grid(ig) = igrid
        ipart = headp(igrid)
        ! Loop over particles
        do jpart = 1, npart1
           ! Save next particle   <--- Very important !!!
           next_part = nextp(ipart)
           ! Select particles
           if (star) then
              if ((.not. static_dm .and. tp(ipart) == 0) .or. (.not. static_stars .and. tp(ipart) /= 0)) then
                 if (ig == 0) then
                    ig = 1
                    ind_grid(ig) = igrid
                 end if
                 ip = ip + 1
                 ind_part(ip) = ipart
                 ind_grid_part(ip) = ig
              endif
           else
              if(.not. static_dm) then
                 if (ig == 0) then
                    ig = 1
                    ind_grid(ig) = igrid
                 end if
                 ip = ip + 1
                 ind_part(ip) = ipart
                 ind_grid_part(ip) = ig
              endif
           endif
           if (ip == nvector) then
              call move1(ind_grid, ind_part, ind_grid_part, ig, ip, ilevel, flux)
              ip = 0
              ig = 0
           end if
           ipart = next_part  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid = next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if (ip > 0) call move1(ind_grid, ind_part, ind_grid_part, ig, ip, ilevel, flux)

111 format('   Entering move_fine for level ', I2)

end subroutine move_fine_static
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine move1(ind_grid, ind_part, ind_grid_part, ng, np, ilevel, flux)
  use amr_commons
  use pm_commons
  use poisson_commons
  use amr_parameters
  use hydro_commons, only: uold, smallr
  use random

  implicit none
  integer :: ng, np, ilevel
  integer, intent(in), dimension(1:nvector) :: ind_grid
  integer, intent(in), dimension(1:nvector) :: ind_grid_part, ind_part
  real(dp), intent(in)                      :: flux(:, :, :, :, :, :)
  !------------------------------------------------------------
  ! This routine computes the force on each particle by
  ! inverse CIC and computes new positions for all particles.
  ! If particle sits entirely in fine level, then CIC is performed
  ! at level ilevel. Otherwise, it is performed at level ilevel - 1.
  ! This routine is called by move_fine.
  !------------------------------------------------------------
  logical :: error
  integer :: i, j, ind, idim, idim2, nx_loc, isink
  real(dp) :: dx, length, dx_loc, scale, vol_loc, r2
  ! Grid - based arrays
  integer , dimension(1:nvector), save :: father_cell
  real(dp), dimension(1:nvector, 1:ndim), save :: x0
  integer , dimension(1:nvector, 1:threetondim), save :: nbors_father_cells
  integer , dimension(1:nvector, 1:twotondim), save :: nbors_father_grids
  ! Particle - based arrays
  logical , dimension(1:nvector), save :: ok
  real(dp), dimension(1:nvector, 1:ndim), save :: x, ff, new_xp, new_vp, dd, dg
  integer , dimension(1:nvector, 1:ndim), save :: ig, id, igg, igd, icg, icd
  real(dp), dimension(1:nvector, 1:twotondim), save :: vol
  integer , dimension(1:nvector, 1:twotondim), save :: igrid, icell, indp, kg
  real(dp), dimension(1:3) :: skip_loc

  ! Loop indexes
  integer :: ii1, ii2, jj1, jj2, kk1, kk2, ii, jj, kk
  real(dp), dimension(1:nvector) :: outflux, mass_of_cell
  real(dp) :: proba, rand
  ! Mesh spacing in that level
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max - icoarse_min + 1)
  skip_loc = (/0.0d0, 0.0d0, 0.0d0/)
  if (ndim>0) skip_loc(1) = dble(icoarse_min)
  if (ndim>1) skip_loc(2) = dble(jcoarse_min)
  if (ndim>2) skip_loc(3) = dble(kcoarse_min)
  scale = boxlen/dble(nx_loc)
  dx_loc = dx*scale
  vol_loc = dx_loc**3

  ! Lower left corner of 3x3x3 grid - cube
  do idim = 1, ndim
     do i = 1, ng
        x0(i, idim) = xg(ind_grid(i), idim) - 3.0D0*dx
     end do
  end do

  ! Gather neighboring father cells (should be present anytime !)
  do i = 1, ng
     father_cell(i) = father(ind_grid(i))
  end do
  call get3cubefather(father_cell, nbors_father_cells, nbors_father_grids, &
       & ng, ilevel)

  ! Rescale particle position at level ilevel
  do idim = 1, ndim
     do j = 1, np
        x(j, idim) = xp(ind_part(j), idim)/scale + skip_loc(idim)
     end do
  end do
  do idim = 1, ndim
     do j = 1, np
        x(j, idim) = x(j, idim) - x0(ind_grid_part(j), idim)
     end do
  end do
  do idim = 1, ndim
     do j = 1, np
        x(j, idim) = x(j, idim)/dx
     end do
  end do

  ! Check for illegal moves
  error = .false.
  do idim = 1, ndim
     do j = 1, np
        if(x(j, idim)<0.5D0 .or. x(j, idim)>5.5D0)error = .true.
     end do
  end do
  if (error) then
     write(*, *)'problem in move'
     do idim = 1, ndim
        do j = 1, np
           if(x(j, idim)<0.5D0.or.x(j, idim)>5.5D0)then
              write(*, *)x(j, 1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim = 1, ndim
     do j = 1, np
        dd(j, idim) = x(j, idim) + 0.5D0
        id(j, idim) = dd(j, idim)
        dd(j, idim) = dd(j, idim) - id(j, idim)
        dg(j, idim) = 1.0D0 - dd(j, idim)
        ig(j, idim) = id(j, idim) - 1
     end do
  end do

   ! Compute parent grids
  do idim = 1, ndim
     do j = 1, np
        igg(j, idim) = ig(j, idim)/2
        igd(j, idim) = id(j, idim)/2
     end do
  end do
#if NDIM ==1
  do j = 1, np
     kg(j, 1) = 1 + igg(j, 1)
     kg(j, 2) = 1 + igd(j, 1)
  end do
#endif
#if NDIM ==2
  do j = 1, np
     kg(j, 1) = 1 + igg(j, 1) + 3*igg(j, 2)
     kg(j, 2) = 1 + igd(j, 1) + 3*igg(j, 2)
     kg(j, 3) = 1 + igg(j, 1) + 3*igd(j, 2)
     kg(j, 4) = 1 + igd(j, 1) + 3*igd(j, 2)
  end do
#endif
#if NDIM ==3
  do j = 1, np
     kg(j, 1) = 1 + igg(j, 1) + 3*igg(j, 2) + 9*igg(j, 3)
     kg(j, 2) = 1 + igd(j, 1) + 3*igg(j, 2) + 9*igg(j, 3)
     kg(j, 3) = 1 + igg(j, 1) + 3*igd(j, 2) + 9*igg(j, 3)
     kg(j, 4) = 1 + igd(j, 1) + 3*igd(j, 2) + 9*igg(j, 3)
     kg(j, 5) = 1 + igg(j, 1) + 3*igg(j, 2) + 9*igd(j, 3)
     kg(j, 6) = 1 + igd(j, 1) + 3*igg(j, 2) + 9*igd(j, 3)
     kg(j, 7) = 1 + igg(j, 1) + 3*igd(j, 2) + 9*igd(j, 3)
     kg(j, 8) = 1 + igd(j, 1) + 3*igd(j, 2) + 9*igd(j, 3)
  end do
#endif
  do ind = 1, twotondim
     do j = 1, np
        igrid(j, ind) = son(nbors_father_cells(ind_grid_part(j), kg(j, ind)))
     end do
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np) = .true.
  do ind = 1, twotondim
     do j = 1, np
        ok(j) = ok(j).and.igrid(j, ind)>0
     end do
  end do

  ! If not, rescale position at level ilevel - 1
  do idim = 1, ndim
     do j = 1, np
        if(.not.ok(j))then
           x(j, idim) = x(j, idim)/2.0D0
        end if
     end do
  end do
  ! If not, redo CIC at level ilevel - 1
  do idim = 1, ndim
     do j = 1, np
        if(.not.ok(j))then
           dd(j, idim) = x(j, idim) + 0.5D0
           id(j, idim) = dd(j, idim)
           dd(j, idim) = dd(j, idim) - id(j, idim)
           dg(j, idim) = 1.0D0 - dd(j, idim)
           ig(j, idim) = id(j, idim) - 1
        end if
     end do
  end do

 ! Compute parent cell position
  do idim = 1, ndim
     do j = 1, np
        if(ok(j))then
           icg(j, idim) = ig(j, idim) - 2*igg(j, idim)
           icd(j, idim) = id(j, idim) - 2*igd(j, idim)
        else
           icg(j, idim) = ig(j, idim)
           icd(j, idim) = id(j, idim)
        end if
     end do
  end do
#if NDIM == 1
  do j = 1, np
     icell(j, 1) = 1 + icg(j, 1)
     icell(j, 2) = 1 + icd(j, 1)
  end do
#endif
#if NDIM == 2
  do j = 1, np
     if(ok(j))then
        icell(j, 1) = 1 + icg(j, 1) + 2*icg(j, 2)
        icell(j, 2) = 1 + icd(j, 1) + 2*icg(j, 2)
        icell(j, 3) = 1 + icg(j, 1) + 2*icd(j, 2)
        icell(j, 4) = 1 + icd(j, 1) + 2*icd(j, 2)
     else
        icell(j, 1) = 1 + icg(j, 1) + 3*icg(j, 2)
        icell(j, 2) = 1 + icd(j, 1) + 3*icg(j, 2)
        icell(j, 3) = 1 + icg(j, 1) + 3*icd(j, 2)
        icell(j, 4) = 1 + icd(j, 1) + 3*icd(j, 2)
     end if
  end do
#endif
#if NDIM == 3
  do j = 1, np
     if(ok(j))then
        icell(j, 1) = 1 + icg(j, 1) + 2*icg(j, 2) + 4*icg(j, 3)
        icell(j, 2) = 1 + icd(j, 1) + 2*icg(j, 2) + 4*icg(j, 3)
        icell(j, 3) = 1 + icg(j, 1) + 2*icd(j, 2) + 4*icg(j, 3)
        icell(j, 4) = 1 + icd(j, 1) + 2*icd(j, 2) + 4*icg(j, 3)
        icell(j, 5) = 1 + icg(j, 1) + 2*icg(j, 2) + 4*icd(j, 3)
        icell(j, 6) = 1 + icd(j, 1) + 2*icg(j, 2) + 4*icd(j, 3)
        icell(j, 7) = 1 + icg(j, 1) + 2*icd(j, 2) + 4*icd(j, 3)
        icell(j, 8) = 1 + icd(j, 1) + 2*icd(j, 2) + 4*icd(j, 3)
     else
        icell(j, 1) = 1 + icg(j, 1) + 3*icg(j, 2) + 9*icg(j, 3)
        icell(j, 2) = 1 + icd(j, 1) + 3*icg(j, 2) + 9*icg(j, 3)
        icell(j, 3) = 1 + icg(j, 1) + 3*icd(j, 2) + 9*icg(j, 3)
        icell(j, 4) = 1 + icd(j, 1) + 3*icd(j, 2) + 9*icg(j, 3)
        icell(j, 5) = 1 + icg(j, 1) + 3*icg(j, 2) + 9*icd(j, 3)
        icell(j, 6) = 1 + icd(j, 1) + 3*icg(j, 2) + 9*icd(j, 3)
        icell(j, 7) = 1 + icg(j, 1) + 3*icd(j, 2) + 9*icd(j, 3)
        icell(j, 8) = 1 + icd(j, 1) + 3*icd(j, 2) + 9*icd(j, 3)
     end if
  end do
#endif

  ! Compute parent cell adresses
  do ind = 1, twotondim
     do j = 1, np
        if(ok(j))then
           indp(j, ind) = ncoarse + (icell(j, ind) - 1)*ngridmax + igrid(j, ind)
        else
           indp(j, ind) = nbors_father_cells(ind_grid_part(j), icell(j, ind))
        end if
     end do
  end do

  ! Compute cloud volumes
#if NDIM == 1
  do j = 1, np
     vol(j, 1) = dg(j, 1)
     vol(j, 2) = dd(j, 1)
  end do
#endif
#if NDIM == 2
  do j = 1, np
     vol(j, 1) = dg(j, 1)*dg(j, 2)
     vol(j, 2) = dd(j, 1)*dg(j, 2)
     vol(j, 3) = dg(j, 1)*dd(j, 2)
     vol(j, 4) = dd(j, 1)*dd(j, 2)
  end do
#endif
#if NDIM == 3
  do j = 1, np
     vol(j, 1) = dg(j, 1)*dg(j, 2)*dg(j, 3)
     vol(j, 2) = dd(j, 1)*dg(j, 2)*dg(j, 3)
     vol(j, 3) = dg(j, 1)*dd(j, 2)*dg(j, 3)
     vol(j, 4) = dd(j, 1)*dd(j, 2)*dg(j, 3)
     vol(j, 5) = dg(j, 1)*dg(j, 2)*dd(j, 3)
     vol(j, 6) = dd(j, 1)*dg(j, 2)*dd(j, 3)
     vol(j, 7) = dg(j, 1)*dd(j, 2)*dd(j, 3)
     vol(j, 8) = dd(j, 1)*dd(j, 2)*dd(j, 3)
  end do
#endif

  ! Gather 3-force
  ff(1:np, 1:ndim) = 0.0D0
  if (tracer .and. hydro) then
     do ind = 1, twotondim
        do idim = 1, ndim
           do j = 1, np
              ff(j, idim) = ff(j, idim) + uold(indp(j, ind), idim + 1) / max(uold(indp(j, ind), 1), smallr)*vol(j, ind)
           end do
        end do
     end do
  endif
  if (poisson) then
     do ind = 1, twotondim
        do idim = 1, ndim
           do j = 1, np
              ff(j, idim) = ff(j, idim) + f(indp(j, ind), idim)*vol(j, ind)
           end do
        end do
#ifdef OUTPUT_PARTICLE_POTENTIAL
        do j = 1, np
           ptcl_phi(ind_part(j)) = phi(indp(j, ind))
        end do
#endif
     end do
  endif

  ! Update velocity
  do idim = 1, ndim
     if (static .or. tracer) then
        do j = 1, np
           new_vp(j, idim) = ff(j, idim)
        end do
     else
        do j = 1, np
           new_vp(j, idim) = vp(ind_part(j), idim) + ff(j, idim)*0.5D0*dtnew(ilevel)
        end do
     endif
  end do

  ! For sink cloud particle only
  if(sink) then
     ! Overwrite cloud particle velocity with sink velocity
     do idim = 1, ndim
        do j = 1, np
           isink = -idp(ind_part(j))
           if(isink>0)then
              new_vp(j, idim) = vsnew(isink, idim, ilevel)
           end if
        end do
     end do
  end if

  ! Store velocity
  do idim = 1, ndim
     do j = 1, np
        vp(ind_part(j), idim) = new_vp(j, idim)
     end do
  end do

  ! Update position for MC tracers
  mass_of_cell = 0d0
  if (MC_tracer) then
     do j = 1, np
        if (mp(ind_part(j)) > 0d0) then
           cycle
        end if
        ! Move tracer part following the flux
        ii1=0; jj1=0; kk1=0 ! TODO define them depending on whether we're on the right side of the grid or left side
        ii2=0; jj2=0; kk2=0
        if (ndim > 0) ii2 = 2
        if (ndim > 1) jj2 = 2
        if (ndim > 2) kk2 = 2

        outflux(j) = 0d0
        ! Compute the total outgoing flux
        do kk = kk1, kk2
        do jj = jj1, jj2
        do ii = ii1, ii2
           do idim = 1, ndim
              ! TODO: don't iterate over kk1,kk2, jj1... but only on the face of
              ! the particles' cell
              if (flux(ind_grid_part(j), ii, jj, kk, 1, idim) > 0) &
                   outflux(j) = outflux(j) + flux(ind_grid_part(j), ii, jj, kk, 1, idim)
           end do
        end do
        end do
        end do

        ! Break if no outflux in this cell
        if (outflux(j) == 0d0) continue

        proba = outflux(j) * dtnew(ilevel) / uold(ind_grid(ind_grid_part(j)), 1)

        ! Decide randomly if we move the particle
        call ranf(localseed, rand)
        if (rand < proba) then
           ! Ok, I'm moving it. Which direction ?
           call ranf(localseed, rand)
           break_now = .false.
           do kk = kk1, kk2
              do jj = jj1, jj2
                 do ii = ii1, ii2
                    if (break_now) continue
                    ! TODO: don't iterate over kk1, kk2, jj1... but only on the face of
                    ! the particles' cell
                    do idim = 1, ndim
                       proba = flux(ind_grid_part(j), ii, jj, kk, 1, idim) * dtnew(ilevel) &
                            / uold(ind_grid(ind_grid_part(j)), 1)
                       if (rand < proba) then
                          ! now find the center of neighbour cell
                          igrid_nbor = son(nbors_father_cells(j, ind_father))

                          if (igrid_nbor == 0) then
                             new_xp(j, :) = (/ 0d0, 0d0, 0d0 /) ! TODO: location of the center of neighbour cell
                          else
                             ! TODO: descend the tree and pick a cell to put the particle into
                          end if
                          break_now = .true.
                          exit
                       else
                          ! decrease the random number, because P(switch to any) = 1
                          rand = rand - proba
                       end if
                    end do
                 end do
                 if (break_now) exit
              end do
              if (break_now) exit
           end do
        end if
     end do
  end if

  do idim = 1, ndim
     if(static)then
        do j = 1, np
           new_xp(j, idim) = xp(ind_part(j), idim)
        end do
     else
        do j = 1, np
           if (MC_tracer .and. mp(ind_part(j)) == 0d0) then
              ! do nothing, already moved
           else
              new_xp(j, idim) = xp(ind_part(j), idim) + new_vp(j, idim)*dtnew(ilevel)
           end if
        end do
     endif
  end do


  do idim = 1, ndim
     do j = 1, np
        xp(ind_part(j), idim) = new_xp(j, idim)
     end do
  end do

end subroutine move1
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
