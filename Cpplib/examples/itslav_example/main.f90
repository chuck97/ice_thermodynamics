program itslav_example
 
    ! itslav module
    use itslav, only : InitThermodynamics, &
                       FinalizeThermodynamics, &
                       UpdateAtmFlux, &
                       UpdateSwRadiation, &
                       UpdateLatentHeatFlux, &
                       Evaluate, &
                       GetIceSurfaceTemperature, &
                       GetIceThickness, &
                       GetSurfaceConductiveFlux, &
                       StoreState, &
                       RestoreState

    
    ! iso_c_binding module for types
    use, intrinsic :: iso_c_binding, only : c_int, c_double
 
    ! variables that are passed to the library should be c_type!
    real(c_double) :: time_step
    real(c_double), dimension(:,:), allocatable :: init_surf_temp, init_base_temp, init_thick
    real(c_double), dimension(:,:), allocatable :: surf_temp, ice_thick, cond_flux
    real(c_double), dimension(:,:), allocatable :: atm_flux, sw_radiation, lh_flux

    ! general variables
    real(kind=8) :: current_time
    integer      :: i, j, iter
    real(c_double), parameter :: rho_a = 1.28
    real(c_double), parameter :: cp_a = 1.1d3
    real(c_double), parameter :: C_sh = 1d-3
    real(c_double), parameter :: u_a = 15.0
    integer(c_int), parameter :: nlon = 12
    integer(c_int), parameter :: nlat = 12

    ! allocate memory
    allocate(init_surf_temp(1:nlon, 1:nlat))
    allocate(init_base_temp(1:nlon, 1:nlat))
    allocate(init_thick(1:nlon, 1:nlat))

    allocate(surf_temp(1:nlon, 1:nlat))
    allocate(ice_thick(1:nlon, 1:nlat))
    allocate(cond_flux(1:nlon, 1:nlat))

    allocate(atm_flux(1:nlon, 1:nlat))
    allocate(sw_radiation(1:nlon, 1:nlat))
    allocate(lh_flux(1:nlon, 1:nlat))

    ! initialization of arrays
    init_surf_temp = -1.0
    init_base_temp = 0.0

    ! ice thickness is zero outside of the circle
    do i = 1, nlon
        do j = 1, nlat
            if ((1.0*i - nlon/2)*(1.0*i - nlon/2) + (1.0*j - nlat/2)*(1.0*j - nlat/2) < 9) then
                init_thick(i, j) = 1.0
            else
                init_thick(i, j) = 0.0
            end if
        end do
    end do

    ! setup time step (seconds)
    time_step = 3600.0

    ! initialization of thermodynamics solver
    call InitThermodynamics(time_step = time_step, &              ! time step (seconds)
                            num_ice_cells = 10, &                 ! number of ice cells in vertical mesh
                            min_ice_thick = 1d-2, &               ! minimal ice thickness (meters)
                            min_lon_ind = 1, &                    ! minimal longitude index
                            max_lon_ind = nlon, &                  ! maximal longitude index
                            min_lat_ind = 1, &                    ! minimal latitude index
                            max_lat_ind = nlat, &                   ! maximal latitude index
                            init_base_temp = init_base_temp, &    ! 2D-array of initial base temp (deg Cel) 
                            init_surf_temp = init_surf_temp, &    ! 2D-array of initial surface temp (deg Cel)
                            init_ice_thick = init_thick)          ! 2D-array of initial ice thickness (meters) 

    

    ! store initial mesh
    call StoreState(min_lon_ind = 1, &
                    max_lon_ind = nlon, &
                    min_lat_ind = 1, &
                    max_lat_ind = nlat)
    
    ! restore initial state for check
    call RestoreState(min_lon_ind = 1, &
                      max_lon_ind = nlon, &
                      min_lat_ind = 1, &
                      max_lat_ind = nlat)

    call GetIceSurfaceTemperature(surf_temp(1:nlon, 1:nlat), &
                                  min_lon_ind = 1, &
                                  max_lon_ind = nlon, &
                                  min_lat_ind = 1, &
                                  max_lat_ind = nlat)

    print *,  "before computations,  Surface temperatures ="
    call print_slice(surf_temp(1:nlon, 1:nlat), nlon, nlat, &
                               2, nlon-2, &
                               2, nlat-2)
    
    ! time stepping
    current_time = 0.0

    do iter = 1, 200

        ! update current time
        current_time = current_time + time_step

        ! update total atm flux
        do i = 4, nlon-4
            do j = 4, nlat-4
                atm_flux(i, j) = rho_a*cp_a*C_sh*u_a*(-5.0)
            end do
        end do

        call UpdateAtmFlux(atm_flux(4:nlon-4, 4:nlat-4), &
                           min_lon_ind = 4, &
                           max_lon_ind = nlon-4, &
                           min_lat_ind = 4, &
                           max_lat_ind = nlat-4)
        
        ! update short-wave radiation
        sw_radiation = 10.0
        call UpdateSwRadiation(sw_radiation(3:nlon-3, 3:nlat-3), &
                               min_lon_ind = 2, &
                               max_lon_ind = nlon-4, &
                               min_lat_ind = 2, &
                               max_lat_ind = nlat-4)
        
        ! update latent heat flux
        lh_flux = 10.0
        call UpdateLatentHeatFlux(lh_flux(2:nlon-2, 2:nlat-2), &
                                  min_lon_ind = 4, &
                                  max_lon_ind = nlon-1, &
                                  min_lat_ind = 4, &
                                  max_lat_ind = nlat-1)


        ! evaluation of thermodynamics solver
        call Evaluate(min_lon_ind = 1, &
                      max_lon_ind = nlon, &
                      min_lat_ind = 1, &
                      max_lat_ind = nlat)

        ! save ice surface temperature
        call GetIceSurfaceTemperature(surf_temp(1:nlon, 1:nlat), &
                                      min_lon_ind = 1, &
                                      max_lon_ind = nlon, &
                                      min_lat_ind = 1, &
                                      max_lat_ind = nlat)
        
        ! save ice thickness
        call GetIceThickness(ice_thick(1:nlon, 1:nlat), &
                             min_lon_ind = 1, &
                             max_lon_ind = nlon, &
                             min_lat_ind = 1, &
                             max_lat_ind = nlat)

        ! save ice surface conductive heat flux
        call GetSurfaceConductiveFlux(cond_flux(1:nlon, 1:nlat), &
                                      min_lon_ind = 1, &
                                      max_lon_ind = nlon, &
                                      min_lat_ind = 1, &
                                      max_lat_ind = nlat)

    
        ! print surface temperature field to console
        if (MOD (iter,10) == 0) then
            print *,  "Iteration ", iter, "Surface temperatures ="
            call print_slice(surf_temp(1:nlon, 1:nlat), nlon, nlat, &
                                       2, nlon-2, &
                                       2, nlat-2)
        end if
    end do

    ! restore initial state and print temps
    call RestoreState(min_lon_ind = 1, &
                      max_lon_ind = nlon, &
                      min_lat_ind = 1, &
                      max_lat_ind = nlat)
    
    call GetIceSurfaceTemperature(surf_temp(1:nlon, 1:nlat), &
                                  min_lon_ind = 1, &
                                  max_lon_ind = nlon, &
                                  min_lat_ind = 1, &
                                  max_lat_ind = nlat)

    print *,  "after computations and restoring,  Surface temperatures ="
    call print_slice(surf_temp(1:nlon, 1:nlat), nlon, nlat, &
                               2, nlon-2, &
                               2, nlat-2)
 
    ! finalization of thermodynamics solver at the end of the program (freeing the memory)
    call FinalizeThermodynamics()

    contains

    ! subroutine for writing slice of 2D array to console
    subroutine print_slice(arr, dim_x, dim_y, &
                           slice_start_x, slice_end_x, &
                           slice_start_y, slice_end_y)

        implicit none
        integer(c_int), intent(in) ::  dim_x, dim_y
        real(c_double), intent(in) ::  arr(dim_x, dim_y)
        integer(c_int), intent(in) ::  slice_start_x, slice_end_x
        integer(c_int), intent(in) ::  slice_start_y, slice_end_y
        integer :: i, j

        do j = slice_start_y, slice_end_y
            do i = slice_start_x, slice_end_x
                write(*, fmt='((F8.4), 1X)', advance="no") arr(i, j)
            end do
            print *
        end do
        
    end subroutine
 
 end program