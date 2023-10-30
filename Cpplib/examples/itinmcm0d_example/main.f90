program itinmcm0d_example
 
    ! itslav module
    use itinmcm0d, only : InitThermodynamics, &
                          FinalizeThermodynamics, &
                          GetSurfaceTemperature, &
                          GetIceThickness, &
                          GetSnowThickness, &
                          GetIsIce, &
                          GetIsSnow

    ! iso_c_binding module for types
    use, intrinsic :: iso_c_binding, only : c_int, c_double, c_bool
 
    ! variables that are passed to the library should be c_type!
    real(c_double), parameter :: time_step = 3600.0
    integer(c_int), parameter :: nlon = 12
    integer(c_int), parameter :: nlat = 12
    logical(c_bool), parameter :: is_verbose = .false.
    real(c_double), dimension(:,:), allocatable :: init_ice_surf_temp, init_ice_base_temp, init_snow_surf_temp
    real(c_double), dimension(:,:), allocatable :: init_ice_thick, init_snow_thick
    logical(c_bool), dimension(:,:), allocatable :: water_marker

    ! variables for output from library
    real(c_double), dimension(:,:), allocatable :: surf_temp
    real(c_double), dimension(:,:), allocatable :: ice_thick, snow_thick
    logical(c_bool), dimension(:,:), allocatable :: ice_presence, snow_presence
    
    ! allocate memory
    allocate(init_ice_surf_temp(1:nlon, 1:nlat))
    allocate(init_ice_base_temp(1:nlon, 1:nlat))
    allocate(init_snow_surf_temp(1:nlon, 1:nlat))
    allocate(init_ice_thick(1:nlon, 1:nlat))
    allocate(init_snow_thick(1:nlon, 1:nlat))
    allocate(water_marker(1:nlon, 1:nlat))
    allocate(surf_temp(1:nlon, 1:nlat))
    allocate(ice_thick(1:nlon, 1:nlat))
    allocate(snow_thick(1:nlon, 1:nlat))
    allocate(ice_presence(1:nlon, 1:nlat))
    allocate(snow_presence(1:nlon, 1:nlat))

    ! initialization of arrays
    init_ice_surf_temp = -1.0
    init_ice_base_temp = 0.0
    init_snow_surf_temp = -5.0
    init_ice_thick = 1.0
    init_snow_thick = 0.1
    water_marker = .true.

    print *,  "Initialization of arrays done!" 

    ! initialization of thermodynamics solver
    call InitThermodynamics(time_step = time_step, &              
                            min_ice_thick = 1d-2, &               
                            min_lon_ind = 1, &                    
                            max_lon_ind = nlon, &                 
                            min_lat_ind = 1, &                    
                            max_lat_ind = nlat, &                 
                            init_ice_base_temp = init_ice_base_temp, &    
                            init_ice_surf_temp = init_ice_surf_temp, &   
                            init_snow_surf_temp = init_snow_surf_temp, &
                            init_ice_thick = init_ice_thick, &
                            init_snow_thick = init_snow_thick, &
                            water_marker = water_marker, &
                            is_verbose = is_verbose)       
    
    print *,  "Initialization of library done!"  
    
    ! get surface temperature
    call GetSurfaceTemperature(array = surf_temp, &
                               min_lon_ind = 1, &
                               max_lon_ind = nlon, &
                               min_lat_ind = 1, &                    
                               max_lat_ind = nlat)
    
    print *, "Output surface temperature: "
    call print_slice(surf_temp, nlon, nlat, 5, 10, 5, 10)

    ! get ice thickness
    call GetIceThickness(array = ice_thick, &
                         min_lon_ind = 1, &
                         max_lon_ind = nlon, &
                         min_lat_ind = 1, &                    
                         max_lat_ind = nlat)
    
    print *, "Output ice thickness: "
    call print_slice(ice_thick, nlon, nlat, 5, 10, 5, 10)

    ! get snow thickness
    call GetSnowThickness(array = snow_thick, &
                          min_lon_ind = 1, &
                          max_lon_ind = nlon, &
                          min_lat_ind = 1, &                    
                          max_lat_ind = nlat)
    
    print *, "Output snow thickness: "
    call print_slice(snow_thick, nlon, nlat, 5, 10, 5, 10)

    ! get is_ice array
    call GetIsIce(array = ice_presence, &
                  min_lon_ind = 1, &
                  max_lon_ind = nlon, &
                  min_lat_ind = 1, &                    
                  max_lat_ind = nlat)
    
    print *, "Output ice presence: "

    do j = 5, 10
        do i = 5, 10
            write(*, fmt='(L, 1X)', advance="no") ice_presence(i, j)
        end do
        print *
    end do

    ! get is_snow array
    call GetIsSnow(array = snow_presence, &
                   min_lon_ind = 1, &
                   max_lon_ind = nlon, &
                   min_lat_ind = 1, &                    
                   max_lat_ind = nlat)
    
    print *, "Output snow presence: "

    do j = 5, 10
        do i = 5, 10
            write(*, fmt='(L, 1X)', advance="no") snow_presence(i, j)
        end do
        print *
    end do

 
    ! finalization of thermodynamics solver at the end of the program (freeing the memory)
    call FinalizeThermodynamics()
    print *,  "Finalization of library done!"

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