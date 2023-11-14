program itinmcm0d_example
 
    ! itslav module
    use itinmcm0d, only : InitThermodynamics, &
                          FinalizeThermodynamics, &
                          UpdateAirTemperature, &
                          UpdateAirPressure, &
                          UpdatePrecipitationRate, &
                          UpdateAirSpecificHumidity, &
                          UpdateAbsWindSpeed, &
                          UpdateSwRadiation, &
                          UpdateLwRadiation, &
                          UpdateShCoeff, &
                          UpdateLhCoeff, &
                          AssembleTotalAtmFlux, &
                          UpdateOceanSalinity, &
                          UpdateOceanFlux, &
                          GetSurfaceTemperature, &
                          GetIceThickness, &
                          GetSnowThickness, &
                          GetIsIce, &
                          GetIsSnow, &
                          Evaluate, &
                          SetComputationMarker

    ! iso_c_binding module for types
    use, intrinsic :: iso_c_binding, only : c_int, c_double, c_bool
 
    ! variables that are passed to the library should be c_type!
    real(c_double), parameter :: time_step = 1800.0
    integer(c_int), parameter :: nlon = 12
    integer(c_int), parameter :: nlat = 12
    logical(c_bool), parameter :: is_verbose = .false.
    real(c_double), dimension(:,:), allocatable :: init_ice_surf_temp, init_ice_base_temp, init_snow_surf_temp
    real(c_double), dimension(:,:), allocatable :: init_ice_thick, init_snow_thick
    logical(c_bool), dimension(:,:), allocatable :: water_marker, do_compute
    real(c_double), dimension(:,:), allocatable :: atm_temp, atm_press, prec_rate, atm_humid, wind_speed, &
                                                   sw_rad, lw_rad, sh_coeff, lh_coeff, ocean_flux, ocean_sal

    ! variables for output from library
    real(c_double), dimension(:,:), allocatable :: surf_temp
    real(c_double), dimension(:,:), allocatable :: ice_thick, snow_thick
    logical(c_bool), dimension(:,:), allocatable :: ice_presence, snow_presence

    real(c_double), dimension(:,:), allocatable :: a

    ! parameters for experiment
    real(c_double), parameter :: prec_rate_value =  5.235659610960671E-009
    real(c_double), parameter :: ocean_sal_value = 30.0382061004639
    real(c_double), parameter :: sw_value =  0.0
    real(c_double), parameter :: lw_value =  193.323524807374
    real(c_double), parameter :: atm_temp_value =  -36.7469253540039
    real(c_double), parameter :: spec_humid_value =  1.995431739487685E-002
    real(c_double), parameter :: sh_coeff_value =  1.189734437502921E-003
    real(c_double), parameter :: lh_coeff_value =  1.189734437502921E-003
    real(c_double), parameter :: abs_wind_speed_value = 2.78687097042342
    real(c_double), parameter :: atm_pressure_value = 103278.687500000
    real(c_double), parameter :: ocean_flux_value =  0.0

    
    ! allocate memory
    allocate(init_ice_surf_temp(1:nlon, 1:nlat))
    allocate(init_ice_base_temp(1:nlon, 1:nlat))
    allocate(init_snow_surf_temp(1:nlon, 1:nlat))
    allocate(init_ice_thick(1:nlon, 1:nlat))
    allocate(init_snow_thick(1:nlon, 1:nlat))
    allocate(water_marker(1:nlon, 1:nlat))
    allocate(do_compute(1:nlon, 1:nlat))
    allocate(ice_presence(1:nlon, 1:nlat))
    allocate(snow_presence(1:nlon, 1:nlat))

    ! initialization of arrays
    init_ice_surf_temp = -5.0
    init_ice_base_temp = 0.0
    init_snow_surf_temp = -10.0
    init_ice_thick = 4.24790578884097
    init_snow_thick = 0.196908870110462
    water_marker = .true.
    do_compute = .true.


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
    print * 
    
    ! update atmosphere temperature
    call SetComputationMarker(do_compute, 1, nlon, 1, nlat)

    print *,  "Set computation marker done!" 
    print * 

    ! update atmosphere temperature
    allocate(atm_temp(1:nlon, 1:nlat))
    atm_temp = atm_temp_value
    call UpdateAirTemperature(atm_temp, 1, nlon, 1, nlat)

    print *,  "Update of atm temperature done!" 
    print * 

    ! update atmosphere pressure
    allocate(atm_press(1:nlon, 1:nlat))
    atm_press = atm_pressure_value
    call UpdateAirPressure(atm_press, 1, nlon, 1, nlat)

    print *,  "Update of atm pressure done!" 
    print * 

    ! update ocean salinity
    allocate(ocean_sal(1:nlon, 1:nlat))
    ocean_sal = ocean_sal_value
    call UpdateOceanSalinity(ocean_sal, 1, nlon, 1, nlat)

    print *,  "Update of ocean salinity done!" 
    print * 

    ! update precipitation rate
    allocate(prec_rate(1:nlon, 1:nlat))
    prec_rate = prec_rate_value
    call UpdatePrecipitationRate(prec_rate, 1, nlon, 1, nlat)

    print *,  "Update of precipitation rate done!" 
    print * 

    ! update atmosphere specific humidity
    allocate(atm_humid(1:nlon, 1:nlat))
    atm_humid = spec_humid_value
    call UpdateAirSpecificHumidity(atm_humid, 1, nlon, 1, nlat)

    print *,  "Update of air specific humidity done!" 
    print * 

    ! update absolute value of wind speed
    allocate(wind_speed(1:nlon, 1:nlat))
    wind_speed = abs_wind_speed_value
    call UpdateAbsWindSpeed(wind_speed, 1, nlon, 1, nlat)

    print *,  "Update of wind speed done!" 
    print * 

    ! update long-wave radiation from atmosphere
    allocate(lw_rad(1:nlon, 1:nlat))
    lw_rad = lw_value
    call UpdateLwRadiation(lw_rad, 1, nlon, 1, nlat)

    print *,  "Update of long-wave radiation done!" 
    print * 

    ! update short-wave radiation from atmosphere
    allocate(sw_rad(1:nlon, 1:nlat))
    sw_rad = sw_value
    call UpdateSwRadiation(sw_rad, 1, nlon, 1, nlat)

    print *,  "Update of short-wave radiation done!" 
    print * 

    ! update sensible heat transfer coefficient
    allocate(sh_coeff(1:nlon, 1:nlat))
    sh_coeff = sh_coeff_value
    call UpdateShCoeff(sh_coeff, 1, nlon, 1, nlat)

    print *,  "Update of sensible heat coefficient done!" 
    print * 

    ! update latent heat transfer coefficient
    allocate(lh_coeff(1:nlon, 1:nlat))
    lh_coeff = lh_coeff_value
    call UpdateLhCoeff(lh_coeff, 1, nlon, 1, nlat)

    print *,  "Update of latent heat coefficient done!" 
    print * 

    ! assemble total atmosphere flux
    call AssembleTotalAtmFlux(1, nlon, 1, nlat)

    print *,  "Assembling of total atmosphere flux done!" 
    print * 

    ! update total flux from ocean
    allocate(ocean_flux(1:nlon, 1:nlat))
    ocean_flux = ocean_flux_value
    call UpdateOceanFlux(ocean_flux, 1, nlon, 1, nlat)

    print *,  "Update of total ocean flux done!" 
    print *

    ! print prognostic values before EVALUATE
    
    ! surface temperature
    allocate(surf_temp(1:nlon, 1:nlat))
    call GetSurfaceTemperature(array = surf_temp, &
                               min_lon_ind = 1, &
                               max_lon_ind = nlon, &
                               min_lat_ind = 1, &                    
                               max_lat_ind = nlat)
    
    print *, "SURFACE TEMP before EVALUATE: "
    call print_slice(surf_temp, nlon, nlat, 5, 10, 5, 10)
    print *

    ! ice thickness
    allocate(ice_thick(1:nlon, 1:nlat))
    call GetIceThickness(array = ice_thick, &
                         min_lon_ind = 1, &
                         max_lon_ind = nlon, &
                         min_lat_ind = 1, &                    
                         max_lat_ind = nlat)
    
    print *, "ICE THICK before EVALUATE: "
    call print_slice(ice_thick, nlon, nlat, 5, 10, 5, 10)
    print *

    ! ice thickness
    allocate(snow_thick(1:nlon, 1:nlat))
    call GetSnowThickness(array = snow_thick, &
                          min_lon_ind = 1, &
                          max_lon_ind = nlon, &
                          min_lat_ind = 1, &                    
                          max_lat_ind = nlat)
    
    print *, "SNOW THICK before EVALUATE: "
    call print_slice(snow_thick, nlon, nlat, 5, 10, 5, 10)
    print *


    ! single-step evaluation
    call Evaluate(1, nlon, 1, nlat)

    print *,  "Single step evaluation done!" 
    print *

    ! print prognostic values after EVALUATE
    
    ! surface temperature
    call GetSurfaceTemperature(array = surf_temp, &
                               min_lon_ind = 1, &
                               max_lon_ind = nlon, &
                               min_lat_ind = 1, &                    
                               max_lat_ind = nlat)
    
    print *, "SURFACE TEMP after EVALUATE: "
    call print_slice(surf_temp, nlon, nlat, 5, 10, 5, 10)
    print *

    ! ice thickness
    call GetIceThickness(array = ice_thick, &
                         min_lon_ind = 1, &
                         max_lon_ind = nlon, &
                         min_lat_ind = 1, &                    
                         max_lat_ind = nlat)
    
    print *, "ICE THICK after EVALUATE: "
    call print_slice(ice_thick, nlon, nlat, 5, 10, 5, 10)
    print *

    ! ice thickness
    call GetSnowThickness(array = snow_thick, &
                          min_lon_ind = 1, &
                          max_lon_ind = nlon, &
                          min_lat_ind = 1, &                    
                          max_lat_ind = nlat)
    
    print *, "SNOW THICK after EVALUATE: "
    call print_slice(snow_thick, nlon, nlat, 5, 10, 5, 10)
    print *

    ! get is_ice array
    !call GetIsIce(array = ice_presence, &
    !              min_lon_ind = 1, &
    !              max_lon_ind = nlon, &
    !              min_lat_ind = 1, &                    
    !              max_lat_ind = nlat)
    !
    !print *, "Output ice presence: "
!
    !do j = 5, 10
    !    do i = 5, 10
    !        write(*, fmt='(L, 1X)', advance="no") ice_presence(i, j)
    !    end do
    !    print *
    !end do

    ! get is_snow array
    !call GetIsSnow(array = snow_presence, &
    !               min_lon_ind = 1, &
    !               max_lon_ind = nlon, &
    !               min_lat_ind = 1, &                    
    !               max_lat_ind = nlat)
    !
    !print *, "Output snow presence: "
!
    !do j = 5, 10
    !    do i = 5, 10
    !        write(*, fmt='(L, 1X)', advance="no") snow_presence(i, j)
    !    end do
    !    print *
    !end do

 
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
                write(*, fmt='((F10.6))', advance="no") arr(i, j)
            end do
            print *
        end do
        
    end subroutine 
 
 end program