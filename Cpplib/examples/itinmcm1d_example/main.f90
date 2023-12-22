program itinmcm1d_example
 
    ! itslav module
    use itinmcm, only : InitThermodynamics1d, &
                        FinalizeThermodynamics1d, &
                        UpdateAirTemperature, &
                        UpdateAirPressure, &
                        UpdatePrecipitationRate, &
                        UpdateAirSpecificHumidity, &
                        UpdateAbsWindSpeed, &
                        UpdateSwRadiation, &
                        UpdateLwRadiation, &
                        UpdateShCoeff, &
                        UpdateLhCoeff, &
                        UpdateSnowThickness, &
                        UpdateIceThickness, &
                        AssembleTotalAtmFlux, &
                        UpdateOceanSalinity, &
                        UpdateOceanFlux, &
                        GetSurfaceTemperature, &
                        GetIceThickness, &
                        GetSnowThickness, &
                        GetIsIce, &
                        GetIsSnow, &
                        AddPrecipitation, &
                        Evaluate, &
                        SetComputationMarker, &
                        UpdateTemperatureProfile, &
                        GetTemperatureProfile

    ! iso_c_binding module for types
    use, intrinsic :: iso_c_binding, only : c_int, c_double, c_bool
 
    ! variables that are passed to the library should be c_type!
    real(c_double), parameter :: time_step = 1800.0
    integer(c_int), parameter :: nlon = 5
    integer(c_int), parameter :: nlat = 5
    integer(c_int), parameter :: n_ice_layers = 5
    logical(c_bool), parameter :: is_verbose = .false.
    real(c_double), dimension(:,:), allocatable :: init_ice_surf_temp, init_ice_base_temp, init_snow_surf_temp
    real(c_double), dimension(:,:), allocatable :: init_ice_thick, init_snow_thick
    logical(c_bool), dimension(:,:), allocatable :: water_marker, do_compute
    real(c_double), dimension(:,:), allocatable :: atm_temp, atm_press, prec_rate, atm_humid, wind_speed, &
                                                   sw_rad, lw_rad, sh_coeff, lh_coeff, ocean_flux, ocean_sal, &
                                                   snoww_thick, icee_thick, base_sal, surf_sal
    real(c_double), dimension(:,:,:), allocatable :: test_temp_profile

    ! variables for output from library
    real(c_double), dimension(:,:), allocatable :: surf_temp
    real(c_double), dimension(:,:), allocatable :: ice_thick, snow_thick
    logical(c_bool), dimension(:,:), allocatable :: ice_presence, snow_presence

    real(c_double), dimension(:,:), allocatable :: a

    ! parameters for experiment
    real(c_double), parameter :: prec_rate_value = 5.429921990575380E-009
    real(c_double), parameter :: ocean_sal_value = 9.994506835937500E-004
    real(c_double), parameter :: sw_value =  34.3713805388029
    real(c_double), parameter :: lw_value =  205.645509767637
    real(c_double), parameter :: atm_temp_value =  -13.1182537078857
    real(c_double), parameter :: spec_humid_value =  1.43520918209106     
    real(c_double), parameter :: sh_coeff_value =  1.144087756983936E-003
    real(c_double), parameter :: lh_coeff_value =  1.144087756983936E-003
    real(c_double), parameter :: abs_wind_speed_value = 3.23311546728502     
    real(c_double), parameter :: atm_pressure_value = 103296.921875000          
    real(c_double), parameter :: ocean_flux_value =  0.0
    real(c_double), parameter :: ssnow_thickness =  6.393566942358753E-002
    real(c_double), parameter :: iice_thickness =  0.913329218477102
    real(c_double), parameter :: bbase_sal = 4.0
    real(c_double), parameter :: ssurf_sal = 1.0


    
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
    allocate(snoww_thick(1:nlon, 1:nlat))
    allocate(icee_thick(1:nlon, 1:nlat))
    allocate(base_sal(1:nlon, 1:nlat))
    allocate(surf_sal(1:nlon, 1:nlat))
    allocate(test_temp_profile(1:nlon, 1:nlat, 1:(n_ice_layers+1)))

    ! initialization of arrays
    init_ice_surf_temp = -5.0
    init_ice_base_temp = -1.0
    init_snow_surf_temp = -10.0
    init_ice_thick =  1.0
    init_snow_thick =  1.0
    water_marker = .true.
    do_compute = .true.
    snoww_thick = ssnow_thickness
    icee_thick = iice_thickness
    base_sal = bbase_sal
    surf_sal = ssurf_sal

    do i = 1, nlon
        do j = 1, nlat
            do k = 1, (n_ice_layers + 1)
                test_temp_profile(i, j, k) = i + j
            end do
        end do
        print *
    end do


    print *,  "Initialization of arrays done!" 

    ! initialization of thermodynamics solver
    call InitThermodynamics1d(time_step = time_step, &  
                              num_ice_layers = n_ice_layers, &   
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
                              init_ice_base_sal = base_sal, &
                              init_ice_surface_sal = surf_sal, &
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

    ! update snow thickness
    call UpdateSnowThickness(snoww_thick, 1, nlon, 1, nlat)
    print *,  "Update of snow thickness done!" 
    print * 

    ! update ice thickness
    call UpdateIceThickness(icee_thick, 1, nlon, 1, nlat)
    print *,  "Update of ice thickness done!" 
    print * 

    ! add snow precipitations
    ! call AddPrecipitation(1, nlon, 1, nlat)
    ! print *,  "Adding precipitations to snow thickness done!" 
    ! print * 

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
    call print_slice(surf_temp, nlon, nlat, 1, nlon, 1, nlat)
    print *

    ! ice thickness
    allocate(ice_thick(1:nlon, 1:nlat))
    call GetIceThickness(array = ice_thick, &
                         min_lon_ind = 1, &
                         max_lon_ind = nlon, &
                         min_lat_ind = 1, &                    
                         max_lat_ind = nlat)
    
    print *, "ICE THICK before EVALUATE: "
    call print_slice(ice_thick, nlon, nlat, 1, nlon, 1, nlat)
    print *

    ! ice thickness
    allocate(snow_thick(1:nlon, 1:nlat))
    call GetSnowThickness(array = snow_thick, &
                          min_lon_ind = 1, &
                          max_lon_ind = nlon, &
                          min_lat_ind = 1, &                    
                          max_lat_ind = nlat)
    
    print *, "SNOW THICK before EVALUATE: "
    call print_slice(snow_thick, nlon, nlat, 1, nlon, 1, nlat)
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
    call print_slice(surf_temp, nlon, nlat, 1, nlon, 1, nlat)
    print *

    ! ice thickness
    call GetIceThickness(array = ice_thick, &
                         min_lon_ind = 1, &
                         max_lon_ind = nlon, &
                         min_lat_ind = 1, &                    
                         max_lat_ind = nlat)
    
    print *, "ICE THICK after EVALUATE: "
    call print_slice(ice_thick, nlon, nlat, 1, nlon, 1, nlat)
    print *

    ! ice thickness
    call GetSnowThickness(array = snow_thick, &
                          min_lon_ind = 1, &
                          max_lon_ind = nlon, &
                          min_lat_ind = 1, &                    
                          max_lat_ind = nlat)
    
    print *, "SNOW THICK after EVALUATE: "
    call print_slice(snow_thick, nlon, nlat, 1, nlon, 1, nlat)
    print *

    ! get is_ice array
    call GetIsIce(array = ice_presence, &
                  min_lon_ind = 1, &
                  max_lon_ind = nlon, &
                  min_lat_ind = 1, &                    
                  max_lat_ind = nlat)
    
    print *, "Output ice presence: "

    do j = 1, nlat
        do i = 1, nlon
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
    do j = 1, nlat
        do i = 1, nlon
            write(*, fmt='(L, 1X)', advance="no") snow_presence(i, j)
        end do
        print *
    end do

    ! call UpdateTemperatureProfile(test_temp_profile, 1, nlon, 1, nlat)
    call GetTemperatureProfile(test_temp_profile, 1, nlon, 1, nlat)

    print *, "TEMP PROFILE after UPDATE PROFILE: "

    !do i = 1, nlon
    !    do j = 1, nlat
    !        print *, "lon = ", i, ", lat = ", j 
    !        do k = 1, (n_ice_layers + 1)
    !            write(*, fmt='((F15.8))', advance="no") test_temp_profile(i, j, k)
    !        end do
    !        print*
    !    end do
    !end do


    call GetSurfaceTemperature(array = surf_temp, &
                               min_lon_ind = 1, &
                               max_lon_ind = nlon, &
                               min_lat_ind = 1, &                    
                               max_lat_ind = nlat)
    
    print *, "SURFACE TEMP after UPDATE PROFILE: "
    call print_slice(surf_temp, nlon, nlat, 1, nlon, 1, nlat)
    print *
 
    ! finalization of thermodynamics solver at the end of the program (freeing the memory)
    call FinalizeThermodynamics1d()
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
                write(*, fmt='((F15.8))', advance="no") arr(i, j)
            end do
            print *
        end do
        
    end subroutine 
 
 end program
