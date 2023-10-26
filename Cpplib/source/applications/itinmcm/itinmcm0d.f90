module itinmcm0d

    use, intrinsic :: iso_c_binding, only : c_int, c_double, c_bool, c_ptr, c_null_ptr, c_associated
 
    implicit none
 
    private
 
    ! initialization of thermodynamics solver 
    interface
        function InitThermodynamics_(time_step,           &  ! time step (seconds)
                                     min_ice_thick,       &  ! minimal ice thickness (meters)
                                     min_lon_ind,         &  ! minimal longitude index
                                     max_lon_ind,         &  ! maximal longitude index
                                     min_lat_ind,         &  ! minimal latitude index
                                     max_lat_ind,         &  ! maximal latitude index
                                     init_ice_base_temp,  &  ! 2D-array of initial ice base temp (deg Cel) 
                                     init_ice_surf_temp,  &  ! 2D-array of initial ice surface temp (deg Cel)
                                     init_snow_surf_temp, &  ! 2D-array of initial snow surface temp (deg Cel)
                                     init_ice_thick,      &  ! 2D-array of initial ice thickness (meters)
                                     init_snow_thick,     &  ! 2D-array of initial snow thickness (meters)
                                     water_marker,        &  ! 2D-boolean-array for water marker
                                     is_verbose)          &  ! single bool value for detailed output     
                                     result(obj)          &  ! pointer to allocated Cpp class
                                     bind(C, name="InitThermodynamics")
             
            import :: c_int, c_double, c_ptr, c_bool       
            implicit none   

            ! Argument list
            real(c_double), intent(in), value :: time_step
            real(c_double), intent(in), value :: min_ice_thick
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            real(c_double), intent(in), dimension(*) :: init_ice_base_temp
            real(c_double), intent(in), dimension(*) :: init_ice_surf_temp
            real(c_double), intent(in), dimension(*) :: init_snow_surf_temp
            real(c_double), intent(in), dimension(*) :: init_ice_thick
            real(c_double), intent(in), dimension(*) :: init_snow_thick
            logical(c_bool), intent(in), dimension(*) :: water_marker
            logical(c_bool), intent(in), value :: is_verbose
            type(c_ptr) :: obj

        end function
    end interface
       
    ! finalization of thermodynamics solver
    interface
       subroutine FinalizeThermodynamics_(obj) & ! pointer to allocated Cpp class
                                          bind(C, name="FinalizeThermodynamics")
        
        import :: c_ptr
        implicit none

        ! Argument list
        type(c_ptr), intent(in), value :: obj

       end subroutine
    end interface

    ! set computation marker
    interface 
       subroutine SetComputationMarker_(obj ,        & ! pointer to allocated Cpp class
                                        do_compute,  & ! 2D-array of boolean computation marker
                                        min_lon_ind, & ! minimal longitude index
                                        max_lon_ind, & ! maximal longitude index
                                        min_lat_ind, & ! minimal latitude index
                                        max_lat_ind) & ! maximal latitude index
                                        bind(C, name="SetComputationMarker")
          import:: c_int, c_bool, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          logical(c_bool), intent(in), dimension(*) :: do_compute
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update descending short-wave radiation
    interface 
       subroutine UpdateSwRadiation_(obj ,        & ! pointer to allocated Cpp class
                                     sw_values,   & ! 2D-array of sw radiation
                                     min_lon_ind, & ! minimal longitude index
                                     max_lon_ind, & ! maximal longitude index
                                     min_lat_ind, & ! minimal latitude index
                                     max_lat_ind) & ! maximal latitude index
                                     bind(C, name="UpdateSwRadiation")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: sw_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update descending long-wave radiation
    interface 
       subroutine UpdateLwRadiation_(obj ,        & ! pointer to allocated Cpp class
                                     lw_values,   & ! 2D-array of lw radiation
                                     min_lon_ind, & ! minimal longitude index
                                     max_lon_ind, & ! maximal longitude index
                                     min_lat_ind, & ! minimal latitude index
                                     max_lat_ind) & ! maximal latitude index
                                     bind(C, name="UpdateLwRadiation")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: lw_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update air temperature
    interface 
       subroutine UpdateAirTemperature_(obj ,        & ! pointer to allocated Cpp class
                                        ta_values,   & ! 2D-array of air temperature
                                        min_lon_ind, & ! minimal longitude index
                                        max_lon_ind, & ! maximal longitude index
                                        min_lat_ind, & ! minimal latitude index
                                        max_lat_ind) & ! maximal latitude index
                                        bind(C, name="UpdateAirTemperature")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: ta_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update precipitation rate
    interface 
       subroutine UpdatePrecipitationRate_(obj ,        & ! pointer to allocated Cpp class
                                           pr_values,   & ! 2D-array of precipitation values
                                           min_lon_ind, & ! minimal longitude index
                                           max_lon_ind, & ! maximal longitude index
                                           min_lat_ind, & ! minimal latitude index
                                           max_lat_ind) & ! maximal latitude index
                                           bind(C, name="UpdatePrecipitationRate")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: pr_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update air pressure
    interface 
       subroutine UpdateAirPressure_(obj ,        & ! pointer to allocated Cpp class
                                     ap_values,   & ! 2D-array of air pressure
                                     min_lon_ind, & ! minimal longitude index
                                     max_lon_ind, & ! maximal longitude index
                                     min_lat_ind, & ! minimal latitude index
                                     max_lat_ind) & ! maximal latitude index
                                     bind(C, name="UpdateAirPressure")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: ap_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update air specific humidity
    interface 
       subroutine UpdateAirSpecificHumidity_(obj ,        & ! pointer to allocated Cpp class
                                             sh_values,   & ! 2D-array of air specific humidity
                                             min_lon_ind, & ! minimal longitude index
                                             max_lon_ind, & ! maximal longitude index
                                             min_lat_ind, & ! minimal latitude index
                                             max_lat_ind) & ! maximal latitude index
                                             bind(C, name="UpdateAirSpecificHumidity")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: sh_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update wind speed
    interface 
       subroutine UpdateAbsWindSpeed_(obj ,        & ! pointer to allocated Cpp class
                                      ws_values,   & ! 2D-array of wind speed
                                      min_lon_ind, & ! minimal longitude index
                                      max_lon_ind, & ! maximal longitude index
                                      min_lat_ind, & ! minimal latitude index
                                      max_lat_ind) & ! maximal latitude index
                                      bind(C, name="UpdateAbsWindSpeed")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: ws_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update sensible heat coefficient
    interface 
       subroutine UpdateShCoeff_(obj ,        & ! pointer to allocated Cpp class
                                 shc_values,  & ! 2D-array of sensible heat coeffs
                                 min_lon_ind, & ! minimal longitude index
                                 max_lon_ind, & ! maximal longitude index
                                 min_lat_ind, & ! minimal latitude index
                                 max_lat_ind) & ! maximal latitude index
                                 bind(C, name="UpdateShCoeff")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: shc_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update latent heat coefficient
    interface 
       subroutine UpdateLhCoeff_(obj ,        & ! pointer to allocated Cpp class
                                 lhc_values,  & ! 2D-array of latent heat coeffs
                                 min_lon_ind, & ! minimal longitude index
                                 max_lon_ind, & ! maximal longitude index
                                 min_lat_ind, & ! minimal latitude index
                                 max_lat_ind) & ! maximal latitude index
                                 bind(C, name="UpdateLhCoeff")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: lhc_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update ocean salinity
    interface 
       subroutine UpdateOceanSalinity_(obj ,        & ! pointer to allocated Cpp class
                                       os_values,   & ! 2D-array of ocean salinity
                                       min_lon_ind, & ! minimal longitude index
                                       max_lon_ind, & ! maximal longitude index
                                       min_lat_ind, & ! minimal latitude index
                                       max_lat_ind) & ! maximal latitude index
                                       bind(C, name="UpdateOceanSalinity")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: os_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update ocean flux
    interface 
       subroutine UpdateOceanFlux_(obj ,        & ! pointer to allocated Cpp class
                                   of_values,   & ! 2D-array of ocean flux
                                   min_lon_ind, & ! minimal longitude index
                                   max_lon_ind, & ! maximal longitude index
                                   min_lat_ind, & ! minimal latitude index
                                   max_lat_ind) & ! maximal latitude index
                                   bind(C, name="UpdateOceanFlux")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: of_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update ice thickness
    interface 
       subroutine UpdateIceThickness_(obj ,         & ! pointer to allocated Cpp class
                                      thick_values, & ! 2D-array of ice thickness
                                      min_lon_ind,  & ! minimal longitude index
                                      max_lon_ind,  & ! maximal longitude index
                                      min_lat_ind,  & ! minimal latitude index
                                      max_lat_ind)  & ! maximal latitude index
                                      bind(C, name="UpdateIceThickness")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: thick_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! update snow thickness
    interface 
       subroutine UpdateSnowThickness_(obj ,         & ! pointer to allocated Cpp class
                                       thick_values, & ! 2D-array of snow thickness
                                       min_lon_ind,  & ! minimal longitude index
                                       max_lon_ind,  & ! maximal longitude index
                                       min_lat_ind,  & ! minimal latitude index
                                       max_lat_ind)  & ! maximal latitude index
                                       bind(C, name="UpdateSnowThickness")
          import:: c_int, c_double, c_ptr
          implicit none
          
          ! Argument list
          type(c_ptr), intent(in), value :: obj
          real(c_double), intent(in), dimension(*) :: thick_values
          integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
          integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

       end subroutine
    end interface

    ! single-step evaluation
    interface
        subroutine Evaluate_(obj, &           ! pointer to allocated Cpp class
                             min_lon_ind, &   ! minimal longitude index
                             max_lon_ind, &   ! maximal longitude index
                             min_lat_ind, &   ! minimal latitude index
                             max_lat_ind) &   ! maximal latitude index
                             bind(C, name="Evaluate")
            
            import:: c_int, c_ptr
            implicit none

            ! Argument list
            type(c_ptr), intent(in), value :: obj
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            
        end subroutine
    end interface

    ! recieve surface temperatures
    interface
        subroutine GetSurfaceTemperature_(obj, &           ! pointer to allocated Cpp class
                                          array, &         ! 2D-array of surface temperature (deg Cel) - output
                                          min_lon_ind, &   ! minimal longitude index
                                          max_lon_ind, &   ! maximal longitude index
                                          min_lat_ind, &   ! minimal latitude index
                                          max_lat_ind) &   ! maximal latitude index
                                          bind(C, name="GetSurfaceTemperature")
            
            import:: c_int, c_double, c_ptr
            implicit none

            ! Argument list
            type(c_ptr), intent(in), value :: obj
            real(c_double), intent(out), dimension(*) :: array
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            
        end subroutine
    end interface

    ! recieve ice thickness
    interface
        subroutine GetIceThickness_(obj, &           ! pointer to allocated Cpp class
                                    array, &         ! 2D-array of ice thickness (m) - output
                                    min_lon_ind, &   ! minimal longitude index
                                    max_lon_ind, &   ! maximal longitude index
                                    min_lat_ind, &   ! minimal latitude index
                                    max_lat_ind) &   ! maximal latitude index
                                    bind(C, name="GetIceThickness")
            
            import:: c_int, c_double, c_ptr
            implicit none

            ! Argument list
            type(c_ptr), intent(in), value :: obj
            real(c_double), intent(out), dimension(*) :: array
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            
        end subroutine
    end interface

    ! recieve snow thickness
    interface
        subroutine GetSnowThickness_(obj, &           ! pointer to allocated Cpp class
                                     array, &         ! 2D-array of snow thickness (m) - output
                                     min_lon_ind, &   ! minimal longitude index
                                     max_lon_ind, &   ! maximal longitude index
                                     min_lat_ind, &   ! minimal latitude index
                                     max_lat_ind) &   ! maximal latitude index
                                     bind(C, name="GetSnowThickness")
            
            import:: c_int, c_double, c_ptr
            implicit none

            ! Argument list
            type(c_ptr), intent(in), value :: obj
            real(c_double), intent(out), dimension(*) :: array
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            
        end subroutine
    end interface

    ! recieve bool is_snow array
    interface
        subroutine GetIsSnow_(obj, &           ! pointer to allocated Cpp class
                              array, &         ! 2D-boolean-array of snow presence - output
                              min_lon_ind, &   ! minimal longitude index
                              max_lon_ind, &   ! maximal longitude index
                              min_lat_ind, &   ! minimal latitude index
                              max_lat_ind) &   ! maximal latitude index
                              bind(C, name="GetIsSnow")
            
            import:: c_int, c_bool, c_ptr
            implicit none

            ! Argument list
            type(c_ptr), intent(in), value :: obj
            logical(c_bool), intent(out), dimension(*) :: array
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            
        end subroutine
    end interface

    ! recieve bool is_ice array
    interface
        subroutine GetIsIce_(obj, &           ! pointer to allocated Cpp class
                             array, &         ! 2D-boolean-array of ice presence - output
                             min_lon_ind, &   ! minimal longitude index
                             max_lon_ind, &   ! maximal longitude index
                             min_lat_ind, &   ! minimal latitude index
                             max_lat_ind) &   ! maximal latitude index
                             bind(C, name="GetIsIce")
            
            import:: c_int, c_bool, c_ptr
            implicit none

            ! Argument list
            type(c_ptr), intent(in), value :: obj
            logical(c_bool), intent(out), dimension(*) :: array
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            
        end subroutine
    end interface
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! HOLD CPP CLASS !!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type(c_ptr), save :: obj = c_null_ptr
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!  INTERFACE FUNCTIONS  !!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public :: InitThermodynamics
    public :: FinalizeThermodynamics
    public :: SetComputationMarker
    public :: UpdateSwRadiation
    public :: UpdateLwRadiation
    public :: UpdateAirTemperature
    public :: UpdatePrecipitationRate
    public :: UpdateAirPressure
    public :: UpdateAirSpecificHumidity
    public :: UpdateAbsWindSpeed
    public :: UpdateShCoeff
    public :: UpdateLhCoeff
    public :: UpdateOceanSalinity
    public :: UpdateOceanFlux
    public :: UpdateIceThickness
    public :: UpdateSnowThickness
    public :: Evaluate
    public :: GetSnowThickness
    public :: GetIceThickness
    public :: GetIsSnow
    public :: GetIsIce
 
 contains
 
    ! initialization of thermodynamics solver 
    subroutine InitThermodynamics(time_step,           &  ! time step (seconds)
                                  min_ice_thick,       &  ! minimal ice thickness (meters)
                                  min_lon_ind,         &  ! minimal longitude index (inclusively)
                                  max_lon_ind,         &  ! maximal longitude index (inclusively)
                                  min_lat_ind,         &  ! minimal latitude index (inclusively)
                                  max_lat_ind,         &  ! maximal latitude index (inclusively)
                                  init_ice_base_temp,  &  ! 2D-array of initial ice base temp (deg Cel) 
                                  init_ice_surf_temp,  &  ! 2D-array of initial ice surface temp (deg Cel)
                                  init_snow_surf_temp, &  ! 2D-array of initial snow surface temp (deg Cel)
                                  init_ice_thick,      &  ! 2D-array of initial ice thickness (meters)
                                  init_snow_thick,     &  ! 2D-array of initial snow thickness (meters)
                                  water_marker,        &  ! 2D-boolean-array for water marker
                                  is_verbose)             ! single bool value for detailed output     

        implicit none   
        
        ! Argument list
        real(c_double), intent(in), value :: time_step
        real(c_double), intent(in), value :: min_ice_thick
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
        real(c_double), intent(in), dimension(*) :: init_ice_base_temp
        real(c_double), intent(in), dimension(*) :: init_ice_surf_temp
        real(c_double), intent(in), dimension(*) :: init_snow_surf_temp
        real(c_double), intent(in), dimension(*) :: init_ice_thick
        real(c_double), intent(in), dimension(*) :: init_snow_thick
        logical(c_bool), intent(in), dimension(*) :: water_marker
        logical(c_bool), intent(in), value :: is_verbose
        
        ! Body
        if (c_associated(obj)) then
            call FinalizeThermodynamics_(obj)
        end if
   
        obj = InitThermodynamics_(time_step, &                    
                                  min_ice_thick, &         
                                  min_lon_ind, &           
                                  max_lon_ind, &           
                                  min_lat_ind, &           
                                  max_lat_ind, &  
                                  init_ice_base_temp, &
                                  init_ice_surf_temp, &
                                  init_snow_surf_temp, &
                                  init_ice_thick, &
                                  init_snow_thick, &
                                  water_marker, &
                                  is_verbose)
    end subroutine

    ! finalization of thermodynamics solver
    subroutine FinalizeThermodynamics()

        implicit none

        ! Body
        if (c_associated(obj)) then
            call FinalizeThermodynamics_(obj)
            obj = c_null_ptr
        end if

    end subroutine

    ! set 2D computation marker
    subroutine SetComputationMarker(do_compute, &  ! 2D-array of boolean computation marker
                                    min_lon_ind, & ! minimal longitude index (inclusively)
                                    max_lon_ind, & ! maximal longitude index (inclusively)
                                    min_lat_ind, & ! minimal latitude index (inclusively)
                                    max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        logical(c_bool), intent(in), dimension(*) :: do_compute
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call SetComputationMarker_(obj, &
                                       do_compute, &  
                                       min_lon_ind, &      
                                       max_lon_ind, &      
                                       min_lat_ind, &      
                                       max_lat_ind)
        end if

    end subroutine 

    ! update 2D descending short-wave radiation
    subroutine UpdateSwRadiation(sw_values, &   ! 2D-array of sw radiation (W m-2)
                                 min_lon_ind, & ! minimal longitude index (inclusively)
                                 max_lon_ind, & ! maximal longitude index (inclusively)
                                 min_lat_ind, & ! minimal latitude index (inclusively)
                                 max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: sw_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateSwRadiation_(obj, &
                                    sw_values, &  
                                    min_lon_ind, &      
                                    max_lon_ind, &      
                                    min_lat_ind, &      
                                    max_lat_ind)
        end if

    end subroutine 

    ! update 2D descending long-wave radiation
    subroutine UpdateLwRadiation(lw_values, &   ! 2D-array of lw radiation (W m-2)
                                 min_lon_ind, & ! minimal longitude index (inclusively)
                                 max_lon_ind, & ! maximal longitude index (inclusively)
                                 min_lat_ind, & ! minimal latitude index (inclusively)
                                 max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: lw_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateLwRadiation_(obj, &
                                    lw_values, &  
                                    min_lon_ind, &      
                                    max_lon_ind, &      
                                    min_lat_ind, &      
                                    max_lat_ind)
        end if

    end subroutine 

    ! update 2D air temperature
    subroutine UpdateAirTemperature(ta_values, &   ! 2D-array of air temperature (deg Cel)
                                    min_lon_ind, & ! minimal longitude index (inclusively)
                                    max_lon_ind, & ! maximal longitude index (inclusively)
                                    min_lat_ind, & ! minimal latitude index (inclusively)
                                    max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: ta_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateAirTemperature_(obj, &
                                       ta_values, &  
                                       min_lon_ind, &      
                                       max_lon_ind, &      
                                       min_lat_ind, &      
                                       max_lat_ind)
        end if

    end subroutine 

    ! update 2D precipitation rate
    subroutine UpdatePrecipitationRate(pr_values, &   ! 2D-array of precipitation rate (mm s-1 not sure!)
                                       min_lon_ind, & ! minimal longitude index (inclusively)
                                       max_lon_ind, & ! maximal longitude index (inclusively)
                                       min_lat_ind, & ! minimal latitude index (inclusively)
                                       max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: pr_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdatePrecipitationRate_(obj, &
                                          pr_values, &  
                                          min_lon_ind, &      
                                          max_lon_ind, &      
                                          min_lat_ind, &      
                                          max_lat_ind)
        end if

    end subroutine 

    ! update 2D atmosphere pressure
    subroutine UpdateAirPressure(ap_values,   & ! 2D-array of atmosphere pressure (Pascals)
                                 min_lon_ind, & ! minimal longitude index (inclusively)
                                 max_lon_ind, & ! maximal longitude index (inclusively)
                                 min_lat_ind, & ! minimal latitude index (inclusively)
                                 max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: ap_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateAirPressure_(obj, &
                                    ap_values, &  
                                    min_lon_ind, &      
                                    max_lon_ind, &      
                                    min_lat_ind, &      
                                    max_lat_ind)
        end if

    end subroutine 

    ! update 2D air specific humidity
    subroutine UpdateAirSpecificHumidity(sh_values,   & ! 2D-array of air specific humidity (g kg-1)
                                         min_lon_ind, & ! minimal longitude index (inclusively)
                                         max_lon_ind, & ! maximal longitude index (inclusively)
                                         min_lat_ind, & ! minimal latitude index (inclusively)
                                         max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: sh_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateAirSpecificHumidity_(obj, &
                                            sh_values, &  
                                            min_lon_ind, &      
                                            max_lon_ind, &      
                                            min_lat_ind, &      
                                            max_lat_ind)
        end if

    end subroutine 

    ! update 2D absolute wind speed
    subroutine UpdateAbsWindSpeed(ws_values,   & ! 2D-array of air wind speed (m s-1)
                                  min_lon_ind, & ! minimal longitude index (inclusively)
                                  max_lon_ind, & ! maximal longitude index (inclusively)
                                  min_lat_ind, & ! minimal latitude index (inclusively)
                                  max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: ws_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateAbsWindSpeed_(obj, &
                                     ws_values, &  
                                     min_lon_ind, &      
                                     max_lon_ind, &      
                                     min_lat_ind, &      
                                     max_lat_ind)
        end if

    end subroutine 

    ! update 2D sensible heat coefficient
    subroutine UpdateShCoeff(shc_values,  & ! 2D-array of sensible heat coeffs (-)
                             min_lon_ind, & ! minimal longitude index (inclusively)
                             max_lon_ind, & ! maximal longitude index (inclusively)
                             min_lat_ind, & ! minimal latitude index (inclusively)
                             max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: shc_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateShCoeff_(obj, &
                                shc_values, &  
                                min_lon_ind, &      
                                max_lon_ind, &      
                                min_lat_ind, &      
                                max_lat_ind)
        end if

    end subroutine

    ! update 2D latent heat coefficient
    subroutine UpdateLhCoeff(lhc_values,  & ! 2D-array of latent heat coeffs (-)
                             min_lon_ind, & ! minimal longitude index (inclusively)
                             max_lon_ind, & ! maximal longitude index (inclusively)
                             min_lat_ind, & ! minimal latitude index (inclusively)
                             max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: lhc_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateLhCoeff_(obj, &
                                lhc_values, &  
                                min_lon_ind, &      
                                max_lon_ind, &      
                                min_lat_ind, &      
                                max_lat_ind)
        end if

    end subroutine

    ! update 2D ocean salinity
    subroutine UpdateOceanSalinity(os_values,  & ! 2D-array of ocean salinity (psu)
                                   min_lon_ind, & ! minimal longitude index (inclusively)
                                   max_lon_ind, & ! maximal longitude index (inclusively)
                                   min_lat_ind, & ! minimal latitude index (inclusively)
                                   max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: os_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateOceanSalinity_(obj, &
                                      os_values, &  
                                      min_lon_ind, &      
                                      max_lon_ind, &      
                                      min_lat_ind, &      
                                      max_lat_ind)
        end if

    end subroutine

    ! update 2D ocean flux
    subroutine UpdateOceanFlux(of_values,   & ! 2D-array of ocean flux (W m-2)
                               min_lon_ind, & ! minimal longitude index (inclusively)
                               max_lon_ind, & ! maximal longitude index (inclusively)
                               min_lat_ind, & ! minimal latitude index (inclusively)
                               max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: of_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateOceanFlux_(obj, &
                                  of_values, &  
                                  min_lon_ind, &      
                                  max_lon_ind, &      
                                  min_lat_ind, &      
                                  max_lat_ind)
        end if

    end subroutine

    ! update 2D ice thuckness
    subroutine UpdateIceThickness(thick_values,   & ! 2D-array of ice thicknesses (m)
                                  min_lon_ind, & ! minimal longitude index (inclusively)
                                  max_lon_ind, & ! maximal longitude index (inclusively)
                                  min_lat_ind, & ! minimal latitude index (inclusively)
                                  max_lat_ind)   ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: thick_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateIceThickness_(obj, &
                                     thick_values, &  
                                     min_lon_ind, &      
                                     max_lon_ind, &      
                                     min_lat_ind, &      
                                     max_lat_ind)
        end if

    end subroutine

    ! update 2D snow thickness
    subroutine UpdateSnowThickness(thick_values, & ! 2D-array of snow thicknesses (m)
                                   min_lon_ind,  & ! minimal longitude index (inclusively)
                                   max_lon_ind,  & ! maximal longitude index (inclusively)
                                   min_lat_ind,  & ! minimal latitude index (inclusively)
                                   max_lat_ind)    ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: thick_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateSnowThickness_(obj, &
                                      thick_values, &  
                                      min_lon_ind, &      
                                      max_lon_ind, &      
                                      min_lat_ind, &      
                                      max_lat_ind)
        end if

    end subroutine

    ! 2D single-step evaluation
    subroutine Evaluate(min_lon_ind,  & ! minimal longitude index (inclusively)
                        max_lon_ind,  & ! maximal longitude index (inclusively)
                        min_lat_ind,  & ! minimal latitude index (inclusively)
                        max_lat_ind)    ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call Evaluate_(obj, &
                           min_lon_ind, &      
                           max_lon_ind, &      
                           min_lat_ind, &      
                           max_lat_ind)
        end if

    end subroutine

    ! update 2D snow thickness
    subroutine GetSnowThickness(array, & ! 2D output array for snow thickness
                                min_lon_ind,  & ! minimal longitude index (inclusively)
                                max_lon_ind,  & ! maximal longitude index (inclusively)
                                min_lat_ind,  & ! minimal latitude index (inclusively)
                                max_lat_ind)    ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(out), dimension(*) :: array
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call GetSnowThickness_(obj, &
                                   array, &  
                                   min_lon_ind, &      
                                   max_lon_ind, &      
                                   min_lat_ind, &      
                                   max_lat_ind)
        end if

    end subroutine

    ! update 2D ice thickness
    subroutine GetIceThickness(array, &        ! 2D output array for ice thickness
                               min_lon_ind,  & ! minimal longitude index (inclusively)
                               max_lon_ind,  & ! maximal longitude index (inclusively)
                               min_lat_ind,  & ! minimal latitude index (inclusively)
                               max_lat_ind)    ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(out), dimension(*) :: array
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call GetIceThickness_(obj, &
                                  array, &  
                                  min_lon_ind, &      
                                  max_lon_ind, &      
                                  min_lat_ind, &      
                                  max_lat_ind)
        end if

    end subroutine
    
    ! update 2D bool snow presence array
    subroutine GetIsSnow(array, &        ! 2D output bool array for snow presence
                         min_lon_ind,  & ! minimal longitude index (inclusively)
                         max_lon_ind,  & ! maximal longitude index (inclusively)
                         min_lat_ind,  & ! minimal latitude index (inclusively)
                         max_lat_ind)    ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        logical(c_bool), intent(out), dimension(*) :: array
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call GetIsSnow_(obj, &
                            array, &  
                            min_lon_ind, &      
                            max_lon_ind, &      
                            min_lat_ind, &      
                            max_lat_ind)
        end if

    end subroutine

    ! update 2D bool ice presence array
    subroutine GetIsIce(array, &        ! 2D output bool array for ice presence
                        min_lon_ind,  & ! minimal longitude index (inclusively)
                        max_lon_ind,  & ! maximal longitude index (inclusively)
                        min_lat_ind,  & ! minimal latitude index (inclusively)
                        max_lat_ind)    ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        logical(c_bool), intent(out), dimension(*) :: array
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call GetIsIce_(obj, &
                           array, &  
                           min_lon_ind, &      
                           max_lon_ind, &      
                           min_lat_ind, &      
                           max_lat_ind)
        end if

    end subroutine
    
end module itinmcm0d