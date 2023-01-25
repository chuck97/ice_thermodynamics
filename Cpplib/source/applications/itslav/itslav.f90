module itslav

    use, intrinsic :: iso_c_binding, only : c_int, c_double, c_ptr, c_null_ptr, c_associated
 
    implicit none
 
    private
 
    ! initialization of thermodynamics solver 
    interface
        function InitThermodynamics_(time_step, &             ! time step (seconds)
                                     num_ice_cells, &         ! number of ice cells in vertical mesh
                                     min_ice_thick, &         ! minimal ice thickness (meters)
                                     min_lon_ind, &           ! minimal longitude index
                                     max_lon_ind, &           ! maximal longitude index
                                     min_lat_ind, &           ! minimal latitude index
                                     max_lat_ind, &           ! maximal latitude index
                                     init_base_temp, &        ! 2D-array of initial base temp (deg Cel) 
                                     init_surf_temp, &        ! 2D-array of initial surface temp (deg Cel)
                                     init_ice_thick) &        ! 2D-array of initial ice thickness (meters)
                                     result(obj) &            ! pointer to allocated Cpp class
                                     bind(C, name="InitThermodynamics")
             
            import :: c_int, c_double, c_ptr       
            implicit none   

            ! Argument list
            real(c_double), intent(in), value :: time_step
            integer(c_int), intent(in), value :: num_ice_cells
            real(c_double), intent(in), value :: min_ice_thick
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            real(c_double), intent(in), dimension(*) :: init_base_temp
            real(c_double), intent(in), dimension(*) :: init_surf_temp
            real(c_double), intent(in), dimension(*) :: init_ice_thick
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

    ! update atmosphere flux
    interface 
        subroutine UpdateAtmFlux_(obj, &              ! pointer to allocated Cpp class
                                  atm_flux_values, &  ! 2D-array of total atm flux values (W m-2)
                                  min_lon_ind, &      ! minimal longitude index
                                  max_lon_ind, &      ! maximal longitude index
                                  min_lat_ind, &      ! minimal latitude index
                                  max_lat_ind) &      ! maximal latitude index
                                  bind(C, name="UpdateAtmFlux")
            
            import:: c_int, c_double, c_ptr
            implicit none

            ! Argument list
            type(c_ptr), intent(in), value :: obj
            real(c_double), intent(in), dimension(*) :: atm_flux_values
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        end subroutine 
    end interface

    ! update short-wave radiation flux
    interface
        subroutine UpdateSwRadiation_(obj, &           ! pointer to allocated Cpp class
                                      sw_values, &     ! 2D-array of short-wave radiation flux (W m-2)
                                      min_lon_ind, &   ! minimal longitude index
                                      max_lon_ind, &   ! maximal longitude index
                                      min_lat_ind, &   ! minimal latitude index
                                      max_lat_ind) &   ! maximal latitude index
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

    ! update latent heat flux
    interface
        subroutine UpdateLatentHeatFlux_(obj, &           ! pointer to allocated Cpp class
                                         lh_values, &     ! 2D-array of latent heat flux (W m-2)
                                         min_lon_ind, &   ! minimal longitude index
                                         max_lon_ind, &   ! maximal longitude index
                                         min_lat_ind, &   ! minimal latitude index
                                         max_lat_ind) &   ! maximal latitude index
                                         bind(C, name="UpdateLatentHeatFlux")
            
            import:: c_int, c_double, c_ptr
            implicit none

            ! Argument list
            type(c_ptr), intent(in), value :: obj
            real(c_double), intent(in), dimension(*) :: lh_values
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            
        end subroutine
    end interface

    ! one-step evaluation
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

    ! return ice surface temperature
    interface
        subroutine GetIceSurfaceTemperature_(obj, &           ! pointer to allocated Cpp class
                                             array, &         ! 2D-array of ice surface temperature (deg Cel) - output
                                             min_lon_ind, &   ! minimal longitude index
                                             max_lon_ind, &   ! maximal longitude index
                                             min_lat_ind, &   ! minimal latitude index
                                             max_lat_ind) &   ! maximal latitude index
                                             bind(C, name="GetIceSurfaceTemperature")
            
            import:: c_int, c_double, c_ptr
            implicit none

            ! Argument list
            type(c_ptr), intent(in), value :: obj
            real(c_double), intent(out), dimension(*) :: array
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            
        end subroutine
    end interface

    ! return ice thickness 
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

    ! return ice thickness 
    interface
        subroutine GetSurfaceConductiveFlux_(obj, &           ! pointer to allocated Cpp class
                                             array, &         ! 2D-array of ice surface conductive flux (W m-2) - output
                                             min_lon_ind, &   ! minimal longitude index
                                             max_lon_ind, &   ! maximal longitude index
                                             min_lat_ind, &   ! minimal latitude index
                                             max_lat_ind) &   ! maximal latitude index
                                             bind(C, name="GetSurfaceConductiveFlux")
            
            import:: c_int, c_double, c_ptr
            implicit none

            ! Argument list
            type(c_ptr), intent(in), value :: obj
            real(c_double), intent(out), dimension(*) :: array
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            
        end subroutine
    end interface

    ! store state
    interface
        subroutine StoreState_(obj, &           ! pointer to allocated Cpp class
                               min_lon_ind, &   ! minimal longitude index
                               max_lon_ind, &   ! maximal longitude index
                               min_lat_ind, &   ! minimal latitude index
                               max_lat_ind) &   ! maximal latitude index
                               bind(C, name="StoreState")
            
            import:: c_int, c_ptr
            implicit none

            ! Argument list
            type(c_ptr), intent(in), value :: obj
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            
        end subroutine
    end interface

    ! restore state
    interface
        subroutine RestoreState_(obj, &           ! pointer to allocated Cpp class
                                 min_lon_ind, &   ! minimal longitude index
                                 max_lon_ind, &   ! maximal longitude index
                                 min_lat_ind, &   ! minimal latitude index
                                 max_lat_ind) &   ! maximal latitude index
                                 bind(C, name="RestoreState")
            
            import:: c_int, c_ptr
            implicit none

            ! Argument list
            type(c_ptr), intent(in), value :: obj
            integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
            integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
            
        end subroutine
    end interface


    
    ! Hold Cpp class
    type(c_ptr), save :: obj = c_null_ptr
 
    public :: InitThermodynamics
    public :: FinalizeThermodynamics
    public :: UpdateAtmFlux
    public :: UpdateSwRadiation
    public :: UpdateLatentHeatFlux
    public :: Evaluate
    public :: GetIceSurfaceTemperature
    public :: GetIceThickness
    public :: GetSurfaceConductiveFlux
    public :: StoreState
    public :: RestoreState
 
 contains
 
    ! initialization of thermodynamics solver 
    subroutine InitThermodynamics(time_step, &             ! time step (seconds)
                                  num_ice_cells, &         ! number of ice cells in vertical mesh
                                  min_ice_thick, &         ! minimal ice thickness (meters)
                                  min_lon_ind, &           ! minimal longitude index (inclusively)
                                  max_lon_ind, &           ! maximal longitude index (inclusively)
                                  min_lat_ind, &           ! minimal latitude index (inclusively)
                                  max_lat_ind, &           ! maximal latitude index (inclusively)
                                  init_base_temp, &        ! 2D-array of initial base temp (deg Cel) 
                                  init_surf_temp, &        ! 2D-array of initial surface temp (deg Cel)
                                  init_ice_thick)          ! 2D-array of initial ice thickness (meters)

        implicit none   
        
        ! Argument list
        real(c_double), intent(in), value :: time_step
        integer(c_int), intent(in), value :: num_ice_cells
        real(c_double), intent(in), value :: min_ice_thick
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind
        real(c_double), intent(in), dimension(*) :: init_base_temp
        real(c_double), intent(in), dimension(*) :: init_surf_temp
        real(c_double), intent(in), dimension(*) :: init_ice_thick
        
        ! Body
        if (c_associated(obj)) then
            call FinalizeThermodynamics_(obj)
        end if
   
        obj = InitThermodynamics_(time_step, &             
                                  num_ice_cells, &         
                                  min_ice_thick, &         
                                  min_lon_ind, &           
                                  max_lon_ind, &           
                                  min_lat_ind, &           
                                  max_lat_ind, &           
                                  init_base_temp, &        
                                  init_surf_temp, &        
                                  init_ice_thick)
    
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

    ! update atmosphere flux
    subroutine UpdateAtmFlux(atm_flux_values, &  ! 2D-array of total atm flux values (W m-2)
                             min_lon_ind, &      ! minimal longitude index (inclusively)
                             max_lon_ind, &      ! maximal longitude index (inclusively)
                             min_lat_ind, &      ! minimal latitude index (inclusively)
                             max_lat_ind)        ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(in), dimension(*) :: atm_flux_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateAtmFlux_(obj, &
                                atm_flux_values, &  
                                min_lon_ind, &      
                                max_lon_ind, &      
                                min_lat_ind, &      
                                max_lat_ind)
        end if

    end subroutine 

    ! update short-wave radiation flux
    subroutine UpdateSwRadiation(sw_values, &     ! 2D-array of short-wave radiation flux (W m-2)
                                 min_lon_ind, &   ! minimal longitude index (inclusively)
                                 max_lon_ind, &   ! maximal longitude index (inclusively)
                                 min_lat_ind, &   ! minimal latitude index (inclusively)
                                 max_lat_ind)     ! maximal latitude index (inclusively)
        
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

    ! update latent heat flux
    subroutine UpdateLatentHeatFlux(lh_values, &     ! 2D-array of latent heat flux (W m-2)
                                    min_lon_ind, &   ! minimal longitude index (inclusively)
                                    max_lon_ind, &   ! maximal longitude index (inclusively)
                                    min_lat_ind, &   ! minimal latitude index (inclusively)
                                    max_lat_ind)     ! maximal latitude index (inclusively)
        
        implicit none
        
        ! Argument list
        real(c_double), intent(in), dimension(*) :: lh_values
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call UpdateLatentHeatFlux_(obj, &
                                       lh_values, &     
                                       min_lon_ind, &   
                                       max_lon_ind, &   
                                       min_lat_ind, &   
                                       max_lat_ind)     
        end if
        
    end subroutine

    ! one-step evaluation
    subroutine Evaluate(min_lon_ind, &   ! minimal longitude index (inclusively)
                        max_lon_ind, &   ! maximal longitude index (inclusively)
                        min_lat_ind, &   ! minimal latitude index (inclusively)
                        max_lat_ind)     ! maximal latitude index (inclusively)
        
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

    ! return ice surface temperature
    subroutine GetIceSurfaceTemperature(array, &         ! 2D-array of ice surface temperature (deg Cel) - output
                                        min_lon_ind, &   ! minimal longitude index (inclusively)
                                        max_lon_ind, &   ! maximal longitude index (inclusively)
                                        min_lat_ind, &   ! minimal latitude index (inclusively)
                                        max_lat_ind)     ! maximal latitude index (inclusively)
        
        implicit none

        ! Argument list
        real(c_double), intent(out), dimension(*) :: array
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call GetIceSurfaceTemperature_(obj, &
                                           array, &         
                                           min_lon_ind, &   
                                           max_lon_ind, &   
                                           min_lat_ind, &  
                                           max_lat_ind)
        end if
        
    end subroutine

    ! return ice thickness 
    subroutine GetIceThickness(array, &         ! 2D-array of ice thickness (m) - output
                               min_lon_ind, &   ! minimal longitude index (inclusively)
                               max_lon_ind, &   ! maximal longitude index (inclusively)
                               min_lat_ind, &   ! minimal latitude index (inclusively)
                               max_lat_ind)     ! maximal latitude index (inclusively)
            
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

    ! return surface conductive heat flux 
    subroutine GetSurfaceConductiveFlux(array, &         ! 2D-array of ice surface conductive flux (W m-2) - output
                                        min_lon_ind, &   ! minimal longitude index (inclusively)
                                        max_lon_ind, &   ! maximal longitude index (inclusively)
                                        min_lat_ind, &   ! minimal latitude index (inclusively)
                                        max_lat_ind)     ! maximal latitude index (inclusively)
            
        implicit none

        ! Argument list
        real(c_double), intent(out), dimension(*) :: array
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call GetSurfaceConductiveFlux_(obj, &
                                           array, &         
                                           min_lon_ind, &   
                                           max_lon_ind, &   
                                           min_lat_ind, &   
                                           max_lat_ind)
        end if        
    end subroutine

    ! store state
    subroutine StoreState(min_lon_ind, &   ! minimal longitude index (inclusively)
                          max_lon_ind, &   ! maximal longitude index (inclusively)
                          min_lat_ind, &   ! minimal latitude index (inclusively)
                          max_lat_ind)     ! maximal latitude index (inclusively)

        implicit none

        ! Argument list
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call StoreState_(obj, &
                             min_lon_ind, &   
                             max_lon_ind, &   
                             min_lat_ind, &   
                             max_lat_ind)
        end if        
    end subroutine

    ! restore state
    subroutine RestoreState(min_lon_ind, &   ! minimal longitude index (inclusively)
                            max_lon_ind, &   ! maximal longitude index (inclusively)
                            min_lat_ind, &   ! minimal latitude index (inclusively)
                            max_lat_ind)     ! maximal latitude index (inclusively)

        implicit none

        ! Argument list
        integer(c_int), intent(in), value :: min_lon_ind, max_lon_ind
        integer(c_int), intent(in), value :: min_lat_ind, max_lat_ind

        ! Body
        if (c_associated(obj)) then
            call RestoreState_(obj, &
                               min_lon_ind, &   
                               max_lon_ind, &   
                               min_lat_ind, &   
                               max_lat_ind)
        end if        
end subroutine

end module itslav