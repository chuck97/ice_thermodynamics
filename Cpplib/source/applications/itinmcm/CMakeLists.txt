include(FortranCInterface)

FortranCInterface_VERIFY(CXX)

set(SOURCE 
    ${SOURCE}
    ${CMAKE_CURRENT_SOURCE_DIR}/itinmcm0d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/itinmcm0d.f90
    PARENT_SCOPE
)

set(HEADER
    ${HEADER}
    ${CMAKE_CURRENT_SOURCE_DIR}/itinmcm.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/itinmcm0d.hpp
)