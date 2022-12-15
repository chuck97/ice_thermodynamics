include(FortranCInterface)

FortranCInterface_VERIFY(CXX)

set(SOURCE 
    ${SOURCE}
    ${CMAKE_CURRENT_SOURCE_DIR}/itslav.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/itslav.f90
    PARENT_SCOPE
)

set(HEADER
    ${HEADER}
    ${CMAKE_CURRENT_SOURCE_DIR}/itslav.hpp
)