find_package(GSL REQUIRED)
#找到 OpenMP
set(SOURCES 
    integration.cpp 
    InverseCompton.cpp 
    RadData.cpp 
    Radiation.cpp 
    Utils.cpp) 
add_library(RadiationLib SHARED ${SOURCES})
target_link_libraries(RadiationLib PUBLIC Eigen3::Eigen GSL::gsl)
target_include_directories(RadiationLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} /../ include)
add_executable(test test.cxx)
target_link_libraries(test PRIVATE RadiationLib)
