find_package(Python 3.11 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)

set(PYTHON_BINDINGS constants_bind.cpp radiation_bind.cpp)
nanobind_add_module(CRadpy_constants constants_bind.cpp)
nanobind_add_module(CRadpy ${PYTHON_BINDINGS})

target_link_libraries(CRadpy PRIVATE RadiationLib Eigen3::Eigen)

