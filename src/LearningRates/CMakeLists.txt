add_library (LearningRates Constant.cpp Wolfe.cpp Adam.cpp)

if(UNIX)
    message("You are using a UNIX system. Will not check for MKL")
else()
    find_package(MKL REQUIRED)
    include_directories(${MKL_INCLUDE_DIR})
endif()

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIR}\\Armadillo\\include)

target_include_directories (LearningRates PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
