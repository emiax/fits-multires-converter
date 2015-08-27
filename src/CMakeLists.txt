set(OPENSPACE_DIR ~/openspace)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(BOOST_INCLUDE_DIR /usr/include)
set(BOOST_LIBRARY_DIR /usr/lib)

set(CCFITS_INCLUDE_DIR /usr/local/include/CCfits)
set(CCFITS_LIBRARY_DIR /usr/local/lib)

find_package(Boost COMPONENTS system filesystem REQUIRED)

include_directories(${Boost_INCLUDE_DIR})
include_directories(${FitsConverter_SOURCE_DIR}/include)
include_directories(${CCFITS_INCLUDE_DIR})

link_directories(${CCFITS_LIBRARY_DIR})
link_directories(${BOOST_LIBRARY_DIR})



add_executable(FitsConverter
               main.cpp
               fitsconverter.cpp)

target_link_libraries(FitsConverter ${Boost_LIBRARIES} libCCfits.so)






