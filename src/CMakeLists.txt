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
	       fitsconverter.cpp
	       lz4compressor.cpp)

target_link_libraries(FitsConverter Ghoul)
get_property(GHOUL_INCLUDE_DIR TARGET Ghoul PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(FitsConverter PUBLIC ${GHOUL_INCLUDE_DIR})
get_property(GHOUL_DEFINITIONS TARGET Ghoul PROPERTY INTERFACE_COMPILE_DEFINITIONS)
target_compile_definitions(FitsConverter PUBLIC ${GHOUL_DEFINITIONS})
set_property(TARGET Lua PROPERTY FOLDER "External")
set_property(TARGET lz4 PROPERTY FOLDER "External")
set_property(TARGET tinyobjloader PROPERTY FOLDER "External")

target_link_libraries(FitsConverter BasicJ2k)
target_link_libraries(FitsConverter ${Boost_LIBRARIES} libCCfits.so)





