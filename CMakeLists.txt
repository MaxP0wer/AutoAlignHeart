# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.8)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)


# This project is designed to be built outside the Insight source tree.
project(autoAlignHeart)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(autoAlignHeart autoAlignHeart.cxx )

target_link_libraries(autoAlignHeart ${ITK_LIBRARIES})
