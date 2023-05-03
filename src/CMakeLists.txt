aux_source_directory(. SRCS)

# TODO: add auto include
#include_directories(PRIVATE "/usr/local/include")
include_directories(PUBLIC "/usr/local/include")

# add pir library
add_library(pir STATIC ${SRCS})

target_include_directories(pir PUBLIC "/usr/local/include")
target_link_libraries(pir HEXL::hexl)
