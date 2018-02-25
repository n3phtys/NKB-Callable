cmake_minimum_required (VERSION 3.4)
project (NKB-Callable)

add_executable(nas-kernel-benchmark-c ./exec/nas.c)
target_link_libraries(nas-kernel-benchmark-c m pthread)


add_executable(nkb_mxm ./exec/main_mxm.c)
target_link_libraries(nkb_mxm m)

add_executable(nkb_gmtry ./exec/main_gmtry.c)
target_link_libraries(nkb_gmtry m)

add_executable(nkb_emit ./exec/main_emit.c)
target_link_libraries(nkb_emit m)

add_executable(nkb_btrix ./exec/main_btrix.c)
target_link_libraries(nkb_btrix m)

add_executable(nkb_cholsky ./exec/main_cholsky.c)
target_link_libraries(nkb_cholsky m)

add_executable(nkb_vpenta ./exec/main_vpenta.c)
target_link_libraries(nkb_vpenta m)

