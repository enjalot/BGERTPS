
#SET(RTPS_INCLUDE_DIR "/panfs/panasas1/users/idj03/research/rtps/rtps")
SET(RTPS_INCLUDE_DIR $ENV{RTPS_DIR}/rtpslib)
#SET(RTPS_LIB_DIR $ENV{RTPS_DIR}/build/rtpslib)

#hacks, we should either do a FIND_LIBRARY properly
#ideally we move RTPS inside blender eventually
IF(APPLE)
    SET(RTPS_LIB_DIR $ENV{RTPS_DIR}/build/rtpslib/librtps.dylib)
ENDIF(APPLE)
IF(UNIX AND NOT APPLE)
    SET(RTPS_LIB_DIR $ENV{RTPS_DIR}/build/rtpslib/librtps.so)
ENDIF(UNIX AND NOT APPLE)
