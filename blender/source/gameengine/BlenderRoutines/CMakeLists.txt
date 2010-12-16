
set(INC
	.
	../../../source/kernel/gen_system
	../../../intern/string
	../../../intern/guardedalloc
	../../../intern/audaspace/intern
	../../../source/gameengine/Rasterizer/RAS_OpenGLRasterizer 
	../../../source/gameengine/Converter
	../../../source/blender/imbuf
	../../../intern/moto/include
	../../../source/gameengine/Ketsji 
	../../../source/blender/blenlib
	../../../source/blender/blenkernel 
	../../../source/blender/blenfont
	../../../source/blender/editors/include
	../../../source/blender/windowmanager
	../../../source/blender
	../../../source/blender/makesdna
	../../../source/blender/makesrna 
	../../../source/gameengine/Rasterizer 
	../../../source/gameengine/GameLogic
	../../../source/gameengine/Expressions 
	../../../source/gameengine/Network
	../../../source/gameengine/SceneGraph 
	../../../source/gameengine/Physics/common
	../../../source/gameengine/Physics/Bullet
	../../../source/gameengine/Network/LoopBackNetwork
	../../../source/blender/blenloader
	../../../source/blender/gpu
	../../../extern/bullet2/src
	../../../extern/glew/include
)

set(SRC
	BL_KetsjiEmbedStart.cpp
	KX_BlenderCanvas.cpp
	KX_BlenderGL.cpp
	KX_BlenderInputDevice.cpp
	KX_BlenderKeyboardDevice.cpp
	KX_BlenderMouseDevice.cpp
	KX_BlenderRenderTools.cpp
	KX_BlenderSystem.cpp

	KX_BlenderCanvas.h
	KX_BlenderGL.h
	KX_BlenderInputDevice.h
	KX_BlenderKeyboardDevice.h
	KX_BlenderMouseDevice.h
	KX_BlenderRenderTools.h
	KX_BlenderSystem.h
)

<<<<<<< HEAD
IF(WITH_OPENCL)
	LIST(APPEND INC ${OPENCL_INCLUDE_DIR})  
	LIST(APPEND INC ${RTPS_INCLUDE_DIR})  
ENDIF(WITH_OPENCL)

IF(WITH_FFMPEG)
	ADD_DEFINITIONS(-DWITH_FFMPEG)
ENDIF(WITH_FFMPEG)
=======
add_definitions(-DGLEW_STATIC)
>>>>>>> 9b28302d9853b8fee3fc46f398f32a5d8ef01113

if(WITH_CODEC_FFMPEG)
	add_definitions(-DWITH_FFMPEG)
endif()

blenderlib(ge_blen_routines "${SRC}" "${INC}")
