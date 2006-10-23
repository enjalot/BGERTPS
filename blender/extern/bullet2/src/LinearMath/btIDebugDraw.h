/*
Copyright (c) 2005 Gino van den Bergen / Erwin Coumans http://continuousphysics.com

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/


#ifndef IDEBUG_DRAW__H
#define IDEBUG_DRAW__H

#include "LinearMath/btVector3.h"


class	btIDebugDraw
{
	public:

	enum	DebugDrawModes
	{
		DBG_NoDebug=0,
		DBG_DrawWireframe = 1,
		DBG_DrawAabb=2,
		DBG_DrawFeaturesText=4,
		DBG_DrawContactPoints=8,
		DBG_NoDeactivation=16,
		DBG_NoHelpText = 32,
		DBG_DrawText=64,
		DBG_ProfileTimings = 128,
		DBG_EnableSatComparison = 256,
		DBG_DisableBulletLCP = 512,
		DBG_EnableCCD = 1024,
		DBG_MAX_DEBUG_DRAW_MODE
	};

	virtual ~btIDebugDraw() {};

	virtual void	drawLine(const btVector3& from,const btVector3& to,const btVector3& color)=0;

	virtual void	drawContactPoint(const btVector3& PointOnB,const btVector3& normalOnB,float distance,int lifeTime,const btVector3& color)=0;

	virtual void	setDebugMode(int debugMode) =0;
	
	virtual int		getDebugMode() const = 0;


};

#endif //IDEBUG_DRAW__H

