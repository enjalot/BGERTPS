/*
* $Id: MOD_enja.c 28152 2010-5-24 enjalot $
*
* ***** BEGIN GPL LICENSE BLOCK *****
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software  Foundation,
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*
* The Original Code is Copyright (C) 2005 by the Blender Foundation.
* All rights reserved.
*
* Contributor(s): Daniel Dunbar
*                 Ton Roosendaal,
*                 Ben Batt,
*                 Brecht Van Lommel,
*                 Campbell Barton,
*                 Ian Johnson
*                 Myrna Merced Serrano
*
* ***** END GPL LICENSE BLOCK *****
*
*/

#include "BKE_cdderivedmesh.h"
#include "BKE_particle.h"

#include "MOD_modifiertypes.h"
#include "MEM_guardedalloc.h"


static void initData(ModifierData *md) 
{
	RTPSModifierData *rtmd = (RTPSModifierData*) md;

    // why are defaults required? Do not appear to be used
    // at least not if there is no associated UI
	rtmd->system = 2;
    rtmd->max_num = 100;
	                     
    rtmd->sub_intervals= 3;
    rtmd->dt= .001f;

    //sph physics parameters
    rtmd->gravity = -9.8f;
    rtmd->gas_constant = 15.0f;
    rtmd->viscosity = .01f;
    rtmd->velocity_limit = 600.0f;
    rtmd->xsph_factor = .1f;
    //sph simulation parameters
    rtmd->boundary_stiffness = 20000.0f;
    rtmd->boundary_dampening = 256.0f;
 

	// GE: scale of radius used by Andrew for improved rendering
	rtmd->render_radius_scale = 3.;
	rtmd->render_blur_scale = 1.;
	rtmd->render_type = 0;

    rtmd->collision = 0;
    rtmd->glsl = 0;
    rtmd->blending = 0;

    //boids stuff
	rtmd->maxspeed = 1.f;
	rtmd->separationdist = 1.f;
	rtmd->searchradius = 1.f;

	rtmd->color_r=255.0f;
	rtmd->color_g=0.0f;
	rtmd->color_b=0.0f;
	rtmd->color_a=255.0f;

    rtmd->w_sep = 1.f;
    rtmd->w_align = 1.f;
    rtmd->w_coh = 1.f; 
}
/*
static void copyData(ModifierData *md, ModifierData *target)
{
	RTPSModifierData *smd = (RTPSModifierData*) md;
	RTPSModifierData *tsmd = (RTPSModifierData*) target;

	tsmd->system = smd->system;
	tsmd->num = smd->num;
	tsmd->sub_intervals = smd->sub_intervals;
	tsmd->dt = smd->dt;
}
*/

static void deformVerts(ModifierData *md, Object *ob, DerivedMesh *derivedData, float (*vertexCos)[3], int numVerts, int useRenderParams, int isFinalCalc)
{
	//sbObjectStep(md->scene, ob, (float)md->scene->r.cfra, vertexCos, numVerts);
    //printf("deformVerts in MOD_rtps.c\n");
    return;
}

static int dependsOnTime(ModifierData *md)
{
	return 1;
}


ModifierTypeInfo modifierType_RTPS = {
	/* name */              "RTPS",
	/* structName */        "RTPSModifierData",
	/* structSize */        sizeof(RTPSModifierData),
	/* type */              eModifierTypeType_OnlyDeform,
	/* flags */             eModifierTypeFlag_AcceptsMesh,
//                          eModifierTypeFlag_AcceptsCVs
//							| eModifierTypeFlag_RequiresOriginalData
//							| eModifierTypeFlag_Single,

	/* copyData */          NULL,//copyData,
	/* deformVerts */       deformVerts,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     NULL,
	/* applyModifierEM */   NULL,
	/* initData */          initData,
	/* requiredDataMask */  NULL,
	/* freeData */          NULL,
	/* isDisabled */        NULL,
	/* updateDepgraph */    NULL,
	/* dependsOnTime */     NULL,//dependsOnTime,
    /* dependsOnNormals */  NULL,
	/* foreachObjectLink */ NULL,
	/* foreachIDLink */     NULL,
};

