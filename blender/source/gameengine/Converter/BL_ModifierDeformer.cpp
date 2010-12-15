/**
 * $Id$
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
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2001-2002 by NaN Holding BV.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#if defined(WIN32) && !defined(FREE_WINDOWS)
#pragma warning (disable : 4786)
#endif //WIN32

#include "MEM_guardedalloc.h"
#include "BL_ModifierDeformer.h"
#include "GEN_Map.h"
#include "STR_HashedString.h"
#include "RAS_IPolygonMaterial.h"
#include "RAS_MeshObject.h"
#include "PHY_IGraphicController.h"

//#include "BL_ArmatureController.h"
#include "DNA_armature_types.h"
#include "DNA_action_types.h"
#include "DNA_key_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_ipo_types.h"
#include "DNA_curve_types.h"
#include "DNA_modifier_types.h"
#include "DNA_scene_types.h"
#include "BKE_armature.h"
#include "BKE_action.h"
#include "BKE_key.h"
#include "BKE_ipo.h"
#include "MT_Point3.h"

#include "RTPS.h"
//These KX might not belong here...
//this should be addressed when moving away from a modifier
#include "KX_Scene.h"   //for getting the camera
//#include "BKE_scene.h"
#include "KX_Camera.h"
#include "KX_GameObject.h"
#include "KX_PythonInit.h"
#include "SG_BBox.h"

extern "C"{
    #include "BKE_cdderivedmesh.h" //added by IJ for RTPS
	#include "BKE_customdata.h"
	#include "BKE_DerivedMesh.h"
	#include "BKE_lattice.h"
	#include "BKE_modifier.h"
}
 #include "BKE_utildefines.h"

#include "BLI_blenlib.h"
#include "BLI_math.h"

#define __NLA_DEFNORMALS
//#undef __NLA_DEFNORMALS


BL_ModifierDeformer::~BL_ModifierDeformer()
{
	if (m_dm) {
		// deformedOnly is used as a user counter
		if (--m_dm->deformedOnly == 0) {
			m_dm->needsFree = 1;
			m_dm->release(m_dm);
		}
	}
};

RAS_Deformer *BL_ModifierDeformer::GetReplica()
{
	BL_ModifierDeformer *result;

	result = new BL_ModifierDeformer(*this);
	result->ProcessReplica();
	return result;
}

void BL_ModifierDeformer::ProcessReplica()
{
	/* Note! - This is not inherited from PyObjectPlus */
	BL_ShapeDeformer::ProcessReplica();
	if (m_dm)
		// by default try to reuse mesh, deformedOnly is used as a user count
		m_dm->deformedOnly++;
	// this will force an update and if the mesh cannot be reused, a new one will be created
	m_lastModifierUpdate = -1;
}

bool BL_ModifierDeformer::HasCompatibleDeformer(Object *ob)
{
	if (!ob->modifiers.first)
		return false;
	// soft body cannot use mesh modifiers
	if ((ob->gameflag & OB_SOFT_BODY) != 0)
		return false;
	ModifierData* md;
	for (md = (ModifierData*)ob->modifiers.first; md; md = (ModifierData*)md->next) {
		if (modifier_dependsOnTime(md))
			continue;
		if (!(md->mode & eModifierMode_Realtime))
			continue;
		/* armature modifier are handled by SkinDeformer, not ModifierDeformer */
		if (md->type == eModifierType_Armature )
			continue;
		return true;
	}
	return false;
}

bool BL_ModifierDeformer::HasRTPSDeformer(Object *ob)
{
    	printf("IJ: in HasRTPSModifier\n");
   	ModifierData* md;

	for (md = (ModifierData*)ob->modifiers.first; md; md = (ModifierData*)md->next) {
        	if(md->type & eModifierType_RTPS){
            		printf("IJ: found RTPS modifier!\n");
            		return true;
        	}
	}

	return false;
}

bool BL_ModifierDeformer::HasArmatureDeformer(Object *ob)
{
	if (!ob->modifiers.first)
		return false;

	ModifierData* md = (ModifierData*)ob->modifiers.first;
	if(md->type == eModifierType_Armature )
		return true;

	return false;
}

bool BL_ModifierDeformer::Update(void)
{
	bool bShapeUpdate = BL_ShapeDeformer::Update();

	if (bShapeUpdate || m_lastModifierUpdate != m_gameobj->GetLastFrame()) {
		// static derived mesh are not updated
		if (m_dm == NULL || m_bDynamic) {
			/* execute the modifiers */
			Object* blendobj = m_gameobj->GetBlendObject();
			/* hack: the modifiers require that the mesh is attached to the object
			   It may not be the case here because of replace mesh actuator */
			Mesh *oldmesh = (Mesh*)blendobj->data;
			blendobj->data = m_bmesh;
			/* execute the modifiers */		
			DerivedMesh *dm = mesh_create_derived_no_virtual(m_scene, blendobj, m_transverts, CD_MASK_MESH);
			/* restore object data */
			blendobj->data = oldmesh;
			/* free the current derived mesh and replace, (dm should never be NULL) */
			if (m_dm != NULL) {
				// HACK! use deformedOnly as a user counter
				if (--m_dm->deformedOnly == 0) {
					m_dm->needsFree = 1;
					m_dm->release(m_dm);
				}
			}
			m_dm = dm;
			// get rid of temporary data
			m_dm->needsFree = 0;
			m_dm->release(m_dm);
			// HACK! use deformedOnly as a user counter
			m_dm->deformedOnly = 1;
			/* update the graphic controller */
			PHY_IGraphicController *ctrl = m_gameobj->GetGraphicController();
			if (ctrl) {
				float min_r[3], max_r[3];
				INIT_MINMAX(min_r, max_r);
				m_dm->getMinMax(m_dm, min_r, max_r);
				ctrl->setLocalAabb(min_r, max_r);
			}
		}
		m_lastModifierUpdate=m_gameobj->GetLastFrame();
		bShapeUpdate = true;
	}
    
    if(m_bIsRTPS)
    {
        int nmat = m_pMeshObject->NumMaterials();

        for(int imat=0; imat<nmat; imat++)
        {
            RAS_MeshMaterial *mmat = m_pMeshObject->GetMeshMaterial(imat);
            RAS_MeshSlot **slot = mmat->m_slots[(void*)m_gameobj];   
    	    if(!slot || !*slot)
                continue;
            //iterate through the slots until we find one with a particle system
            if((*slot)->m_bRTPS && (*slot)->m_pRTPS)
            {
                rtps::RTPS* rtps = (*slot)->m_pRTPS;

                KX_Scene* kxs = KX_GetActiveScene();
                //handle camera scale for proper render size
                KX_Camera* kxc = kxs->GetActiveCamera();
                if(!kxc)
                {
                    //(*slot)->m_pEnjaParticles->point_scale = 1.0;
                }
                else
                {
                    float scale = kxc->GetScale();
                    //printf("SCALE: %g\n", scale);
                    //(*slot)->m_pEnjaParticles->point_scale = scale;
                }


                //deal with emitters if SPH
                //maybe this should be done with a new modifier?
                if(rtps->settings.system == rtps::RTPSettings::SPH)
                {
                    //loop through the objects looking for objects with collider property
                    CListValue* oblist = kxs->GetObjectList();
                    int num_objects = oblist->GetCount();
                    for(int iob = 0; iob < num_objects; iob++)
                    {
                        KX_GameObject* gobj = (KX_GameObject*)oblist->GetValue(iob);
                        //print out some info about the object
                        //int mesh_count = gobj->GetMeshCount();
                        STR_String name = gobj->GetName();


                        //Check if object is a collider
                        bool collider = false;
                        CBoolValue* boolprop = (CBoolValue*)gobj->GetProperty("collider");
                        if(boolprop)
                        {
                            printf("obj: %s, collider: %d\n", name.Ptr(), boolprop->GetBool());
                            collider = boolprop->GetBool();
                        }


                        using namespace rtps;
                        //Check if object is an emitter
                        //for now we are just doing boxes 
                        CIntValue* intprop = (CIntValue*)gobj->GetProperty("num");
                        if(intprop)
                        {
                            //get the number of particles in this emitter
                            int num = (int)intprop->GetInt();
                            int nn = 0; //the number of particles to emit
                            if (num == 0) { continue;} //out of particles
                            boolprop = (CBoolValue*)gobj->GetProperty("hose");
                            if(boolprop && boolprop->GetBool())
                            {
                                //particles per frame
                                intprop = (CIntValue*)gobj->GetProperty("ppf");
                                int ppf = (int)intprop->GetInt();
                                intprop = (CIntValue*)gobj->GetProperty("freq");
                                int freq = (int)intprop->GetInt();
                                intprop = (CIntValue*)gobj->GetProperty("ioff");
                                int ioff = (int)intprop->GetInt();
                                if(ioff == 0)
                                {
                                    num -= ppf;
                                    if (num < 0)
                                    {
                                        //inject the remaining particles
                                        nn = num + ppf;
                                    }
                                    else
                                    {
                                        //inject the number of particles for this frame
                                        nn = ppf;
                                    }
                                }
                                ioff++;
                                if(ioff >= freq) ioff = 0;
                                printf("ioff: %d\n", ioff);
				                CIntValue * newioffprop = new CIntValue(ioff);            

                                gobj->SetProperty("ioff", newioffprop);
                            }
                            else
                            {
                                //not a hose
                                nn = num;
                            }


                            if(nn <= 0)
                            {
                                continue;
                            }
                            printf("obj: %s, nn: %d\n", name.Ptr(), nn);
                            MT_Point3 bbpts[8];
                            gobj->GetSGNode()->getAABBox(bbpts);
                            float4 min = float4(bbpts[0].x(), bbpts[0].y(), bbpts[0].z(), 0);
                            float4 max = float4(bbpts[7].x(), bbpts[7].y(), bbpts[7].z(), 0);
                            rtps->system->addBox(nn, min, max, false);

                        }
                       
                        //printf("obj: %s mesh_count: %d collider: %d\n", name.Ptr(), mesh_count, collider);
     
                    }
                }
                


                (*slot)->m_pRTPS->update();
	        }
        }//for loop over materials
    }//if(m_bIsRTPS)


	return bShapeUpdate;
}

bool BL_ModifierDeformer::Apply(RAS_IPolyMaterial *mat)
{
	if (!Update())
		return false;

	// drawing is based on derived mesh, must set it in the mesh slots
	int nmat = m_pMeshObject->NumMaterials();
	for (int imat=0; imat<nmat; imat++) {
		RAS_MeshMaterial *mmat = m_pMeshObject->GetMeshMaterial(imat);
		RAS_MeshSlot **slot = mmat->m_slots[(void*)m_gameobj];
		if(!slot || !*slot)
			continue;
		(*slot)->m_pDerivedMesh = m_dm;
        
        if(m_bIsRTPS)
        {
            printf("IJ: set the mesh slot m_bRTPS to true (in ModifierDeformer::Apply()\n");
            (*slot)->m_bRTPS = true;
            //initialize the particle system
            printf("IJ: initialize particle system\n");
            
            ModifierData* md;
            for (md = (ModifierData*)m_objMesh->modifiers.first; md; md = (ModifierData*)md->next) 
            {
                if(md->type & eModifierType_RTPS)
                {
                    RTPSModifierData* rtmd = (RTPSModifierData*)md;
                    //printf("SysType::Simple %d\n", rtps::RTPSettings::Simple);
                    //printf("rtmd->system: %d\n", rtmd->system);

                    rtps::RTPSettings::SysType sys = (rtps::RTPSettings::SysType)rtmd->system;
                    //printf("sys: %d\n", sys);
        
                    //get the bounding box of the object for use as domain bounds
                    KX_GameObject* gobj = (KX_GameObject*)m_gameobj;
                    MT_Point3 bbpts[8];
                    gobj->GetSGNode()->getAABBox(bbpts);
                    MT_Point3 min = bbpts[0];
                    MT_Point3 max = bbpts[7];
                    using namespace rtps;
                    rtps::Domain grid(float4(min.x(), min.y(), min.z(), 0), float4(max.x(), max.y(), max.z(), 0));


                    if (sys == rtps::RTPSettings::SimpleFlock)
                    {
                        float color[3] = {rtmd->color_r, rtmd->color_g, rtmd->color_b};
                        rtps::RTPSettings settings(rtmd->num, rtmd->maxspeed, rtmd->separationdist, rtmd->perceptionrange, color);
                        (*slot)->m_pRTPS = new rtps::RTPS(settings);
                    }
                    else if (sys == rtps::RTPSettings::SPH) 
                    {
                        rtps::RTPSettings settings(sys, rtmd->num, rtmd->dt, grid);
                        (*slot)->m_pRTPS = new rtps::RTPS(settings);

                        rtps::RTPS* rtps = (*slot)->m_pRTPS;
                        KX_Scene* kxs = KX_GetActiveScene();

                        //Loop through objects to look for emitters
                        //this should probably use modifiers/particle systems
                        //instead of game props but lets just get something going
                        CListValue* oblist = kxs->GetObjectList();
                        int num_objects = oblist->GetCount();
                        for(int iob = 0; iob < num_objects; iob++)
                        {
                            KX_GameObject* gobj = (KX_GameObject*)oblist->GetValue(iob);
                            STR_String name = gobj->GetName();

                            //Check if object is an emitter
                            //for now we are just doing boxes 
                            CIntValue* intprop = (CIntValue*)gobj->GetProperty("emitter");
                            if(intprop)
                            {
                                int nn = (int)intprop->GetInt();
                                printf("obj: %s, emitter: %d\n", name.Ptr(), nn);
                                MT_Point3 bbpts[8];
                                gobj->GetSGNode()->getAABBox(bbpts);
                                float4 min = float4(bbpts[0].x(), bbpts[0].y(), bbpts[0].z(), 0);
                                float4 max = float4(bbpts[7].x(), bbpts[7].y(), bbpts[7].z(), 0);
                                rtps->system->addBox(nn, min, max, false);

                            }
                            
                            //printf("obj: %s mesh_count: %d collider: %d\n", name.Ptr(), mesh_count, collider);
         
                        }

                            

                    }
                    else 
                    {
                        rtps::RTPSettings settings(sys, rtmd->num, rtmd->dt, grid);
                        (*slot)->m_pRTPS = new rtps::RTPS(settings);
                    }

                }
            }
        }
	}
	return true;
}
