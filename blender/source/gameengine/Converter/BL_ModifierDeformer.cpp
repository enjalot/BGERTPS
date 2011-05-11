/*
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

/** \file gameengine/Converter/BL_ModifierDeformer.cpp
 *  \ingroup bgeconv
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
#include "BLI_utildefines.h"
#include "BKE_armature.h"
#include "BKE_action.h"
#include "BKE_key.h"
#include "BKE_ipo.h"
#include "MT_Point3.h"

// ==== begin: Included by enjalot
#include "RTPS.h"
//#include "timege.h"
//These KX might not belong here...
//this should be addressed when moving away from a modifier
#include "KX_Scene.h"   //for getting the camera
//#include "BKE_scene.h"
#include "KX_Camera.h"
#include "KX_GameObject.h"
#include "KX_PythonInit.h"
#include "SG_BBox.h"
#include "FloatValue.h"
#include "BLI_path_util.h"
// ==== end included by enjalot

extern "C"{
    #include "BKE_cdderivedmesh.h" //added by enjalot for RTPS
	#include "BKE_customdata.h"
	#include "BKE_DerivedMesh.h"
	#include "BKE_lattice.h"
	#include "BKE_modifier.h"
}
 

#include "BLI_blenlib.h"
#include "BLI_math.h"

#define __NLA_DEFNORMALS
//#undef __NLA_DEFNORMALS

/* helper function declarations used to get data from
 * blender to the RTPS library
 */

using namespace rtps;
void getTriangle(MFace& face, MT_Matrix3x3& grot, MT_Point3& gp, DerivedMesh* dm, Triangle& tri);
int makeEmitter(int num, KX_GameObject* gobj);
int getTriangles(KX_GameObject* gobj, std::vector<Triangle> &triangles);


BL_ModifierDeformer::~BL_ModifierDeformer()
{
	if (m_dm) {
		// deformedOnly is used as a user counter
		if (--m_dm->deformedOnly == 0) {
			m_dm->needsFree = 1;
			m_dm->release(m_dm);
		}
	}
	//RTPS_NOTE seriously wtf can't i destruct any objects from rtps?
	//printf("deleting settings\n");
	//delete m_RTPS_settings;
	printf("is it them?\n");
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

// return a deformed mesh that supports mapping (with a valid CD_ORIGINDEX layer)
struct DerivedMesh* BL_ModifierDeformer::GetPhysicsMesh()
{
	// we need to compute the deformed mesh taking into account the current
	// shape and skin deformers, we cannot just call mesh_create_derived_physics()
	// because that would use the m_transvers already deformed previously by BL_ModifierDeformer::Update(),
	// so restart from scratch by forcing a full update the shape/skin deformers 
	// (will do nothing if there is no such deformer)
	BL_ShapeDeformer::ForceUpdate();
	BL_ShapeDeformer::Update();
	// now apply the modifiers but without those that don't support mapping
	Object* blendobj = m_gameobj->GetBlendObject();
	/* hack: the modifiers require that the mesh is attached to the object
	   It may not be the case here because of replace mesh actuator */
	Mesh *oldmesh = (Mesh*)blendobj->data;
	blendobj->data = m_bmesh;
	DerivedMesh *dm = mesh_create_derived_physics(m_scene, blendobj, m_transverts, CD_MASK_MESH);
	/* restore object data */
	blendobj->data = oldmesh;
	/* m_transverts is correct here (takes into account deform only modifiers) */
	/* the derived mesh returned by this function must be released by the caller !!! */
	return dm;
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
    
	//----------------------------------------
    if(m_bIsRTPS)
    {
        //timers[TI_UPDATE]->start();
        int nmat = m_pMeshObject->NumMaterials();

        for(int imat=0; imat<nmat; imat++)
        {
            RAS_MeshMaterial *mmat = m_pMeshObject->GetMeshMaterial(imat);
            RAS_MeshSlot **slot = mmat->m_slots[(void*)m_gameobj];   
            
            unsigned char color[4];
            mmat->m_bucket->GetPolyMaterial()->GetMaterialRGBAColor(color);
            float4 col(color[0]/255.f,color[1]/255.f,color[2]/255.f,color[3]/255.f);
            
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
                    /*
                    int left = kxc->GetViewportLeft();
                    int right= kxc->GetViewportRight();
                    int bottom = kxc->GetViewportBottom();
                    int top = kxc->GetViewportTop();
                    printf("left %d right %d bottom %d top %d\n", left, right, bottom, top);
                    */
                    //printf("SCALE: %g\n", scale);
                    //(*slot)->m_pEnjaParticles->point_scale = scale;
                }


                //deal with emitters if SPH
                //maybe this should be done with a new modifier?
                if(rtps->settings->system == rtps::RTPSettings::SPH || rtps->settings->system == rtps::RTPSettings::FLOCK)
                {
                    //loop through the objects looking for objects with collider property
                    CListValue* oblist = kxs->GetObjectList();
                    int num_objects = oblist->GetCount();

                    std::vector<Triangle> triangles;

                    //timers[TI_EMIT]->start();
                    
                    for(int iob = 0; iob < num_objects; iob++)
                    {
                        
                        KX_GameObject* gobj = (KX_GameObject*)oblist->GetValue(iob);
                        //print out some info about the object
                        //int mesh_count = gobj->GetMeshCount();
                        STR_String name = gobj->GetName();

                        //Check if object is a collider
                        bool collider = false;
                        CBoolValue* colliderprop = (CBoolValue*)gobj->GetProperty("collider");
                        if(colliderprop)
                        {
                            //printf("obj: %s, collider: %d\n", name.Ptr(), boolprop->GetBool());
                            collider = colliderprop->GetBool();
                            if(collider)
                            {
                                getTriangles(gobj, triangles);
                                //printf("triangles size: %d\n", triangles.size());
                            }
                        }

                        //float4 col((*slot)->m_RGBAcolor.x(),(*slot)->m_RGBAcolor.y(),(*slot)->m_RGBAcolor.z(),(*slot)->m_RGBAcolor.w());
                        //float4 col;
                        //(*slot)->m_RGBAcolor.getValue(&col);
                        

                        //Check if object is an emitter
                        //for now we are just doing boxes 
                        CIntValue* intprop = (CIntValue*)gobj->GetProperty("num");
                        if(intprop)
                        {  
                             //get the number of particles in this emitter
                            int num = (int)intprop->GetInt();
                            if (num == 0) { continue;} //out of particles


                            //hose
                            CBoolValue* hoseprop = (CBoolValue*)gobj->GetProperty("hose");
                            if(hoseprop)
                            {
                                
                                //number of particles for hose to emit
                                CIntValue* numprop = (CIntValue*)gobj->GetProperty("num");
                                int num = (int)numprop->GetInt();

                                CFloatValue* velprop = (CFloatValue*)gobj->GetProperty("speed");
                                float speed = velprop->GetFloat();

                                CFloatValue* radprop = (CFloatValue*)gobj->GetProperty("radius");
                                float radius = radprop->GetFloat();

                                //get center and direction from world transformations
                                MT_Point3 gp = gobj->NodeGetWorldPosition();
                                MT_Matrix3x3 grot = gobj->NodeGetWorldOrientation();
                                MT_Matrix3x3 lrot = gobj->NodeGetLocalOrientation();
                                //we shoot the hose in the object's Y direction
                                MT_Vector3 dir(0., 1., 0.);
                                dir = grot * dir;
                                float4 center(gp[0], gp[1], gp[2], 1.f); 
                                float4 velocity(dir[0], dir[1], dir[2], 0.);
                                velocity = velocity * speed;

                                CIntValue* indexprop = (CIntValue*)gobj->GetProperty("index");
                                if(hoseprop->GetBool())
                                {
                                    int index = rtps->system->addHose(num, center, velocity, radius, col);
                                    CIntValue* hose_index = new CIntValue(index);
                                    gobj->SetProperty("index", hose_index);
                                    //delete hose_index;

                                    CBoolValue setfalse(false);
                                    hoseprop->SetValue(&setfalse);
                                }
                                else
                                {
                                    int index = (int)indexprop->GetInt();
                                    rtps->system->updateHose(index, center, velocity, radius, col);

                                    CIntValue* refillprop = (CIntValue*)gobj->GetProperty("refill");
                                    if(refillprop)
                                    {
                                        int refill = (int)refillprop->GetInt();
                                        if (refill > 0)
                                        {
                                            rtps->system->refillHose(index, refill);
                                            CIntValue* newrefill = new CIntValue(0);
                                            refillprop->SetValue(newrefill);
                                        }
                                    }

                                }
                            }
                            else
                            {
                                int nn = makeEmitter(num, gobj);
                                if( nn == 0) {continue;}
                                //printf("obj: %s, nn: %d\n", name.Ptr(), nn);
                                MT_Point3 bbpts[8];
                                gobj->GetSGNode()->getAABBox(bbpts);
         
                                float4 min = float4(bbpts[0].x(), bbpts[0].y(), bbpts[0].z(), 0);
                                float4 max = float4(bbpts[7].x(), bbpts[7].y(), bbpts[7].z(), 0);
                                rtps->system->addBox(nn, min, max, false, col);
                            }
                           
                        }//if emitters
                        
                        //printf("obj: %s mesh_count: %d collider: %d\n", name.Ptr(), mesh_count, collider);
     
                    }//for loop of objects

                    //timers[TI_EMIT]->end();


                    if(triangles.size() > 0 && rtps->settings->tri_collision)
                    {
                        //printf("about to load triangles\n");
                        //printf("triangles size: %d\n", triangles.size());
                        rtps->system->loadTriangles(triangles);
						//printf("after load triangles\n");
                    }

                    //timers[TI_RTPSUP]->start();
                    (*slot)->m_pRTPS->update();
                    //timers[TI_RTPSUP]->end();

                }
                

   
	        }
        }//for loop over materials
        //timers[TI_UPDATE]->end();

    }//if(m_bIsRTPS)
	//----------------------------------------


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
            printf("RTPS: set the mesh slot m_bRTPS to true (in ModifierDeformer::Apply()\n");
            (*slot)->m_bRTPS = true;
            //initialize the particle system
            printf("RTPS: initialize particle system\n");
            //printf("blender user path system: %s\n", BLI_get_folder_version(BLENDER_RESOURCE_PATH_SYSTEM, BLENDER_VERSION, FALSE));

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
                    MT_Point3 dmin = bbpts[0];
                    MT_Point3 dmax = bbpts[7];
                    //apply the object's transforms so we use global coordinates
                    MT_Matrix3x3 grot = gobj->NodeGetWorldOrientation();
                    MT_Point3 gp = gobj->NodeGetWorldPosition();
                    dmin = grot*dmin + gp;
                    dmax = grot*dmax + gp;

                    printf("RTPS: declaring the path\n");
                    //printf("blender rtps dir: %s\n", BLI_get_folder(BLENDER_RTPS, NULL));
                    std::string rtps_path( BLI_get_folder(BLENDER_RTPS, NULL) );

                    using namespace rtps;
                    printf("RTPS: created the domain\n");
                    rtps::Domain* grid = new Domain(float4(dmin.x(), dmin.y(), dmin.z(), 0), float4(dmax.x(), dmax.y(), dmax.z(), 0));

                    if (sys == rtps::RTPSettings::SPH) 
                    {
                        printf("RTPS: inside SPH system\n");
                        //rtps::RTPSettings settings(sys, grid);
						m_RTPS_settings = new rtps::RTPSettings(sys, rtmd->max_num, rtmd->dt, grid, rtmd->collision);
                        
                        //this path gives the location of the opencl and shader source files
                        printf("rtps path: %s\n", rtps_path.c_str());
                        m_RTPS_settings->SetSetting("rtps_path", rtps_path);


						m_RTPS_settings->setRadiusScale(rtmd->render_radius_scale);
						m_RTPS_settings->setRenderType((rtps::RTPSettings::RenderType)rtmd->render_type);
						m_RTPS_settings->setBlurScale(rtmd->render_blur_scale);
						m_RTPS_settings->setUseGLSL(rtmd->glsl);
						m_RTPS_settings->setUseAlphaBlending(rtmd->blending);


                        m_RTPS_settings->SetSetting("render_texture", "firejet_blast.png");
                        m_RTPS_settings->SetSetting("render_frag_shader", "sprite_tex_frag.glsl");
                        m_RTPS_settings->SetSetting("render_use_alpha", true);
                        //settings.SetSetting("render_use_alpha", false);
                        m_RTPS_settings->SetSetting("render_alpha_function", "add");
                        m_RTPS_settings->SetSetting("lt_increment", -.00);
                        m_RTPS_settings->SetSetting("lt_cl", "lifetime.cl");

                        

                        rtps::RTPS* ps = new rtps::RTPS(m_RTPS_settings);
                        (*slot)->m_pRTPS = ps;

                        //dynamic params
                        ps->settings->SetSetting("Gas Constant", rtmd->gas_constant);
                        ps->settings->SetSetting("Viscosity", rtmd->viscosity);
                        ps->settings->SetSetting("Velocity Limit", rtmd->velocity_limit);
                        ps->settings->SetSetting("XSPH Factor", rtmd->xsph_factor);
                        ps->settings->SetSetting("Gravity", rtmd->gravity); // -9.8 m/sec^2

                        ps->settings->SetSetting("Boundary Stiffness", rtmd->boundary_stiffness);
                        ps->settings->SetSetting("Boundary Dampening", rtmd->boundary_dampening);
                        //settings.SetSetting("Friction Kinetic", rtmd->friction_kinetic);
                        //settings.SetSetting("Friction Static", rtmd->friction_static);
                        ps->settings->SetSetting("sub_intervals", rtmd->sub_intervals);

                        //color hack for now
                        ps->settings->SetSetting("color_r", rtmd->color_r/255.0f);
                        ps->settings->SetSetting("color_g", rtmd->color_g/255.0f);
                        ps->settings->SetSetting("color_b", rtmd->color_b/255.0f);
                        ps->settings->SetSetting("color_a", rtmd->color_a/255.0f);


                    }
                   /* else if (sys == rtps::RTPSettings::FLOCK) 
                    {
						//printf("*** scale radius** = %f\n", rtmd->render_radius_scale);
						//printf("*** dt ** = %f\n", rtmd->dt);
						m_RTPS_settings = new rtps::RTPSettings(sys, rtmd->max_num, rtmd->dt, grid, rtmd->collision);
						//GE should automate with python
						m_RTPS_settings->setRadiusScale(rtmd->render_radius_scale);
						m_RTPS_settings->setRenderType((rtps::RTPSettings::RenderType)rtmd->render_type);
						m_RTPS_settings->setBlurScale(rtmd->render_blur_scale);
						m_RTPS_settings->setUseGLSL(rtmd->glsl);
						m_RTPS_settings->setUseAlphaBlending(rtmd->blending);

                        (*slot)->m_pRTPS = new rtps::RTPS(m_RTPS_settings);
        
                    } */
                    else if (sys == rtps::RTPSettings::FLOCK)
                    {
                        printf("RTPS: inside FLOCK system\n");
                        float color[3] = {rtmd->color_r, rtmd->color_g, rtmd->color_b};
                        //rtps::RTPSettings settings(sys, rtmd->max_num, rtmd->dt, grid, rtmd->maxspeed, rtmd->separationdist, rtmd->searchradius, color, rtmd->w_sep, rtmd->w_align, rtmd->w_coh);
                        //rtps::RTPSettings* settings = new rtps::RTPSettings(sys, rtmd->max_num, rtmd->dt, grid, rtmd->collision);
						m_RTPS_settings = new rtps::RTPSettings(sys, rtmd->max_num, rtmd->dt, grid);
                       
                        printf("RTPS: about to print the path\n");

                        //this path gives the location of the opencl and shader source files
                        printf("rtps path: %s\n", rtps_path.c_str());
                        m_RTPS_settings->SetSetting("rtps_path", rtps_path);

                        m_RTPS_settings->setRenderType((rtps::RTPSettings::RenderType)rtmd->render_type);
						m_RTPS_settings->setRadiusScale(rtmd->render_radius_scale);
						m_RTPS_settings->setBlurScale(rtmd->render_blur_scale);
						//settings.setUseGLSL(rtmd->glsl);
						//settings.setUseAlphaBlending(rtmd->blending);
						
                        //settings.setRadiusScale(1.0);
						//settings.setBlurScale(1.0);
						//settings.setUseGLSL(1);
                        
                        m_RTPS_settings->SetSetting("render_texture", "nemo.png");
                        m_RTPS_settings->SetSetting("render_frag_shader", "boid_tex_frag.glsl");
                        //settings->SetSetting("render_texture", "fsu_seal.jpg");
                        //settings->SetSetting("render_frag_shader", "sprite_tex_frag.glsl");
                        //settings->SetSetting("render_use_alpha", true);
                        m_RTPS_settings->SetSetting("render_use_alpha", false);
                        //settings->SetSetting("render_alpha_function", "add");
                        m_RTPS_settings->SetSetting("lt_increment", -.00);
                        m_RTPS_settings->SetSetting("lt_cl", "lifetime.cl");

                        rtps::RTPS* ps = new rtps::RTPS(m_RTPS_settings);
                        (*slot)->m_pRTPS = ps;

                        ps->settings->SetSetting("Max Speed", rtmd->maxspeed);
                        ps->settings->SetSetting("Min Separation Distance", rtmd->separationdist);
                        ps->settings->SetSetting("Searching Radius", rtmd->searchradius);
                        ps->settings->SetSetting("Separation Weight", rtmd->w_sep);
                        ps->settings->SetSetting("Alignment Weight", rtmd->w_align);
                        ps->settings->SetSetting("Cohesion Weight", rtmd->w_coh);
                    }
                    else 
                    {
						m_RTPS_settings = new rtps::RTPSettings(sys, rtmd->max_num, rtmd->dt, grid);
                        (*slot)->m_pRTPS = new rtps::RTPS(m_RTPS_settings);
                    }

                }
            }
        }
	}
	return true;
}


int makeEmitter(int num, KX_GameObject* gobj)
{
    CIntValue* intprop;
    CBoolValue* boolprop;

    int nn = 0; //the number of particles to emit
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
        CIntValue *newioffprop = new CIntValue(ioff);            
        gobj->SetProperty("ioff", newioffprop);

        CIntValue *numprop = new CIntValue(num);
        gobj->SetProperty("num", numprop);
    }
    else
    {
        //not a hose
        nn = num;
        CIntValue *numprop = new CIntValue(0);
        gobj->SetProperty("num", numprop);
    }

    return nn;

}


int getTriangles(KX_GameObject* gobj, std::vector<Triangle> &triangles)
{
    //get the mesh and the derived mesh so we can get faces/verts
    RAS_MeshObject* meshobj = gobj->GetMesh(0);
    Mesh* mesh = meshobj->GetMesh();
    DerivedMesh *dm = CDDM_from_mesh(mesh, NULL);
    
    //for now we are assuming triangles
    int n_faces = dm->getNumFaces(dm);
    //printf("**** num faces: %d\n", n_faces);
    MFace* face = dm->getFaceArray(dm);
    
    MT_Matrix3x3 grot = gobj->NodeGetWorldOrientation();
    MT_Point3 gp = gobj->NodeGetWorldPosition();

    Triangle tri;
    MT_Vector3 fv, ftv;

    for(int i = 0; i < n_faces; i++)
    {
        getTriangle(face[i], grot, gp, dm, tri);
        triangles.push_back(tri);
    } 
    //printf("in getTriangles: triangles.size(): %zd\n", triangles.size());

    return triangles.size();

}

void getTriangle(MFace& face, MT_Matrix3x3& grot, MT_Point3& gp, DerivedMesh* dm, Triangle& tri)
{
    MT_Vector3 fv, ftv;
    float mv[3];

    dm->getVertCo(dm, face.v1, mv);
    //rotate and translate by global rotation/translation to get global coords
    fv = MT_Vector3(mv[0], mv[1], mv[2]);
    ftv = grot*fv + gp;
    tri.verts[0] = float4(ftv.x(), ftv.y(), ftv.z(), 1);

    dm->getVertCo(dm, face.v2, mv);
    fv = MT_Vector3(mv[0], mv[1], mv[2]);
    ftv = grot*fv + gp;
    tri.verts[1] = float4(ftv.x(), ftv.y(), ftv.z(), 1);

    dm->getVertCo(dm, face.v3, mv);
    fv = MT_Vector3(mv[0], mv[1], mv[2]);
    ftv = grot*fv + gp;
    tri.verts[2] = float4(ftv.x(), ftv.y(), ftv.z(), 1);
    //printf("ftv: %f %f %f\n", ftv.x(), ftv.y(), ftv.z());

    float4 e1 = float4(tri.verts[1].x - tri.verts[0].x, tri.verts[1].y - tri.verts[0].y, tri.verts[1].z - tri.verts[0].z, 0);
    float4 e2 = float4(tri.verts[2].x - tri.verts[0].x, tri.verts[2].y - tri.verts[0].y, tri.verts[2].z - tri.verts[0].z, 0);
    tri.normal = normalize(cross(e1, e2));
}

