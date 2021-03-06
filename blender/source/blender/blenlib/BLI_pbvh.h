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
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef BLI_PBVH_H
#define BLI_PBVH_H

/** \file BLI_pbvh.h
 *  \ingroup bli
 *  \brief A BVH for high poly meshes.
 */

struct MFace;
struct MVert;
struct DMGridAdjacency;
struct DMGridData;
struct PBVH;
struct PBVHNode;
struct ListBase;

typedef struct PBVH PBVH;
typedef struct PBVHNode PBVHNode;

typedef struct {
	float (*co)[3];
} PBVHProxyNode;

/* Callbacks */

/* returns 1 if the search should continue from this node, 0 otherwise */
typedef int (*BLI_pbvh_SearchCallback)(PBVHNode *node, void *data);

typedef void (*BLI_pbvh_HitCallback)(PBVHNode *node, void *data);
typedef void (*BLI_pbvh_HitOccludedCallback)(PBVHNode *node, void *data, float* tmin);

/* Building */

PBVH *BLI_pbvh_new(void);
void BLI_pbvh_build_mesh(PBVH *bvh, struct MFace *faces, struct MVert *verts,
			int totface, int totvert);
void BLI_pbvh_build_grids(PBVH *bvh, struct DMGridData **grids,
	struct DMGridAdjacency *gridadj, int totgrid,
	int gridsize, void **gridfaces);
void BLI_pbvh_free(PBVH *bvh);

/* Hierarchical Search in the BVH, two methods:
   * for each hit calling a callback
   * gather nodes in an array (easy to multithread) */

void BLI_pbvh_search_callback(PBVH *bvh,
	BLI_pbvh_SearchCallback scb, void *search_data,
	BLI_pbvh_HitCallback hcb, void *hit_data);

void BLI_pbvh_search_gather(PBVH *bvh,
	BLI_pbvh_SearchCallback scb, void *search_data,
	PBVHNode ***array, int *tot);

/* Raycast
   the hit callback is called for all leaf nodes intersecting the ray;
   it's up to the callback to find the primitive within the leaves that is
   hit first */

void BLI_pbvh_raycast(PBVH *bvh, BLI_pbvh_HitOccludedCallback cb, void *data,
			  float ray_start[3], float ray_normal[3], int original);
int BLI_pbvh_node_raycast(PBVH *bvh, PBVHNode *node, float (*origco)[3],
	float ray_start[3], float ray_normal[3], float *dist);

/* Drawing */

void BLI_pbvh_node_draw(PBVHNode *node, void *data);
int BLI_pbvh_node_planes_contain_AABB(PBVHNode *node, void *data);
void BLI_pbvh_draw(PBVH *bvh, float (*planes)[4], float (*face_nors)[3], int smooth);

/* Node Access */

typedef enum {
	PBVH_Leaf = 1,

	PBVH_UpdateNormals = 2,
	PBVH_UpdateBB = 4,
	PBVH_UpdateOriginalBB = 8,
	PBVH_UpdateDrawBuffers = 16,
	PBVH_UpdateRedraw = 32
} PBVHNodeFlags;

void BLI_pbvh_node_mark_update(PBVHNode *node);

void BLI_pbvh_node_get_grids(PBVH *bvh, PBVHNode *node,
	int **grid_indices, int *totgrid, int *maxgrid, int *gridsize,
	struct DMGridData ***griddata, struct DMGridAdjacency **gridadj);
void BLI_pbvh_node_num_verts(PBVH *bvh, PBVHNode *node,
	int *uniquevert, int *totvert);
void BLI_pbvh_node_get_verts(PBVH *bvh, PBVHNode *node,
	int **vert_indices, struct MVert **verts);

void BLI_pbvh_node_get_BB(PBVHNode *node, float bb_min[3], float bb_max[3]);
void BLI_pbvh_node_get_original_BB(PBVHNode *node, float bb_min[3], float bb_max[3]);

float BLI_pbvh_node_get_tmin(PBVHNode* node);

/* Update Normals/Bounding Box/Draw Buffers/Redraw and clear flags */

void BLI_pbvh_update(PBVH *bvh, int flags, float (*face_nors)[3]);
void BLI_pbvh_redraw_BB(PBVH *bvh, float bb_min[3], float bb_max[3]);
void BLI_pbvh_get_grid_updates(PBVH *bvh, int clear, void ***gridfaces, int *totface);
void BLI_pbvh_grids_update(PBVH *bvh, struct DMGridData **grids,
	struct DMGridAdjacency *gridadj, void **gridfaces);

/* vertex deformer */
float (*BLI_pbvh_get_vertCos(struct PBVH *pbvh))[3];
void BLI_pbvh_apply_vertCos(struct PBVH *pbvh, float (*vertCos)[3]);
int BLI_pbvh_isDeformed(struct PBVH *pbvh);


/* Vertex Iterator */

/* this iterator has quite a lot of code, but it's designed to:
   - allow the compiler to eliminate dead code and variables
   - spend most of the time in the relatively simple inner loop */

#define PBVH_ITER_ALL		0
#define PBVH_ITER_UNIQUE	1

typedef struct PBVHVertexIter {
	/* iteration */
	int g;
	int width;
	int height;
	int skip;
	int gx;
	int gy;
	int i;

	/* grid */
	struct DMGridData **grids;
	struct DMGridData *grid;
	int *grid_indices;
	int totgrid;
	int gridsize;

	/* mesh */
	struct MVert *mverts;
	int totvert;
	int *vert_indices;

	/* result: these are all computed in the macro, but we assume
	   that compiler optimizations will skip the ones we don't use */
	struct MVert *mvert;
	float *co;
	short *no;
	float *fno;
} PBVHVertexIter;

#ifdef _MSC_VER
#pragma warning (disable:4127) // conditional expression is constant
#endif

#define BLI_pbvh_vertex_iter_begin(bvh, node, vi, mode) \
	{ \
		struct DMGridData **grids; \
		struct MVert *verts; \
		int *grid_indices, totgrid, gridsize, *vert_indices, uniq_verts, totvert; \
		\
		vi.grid= 0; \
		vi.no= 0; \
		vi.fno= 0; \
		vi.mvert= 0; \
		vi.skip= 0; \
		\
		BLI_pbvh_node_get_grids(bvh, node, &grid_indices, &totgrid, NULL, &gridsize, &grids, NULL); \
		BLI_pbvh_node_num_verts(bvh, node, &uniq_verts, &totvert); \
		BLI_pbvh_node_get_verts(bvh, node, &vert_indices, &verts); \
		\
		vi.grids= grids; \
		vi.grid_indices= grid_indices; \
		vi.totgrid= (grids)? totgrid: 1; \
		vi.gridsize= gridsize; \
		\
		if(mode == PBVH_ITER_ALL) \
			vi.totvert = totvert; \
		else \
			vi.totvert= uniq_verts; \
		vi.vert_indices= vert_indices; \
		vi.mverts= verts; \
	}\
	\
	for(vi.i=0, vi.g=0; vi.g<vi.totgrid; vi.g++) { \
		if(vi.grids) { \
			vi.width= vi.gridsize; \
			vi.height= vi.gridsize; \
			vi.grid= vi.grids[vi.grid_indices[vi.g]]; \
			vi.skip= 0; \
			 \
			/*if(mode == PVBH_ITER_UNIQUE) { \
				vi.grid += subm->grid.offset; \
				vi.skip= subm->grid.skip; \
				vi.grid -= skip; \
			}*/ \
		} \
		else { \
			vi.width= vi.totvert; \
			vi.height= 1; \
		} \
		 \
		for(vi.gy=0; vi.gy<vi.height; vi.gy++) { \
			if(vi.grid) vi.grid += vi.skip; \
			\
			for(vi.gx=0; vi.gx<vi.width; vi.gx++, vi.i++) { \
				if(vi.grid) { \
					vi.co= vi.grid->co; \
					vi.fno= vi.grid->no; \
					vi.grid++; \
				} \
				else { \
					vi.mvert= &vi.mverts[vi.vert_indices[vi.gx]]; \
					vi.co= vi.mvert->co; \
					vi.no= vi.mvert->no; \
				} \

#define BLI_pbvh_vertex_iter_end \
			} \
		} \
	}

void BLI_pbvh_node_get_proxies(PBVHNode* node, PBVHProxyNode** proxies, int* proxy_count);
void BLI_pbvh_node_free_proxies(PBVHNode* node);
PBVHProxyNode* BLI_pbvh_node_add_proxy(PBVH* bvh, PBVHNode* node);
void BLI_pbvh_gather_proxies(PBVH* pbvh, PBVHNode*** nodes,  int* totnode);

//void BLI_pbvh_node_BB_reset(PBVHNode* node);
//void BLI_pbvh_node_BB_expand(PBVHNode* node, float co[3]);

#endif /* BLI_PBVH_H */

