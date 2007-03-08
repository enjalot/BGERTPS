/* 
 * $Id$
 *
 * ***** BEGIN GPL/BL DUAL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version. The Blender
 * Foundation also sells licenses for use in proprietary software under
 * the Blender License.  See http://www.blender.org/BL/ for information
 * about this.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA	02111-1307, USA.
 *
 * The Original Code is Copyright (C) 2001-2002 by NaN Holding BV.
 * All rights reserved.
 *
 * This is a new part of Blender.
 *
 * Contributor(s): Michel Selten, Willian P. Germano, Alex Mole, Joseph Gilbert
 *
 * ***** END GPL/BL DUAL LICENSE BLOCK *****
*/

#ifndef EXPP_gen_utils_h
#define EXPP_gen_utils_h

#include <Python.h>
#include "DNA_ID.h"
#include "DNA_scriptlink_types.h"
#include "DNA_listBase.h"

#include "constant.h"

#define Py_PI  3.14159265358979323846
#define Py_WRAP 1024
#define Py_NEW  2048

/* 
   Py_RETURN_NONE
   Python 2.4 macro.  
   defined here until we switch to 2.4
*/
#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE	return Py_BuildValue("O", Py_None)
#endif
#ifndef Py_RETURN_FALSE
#define Py_RETURN_FALSE  return PyBool_FromLong(0) 
#endif
#ifndef Py_RETURN_TRUE
#define Py_RETURN_TRUE  return PyBool_FromLong(1)
#endif

/* name of list of Armature weak refs built into __main__ */
#define ARM_WEAKREF_LIST_NAME "__arm_weakrefs"

int EXPP_FloatsAreEqual(float A, float B, int floatSteps);
int EXPP_VectorsAreEqual(float *vecA, float *vecB, int size, int floatSteps);

PyObject *EXPP_GetModuleConstant(char *module, char *constant);

int StringEqual( const char *string1, const char *string2 );
char *GetIdName( ID * id );
int SetIdFakeUser( ID * id, PyObject *value);

ID *GetIdFromList( ListBase * list, char *name );

PyObject *PythonReturnErrorObject( PyObject * type, char *error_msg );
PyObject *PythonIncRef( PyObject * object );

char *event_to_name( short event );

float EXPP_ClampFloat( float value, float min, float max );
int EXPP_ClampInt( int value, int min, int max );

void EXPP_incr2( PyObject * ob1, PyObject * ob2 );
void EXPP_incr3( PyObject * ob1, PyObject * ob2, PyObject * ob3 );
void EXPP_decr2( PyObject * ob1, PyObject * ob2 );
void EXPP_decr3( PyObject * ob1, PyObject * ob2, PyObject * ob3 );
PyObject *EXPP_incr_ret( PyObject * object );
PyObject *EXPP_incr_ret_True(void);
PyObject *EXPP_incr_ret_False(void);
PyObject *EXPP_ReturnPyObjError( PyObject * type, char *error_msg );
int EXPP_ReturnIntError( PyObject * type, char *error_msg );

PyObject *EXPP_objError(PyObject *type, const char *format, ...);
int EXPP_intError(PyObject *type, const char *format, ...);

int EXPP_check_sequence_consistency( PyObject * seq, PyTypeObject * against );
PyObject *EXPP_tuple_repr( PyObject * self, int size );

/* mapping utilities - see Texture.c for an example of how to use these */
typedef struct {
	const char *sval;
	int ival;
} EXPP_map_pair;

/* maps must end with a pair that has NULL as sval */
int EXPP_map_getIntVal( const EXPP_map_pair * map,
			const char *sval, int *ival );
int EXPP_map_case_getIntVal( const EXPP_map_pair * map,
			     const char *sval, int *ival );
int EXPP_map_getShortVal( const EXPP_map_pair * map,
			  const char *sval, short *ival );
int EXPP_map_getStrVal( const EXPP_map_pair * map,
			int ival, const char **sval );

/* clamping and range-checking utilities */

int EXPP_setIValueClamped( PyObject *value, void *param,
		int min, int max, char type );
int EXPP_setFloatClamped ( PyObject *value, float *param,
			float min, float max);
int EXPP_setVec3Clamped ( PyObject *value, float *param,
			float min, float max);
int EXPP_setIValueRange( PyObject *value, void *param,
		int min, int max, char type );
int EXPP_setFloatRange ( PyObject *value, float *param,
			float min, float max);

/* utility routine for PyType attributes setters with module constant */

int EXPP_setModuleConstant ( BPy_constant *constant, void *param,
			char type );

/* utilities to get/set bits in bitfields */

PyObject *EXPP_getBitfield( void *param, int setting, char type );
int EXPP_setBitfield( PyObject * value, void *param, int setting, char type );

/*
 * Procedures to handle older setStuff() methods, which now access 
 * a PyType's setter attributes using the tp_getset mechanism.
 */

PyObject *EXPP_setterWrapper ( PyObject * self, PyObject * args,
				setter func);

PyObject *EXPP_setterWrapperTuple ( PyObject * self, PyObject * args,
				setter func);

/* scriplinks-related: */
PyObject *EXPP_getScriptLinks(ScriptLink *slink, PyObject *args, int is_scene);
PyObject *EXPP_addScriptLink(ScriptLink *slink, PyObject *args, int is_scene);
PyObject *EXPP_clearScriptLinks(ScriptLink *slink, PyObject *args);

/* this queues redraws if we're not in background mode: */
void EXPP_allqueue(unsigned short event, short val);

/* helper to keep dictionaries from causing memory leaks */
int EXPP_dict_set_item_str( PyObject *dict, char *key, PyObject *value);





/* Dummy struct for getting the ID from a libdata BPyObject */
typedef struct {
	PyObject_HEAD		/* required python macro */
	ID *id;
} BPy_GenericLib;


/* ID functions for all libdata */
#define	GENERIC_LIB_GETSETATTR \
	{"name",\
	 (getter)GenericLib_getName, (setter)GenericLib_setName,\
	 "name",\
	 NULL},\
	{"lib",\
	 (getter)GenericLib_getLib, (setter)NULL,\
	 "external library path",\
	 NULL},\
	{"users",\
	 (getter)GenericLib_getUsers, (setter)NULL,\
	 "user count",\
	 NULL},\
	{"fakeUser",\
	 (getter)GenericLib_getFakeUser, (setter)GenericLib_setFakeUser,\
	 "fake user state",\
	 NULL},\
	{"properties",\
	 (getter)GenericLib_getProperties, (setter)NULL,\
	 "properties",\
	 NULL}


int GenericLib_setName( void *self, PyObject *value );
PyObject *GenericLib_getName( void *self );
PyObject *GenericLib_getFakeUser( void *self );
int GenericLib_setFakeUser( void *self, PyObject *value );
PyObject *GenericLib_getLib( void *self );
PyObject *GenericLib_getUsers( void *self );
PyObject *GenericLib_getProperties( void *self );

/* use this for oldstyle somedata.getName("name") */
PyObject * GenericLib_setName_with_method( void *self, PyObject *value ); 

int GenericLib_assignData(PyObject *value, void **data, void **ndata, short refcount, short type, short subtype);
short GenericLib_getType(PyObject * pydata);

#endif				/* EXPP_gen_utils_h */
