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
* Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
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

#include <stdio.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "DNA_listBase.h"
#include "DNA_userdef_types.h" 

#include "BIF_cursors.h"
#include "BIF_resources.h"
#include "BIF_graphics.h"
#include "BIF_screen.h"

#include "winlay.h"


/* ****************************************************************** 
Cursor Description:

Each bit represents a pixel, so 1 byte = 8 pixels, 
the bytes go Left to Right. Top to bottom
the bits in a byte go right to left
(ie;  0x01, 0x80  represents a line of 16 pix with the first and last pix set.) 

A 0 in the bitmap = bg_color, a 1 fg_color
a 0 in the mask   = transparent pix.

Until 32x32 cursors are supported on all platforms, the size of the 
small cursors MUST be 16x16.

Large cursors have a MAXSIZE of 32x32.

Other than that, the specified size of the cursors is just a guideline, 
However, the char array that defines the BM and MASK must be byte aligned.
ie a 17x17 cursor needs 3 bytes (cols) * 17 bytes (rows) 
(3 bytes = 17 bits rounded up to nearest whole byte).  Pad extra bits
in mask with 0's.

Setting big_bm=NULL disables the large version of the cursor.

******************************************************************* 

There is a nice Python GUI utility that can be used for drawing cursors in
this format in the Blender source distribution, in 
blender/source/tools/MakeCursor.py . Start it with $ python MakeCursor.py
It will copy its output to the console when you press 'Do it'.

*/

/* Because defining a cursor mixes declarations and executable code
   each cursor needs it's own scoping block or it would be split up 
   over several hundred lines of code.  To enforce/document this better
   I define 2 pretty braindead macros so it's obvious what the extra "[]"
   are for */

#define BEGIN_CURSOR_BLOCK {
#define END_CURSOR_BLOCK   }

/* Cursor Globals */
static BCursor *BlenderCursor[BC_NUMCURSORS]; /*Points to static BCursor Structs */
static short CurrentCursor=-1, LastCursor=-1;

short GetBlenderCursor(void) {
	return CurrentCursor;
}

void SetBlenderCursor(short curs){
	Window *win;

	if ((curs<LASTCURSOR)||(curs>=BC_NUMCURSORS)) return;

	win=winlay_get_active_window();

	if (win==NULL) return; /* Can't set custom cursor before Window init */ 

	LastCursor=CurrentCursor;
	CurrentCursor=curs;

	if (curs==LASTCURSOR) curs=LastCursor;

	if (curs==SYSCURSOR) {  /* System default Cursor */
		window_set_cursor(win, CURSOR_STD);
		//set_cursor(CURSOR_STD);
	}
	else if ( (U.curssize==0) || (BlenderCursor[curs]->big_bm == NULL) ) {
		window_set_custom_cursor_ex(win, BlenderCursor[curs], 0);
	}
	else {
		window_set_custom_cursor_ex(win, BlenderCursor[curs], 1);
	}
}

/* unused no prototypes 
short GetCurrentCursor(void){
	return(CurrentCursor);
}
*/

void InitCursorData(void){

	/********************** NW_ARROW Cursor **************************/
BEGIN_CURSOR_BLOCK
		static char nw_sbm[]={
			0x03,  0x00,  0x05,  0x00,  0x09,  0x00,  0x11,  0x00,
				0x21,  0x00,  0x41,  0x00,  0x81,  0x00,  0x01,  0x01,
				0x01,  0x02,  0xc1,  0x03,  0x49,  0x00,  0x8d,  0x00,
				0x8b,  0x00,  0x10,  0x01,  0x90,  0x01,  0x60,  0x00,
		};

		static char nw_smsk[]={
			0x03,  0x00,  0x07,  0x00,  0x0f,  0x00,  0x1f,  0x00,
				0x3f,  0x00,  0x7f,  0x00,  0xff,  0x00,  0xff,  0x01,
				0xff,  0x03,  0xff,  0x03,  0x7f,  0x00,  0xff,  0x00,
				0xfb,  0x00,  0xf0,  0x01,  0xf0,  0x01,  0x60,  0x00,
		};

		static BCursor NWArrowCursor = {
			/*small*/
			nw_sbm, nw_smsk,
				16, 16, 
				6,  7,
				/*big*/
				NULL, NULL,
				32,32, 
				15, 15,
				/*color*/
				BC_BLACK, BC_WHITE
		};

		BlenderCursor[BC_NW_ARROWCURSOR]=&NWArrowCursor;
END_CURSOR_BLOCK

	///********************** NS_ARROW Cursor *************************/
BEGIN_CURSOR_BLOCK
		static char ns_sbm[]={
			0x40,  0x01,  0x20,  0x02,  0x10,  0x04,  0x08,  0x08,
				0x04,  0x10,  0x3c,  0x1e,  0x20,  0x02,  0x20,  0x02,
				0x20,  0x02,  0x20,  0x02,  0x3c,  0x1e,  0x04,  0x10,
				0x08,  0x08,  0x10,  0x04,  0x20,  0x02,  0x40,  0x01
		};

		static char ns_smsk[]={
			0xc0,  0x01,  0xe0,  0x03,  0xf0,  0x07,  0xf8,  0x0f,
				0xfc,  0x1f,  0xfc,  0x1f,  0xe0,  0x03,  0xe0,  0x03,
				0xe0,  0x03,  0xe0,  0x03,  0xfc,  0x1f,  0xfc,  0x1f,
				0xf8,  0x0f,  0xf0,  0x07,  0xe0,  0x03,  0xc0,  0x01
		};

		static BCursor NSArrowCursor = {
			/*small*/
			ns_sbm, ns_smsk,
				16, 16, 
				6,  7,
				/*big*/
				NULL, NULL,
				32,32, 
				15, 15,
				/*color*/
				BC_BLACK, BC_WHITE
		};

		BlenderCursor[BC_NS_ARROWCURSOR]=&NSArrowCursor;
		
END_CURSOR_BLOCK
	/********************** EW_ARROW Cursor *************************/
BEGIN_CURSOR_BLOCK
	static char ew_sbm[]={
		0x00,  0x00,  0x00,  0x00,  0x10,  0x08,  0x38,  0x1c,
		0x2c,  0x34,  0xe6,  0x67,  0x03,  0xc0,  0x01,  0x80,
		0x03,  0xc0,  0xe6,  0x67,  0x2c,  0x34,  0x38,  0x1c,
		0x10,  0x08,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
	};

	static char ew_smsk[]={
		0x00,  0x00,  0x00,  0x00,  0x10,  0x08,  0x38,  0x1c,
		0x3c,  0x3c,  0xfe,  0x7f,  0xff,  0xff,  0x3f,  0xfc,
		0xff,  0xff,  0xfe,  0x7f,  0x3c,  0x3c,  0x38,  0x1c,
		0x10,  0x08,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
	};

	static BCursor EWArrowCursor = {
		/*small*/
		ew_sbm, ew_smsk,
		16, 16, 
		7,  6,
		/*big*/
		NULL, NULL,
		32,32, 
		15, 15,
		/*color*/
		BC_BLACK, BC_WHITE
	};

	BlenderCursor[BC_EW_ARROWCURSOR]=&EWArrowCursor;
END_CURSOR_BLOCK

	/********************** Wait Cursor *****************************/
BEGIN_CURSOR_BLOCK
	static char wait_sbm[]={
		0xfe,  0x7f,  0x02,  0x40,  0x02,  0x40,  0x84,  0x21,
		0xc8,  0x13,  0xd0,  0x0b,  0xa0,  0x04,  0x20,  0x05,
		0xa0,  0x04,  0x10,  0x09,  0x88,  0x11,  0xc4,  0x23,
		0xe2,  0x47,  0xfa,  0x5f,  0x02,  0x40,  0xfe,  0x7f,
	};

	static char wait_smsk[]={
		0xfe,  0x7f,  0xfe,  0x7f,  0x06,  0x60,  0x8c,  0x31,
		0xd8,  0x1b,  0xf0,  0x0f,  0xe0,  0x06,  0x60,  0x07,
		0xe0,  0x06,  0x30,  0x0d,  0x98,  0x19,  0xcc,  0x33,
		0xe6,  0x67,  0xfe,  0x7f,  0xfe,  0x7f,  0xfe,  0x7f,
	};

	static char wait_lbm[]={
		0xfc,  0xff,  0xff,  0x3f,  0xfc,  0xff,  0xff,  0x3f,
		0x0c,  0x00,  0x00,  0x30,  0x0c,  0x00,  0x00,  0x30,
		0x0c,  0x00,  0x00,  0x30,  0x0c,  0x00,  0x00,  0x18,
		0x18,  0xc0,  0x03,  0x0c,  0x30,  0x20,  0x07,  0x06,
		0x60,  0xf0,  0x0f,  0x03,  0xc0,  0xd0,  0x8d,  0x01,
		0x80,  0x79,  0xcf,  0x00,  0x00,  0xf3,  0x67,  0x00,
		0x00,  0x66,  0x37,  0x00,  0x00,  0x8c,  0x33,  0x00,
		0x00,  0x0c,  0x32,  0x00,  0x00,  0xcc,  0x33,  0x00,
		0x00,  0x8c,  0x30,  0x00,  0x00,  0x46,  0x61,  0x00,
		0x00,  0x03,  0xc3,  0x00,  0x80,  0x01,  0x83,  0x01,
		0xc0,  0xc0,  0x03,  0x03,  0x60,  0xa0,  0x05,  0x06,
		0x30,  0xf0,  0x0f,  0x0c,  0x18,  0xf8,  0x1d,  0x18,
		0x0c,  0x5c,  0x3f,  0x30,  0x0c,  0xff,  0x5f,  0x30,
		0x0c,  0xf7,  0xfe,  0x31,  0xcc,  0xfb,  0x9f,  0x33,
		0x0c,  0x00,  0x00,  0x30,  0x0c,  0x00,  0x00,  0x30,
		0xfc,  0xff,  0xff,  0x3f,  0xfc,  0xff,  0xff,  0x3f,
	};

	static char wait_lmsk[]={
		0xfc,  0xff,  0xff,  0x3f,  0xfc,  0xff,  0xff,  0x3f,
		0xfc,  0xff,  0xff,  0x3f,  0xfc,  0xff,  0xff,  0x3f,
		0x3c,  0x00,  0x00,  0x3c,  0x3c,  0x00,  0x00,  0x1e,
		0x78,  0xc0,  0x03,  0x0f,  0xf0,  0xa0,  0x87,  0x07,
		0xe0,  0xf1,  0xcf,  0x03,  0xc0,  0xf3,  0xef,  0x01,
		0x80,  0xff,  0xff,  0x00,  0x00,  0xff,  0x7f,  0x00,
		0x00,  0xfe,  0x3f,  0x00,  0x00,  0xfc,  0x3f,  0x00,
		0x00,  0x3c,  0x3f,  0x00,  0x00,  0xfc,  0x3f,  0x00,
		0x00,  0xbc,  0x3c,  0x00,  0x00,  0xde,  0x79,  0x00,
		0x00,  0x0f,  0xf3,  0x00,  0x80,  0x07,  0xe3,  0x01,
		0xc0,  0xc3,  0xc3,  0x03,  0xe0,  0xe1,  0x87,  0x07,
		0xf0,  0xf0,  0x0f,  0x0f,  0x78,  0xf8,  0x1f,  0x1e,
		0x3c,  0x7c,  0x3f,  0x3c,  0x3c,  0xff,  0x7f,  0x3c,
		0xbc,  0xff,  0xff,  0x3d,  0xfc,  0xfb,  0xbf,  0x3f,
		0xfc,  0xff,  0xff,  0x3f,  0xfc,  0xff,  0xff,  0x3f,
		0xfc,  0xff,  0xff,  0x3f,  0xfc,  0xff,  0xff,  0x3f,
	};

	static BCursor WaitCursor = {
		/*small*/
	wait_sbm, wait_smsk,
		16, 16, 
		7,  7,
		/*big*/
		wait_lbm, wait_lmsk,
		32,32, 
		15, 15,
		/*color*/
		BC_BLACK, BC_WHITE
	};

	BlenderCursor[BC_WAITCURSOR]=&WaitCursor;
END_CURSOR_BLOCK

	/********************** Cross Cursor ***************************/
BEGIN_CURSOR_BLOCK
	static char cross_sbm[]={
		0x00,  0x00,  0x80,  0x01,  0x80,  0x01,  0x80,  0x01,
		0x80,  0x01,  0x80,  0x01,  0x80,  0x01,  0x7e,  0x7e,
		0x7e,  0x7e,  0x80,  0x01,  0x80,  0x01,  0x80,  0x01,
		0x80,  0x01,  0x80,  0x01,  0x80,  0x01,  0x00,  0x00,
	};

	static char cross_smsk[]={
		0x80,  0x01,  0x80,  0x01,  0x80,  0x01,  0x80,  0x01,
		0x80,  0x01,  0x80,  0x01,  0xc0,  0x03,  0x7f,  0xfe,
		0x7f,  0xfe,  0xc0,  0x03,  0x80,  0x01,  0x80,  0x01,
		0x80,  0x01,  0x80,  0x01,  0x80,  0x01,  0x80,  0x01,
	};
	static char cross_lbm[]={
		0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0xc0,  0x03,  0x00,
		0x00,  0xc0,  0x03,  0x00,  0x00,  0x40,  0x02,  0x00,
		0x00,  0x78,  0x1e,  0x00,  0xfc,  0x1f,  0xf8,  0x3f,
		0xfc,  0x1f,  0xf8,  0x3f,  0x00,  0x78,  0x1e,  0x00,
		0x00,  0x40,  0x02,  0x00,  0x00,  0xc0,  0x03,  0x00,
		0x00,  0xc0,  0x03,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
	};

	static char cross_lmsk[]={
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0xc0,  0x03,  0x00,
		0x00,  0xe0,  0x07,  0x00,  0x00,  0x70,  0x0e,  0x00,
		0x00,  0x78,  0x1e,  0x00,  0xff,  0x1f,  0xf8,  0xff,
		0xff,  0x1f,  0xf8,  0xff,  0x00,  0x78,  0x1e,  0x00,
		0x00,  0x70,  0x0e,  0x00,  0x00,  0xe0,  0x07,  0x00,
		0x00,  0xc0,  0x03,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
		0x00,  0x80,  0x01,  0x00,  0x00,  0x80,  0x01,  0x00,
	};

	static BCursor CrossCursor = {
		/*small*/
		cross_sbm, cross_smsk,
		16, 16, 
		7,  7,
		/*big*/
		cross_lbm, cross_lmsk,
		32,32, 
		15, 15,
		/*color*/
		BC_BLACK, BC_WHITE
	};

	BlenderCursor[BC_CROSSCURSOR]=&CrossCursor;
END_CURSOR_BLOCK

	/********************** EditCross Cursor ***********************/	
BEGIN_CURSOR_BLOCK
	static char editcross_sbm[]={
		0x0e,  0x00,  0x11,  0x00,  0x1d,  0x00,  0x19,  0x03,
		0x1d,  0x03,  0x11,  0x03,  0x0e,  0x03,  0x00,  0x03,
		0xf8,  0x7c,  0xf8,  0x7c,  0x00,  0x03,  0x00,  0x03,
		0x00,  0x03,  0x00,  0x03,  0x00,  0x03,  0x00,  0x00,
	};

	static char editcross_smsk[]={
		0x0e,  0x00,  0x1f,  0x00,  0x1f,  0x03,  0x1f,  0x03,
		0x1f,  0x03,  0x1f,  0x03,  0x0e,  0x03,  0x80,  0x07,
		0xfc,  0xfc,  0xfc,  0xfc,  0x80,  0x07,  0x00,  0x03,
		0x00,  0x03,  0x00,  0x03,  0x00,  0x03,  0x00,  0x03,
	};

	static BCursor EditCrossCursor = {
		/*small*/
		editcross_sbm, editcross_smsk,
		16, 16, 
		9,  8,
		/*big*/
		NULL, NULL,
		32,32, 
		15, 15,
		/*color*/
		BC_BLACK, BC_WHITE
	};

	BlenderCursor[BC_EDITCROSSCURSOR]=&EditCrossCursor;
END_CURSOR_BLOCK

	/********************** Box Select *************************/
BEGIN_CURSOR_BLOCK
	static char box_sbm[32]={
	0x7f,  0x00,  0x41,  0x00,  0x41,  0x00,  0x41,  0x06,
		0x41,  0x06,  0x41,  0x06,  0x7f,  0x06,  0x00,  0x06,
		0xe0,  0x79,  0xe0,  0x79,  0x00,  0x06,  0x00,  0x06,
		0x00,  0x06,  0x00,  0x06,  0x00,  0x06,  0x00,  0x00,
	};

	static char box_smsk[32]={
	0x7f,  0x00,  0x7f,  0x00,  0x63,  0x06,  0x63,  0x06,
		0x63,  0x06,  0x7f,  0x06,  0x7f,  0x06,  0x00,  0x0f,
		0xf0,  0xf9,  0xf0,  0xf9,  0x00,  0x0f,  0x00,  0x06,
		0x00,  0x06,  0x00,  0x06,  0x00,  0x06,  0x00,  0x06,

	};

	static BCursor BoxSelCursor = {
		/*small*/
		box_sbm, box_smsk,
		16, 16, 
		9,  8,
		/*big*/
		NULL, NULL,
		32,32, 
		15, 15,
		/*color*/
		BC_BLACK, BC_WHITE
	};

	BlenderCursor[BC_BOXSELCURSOR]=&BoxSelCursor;

END_CURSOR_BLOCK
	/********************** Knife Cursor ***********************/
BEGIN_CURSOR_BLOCK
	static char knife_sbm[]={
		0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x2c,
		0x00, 0x5a, 0x00, 0x34, 0x00, 0x2a, 0x00, 0x17,
		0x80, 0x06, 0x40, 0x03, 0xa0, 0x03, 0xd0, 0x01,
		0x68, 0x00, 0x1c, 0x00, 0x06, 0x00, 0x00, 0x00
	};

	static char knife_smsk[]={
		0x00, 0x60, 0x00, 0xf0, 0x00, 0xfc, 0x00, 0xfe,
		0x00, 0xfe, 0x00, 0x7e, 0x00, 0x7f, 0x80, 0x3f,
		0xc0, 0x0e, 0x60, 0x07, 0xb0, 0x07, 0xd8, 0x03,
		0xec, 0x01, 0x7e, 0x00, 0x1f, 0x00, 0x07, 0x00
	};

	static char knife_lbm[]={
		0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
		0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
		0x00,  0x00,  0x00,  0x08,  0x00,  0x00,  0x00,  0x1c,
		0x00,  0x00,  0x00,  0x3e,  0x00,  0x00,  0x00,  0x7f,
		0x00,  0x00,  0x80,  0xbf,  0x00,  0x00,  0xc0,  0x5f,
		0x00,  0x00,  0xc0,  0x6f,  0x00,  0x00,  0xc0,  0x37,
		0x00,  0x00,  0xa8,  0x1b,  0x00,  0x00,  0x54,  0x0d,
		0x00,  0x00,  0xa8,  0x00,  0x00,  0x00,  0x54,  0x00,
		0x00,  0x00,  0xa8,  0x00,  0x00,  0x00,  0x53,  0x00,
		0x00,  0xc0,  0x07,  0x00,  0x00,  0xe0,  0x0f,  0x00,
		0x00,  0xd0,  0x0f,  0x00,  0x00,  0xe8,  0x07,  0x00,
		0x00,  0xf4,  0x07,  0x00,  0x00,  0xfa,  0x00,  0x00,
		0x00,  0x3d,  0x00,  0x00,  0x80,  0x0e,  0x00,  0x00,
		0xc0,  0x03,  0x00,  0x00,  0xe0,  0x00,  0x00,  0x00,
		0x30,  0x00,  0x00,  0x00,  0x08,  0x00,  0x00,  0x00,
		0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,

	};

	static char knife_lmsk[]={
		0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
		0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x18,
		0x00,  0x00,  0x00,  0x3c,  0x00,  0x00,  0x00,  0x7e,
		0x00,  0x00,  0x00,  0xff,  0x00,  0x00,  0x80,  0xff,
		0x00,  0x00,  0xc0,  0xbf,  0x00,  0x00,  0xe0,  0xdf,
		0x00,  0x00,  0xe0,  0xef,  0x00,  0x00,  0xf8,  0x77,
		0x00,  0x00,  0xfc,  0x3b,  0x00,  0x00,  0xfe,  0x1d,
		0x00,  0x00,  0xfe,  0x0f,  0x00,  0x00,  0xfe,  0x01,
		0x00,  0x00,  0xff,  0x01,  0x00,  0xc0,  0xff,  0x00,
		0x00,  0xe0,  0x7f,  0x00,  0x00,  0xf0,  0x1f,  0x00,
		0x00,  0xd8,  0x1f,  0x00,  0x00,  0xec,  0x0f,  0x00,
		0x00,  0xf6,  0x0f,  0x00,  0x00,  0xfb,  0x06,  0x00,
		0x80,  0xbd,  0x01,  0x00,  0xc0,  0x6e,  0x00,  0x00,
		0xe0,  0x1b,  0x00,  0x00,  0xf0,  0x06,  0x00,  0x00,
		0xb8,  0x01,  0x00,  0x00,  0x6c,  0x00,  0x00,  0x00,
		0x1c,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,

	};

	static BCursor KnifeCursor = {
		/*small*/
	knife_sbm, knife_smsk,
		16, 16, 
		0,  15,
		/*big*/
		knife_lbm, knife_lmsk,
		32,32, 
		0, 31,
		/*color*/
		BC_BLACK, BC_WHITE
	};

	BlenderCursor[BC_KNIFECURSOR]=&KnifeCursor;

END_CURSOR_BLOCK
	
	/********************** Loop Select Cursor ***********************/
BEGIN_CURSOR_BLOCK

static char vloop_sbm[]={
        0x00,  0x00,  0x7e,  0x00,  0x3e,  0x00,  0x1e,  0x00,
        0x0e,  0x00,  0x66,  0x60,  0x62,  0x6f,  0x00,  0x00,
        0x20,  0x20,  0x20,  0x20,  0x20,  0x20,  0x20,  0x20,
        0x00,  0x00,  0x60,  0x60,  0x60,  0x6f,  0x00,  0x00,
};

static char vloop_smsk[]={
        0xff,  0x01,  0xff,  0x00,  0x7f,  0x00,  0x3f,  0x00,
        0xff,  0xf0,  0xff,  0xff,  0xf7,  0xff,  0xf3,  0xf0,
        0x61,  0x60,  0x60,  0x60,  0x60,  0x60,  0x60,  0x60,
        0xf0,  0xf0,  0xf0,  0xff,  0xf0,  0xff,  0xf0,  0xf0,
};



static char vloop_lbm[]={
        0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
        0xfc,  0x3f,  0x00,  0x00,  0xfc,  0x3f,  0x00,  0x00,
        0xfc,  0x0f,  0x00,  0x00,  0xfc,  0x0f,  0x00,  0x00,
        0xfc,  0x03,  0x00,  0x00,  0xfc,  0x03,  0x00,  0x00,
        0xfc,  0x00,  0x00,  0x00,  0xfc,  0x00,  0x00,  0x00,
        0x3c,  0x3c,  0x00,  0x3c,  0x3c,  0x3c,  0x00,  0x3c,
        0x0c,  0x3c,  0xff,  0x3c,  0x0c,  0x3c,  0xff,  0x3c,
        0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
        0x00,  0x0c,  0x00,  0x0c,  0x00,  0x0c,  0x00,  0x0c,
        0x00,  0x0c,  0x00,  0x0c,  0x00,  0x0c,  0x00,  0x0c,
        0x00,  0x0c,  0x00,  0x0c,  0x00,  0x0c,  0x00,  0x0c,
        0x00,  0x0c,  0x00,  0x0c,  0x00,  0x0c,  0x00,  0x0c,
        0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
        0x00,  0x3c,  0x00,  0x3c,  0x00,  0x3c,  0x00,  0x3c,
        0x00,  0x3c,  0xff,  0x3c,  0x00,  0x3c,  0xff,  0x3c,
        0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
};

static char vloop_lmsk[]={
        0xff,  0xff,  0x03,  0x00,  0xff,  0xff,  0x03,  0x00,
        0xff,  0xff,  0x00,  0x00,  0xff,  0xff,  0x00,  0x00,
        0xff,  0x3f,  0x00,  0x00,  0xff,  0x3f,  0x00,  0x00,
        0xff,  0x0f,  0x00,  0x00,  0xff,  0x0f,  0x00,  0x00,
        0xff,  0xff,  0x00,  0xff,  0xff,  0xff,  0x00,  0xff,
        0xff,  0xff,  0xff,  0xff,  0xff,  0xff,  0xff,  0xff,
        0x3f,  0xff,  0xff,  0xff,  0x3f,  0xff,  0xff,  0xff,
        0x0f,  0xff,  0x00,  0xff,  0x0f,  0xff,  0x00,  0xff,
        0x03,  0x3c,  0x00,  0x3c,  0x03,  0x3c,  0x00,  0x3c,
        0x00,  0x3c,  0x00,  0x3c,  0x00,  0x3c,  0x00,  0x3c,
        0x00,  0x3c,  0x00,  0x3c,  0x00,  0x3c,  0x00,  0x3c,
        0x00,  0x3c,  0x00,  0x3c,  0x00,  0x3c,  0x00,  0x3c,
        0x00,  0xff,  0x00,  0xff,  0x00,  0xff,  0x00,  0xff,
        0x00,  0xff,  0xff,  0xff,  0x00,  0xff,  0xff,  0xff,
        0x00,  0xff,  0xff,  0xff,  0x00,  0xff,  0xff,  0xff,
        0x00,  0xff,  0x00,  0xff,  0x00,  0xff,  0x00,  0xff,
};



	static BCursor VLoopCursor = {
		/*small*/
	vloop_sbm, vloop_smsk,
		16, 16, 
		0,  0,
		/*big*/
		vloop_lbm, vloop_lmsk,
		32,32, 
		0, 0,
		/*color*/
		BC_BLACK, BC_WHITE
	};

	BlenderCursor[BC_VLOOPCURSOR]=&VLoopCursor;

END_CURSOR_BLOCK	
	

	/********************** TextEdit Cursor ***********************/	
BEGIN_CURSOR_BLOCK
	static char textedit_sbm[]={
		0xe0,  0x03,  0x10,  0x04,  0x60,  0x03,  0x40,  0x01,
		0x40,  0x01,  0x40,  0x01,  0x40,  0x01,  0x40,  0x01,
		0x40,  0x01,  0x40,  0x01,  0x40,  0x01,  0x40,  0x01,
		0x40,  0x01,  0x60,  0x03,  0x10,  0x04,  0xe0,  0x03,
	};

	static char textedit_smsk[]={
		0xe0,  0x03,  0xf0,  0x07,  0xe0,  0x03,  0xc0,  0x01,
		0xc0,  0x01,  0xc0,  0x01,  0xc0,  0x01,  0xc0,  0x01,
		0xc0,  0x01,  0xc0,  0x01,  0xc0,  0x01,  0xc0,  0x01,
		0xc0,  0x01,  0xe0,  0x03,  0xf0,  0x07,  0xe0,  0x03,
	};

	static BCursor TextEditCursor = {
		/*small*/
		textedit_sbm, textedit_smsk,
		16, 16, 
		9,  8,
		/*big*/
		NULL, NULL,
		32,32, 
		15, 15,
		/*color*/
		BC_BLACK, BC_WHITE
	};

	BlenderCursor[BC_TEXTEDITCURSOR]=&TextEditCursor;
END_CURSOR_BLOCK


	/********************** Paintbrush Cursor ***********************/	
BEGIN_CURSOR_BLOCK
	static char paintbrush_sbm[]={

		0x00,  0xe0,  0x00,  0x98,  0x00,  0x44,  0x00,  0x42,
		0x00,  0x21,  0x80,  0x20,  0x40,  0x13,  0x40,  0x17,
		0xa0,  0x0b,  0x98,  0x05,  0x04,  0x02,  0x02,  0x01,
		0x02,  0x01,  0x02,  0x01,  0x81,  0x00,  0x7f,  0x00,



	};

	static char paintbrush_smsk[]={
		0x00,  0xe0,  0x00,  0xf8,  0x00,  0x7c,  0x00,  0x7e,
		0x00,  0x3f,  0x80,  0x3f,  0xc0,  0x1f,  0xc0,  0x1f,
		0xe0,  0x0f,  0xf8,  0x07,  0xfc,  0x03,  0xfe,  0x01,
		0xfe,  0x01,  0xfe,  0x01,  0xff,  0x00,  0x7f,  0x00,


	};

	static BCursor PaintBrushCursor = {
		/*small*/
		paintbrush_sbm, paintbrush_smsk,
		16, 16, 
		0,  15,
		/*big*/
		NULL, NULL,
		32,32, 
		15, 15,
		/*color*/
		BC_BLACK, BC_WHITE
	};

	BlenderCursor[BC_PAINTBRUSHCURSOR]=&PaintBrushCursor;
END_CURSOR_BLOCK


/********************** Hand Cursor ***********************/
BEGIN_CURSOR_BLOCK

static char hand_sbm[]={ 
	0x00,  0x00,  0x00,  0x00,  0x80,  0x01,  0x80,  0x0d,  
	0x98,  0x6d,  0x98,  0x6d,  0xb0,  0x6d,  0xb0,  0x6d,  
	0xe0,  0x6f,  0xe6,  0x7f,  0xee,  0x7f,  0xfc,  0x3f,  
	0xf8,  0x3f,  0xf0,  0x1f,  0xc0,  0x1f,  0xc0,  0x1f,  
};

static char hand_smsk[]={ 
	0x00,  0x00,  0x80,  0x01,  0xc0,  0x0f,  0xd8,  0x7f,  
	0xfc,  0xff,  0xfc,  0xff,  0xf8,  0xff,  0xf8,  0xff,  
	0xf6,  0xff,  0xff,  0xff,  0xff,  0xff,  0xfe,  0x7f,  
	0xfc,  0x7f,  0xf8,  0x3f,  0xf0,  0x3f,  0xe0,  0x3f,  
};


static BCursor HandCursor = {
	/*small*/
	hand_sbm, hand_smsk,
	16, 16, 
	8,  8,
	/*big*/
	NULL, NULL,
	32,32, 
	15, 15,
	/*color*/
	BC_BLACK, BC_WHITE
};

BlenderCursor[BC_HANDCURSOR]=&HandCursor;

END_CURSOR_BLOCK

/********************** NSEW Scroll Cursor ***********************/
BEGIN_CURSOR_BLOCK

static char nsewscroll_sbm[]={ 
	0x00,  0x00,  0x80,  0x01,  0xc0,  0x03,  0xc0,  0x03,  
	0x00,  0x00,  0x00,  0x00,  0x0c,  0x30,  0x0e,  0x70,  
	0x0e,  0x70,  0x0c,  0x30,  0x00,  0x00,  0x00,  0x00,  
	0xc0,  0x03,  0xc0,  0x03,  0x80,  0x01,  0x00,  0x00, 
};

static char nsewscroll_smsk[]={ 
	0x80,  0x01,  0xc0,  0x03,  0xe0,  0x07,  0xe0,  0x07,  
	0xc0,  0x03,  0x0c,  0x30,  0x1e,  0x78,  0x1f,  0xf8,  
	0x1f,  0xf8,  0x1e,  0x78,  0x0c,  0x30,  0xc0,  0x03,  
	0xe0,  0x07,  0xe0,  0x07,  0xc0,  0x03,  0x80,  0x01, 
};


static BCursor NSEWScrollCursor = {
	/*small*/
	nsewscroll_sbm, nsewscroll_smsk,
	16, 16, 
	8, 8,
	/*big*/
	NULL, NULL,
	32,32, 
	15, 15,
	/*color*/
	BC_BLACK, BC_WHITE
};

BlenderCursor[BC_NSEW_SCROLLCURSOR]=&NSEWScrollCursor;

END_CURSOR_BLOCK


/********************** NS Scroll Cursor ***********************/
BEGIN_CURSOR_BLOCK

static char nsscroll_sbm[]={ 
	0x00,  0x00,  0x80,  0x01,  0xc0,  0x03,  0xc0,  0x03,  
	0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  
	0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  
	0xc0,  0x03,  0xc0,  0x03,  0x80,  0x01,  0x00,  0x00,
};

static char nsscroll_smsk[]={ 
	0x80,  0x01,  0xc0,  0x03,  0xe0,  0x07,  0xe0,  0x07,  
	0xc0,  0x03,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  
	0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0xc0,  0x03,  
	0xe0,  0x07,  0xe0,  0x07,  0xc0,  0x03,  0x80,  0x01,
};


static BCursor NSScrollCursor = {
	/*small*/
	nsscroll_sbm, nsscroll_smsk,
	16, 16, 
	8, 8,
	/*big*/
	NULL, NULL,
	32,32, 
	15, 15,
	/*color*/
	BC_BLACK, BC_WHITE
};

BlenderCursor[BC_NS_SCROLLCURSOR]=&NSScrollCursor;

END_CURSOR_BLOCK


/********************** EW Scroll Cursor ***********************/
BEGIN_CURSOR_BLOCK

static char ewscroll_sbm[]={ 
	0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  
	0x00,  0x00,  0x00,  0x00,  0x0c,  0x30,  0x0e,  0x70,  
	0x0e,  0x70,  0x0c,  0x30,  0x00,  0x00,  0x00,  0x00,  
	0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
};

static char ewscroll_smsk[]={ 
	0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  
	0x00,  0x00,  0x0c,  0x30,  0x1e,  0x78,  0x1f,  0xf8,  
	0x1f,  0xf8,  0x1e,  0x78,  0x0c,  0x30,  0x00,  0x00,  
	0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,  0x00,
};


static BCursor EWScrollCursor = {
	/*small*/
	ewscroll_sbm, ewscroll_smsk,
	16, 16, 
	8, 8,
	/*big*/
	NULL, NULL,
	32,32, 
	15, 15,
	/*color*/
	BC_BLACK, BC_WHITE
};

BlenderCursor[BC_EW_SCROLLCURSOR]=&EWScrollCursor;

END_CURSOR_BLOCK

/********************** Eyedropper Cursor ***********************/
BEGIN_CURSOR_BLOCK

static char eyedropper_sbm[]={ 
	0x00,  0x30,  0x00,  0x48,  0x00,  0x85,  0x80,  0x82,  
	0x40,  0x40,  0x80,  0x20,  0x40,  0x11,  0xa0,  0x23,  
	0xd0,  0x15,  0xe8,  0x0a,  0x74,  0x01,  0xb4,  0x00,  
	0x4a,  0x00,  0x35,  0x00,  0x08,  0x00,  0x04,  0x00,
};

static char eyedropper_smsk[]={ 
	0x00,  0x30,  0x00,  0x78,  0x00,  0xfd,  0x80,  0xff,  
	0xc0,  0x7f,  0x80,  0x3f,  0xc0,  0x1f,  0xe0,  0x3f,  
	0xf0,  0x1f,  0xf8,  0x0b,  0xfc,  0x01,  0xfc,  0x00,  
	0x7e,  0x00,  0x3f,  0x00,  0x0c,  0x00,  0x04,  0x00, 
};


static BCursor EyedropperCursor = {
	/*small*/
	eyedropper_sbm, eyedropper_smsk,
	16, 16, 
	1, 15,
	/*big*/
	NULL, NULL,
	32,32, 
	15, 15,
	/*color*/
	BC_BLACK, BC_WHITE
};

BlenderCursor[BC_EYEDROPPER_CURSOR]=&EyedropperCursor;

END_CURSOR_BLOCK

/********************** Put the cursors in the array ***********************/
	


}



