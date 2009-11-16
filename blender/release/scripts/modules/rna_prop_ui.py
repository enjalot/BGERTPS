# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>

import bpy

EVIL_PROP = "act_property"
EVIL_PROP_PATH = EVIL_PROP + '_path'
EVIL_PROP_VALUE = EVIL_PROP + '_value'
EVIL_PROP_PROP = EVIL_PROP + '_prop'
EVIL_PROP_PROP_ORIG = EVIL_PROP + '_prop_orig'

# nasty!, use a scene property to store the active edit item
def evil_prop_init():
    Scene = bpy.types.Scene
    if EVIL_PROP_PROP_ORIG not in Scene.bl_rna.properties:
        Scene.StringProperty(attr=EVIL_PROP_PATH)
        Scene.StringProperty(attr=EVIL_PROP_VALUE)
        Scene.StringProperty(attr=EVIL_PROP_PROP)
        Scene.StringProperty(attr=EVIL_PROP_PROP_ORIG)


def draw(layout, context, context_member):
    
    def assign_props(prop, val, key):
        prop.path = context_member
        prop.property = key
        
        try:
            prop.value = str(val)
        except:
            pass
    
    rna_item = eval("context." + context_member)
    
    evil_prop_init()

    scene = context.scene
    
    global_path = getattr(scene, EVIL_PROP_PATH)
    global_value = getattr(scene, EVIL_PROP_VALUE)
    global_prop = getattr(scene, EVIL_PROP_PROP)
    global_prop_orig = getattr(scene, EVIL_PROP_PROP_ORIG)
    
    # print((global_path, global_value, global_prop, global_prop_orig))

    items = rna_item.items()
    items.sort()
    
    row = layout.row()
    props = row.itemO("wm.properties_add", properties=True, text="Add")
    props.path = context_member
    
    for key, val in items:
        row = layout.row()
        
        try:
            val_draw = val = val.convert_to_pyobject()
        except:
            if type(val)==str:
                val_draw = '"' + val + '"'
            else:
                val_draw = val
        
        box = row.box()
        
        
        if key == global_prop_orig and context_member == global_path:
            split = box.split(percentage=0.7)
            
            col = split.row()
            col.itemR(scene, EVIL_PROP_PROP)
            col.itemR(scene, EVIL_PROP_VALUE)
            
            col = split.column()
            prop = col.itemO("wm.properties_edit_end", properties=True, text="done")
            assign_props(prop, val, key)
            
        else:
            split = box.split(percentage=0.5)
            col = split.column()
            col.itemL(text="%s = %s" % (key, str(val_draw)))
            
            col = split.column()
            prop = col.itemO("wm.properties_edit_begin", properties=True, text="edit")
            assign_props(prop, val, key)
            
            col = split.column()
            prop = col.itemO("wm.properties_remove", properties=True, text="del")
            assign_props(prop, val, key)
    

from bpy.props import *


rna_path = StringProperty(name="Property Edit",
    description="Property path edit", maxlen=1024, default="")

rna_value = StringProperty(name="Property Value",
    description="Property value edit", maxlen=1024, default="")

rna_property = StringProperty(name="Property Name",
    description="Property name edit", maxlen=1024, default="")

class WM_OT_properties_edit_begin(bpy.types.Operator):
    '''Internal use (edit a property path)'''
    bl_idname = "wm.properties_edit_begin"
    bl_label = "Edit Property"

    path = rna_path
    value = rna_value
    property = rna_property

    def execute(self, context):
        scene = context.scene
        
        setattr(scene, EVIL_PROP_PATH, self.path)
        setattr(scene, EVIL_PROP_VALUE, self.value)
        setattr(scene, EVIL_PROP_PROP, self.property)
        setattr(scene, EVIL_PROP_PROP_ORIG, self.property)
        
        return ('FINISHED',)


class WM_OT_properties_edit_end(bpy.types.Operator):
    '''Internal use (edit a property path)'''
    bl_idname = "wm.properties_edit_end"
    bl_label = "Edit Property"

    path = rna_path
    value = rna_value
    property = rna_property

    def execute(self, context):
        
        scene = context.scene
        global_path = getattr(scene, EVIL_PROP_PATH)
        global_value = getattr(scene, EVIL_PROP_VALUE)
        global_prop = getattr(scene, EVIL_PROP_PROP)
        
        setattr(scene, EVIL_PROP_PATH, "")
        setattr(scene, EVIL_PROP_VALUE, "")
        setattr(scene, EVIL_PROP_PROP, "")
        setattr(scene, EVIL_PROP_PROP_ORIG, "")
        
        try:
            value = eval(global_value)
        except:
            value = '"' + global_value + '"' # keep as a string
        
        
        # First remove
        exec_str = "del context.%s['%s']" % (global_path, self.property)
        # print(exec_str)
        exec(exec_str)
        
        
        # Reassign
        exec_str = "context.%s['%s'] = %s" % (global_path, global_prop, value)
        # print(exec_str)
        exec(exec_str)
        
        return ('FINISHED',)


class WM_OT_properties_add(bpy.types.Operator):
    '''Internal use (edit a property path)'''
    bl_idname = "wm.properties_add"
    bl_label = "Add Property"

    path = rna_path

    def execute(self, context):
        item = eval("context.%s" % self.path)
        
        def unique_name(names):
            prop = 'prop'
            prop_new = prop
            i = 1
            while prop_new in names:
                prop_new = prop + str(i)
                i+=1
        
            return prop_new
        
        property = unique_name(item.keys())
        
        item[property] = 1.0
        return ('FINISHED',)
        
class WM_OT_properties_remove(bpy.types.Operator):
    '''Internal use (edit a property path)'''
    bl_idname = "wm.properties_remove"
    bl_label = "Add Property"

    path = rna_path
    property = rna_property

    def execute(self, context):
        item = eval("context.%s" % self.path)
        del item[self.property]
        return ('FINISHED',)        
