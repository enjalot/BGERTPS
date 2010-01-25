# Blender.Mathutils module and its subtypes

"""
The Blender.Mathutils submodule.

Mathutils
=========
(when accessing it from the Game Engine use Mathutils instead of Blender.Mathutils)

This module provides access to matrices, eulers, quaternions and vectors.

Example::
  import Blender
  from Blender import Mathutils
  from Blender.Mathutils import *

  vec = Vector([1,2,3])
  mat = RotationMatrix(90, 4, 'x')
  matT = TranslationMatrix(vec)

  matTotal = mat * matT
  matTotal.invert()

  mat3 = matTotal.rotationPart
  quat1 = mat.to_quat()
  quat2 = mat3.to_quat()

  angle = DifferenceQuats(quat1, quat2)
  print angle  
"""

class Vector:
  """
  The Vector object
  =================
    This object gives access to Vectors in Blender.
  @group Axises: x, y, z, w
  @ivar x: The x value.
  @ivar y: The y value.
  @ivar z: The z value (if any).
  @ivar w: The w value (if any).
  @ivar length: The magnitude of the vector.
  @ivar magnitude: This is a synonym for length.
  @ivar wrapped: Whether or not this item is wrapped data
  @note: Comparison operators can be done on Vector classes:
      - >, >=, <, <= test the vector magnitude
      - ==, != test vector values e.g. 1,2,3 != 1,2,4 even if they are the same length
  @note: Math can be performed on Vector classes
      - vec + vec
      - vec - vec
      - vec * float/int
      - vec * matrix
      - vec * vec
      - vec * quat
      - -vec
  @note: You can access a vector object like a sequence
      - x = vector[0]
      - vec_a[:] vec_b
      - vec2d[:] vec3d[:2]
  @note: Vectors support 'swizzle' operations
      - vec.xyz = vec.zyx
      - vec.xy = vec.zw
      - vec.xxy = vec.wzz
      - vec.yzyz = vec.yxyx

      See U{http://en.wikipedia.org/wiki/Swizzling_(computer_graphics)}
  
  @attention: Vector data can be wrapped or non-wrapped. When a object is wrapped it
  means that the object will give you direct access to the data inside of blender. Modification
  of this object will directly change the data inside of blender. To copy a wrapped object
  you need to use the object's constructor. If you copy and object by assignment you will not get
  a second copy but a second reference to the same data. Only certain functions will return 
  wrapped data. This will be indicated in the method description.
  Example::
      wrappedObject = Object.getAttribute() #this is wrapped data
      print wrappedObject.wrapped #prints 'True'
      copyOfObject = wrappedObject.copy() #creates a copy of the object
      secondPointer = wrappedObject #creates a second pointer to the same data
      print wrappedObject.attribute #prints '5'
      secondPointer.attribute = 10
      print wrappedObject.attribute #prints '10'
      print copyOfObject.attribute #prints '5'
  """

  def __init__(list = None):
    """
    Create a new 2d, 3d, or 4d Vector object from a list of floating point numbers.
    @note: that python uses higher precission floating point numbers, so values assigned to a vector may have some rounding error.
    

    Example::
      v = Vector(1,0,0)
      v = Vector(myVec)
      v = Vector(list)
    @type list: PyList of float or int
    @param list: The list of values for the Vector object. Can be a sequence or raw numbers.
    Must be 2, 3, or 4 values. The list is mapped to the parameters as [x,y,z,w].
    @rtype: Vector object.
    @return: It depends wheter a parameter was passed:
        - (list): Vector object initialized with the given values;
        - ():     An empty 3 dimensional vector.
    """

class Euler:
  """
  The Euler object
  ================
    This object gives access to Eulers in Blender.
  @group Axises: x, y, z
  @ivar x: The heading value in degrees.
  @ivar y: The pitch value in degrees.
  @ivar z: The roll value in degrees.
  @ivar wrapped: Whether or not this object is wrapping data directly
  @note: You can access a euler object like a sequence
      - x = euler[0]
  @note: Comparison operators can be done:
      - ==, != test numeric values within epsilon
  @attention: Euler data can be wrapped or non-wrapped. When a object is wrapped it
  means that the object will give you direct access to the data inside of blender. Modification
  of this object will directly change the data inside of blender. To copy a wrapped object
  you need to use the object's constructor. If you copy and object by assignment you will not get
  a second copy but a second reference to the same data. Only certain functions will return 
  wrapped data. This will be indicated in the method description.
  Example::
      wrappedObject = Object.getAttribute() #this is wrapped data
      print wrappedObject.wrapped #prints 'True'
      copyOfObject = wrappedObject.copy() #creates a copy of the object
      secondPointer = wrappedObject #creates a second pointer to the same data
      print wrappedObject.attribute #prints '5'
      secondPointer.attribute = 10
      print wrappedObject.attribute #prints '10'
      print copyOfObject.attribute #prints '5'
  """

  def __init__(list = None):
    """
    Create a new euler object.

    Example::
      euler = Euler(45,0,0)
      euler = Euler(myEuler)
      euler = Euler(sequence)
    @type list: PyList of float/int
    @param list: 3d list to initialize euler
    @rtype: Euler object
    @return: Euler representing heading, pitch, bank.
    @note: Values are in degrees.
    """

class Quaternion:
  """
  The Quaternion object
  =====================
    This object gives access to Quaternions in Blender.
  @group Axises: x, y, z, w
  @ivar w: The w value.
  @ivar x: The x value.
  @ivar y: The y value.
  @ivar z: The z value.
  @ivar wrapped: Wether or not this object wraps data directly
  @ivar magnitude: The magnitude of the quaternion.
  @ivar axis: Vector representing the axis of rotation.
  @ivar angle: A scalar representing the amount of rotation
  in degrees.
  @note: Comparison operators can be done:
      - ==, != test numeric values within epsilon
  @note: Math can be performed on Quaternion classes
      - quat + quat
      - quat - quat 
      - quat * float/int
      - quat * vec
      - quat * quat
  @note: You can access a quaternion object like a sequence
      - x = quat[0]
  @attention: Quaternion data can be wrapped or non-wrapped. When a object is wrapped it
  means that the object will give you direct access to the data inside of blender. Modification
  of this object will directly change the data inside of blender. To copy a wrapped object
  you need to use the object's constructor. If you copy and object by assignment you will not get
  a second copy but a second reference to the same data. Only certain functions will return 
  wrapped data. This will be indicated in the method description.
  Example::
      wrappedObject = Object.getAttribute() #this is wrapped data
      print wrappedObject.wrapped #prints 'True'
      copyOfObject = wrappedObject.copy() #creates a copy of the object
      secondPointer = wrappedObject #creates a second pointer to the same data
      print wrappedObject.attribute #prints '5'
      secondPointer.attribute = 10
      print wrappedObject.attribute #prints '10'
      print copyOfObject.attribute #prints '5'
  """

  def __init__(list, angle = None):
    """  
    Create a new quaternion object from initialized values.

    Example::
      quat = Quaternion(1,2,3,4)
      quat = Quaternion(axis, angle)
    quat = Quaternion()
    quat = Quaternion(180, list)

    @type list: PyList of int/float
    @param list: A 3d or 4d list to initialize quaternion.
        4d if intializing [w,x,y,z], 3d if used as an axis of rotation.
    @type angle: float (optional)
    @param angle: An arbitrary rotation amount around 'list'.
        List is used as an axis of rotation in this case.
    @rtype: New quaternion object.
    @return: It depends wheter a parameter was passed:
        - (list/angle): Quaternion object initialized with the given values;
        - ():     An identity 4 dimensional quaternion.
    """

class Matrix:
  """
  The Matrix Object
  =================
    This object gives access to Matrices in Blender.
  @ivar rowSize: The row size of the matrix.
  @ivar colSize: The column size of the matrix.
  @ivar wrapped: Whether or not this object wrapps internal data
  @note: Math can be performed on Matrix classes
      - mat + mat 
      - mat - mat 
      - mat * float/int
      - mat * vec
      - mat * mat 
  @note: Comparison operators can be done:
      - ==, != test numeric values within epsilon
  @note: You can access a quaternion object like a 2d sequence
      - x = matrix[0][1]
      - vector = matrix[2]
  @attention: Quaternion data can be wrapped or non-wrapped. When a object is wrapped it
  means that the object will give you direct access to the data inside of blender. Modification
  of this object will directly change the data inside of blender. To copy a wrapped object
  you need to use the object's constructor. If you copy and object by assignment you will not get
  a second copy but a second reference to the same data. Only certain functions will return 
  wrapped data. This will be indicated in the method description.
  Example::
      wrappedObject = Object.getAttribute() #this is wrapped data
      print wrappedObject.wrapped #prints 'True'
      copyOfObject = wrappedObject.copy() #creates a copy of the object
      secondPointer = wrappedObject #creates a second pointer to the same data
      print wrappedObject.attribute #prints '5'
      secondPointer.attribute = 10
      print wrappedObject.attribute #prints '10'
      print copyOfObject.attribute #prints '5'
  """

  def __init__(list1 = None, list2 = None, list3 = None, list4 = None):
    """  
    Create a new matrix object from initialized values.

    Example::
      matrix = Matrix([1,1,1],[0,1,0],[1,0,0])
      matrix = Matrix(mat)
      matrix = Matrix(seq1, seq2, vector)

    @type list1: PyList of int/float
    @param list1: A 2d,3d or 4d list.
    @type list2: PyList of int/float
    @param list2: A 2d,3d or 4d list.
    @type list3: PyList of int/float
    @param list3: A 2d,3d or 4d list.
    @type list4: PyList of int/float
    @param list4: A 2d,3d or 4d list.
    @rtype: New matrix object.
    @return: It depends wheter a parameter was passed:
        - (list1, etc.): Matrix object initialized with the given values;
        - ():     An empty 3 dimensional matrix.
    """