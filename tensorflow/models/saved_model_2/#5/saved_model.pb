��	
�!�!
W
AddN
inputs"T*N
sum"T"
Nint(0"!
Ttype:
2	��
�
ArgMax

input"T
	dimension"Tidx
output"output_type" 
Ttype:
2	"
Tidxtype0:
2	"
output_typetype0	:
2	
�
AsString

input"T

output"
Ttype:
2		
"
	precisionint���������"

scientificbool( "
shortestbool( "
widthint���������"
fillstring 
x
Assign
ref"T�

value"T

output_ref"T�"	
Ttype"
validate_shapebool("
use_lockingbool(�
B
AssignVariableOp
resource
value"dtype"
dtypetype�
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
N
Cast	
x"SrcT	
y"DstT"
SrcTtype"
DstTtype"
Truncatebool( 
h
ConcatV2
values"T*N
axis"Tidx
output"T"
Nint(0"	
Ttype"
Tidxtype0:
2	
8
Const
output"dtype"
valuetensor"
dtypetype
B
Equal
x"T
y"T
z
"
Ttype:
2	
�
W

ExpandDims

input"T
dim"Tdim
output"T"	
Ttype"
Tdimtype0:
2	
.
Identity

input"T
output"T"	
Ttype
?
	LessEqual
x"T
y"T
z
"
Ttype:
2	
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
N
Merge
inputs"T*N
output"T
value_index"	
Ttype"
Nint(0
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(�
=
Mul
x"T
y"T
z"T"
Ttype:
2	�

NoOp
E
NotEqual
x"T
y"T
z
"
Ttype:
2	
�
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
�
ParseExample

serialized	
names
sparse_keys*Nsparse

dense_keys*Ndense
dense_defaults2Tdense
sparse_indices	*Nsparse
sparse_values2sparse_types
sparse_shapes	*Nsparse
dense_values2Tdense"
Nsparseint("
Ndenseint("%
sparse_types
list(type)(:
2	"
Tdense
list(type)(:
2	"
dense_shapeslist(shape)(
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
X
PlaceholderWithDefault
input"dtype
output"dtype"
dtypetype"
shapeshape
a
Range
start"Tidx
limit"Tidx
delta"Tidx
output"Tidx"
Tidxtype0:	
2	
@
ReadVariableOp
resource
value"dtype"
dtypetype�
>
RealDiv
x"T
y"T
z"T"
Ttype:
2	
[
Reshape
tensor"T
shape"Tshape
output"T"	
Ttype"
Tshapetype0:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
P
ScalarSummary
tags
values"T
summary"
Ttype:
2	
P
Shape

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
H
ShardedFilename
basename	
shard

num_shards
filename
0
Sigmoid
x"T
y"T"
Ttype:

2
9
Softmax
logits"T
softmax"T"
Ttype:
2
�
StridedSlice

input"T
begin"Index
end"Index
strides"Index
output"T"	
Ttype"
Indextype:
2	"

begin_maskint "
end_maskint "
ellipsis_maskint "
new_axis_maskint "
shrink_axis_maskint 
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
:
Sub
x"T
y"T
z"T"
Ttype:
2	
�
Sum

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
M
Switch	
data"T
pred

output_false"T
output_true"T"	
Ttype
c
Tile

input"T
	multiples"
Tmultiples
output"T"	
Ttype"

Tmultiplestype0:
2	
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape�
9
VarIsInitializedOp
resource
is_initialized
�
s

VariableV2
ref"dtype�"
shapeshape"
dtypetype"
	containerstring "
shared_namestring �
&
	ZerosLike
x"T
y"T"	
Ttype"serve*1.14.02unknown8ɧ

global_step/Initializer/zerosConst*
value	B	 R *
_class
loc:@global_step*
dtype0	*
_output_shapes
: 
k
global_step
VariableV2*
shape: *
_class
loc:@global_step*
dtype0	*
_output_shapes
: 
�
global_step/AssignAssignglobal_stepglobal_step/Initializer/zeros*
T0	*
_class
loc:@global_step*
_output_shapes
: 
j
global_step/readIdentityglobal_step*
T0	*
_class
loc:@global_step*
_output_shapes
: 
o
input_example_tensorPlaceholder*
shape:���������*
dtype0*#
_output_shapes
:���������
U
ParseExample/ConstConst*
valueB *
dtype0*
_output_shapes
: 
W
ParseExample/Const_1Const*
valueB *
dtype0*
_output_shapes
: 
W
ParseExample/Const_2Const*
valueB *
dtype0*
_output_shapes
: 
W
ParseExample/Const_3Const*
valueB *
dtype0*
_output_shapes
: 
W
ParseExample/Const_4Const*
valueB *
dtype0*
_output_shapes
: 
W
ParseExample/Const_5Const*
valueB *
dtype0*
_output_shapes
: 
b
ParseExample/ParseExample/namesConst*
valueB *
dtype0*
_output_shapes
: 
r
&ParseExample/ParseExample/dense_keys_0Const*
valueB Bdirection.x*
dtype0*
_output_shapes
: 
r
&ParseExample/ParseExample/dense_keys_1Const*
valueB Bdirection.y*
dtype0*
_output_shapes
: 
r
&ParseExample/ParseExample/dense_keys_2Const*
valueB Bdirection.z*
dtype0*
_output_shapes
: 
o
&ParseExample/ParseExample/dense_keys_3Const*
valueB Borigin.x*
dtype0*
_output_shapes
: 
o
&ParseExample/ParseExample/dense_keys_4Const*
valueB Borigin.y*
dtype0*
_output_shapes
: 
o
&ParseExample/ParseExample/dense_keys_5Const*
valueB Borigin.z*
dtype0*
_output_shapes
: 
�
ParseExample/ParseExampleParseExampleinput_example_tensorParseExample/ParseExample/names&ParseExample/ParseExample/dense_keys_0&ParseExample/ParseExample/dense_keys_1&ParseExample/ParseExample/dense_keys_2&ParseExample/ParseExample/dense_keys_3&ParseExample/ParseExample/dense_keys_4&ParseExample/ParseExample/dense_keys_5ParseExample/ConstParseExample/Const_1ParseExample/Const_2ParseExample/Const_3ParseExample/Const_4ParseExample/Const_5*
Nsparse *6
dense_shapes&
$::::::*
sparse_types
 *
Tdense

2*
Ndense*�
_output_shapest
r:���������:���������:���������:���������:���������:���������
�
@linear/linear_model/direction.x/weights/part_0/Initializer/zerosConst*
valueB*    *A
_class7
53loc:@linear/linear_model/direction.x/weights/part_0*
dtype0*
_output_shapes

:
�
.linear/linear_model/direction.x/weights/part_0VarHandleOp*
shape
:*?
shared_name0.linear/linear_model/direction.x/weights/part_0*A
_class7
53loc:@linear/linear_model/direction.x/weights/part_0*
dtype0*
_output_shapes
: 
�
Olinear/linear_model/direction.x/weights/part_0/IsInitialized/VarIsInitializedOpVarIsInitializedOp.linear/linear_model/direction.x/weights/part_0*
_output_shapes
: 
�
5linear/linear_model/direction.x/weights/part_0/AssignAssignVariableOp.linear/linear_model/direction.x/weights/part_0@linear/linear_model/direction.x/weights/part_0/Initializer/zeros*A
_class7
53loc:@linear/linear_model/direction.x/weights/part_0*
dtype0
�
Blinear/linear_model/direction.x/weights/part_0/Read/ReadVariableOpReadVariableOp.linear/linear_model/direction.x/weights/part_0*A
_class7
53loc:@linear/linear_model/direction.x/weights/part_0*
dtype0*
_output_shapes

:
�
@linear/linear_model/direction.y/weights/part_0/Initializer/zerosConst*
valueB*    *A
_class7
53loc:@linear/linear_model/direction.y/weights/part_0*
dtype0*
_output_shapes

:
�
.linear/linear_model/direction.y/weights/part_0VarHandleOp*
shape
:*?
shared_name0.linear/linear_model/direction.y/weights/part_0*A
_class7
53loc:@linear/linear_model/direction.y/weights/part_0*
dtype0*
_output_shapes
: 
�
Olinear/linear_model/direction.y/weights/part_0/IsInitialized/VarIsInitializedOpVarIsInitializedOp.linear/linear_model/direction.y/weights/part_0*
_output_shapes
: 
�
5linear/linear_model/direction.y/weights/part_0/AssignAssignVariableOp.linear/linear_model/direction.y/weights/part_0@linear/linear_model/direction.y/weights/part_0/Initializer/zeros*A
_class7
53loc:@linear/linear_model/direction.y/weights/part_0*
dtype0
�
Blinear/linear_model/direction.y/weights/part_0/Read/ReadVariableOpReadVariableOp.linear/linear_model/direction.y/weights/part_0*A
_class7
53loc:@linear/linear_model/direction.y/weights/part_0*
dtype0*
_output_shapes

:
�
@linear/linear_model/direction.z/weights/part_0/Initializer/zerosConst*
valueB*    *A
_class7
53loc:@linear/linear_model/direction.z/weights/part_0*
dtype0*
_output_shapes

:
�
.linear/linear_model/direction.z/weights/part_0VarHandleOp*
shape
:*?
shared_name0.linear/linear_model/direction.z/weights/part_0*A
_class7
53loc:@linear/linear_model/direction.z/weights/part_0*
dtype0*
_output_shapes
: 
�
Olinear/linear_model/direction.z/weights/part_0/IsInitialized/VarIsInitializedOpVarIsInitializedOp.linear/linear_model/direction.z/weights/part_0*
_output_shapes
: 
�
5linear/linear_model/direction.z/weights/part_0/AssignAssignVariableOp.linear/linear_model/direction.z/weights/part_0@linear/linear_model/direction.z/weights/part_0/Initializer/zeros*A
_class7
53loc:@linear/linear_model/direction.z/weights/part_0*
dtype0
�
Blinear/linear_model/direction.z/weights/part_0/Read/ReadVariableOpReadVariableOp.linear/linear_model/direction.z/weights/part_0*A
_class7
53loc:@linear/linear_model/direction.z/weights/part_0*
dtype0*
_output_shapes

:
�
=linear/linear_model/origin.x/weights/part_0/Initializer/zerosConst*
valueB*    *>
_class4
20loc:@linear/linear_model/origin.x/weights/part_0*
dtype0*
_output_shapes

:
�
+linear/linear_model/origin.x/weights/part_0VarHandleOp*
shape
:*<
shared_name-+linear/linear_model/origin.x/weights/part_0*>
_class4
20loc:@linear/linear_model/origin.x/weights/part_0*
dtype0*
_output_shapes
: 
�
Llinear/linear_model/origin.x/weights/part_0/IsInitialized/VarIsInitializedOpVarIsInitializedOp+linear/linear_model/origin.x/weights/part_0*
_output_shapes
: 
�
2linear/linear_model/origin.x/weights/part_0/AssignAssignVariableOp+linear/linear_model/origin.x/weights/part_0=linear/linear_model/origin.x/weights/part_0/Initializer/zeros*>
_class4
20loc:@linear/linear_model/origin.x/weights/part_0*
dtype0
�
?linear/linear_model/origin.x/weights/part_0/Read/ReadVariableOpReadVariableOp+linear/linear_model/origin.x/weights/part_0*>
_class4
20loc:@linear/linear_model/origin.x/weights/part_0*
dtype0*
_output_shapes

:
�
=linear/linear_model/origin.y/weights/part_0/Initializer/zerosConst*
valueB*    *>
_class4
20loc:@linear/linear_model/origin.y/weights/part_0*
dtype0*
_output_shapes

:
�
+linear/linear_model/origin.y/weights/part_0VarHandleOp*
shape
:*<
shared_name-+linear/linear_model/origin.y/weights/part_0*>
_class4
20loc:@linear/linear_model/origin.y/weights/part_0*
dtype0*
_output_shapes
: 
�
Llinear/linear_model/origin.y/weights/part_0/IsInitialized/VarIsInitializedOpVarIsInitializedOp+linear/linear_model/origin.y/weights/part_0*
_output_shapes
: 
�
2linear/linear_model/origin.y/weights/part_0/AssignAssignVariableOp+linear/linear_model/origin.y/weights/part_0=linear/linear_model/origin.y/weights/part_0/Initializer/zeros*>
_class4
20loc:@linear/linear_model/origin.y/weights/part_0*
dtype0
�
?linear/linear_model/origin.y/weights/part_0/Read/ReadVariableOpReadVariableOp+linear/linear_model/origin.y/weights/part_0*>
_class4
20loc:@linear/linear_model/origin.y/weights/part_0*
dtype0*
_output_shapes

:
�
=linear/linear_model/origin.z/weights/part_0/Initializer/zerosConst*
valueB*    *>
_class4
20loc:@linear/linear_model/origin.z/weights/part_0*
dtype0*
_output_shapes

:
�
+linear/linear_model/origin.z/weights/part_0VarHandleOp*
shape
:*<
shared_name-+linear/linear_model/origin.z/weights/part_0*>
_class4
20loc:@linear/linear_model/origin.z/weights/part_0*
dtype0*
_output_shapes
: 
�
Llinear/linear_model/origin.z/weights/part_0/IsInitialized/VarIsInitializedOpVarIsInitializedOp+linear/linear_model/origin.z/weights/part_0*
_output_shapes
: 
�
2linear/linear_model/origin.z/weights/part_0/AssignAssignVariableOp+linear/linear_model/origin.z/weights/part_0=linear/linear_model/origin.z/weights/part_0/Initializer/zeros*>
_class4
20loc:@linear/linear_model/origin.z/weights/part_0*
dtype0
�
?linear/linear_model/origin.z/weights/part_0/Read/ReadVariableOpReadVariableOp+linear/linear_model/origin.z/weights/part_0*>
_class4
20loc:@linear/linear_model/origin.z/weights/part_0*
dtype0*
_output_shapes

:
�
9linear/linear_model/bias_weights/part_0/Initializer/zerosConst*
valueB*    *:
_class0
.,loc:@linear/linear_model/bias_weights/part_0*
dtype0*
_output_shapes
:
�
'linear/linear_model/bias_weights/part_0VarHandleOp*
shape:*8
shared_name)'linear/linear_model/bias_weights/part_0*:
_class0
.,loc:@linear/linear_model/bias_weights/part_0*
dtype0*
_output_shapes
: 
�
Hlinear/linear_model/bias_weights/part_0/IsInitialized/VarIsInitializedOpVarIsInitializedOp'linear/linear_model/bias_weights/part_0*
_output_shapes
: 
�
.linear/linear_model/bias_weights/part_0/AssignAssignVariableOp'linear/linear_model/bias_weights/part_09linear/linear_model/bias_weights/part_0/Initializer/zeros*:
_class0
.,loc:@linear/linear_model/bias_weights/part_0*
dtype0
�
;linear/linear_model/bias_weights/part_0/Read/ReadVariableOpReadVariableOp'linear/linear_model/bias_weights/part_0*:
_class0
.,loc:@linear/linear_model/bias_weights/part_0*
dtype0*
_output_shapes
:
�
?linear/linear_model/linear_model/linear_model/direction.x/ShapeShapeParseExample/ParseExample*
T0*
_output_shapes
:
�
Mlinear/linear_model/linear_model/linear_model/direction.x/strided_slice/stackConst*
valueB: *
dtype0*
_output_shapes
:
�
Olinear/linear_model/linear_model/linear_model/direction.x/strided_slice/stack_1Const*
valueB:*
dtype0*
_output_shapes
:
�
Olinear/linear_model/linear_model/linear_model/direction.x/strided_slice/stack_2Const*
valueB:*
dtype0*
_output_shapes
:
�
Glinear/linear_model/linear_model/linear_model/direction.x/strided_sliceStridedSlice?linear/linear_model/linear_model/linear_model/direction.x/ShapeMlinear/linear_model/linear_model/linear_model/direction.x/strided_slice/stackOlinear/linear_model/linear_model/linear_model/direction.x/strided_slice/stack_1Olinear/linear_model/linear_model/linear_model/direction.x/strided_slice/stack_2*
shrink_axis_mask*
Index0*
T0*
_output_shapes
: 
�
Ilinear/linear_model/linear_model/linear_model/direction.x/Reshape/shape/1Const*
value	B :*
dtype0*
_output_shapes
: 
�
Glinear/linear_model/linear_model/linear_model/direction.x/Reshape/shapePackGlinear/linear_model/linear_model/linear_model/direction.x/strided_sliceIlinear/linear_model/linear_model/linear_model/direction.x/Reshape/shape/1*
T0*
N*
_output_shapes
:
�
Alinear/linear_model/linear_model/linear_model/direction.x/ReshapeReshapeParseExample/ParseExampleGlinear/linear_model/linear_model/linear_model/direction.x/Reshape/shape*
T0*'
_output_shapes
:���������
�
6linear/linear_model/direction.x/weights/ReadVariableOpReadVariableOp.linear/linear_model/direction.x/weights/part_0*
dtype0*
_output_shapes

:
�
'linear/linear_model/direction.x/weightsIdentity6linear/linear_model/direction.x/weights/ReadVariableOp*
T0*
_output_shapes

:
�
Flinear/linear_model/linear_model/linear_model/direction.x/weighted_sumMatMulAlinear/linear_model/linear_model/linear_model/direction.x/Reshape'linear/linear_model/direction.x/weights*
T0*'
_output_shapes
:���������
�
?linear/linear_model/linear_model/linear_model/direction.y/ShapeShapeParseExample/ParseExample:1*
T0*
_output_shapes
:
�
Mlinear/linear_model/linear_model/linear_model/direction.y/strided_slice/stackConst*
valueB: *
dtype0*
_output_shapes
:
�
Olinear/linear_model/linear_model/linear_model/direction.y/strided_slice/stack_1Const*
valueB:*
dtype0*
_output_shapes
:
�
Olinear/linear_model/linear_model/linear_model/direction.y/strided_slice/stack_2Const*
valueB:*
dtype0*
_output_shapes
:
�
Glinear/linear_model/linear_model/linear_model/direction.y/strided_sliceStridedSlice?linear/linear_model/linear_model/linear_model/direction.y/ShapeMlinear/linear_model/linear_model/linear_model/direction.y/strided_slice/stackOlinear/linear_model/linear_model/linear_model/direction.y/strided_slice/stack_1Olinear/linear_model/linear_model/linear_model/direction.y/strided_slice/stack_2*
shrink_axis_mask*
Index0*
T0*
_output_shapes
: 
�
Ilinear/linear_model/linear_model/linear_model/direction.y/Reshape/shape/1Const*
value	B :*
dtype0*
_output_shapes
: 
�
Glinear/linear_model/linear_model/linear_model/direction.y/Reshape/shapePackGlinear/linear_model/linear_model/linear_model/direction.y/strided_sliceIlinear/linear_model/linear_model/linear_model/direction.y/Reshape/shape/1*
T0*
N*
_output_shapes
:
�
Alinear/linear_model/linear_model/linear_model/direction.y/ReshapeReshapeParseExample/ParseExample:1Glinear/linear_model/linear_model/linear_model/direction.y/Reshape/shape*
T0*'
_output_shapes
:���������
�
6linear/linear_model/direction.y/weights/ReadVariableOpReadVariableOp.linear/linear_model/direction.y/weights/part_0*
dtype0*
_output_shapes

:
�
'linear/linear_model/direction.y/weightsIdentity6linear/linear_model/direction.y/weights/ReadVariableOp*
T0*
_output_shapes

:
�
Flinear/linear_model/linear_model/linear_model/direction.y/weighted_sumMatMulAlinear/linear_model/linear_model/linear_model/direction.y/Reshape'linear/linear_model/direction.y/weights*
T0*'
_output_shapes
:���������
�
?linear/linear_model/linear_model/linear_model/direction.z/ShapeShapeParseExample/ParseExample:2*
T0*
_output_shapes
:
�
Mlinear/linear_model/linear_model/linear_model/direction.z/strided_slice/stackConst*
valueB: *
dtype0*
_output_shapes
:
�
Olinear/linear_model/linear_model/linear_model/direction.z/strided_slice/stack_1Const*
valueB:*
dtype0*
_output_shapes
:
�
Olinear/linear_model/linear_model/linear_model/direction.z/strided_slice/stack_2Const*
valueB:*
dtype0*
_output_shapes
:
�
Glinear/linear_model/linear_model/linear_model/direction.z/strided_sliceStridedSlice?linear/linear_model/linear_model/linear_model/direction.z/ShapeMlinear/linear_model/linear_model/linear_model/direction.z/strided_slice/stackOlinear/linear_model/linear_model/linear_model/direction.z/strided_slice/stack_1Olinear/linear_model/linear_model/linear_model/direction.z/strided_slice/stack_2*
shrink_axis_mask*
Index0*
T0*
_output_shapes
: 
�
Ilinear/linear_model/linear_model/linear_model/direction.z/Reshape/shape/1Const*
value	B :*
dtype0*
_output_shapes
: 
�
Glinear/linear_model/linear_model/linear_model/direction.z/Reshape/shapePackGlinear/linear_model/linear_model/linear_model/direction.z/strided_sliceIlinear/linear_model/linear_model/linear_model/direction.z/Reshape/shape/1*
T0*
N*
_output_shapes
:
�
Alinear/linear_model/linear_model/linear_model/direction.z/ReshapeReshapeParseExample/ParseExample:2Glinear/linear_model/linear_model/linear_model/direction.z/Reshape/shape*
T0*'
_output_shapes
:���������
�
6linear/linear_model/direction.z/weights/ReadVariableOpReadVariableOp.linear/linear_model/direction.z/weights/part_0*
dtype0*
_output_shapes

:
�
'linear/linear_model/direction.z/weightsIdentity6linear/linear_model/direction.z/weights/ReadVariableOp*
T0*
_output_shapes

:
�
Flinear/linear_model/linear_model/linear_model/direction.z/weighted_sumMatMulAlinear/linear_model/linear_model/linear_model/direction.z/Reshape'linear/linear_model/direction.z/weights*
T0*'
_output_shapes
:���������
�
<linear/linear_model/linear_model/linear_model/origin.x/ShapeShapeParseExample/ParseExample:3*
T0*
_output_shapes
:
�
Jlinear/linear_model/linear_model/linear_model/origin.x/strided_slice/stackConst*
valueB: *
dtype0*
_output_shapes
:
�
Llinear/linear_model/linear_model/linear_model/origin.x/strided_slice/stack_1Const*
valueB:*
dtype0*
_output_shapes
:
�
Llinear/linear_model/linear_model/linear_model/origin.x/strided_slice/stack_2Const*
valueB:*
dtype0*
_output_shapes
:
�
Dlinear/linear_model/linear_model/linear_model/origin.x/strided_sliceStridedSlice<linear/linear_model/linear_model/linear_model/origin.x/ShapeJlinear/linear_model/linear_model/linear_model/origin.x/strided_slice/stackLlinear/linear_model/linear_model/linear_model/origin.x/strided_slice/stack_1Llinear/linear_model/linear_model/linear_model/origin.x/strided_slice/stack_2*
shrink_axis_mask*
Index0*
T0*
_output_shapes
: 
�
Flinear/linear_model/linear_model/linear_model/origin.x/Reshape/shape/1Const*
value	B :*
dtype0*
_output_shapes
: 
�
Dlinear/linear_model/linear_model/linear_model/origin.x/Reshape/shapePackDlinear/linear_model/linear_model/linear_model/origin.x/strided_sliceFlinear/linear_model/linear_model/linear_model/origin.x/Reshape/shape/1*
T0*
N*
_output_shapes
:
�
>linear/linear_model/linear_model/linear_model/origin.x/ReshapeReshapeParseExample/ParseExample:3Dlinear/linear_model/linear_model/linear_model/origin.x/Reshape/shape*
T0*'
_output_shapes
:���������
�
3linear/linear_model/origin.x/weights/ReadVariableOpReadVariableOp+linear/linear_model/origin.x/weights/part_0*
dtype0*
_output_shapes

:
�
$linear/linear_model/origin.x/weightsIdentity3linear/linear_model/origin.x/weights/ReadVariableOp*
T0*
_output_shapes

:
�
Clinear/linear_model/linear_model/linear_model/origin.x/weighted_sumMatMul>linear/linear_model/linear_model/linear_model/origin.x/Reshape$linear/linear_model/origin.x/weights*
T0*'
_output_shapes
:���������
�
<linear/linear_model/linear_model/linear_model/origin.y/ShapeShapeParseExample/ParseExample:4*
T0*
_output_shapes
:
�
Jlinear/linear_model/linear_model/linear_model/origin.y/strided_slice/stackConst*
valueB: *
dtype0*
_output_shapes
:
�
Llinear/linear_model/linear_model/linear_model/origin.y/strided_slice/stack_1Const*
valueB:*
dtype0*
_output_shapes
:
�
Llinear/linear_model/linear_model/linear_model/origin.y/strided_slice/stack_2Const*
valueB:*
dtype0*
_output_shapes
:
�
Dlinear/linear_model/linear_model/linear_model/origin.y/strided_sliceStridedSlice<linear/linear_model/linear_model/linear_model/origin.y/ShapeJlinear/linear_model/linear_model/linear_model/origin.y/strided_slice/stackLlinear/linear_model/linear_model/linear_model/origin.y/strided_slice/stack_1Llinear/linear_model/linear_model/linear_model/origin.y/strided_slice/stack_2*
shrink_axis_mask*
Index0*
T0*
_output_shapes
: 
�
Flinear/linear_model/linear_model/linear_model/origin.y/Reshape/shape/1Const*
value	B :*
dtype0*
_output_shapes
: 
�
Dlinear/linear_model/linear_model/linear_model/origin.y/Reshape/shapePackDlinear/linear_model/linear_model/linear_model/origin.y/strided_sliceFlinear/linear_model/linear_model/linear_model/origin.y/Reshape/shape/1*
T0*
N*
_output_shapes
:
�
>linear/linear_model/linear_model/linear_model/origin.y/ReshapeReshapeParseExample/ParseExample:4Dlinear/linear_model/linear_model/linear_model/origin.y/Reshape/shape*
T0*'
_output_shapes
:���������
�
3linear/linear_model/origin.y/weights/ReadVariableOpReadVariableOp+linear/linear_model/origin.y/weights/part_0*
dtype0*
_output_shapes

:
�
$linear/linear_model/origin.y/weightsIdentity3linear/linear_model/origin.y/weights/ReadVariableOp*
T0*
_output_shapes

:
�
Clinear/linear_model/linear_model/linear_model/origin.y/weighted_sumMatMul>linear/linear_model/linear_model/linear_model/origin.y/Reshape$linear/linear_model/origin.y/weights*
T0*'
_output_shapes
:���������
�
<linear/linear_model/linear_model/linear_model/origin.z/ShapeShapeParseExample/ParseExample:5*
T0*
_output_shapes
:
�
Jlinear/linear_model/linear_model/linear_model/origin.z/strided_slice/stackConst*
valueB: *
dtype0*
_output_shapes
:
�
Llinear/linear_model/linear_model/linear_model/origin.z/strided_slice/stack_1Const*
valueB:*
dtype0*
_output_shapes
:
�
Llinear/linear_model/linear_model/linear_model/origin.z/strided_slice/stack_2Const*
valueB:*
dtype0*
_output_shapes
:
�
Dlinear/linear_model/linear_model/linear_model/origin.z/strided_sliceStridedSlice<linear/linear_model/linear_model/linear_model/origin.z/ShapeJlinear/linear_model/linear_model/linear_model/origin.z/strided_slice/stackLlinear/linear_model/linear_model/linear_model/origin.z/strided_slice/stack_1Llinear/linear_model/linear_model/linear_model/origin.z/strided_slice/stack_2*
shrink_axis_mask*
Index0*
T0*
_output_shapes
: 
�
Flinear/linear_model/linear_model/linear_model/origin.z/Reshape/shape/1Const*
value	B :*
dtype0*
_output_shapes
: 
�
Dlinear/linear_model/linear_model/linear_model/origin.z/Reshape/shapePackDlinear/linear_model/linear_model/linear_model/origin.z/strided_sliceFlinear/linear_model/linear_model/linear_model/origin.z/Reshape/shape/1*
T0*
N*
_output_shapes
:
�
>linear/linear_model/linear_model/linear_model/origin.z/ReshapeReshapeParseExample/ParseExample:5Dlinear/linear_model/linear_model/linear_model/origin.z/Reshape/shape*
T0*'
_output_shapes
:���������
�
3linear/linear_model/origin.z/weights/ReadVariableOpReadVariableOp+linear/linear_model/origin.z/weights/part_0*
dtype0*
_output_shapes

:
�
$linear/linear_model/origin.z/weightsIdentity3linear/linear_model/origin.z/weights/ReadVariableOp*
T0*
_output_shapes

:
�
Clinear/linear_model/linear_model/linear_model/origin.z/weighted_sumMatMul>linear/linear_model/linear_model/linear_model/origin.z/Reshape$linear/linear_model/origin.z/weights*
T0*'
_output_shapes
:���������
�
Blinear/linear_model/linear_model/linear_model/weighted_sum_no_biasAddNFlinear/linear_model/linear_model/linear_model/direction.x/weighted_sumFlinear/linear_model/linear_model/linear_model/direction.y/weighted_sumFlinear/linear_model/linear_model/linear_model/direction.z/weighted_sumClinear/linear_model/linear_model/linear_model/origin.x/weighted_sumClinear/linear_model/linear_model/linear_model/origin.y/weighted_sumClinear/linear_model/linear_model/linear_model/origin.z/weighted_sum*
T0*
N*'
_output_shapes
:���������
�
/linear/linear_model/bias_weights/ReadVariableOpReadVariableOp'linear/linear_model/bias_weights/part_0*
dtype0*
_output_shapes
:
�
 linear/linear_model/bias_weightsIdentity/linear/linear_model/bias_weights/ReadVariableOp*
T0*
_output_shapes
:
�
:linear/linear_model/linear_model/linear_model/weighted_sumBiasAddBlinear/linear_model/linear_model/linear_model/weighted_sum_no_bias linear/linear_model/bias_weights*
T0*'
_output_shapes
:���������
y
linear/ReadVariableOpReadVariableOp'linear/linear_model/bias_weights/part_0*
dtype0*
_output_shapes
:
d
linear/strided_slice/stackConst*
valueB: *
dtype0*
_output_shapes
:
f
linear/strided_slice/stack_1Const*
valueB:*
dtype0*
_output_shapes
:
f
linear/strided_slice/stack_2Const*
valueB:*
dtype0*
_output_shapes
:
�
linear/strided_sliceStridedSlicelinear/ReadVariableOplinear/strided_slice/stacklinear/strided_slice/stack_1linear/strided_slice/stack_2*
shrink_axis_mask*
Index0*
T0*
_output_shapes
: 
\
linear/bias/tagsConst*
valueB Blinear/bias*
dtype0*
_output_shapes
: 
e
linear/biasScalarSummarylinear/bias/tagslinear/strided_slice*
T0*
_output_shapes
: 
�
3linear/zero_fraction/total_size/Size/ReadVariableOpReadVariableOp.linear/linear_model/direction.x/weights/part_0*
dtype0*
_output_shapes

:
f
$linear/zero_fraction/total_size/SizeConst*
value	B	 R*
dtype0	*
_output_shapes
: 
�
5linear/zero_fraction/total_size/Size_1/ReadVariableOpReadVariableOp.linear/linear_model/direction.y/weights/part_0*
dtype0*
_output_shapes

:
h
&linear/zero_fraction/total_size/Size_1Const*
value	B	 R*
dtype0	*
_output_shapes
: 
�
5linear/zero_fraction/total_size/Size_2/ReadVariableOpReadVariableOp.linear/linear_model/direction.z/weights/part_0*
dtype0*
_output_shapes

:
h
&linear/zero_fraction/total_size/Size_2Const*
value	B	 R*
dtype0	*
_output_shapes
: 
�
5linear/zero_fraction/total_size/Size_3/ReadVariableOpReadVariableOp+linear/linear_model/origin.x/weights/part_0*
dtype0*
_output_shapes

:
h
&linear/zero_fraction/total_size/Size_3Const*
value	B	 R*
dtype0	*
_output_shapes
: 
�
5linear/zero_fraction/total_size/Size_4/ReadVariableOpReadVariableOp+linear/linear_model/origin.y/weights/part_0*
dtype0*
_output_shapes

:
h
&linear/zero_fraction/total_size/Size_4Const*
value	B	 R*
dtype0	*
_output_shapes
: 
�
5linear/zero_fraction/total_size/Size_5/ReadVariableOpReadVariableOp+linear/linear_model/origin.z/weights/part_0*
dtype0*
_output_shapes

:
h
&linear/zero_fraction/total_size/Size_5Const*
value	B	 R*
dtype0	*
_output_shapes
: 
�
$linear/zero_fraction/total_size/AddNAddN$linear/zero_fraction/total_size/Size&linear/zero_fraction/total_size/Size_1&linear/zero_fraction/total_size/Size_2&linear/zero_fraction/total_size/Size_3&linear/zero_fraction/total_size/Size_4&linear/zero_fraction/total_size/Size_5*
T0	*
N*
_output_shapes
: 
g
%linear/zero_fraction/total_zero/ConstConst*
value	B	 R *
dtype0	*
_output_shapes
: 
�
%linear/zero_fraction/total_zero/EqualEqual$linear/zero_fraction/total_size/Size%linear/zero_fraction/total_zero/Const*
T0	*
_output_shapes
: 
�
1linear/zero_fraction/total_zero/zero_count/SwitchSwitch%linear/zero_fraction/total_zero/Equal%linear/zero_fraction/total_zero/Equal*
T0
*
_output_shapes
: : 
�
3linear/zero_fraction/total_zero/zero_count/switch_tIdentity3linear/zero_fraction/total_zero/zero_count/Switch:1*
T0
*
_output_shapes
: 
�
3linear/zero_fraction/total_zero/zero_count/switch_fIdentity1linear/zero_fraction/total_zero/zero_count/Switch*
T0
*
_output_shapes
: 
�
2linear/zero_fraction/total_zero/zero_count/pred_idIdentity%linear/zero_fraction/total_zero/Equal*
T0
*
_output_shapes
: 
�
0linear/zero_fraction/total_zero/zero_count/ConstConst4^linear/zero_fraction/total_zero/zero_count/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Glinear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOpReadVariableOpNlinear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOp/Switch*
dtype0*
_output_shapes

:
�
Nlinear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOp/SwitchSwitch.linear/linear_model/direction.x/weights/part_02linear/zero_fraction/total_zero/zero_count/pred_id*
T0*A
_class7
53loc:@linear/linear_model/direction.x/weights/part_0*
_output_shapes
: : 
�
=linear/zero_fraction/total_zero/zero_count/zero_fraction/SizeConst4^linear/zero_fraction/total_zero/zero_count/switch_f*
value	B	 R*
dtype0	*
_output_shapes
: 
�
Dlinear/zero_fraction/total_zero/zero_count/zero_fraction/LessEqual/yConst4^linear/zero_fraction/total_zero/zero_count/switch_f*
valueB	 R����*
dtype0	*
_output_shapes
: 
�
Blinear/zero_fraction/total_zero/zero_count/zero_fraction/LessEqual	LessEqual=linear/zero_fraction/total_zero/zero_count/zero_fraction/SizeDlinear/zero_fraction/total_zero/zero_count/zero_fraction/LessEqual/y*
T0	*
_output_shapes
: 
�
Dlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/SwitchSwitchBlinear/zero_fraction/total_zero/zero_count/zero_fraction/LessEqualBlinear/zero_fraction/total_zero/zero_count/zero_fraction/LessEqual*
T0
*
_output_shapes
: : 
�
Flinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/switch_tIdentityFlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/Switch:1*
T0
*
_output_shapes
: 
�
Flinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/switch_fIdentityDlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/Switch*
T0
*
_output_shapes
: 
�
Elinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/pred_idIdentityBlinear/zero_fraction/total_zero/zero_count/zero_fraction/LessEqual*
T0
*
_output_shapes
: 
�
Qlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/zerosConstG^linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Tlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/NotEqualNotEqual]linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/NotEqual/Switch:1Qlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/zeros*
T0*
_output_shapes

:
�
[linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/NotEqual/SwitchSwitchGlinear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOpElinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/pred_id*
T0*Z
_classP
NLloc:@linear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOp*(
_output_shapes
::
�
Plinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/CastCastTlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/NotEqual*

SrcT0
*
_output_shapes

:*

DstT0
�
Qlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/ConstConstG^linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/switch_t*
valueB"       *
dtype0*
_output_shapes
:
�
Ylinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/nonzero_countSumPlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/CastQlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/Const*
T0*
_output_shapes
: 
�
Blinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/CastCastYlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/nonzero_count*

SrcT0*
_output_shapes
: *

DstT0	
�
Slinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/zerosConstG^linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/switch_f*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Vlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/NotEqualNotEqual]linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/NotEqual/SwitchSlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/zeros*
T0*
_output_shapes

:
�
]linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/NotEqual/SwitchSwitchGlinear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOpElinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/pred_id*
T0*Z
_classP
NLloc:@linear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOp*(
_output_shapes
::
�
Rlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/CastCastVlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/NotEqual*

SrcT0
*
_output_shapes

:*

DstT0	
�
Slinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/ConstConstG^linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/switch_f*
valueB"       *
dtype0*
_output_shapes
:
�
[linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/nonzero_countSumRlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/CastSlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/Const*
T0	*
_output_shapes
: 
�
Clinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/MergeMerge[linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/nonzero_countBlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/Cast*
T0	*
N*
_output_shapes
: : 
�
Olinear/zero_fraction/total_zero/zero_count/zero_fraction/counts_to_fraction/subSub=linear/zero_fraction/total_zero/zero_count/zero_fraction/SizeClinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/Merge*
T0	*
_output_shapes
: 
�
Plinear/zero_fraction/total_zero/zero_count/zero_fraction/counts_to_fraction/CastCastOlinear/zero_fraction/total_zero/zero_count/zero_fraction/counts_to_fraction/sub*

SrcT0	*
_output_shapes
: *

DstT0
�
Rlinear/zero_fraction/total_zero/zero_count/zero_fraction/counts_to_fraction/Cast_1Cast=linear/zero_fraction/total_zero/zero_count/zero_fraction/Size*

SrcT0	*
_output_shapes
: *

DstT0
�
Slinear/zero_fraction/total_zero/zero_count/zero_fraction/counts_to_fraction/truedivRealDivPlinear/zero_fraction/total_zero/zero_count/zero_fraction/counts_to_fraction/CastRlinear/zero_fraction/total_zero/zero_count/zero_fraction/counts_to_fraction/Cast_1*
T0*
_output_shapes
: 
�
Alinear/zero_fraction/total_zero/zero_count/zero_fraction/fractionIdentitySlinear/zero_fraction/total_zero/zero_count/zero_fraction/counts_to_fraction/truediv*
T0*
_output_shapes
: 
�
2linear/zero_fraction/total_zero/zero_count/ToFloatCast9linear/zero_fraction/total_zero/zero_count/ToFloat/Switch*

SrcT0	*
_output_shapes
: *

DstT0
�
9linear/zero_fraction/total_zero/zero_count/ToFloat/SwitchSwitch$linear/zero_fraction/total_size/Size2linear/zero_fraction/total_zero/zero_count/pred_id*
T0	*7
_class-
+)loc:@linear/zero_fraction/total_size/Size*
_output_shapes
: : 
�
.linear/zero_fraction/total_zero/zero_count/mulMulAlinear/zero_fraction/total_zero/zero_count/zero_fraction/fraction2linear/zero_fraction/total_zero/zero_count/ToFloat*
T0*
_output_shapes
: 
�
0linear/zero_fraction/total_zero/zero_count/MergeMerge.linear/zero_fraction/total_zero/zero_count/mul0linear/zero_fraction/total_zero/zero_count/Const*
T0*
N*
_output_shapes
: : 
i
'linear/zero_fraction/total_zero/Const_1Const*
value	B	 R *
dtype0	*
_output_shapes
: 
�
'linear/zero_fraction/total_zero/Equal_1Equal&linear/zero_fraction/total_size/Size_1'linear/zero_fraction/total_zero/Const_1*
T0	*
_output_shapes
: 
�
3linear/zero_fraction/total_zero/zero_count_1/SwitchSwitch'linear/zero_fraction/total_zero/Equal_1'linear/zero_fraction/total_zero/Equal_1*
T0
*
_output_shapes
: : 
�
5linear/zero_fraction/total_zero/zero_count_1/switch_tIdentity5linear/zero_fraction/total_zero/zero_count_1/Switch:1*
T0
*
_output_shapes
: 
�
5linear/zero_fraction/total_zero/zero_count_1/switch_fIdentity3linear/zero_fraction/total_zero/zero_count_1/Switch*
T0
*
_output_shapes
: 
�
4linear/zero_fraction/total_zero/zero_count_1/pred_idIdentity'linear/zero_fraction/total_zero/Equal_1*
T0
*
_output_shapes
: 
�
2linear/zero_fraction/total_zero/zero_count_1/ConstConst6^linear/zero_fraction/total_zero/zero_count_1/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Ilinear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOpReadVariableOpPlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOp/Switch*
dtype0*
_output_shapes

:
�
Plinear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOp/SwitchSwitch.linear/linear_model/direction.y/weights/part_04linear/zero_fraction/total_zero/zero_count_1/pred_id*
T0*A
_class7
53loc:@linear/linear_model/direction.y/weights/part_0*
_output_shapes
: : 
�
?linear/zero_fraction/total_zero/zero_count_1/zero_fraction/SizeConst6^linear/zero_fraction/total_zero/zero_count_1/switch_f*
value	B	 R*
dtype0	*
_output_shapes
: 
�
Flinear/zero_fraction/total_zero/zero_count_1/zero_fraction/LessEqual/yConst6^linear/zero_fraction/total_zero/zero_count_1/switch_f*
valueB	 R����*
dtype0	*
_output_shapes
: 
�
Dlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/LessEqual	LessEqual?linear/zero_fraction/total_zero/zero_count_1/zero_fraction/SizeFlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/LessEqual/y*
T0	*
_output_shapes
: 
�
Flinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/SwitchSwitchDlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/LessEqualDlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/LessEqual*
T0
*
_output_shapes
: : 
�
Hlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/switch_tIdentityHlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/Switch:1*
T0
*
_output_shapes
: 
�
Hlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/switch_fIdentityFlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/Switch*
T0
*
_output_shapes
: 
�
Glinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/pred_idIdentityDlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/LessEqual*
T0
*
_output_shapes
: 
�
Slinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/zerosConstI^linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Vlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/NotEqualNotEqual_linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/NotEqual/Switch:1Slinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/zeros*
T0*
_output_shapes

:
�
]linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/NotEqual/SwitchSwitchIlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOpGlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/pred_id*
T0*\
_classR
PNloc:@linear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOp*(
_output_shapes
::
�
Rlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/CastCastVlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/NotEqual*

SrcT0
*
_output_shapes

:*

DstT0
�
Slinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/ConstConstI^linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/switch_t*
valueB"       *
dtype0*
_output_shapes
:
�
[linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/nonzero_countSumRlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/CastSlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/Const*
T0*
_output_shapes
: 
�
Dlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/CastCast[linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/nonzero_count*

SrcT0*
_output_shapes
: *

DstT0	
�
Ulinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/zerosConstI^linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/switch_f*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Xlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/NotEqualNotEqual_linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/NotEqual/SwitchUlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/zeros*
T0*
_output_shapes

:
�
_linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/NotEqual/SwitchSwitchIlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOpGlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/pred_id*
T0*\
_classR
PNloc:@linear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOp*(
_output_shapes
::
�
Tlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/CastCastXlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/NotEqual*

SrcT0
*
_output_shapes

:*

DstT0	
�
Ulinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/ConstConstI^linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/switch_f*
valueB"       *
dtype0*
_output_shapes
:
�
]linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/nonzero_countSumTlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/CastUlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/Const*
T0	*
_output_shapes
: 
�
Elinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/MergeMerge]linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/nonzero_countDlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/Cast*
T0	*
N*
_output_shapes
: : 
�
Qlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/counts_to_fraction/subSub?linear/zero_fraction/total_zero/zero_count_1/zero_fraction/SizeElinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/Merge*
T0	*
_output_shapes
: 
�
Rlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/counts_to_fraction/CastCastQlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/counts_to_fraction/sub*

SrcT0	*
_output_shapes
: *

DstT0
�
Tlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/counts_to_fraction/Cast_1Cast?linear/zero_fraction/total_zero/zero_count_1/zero_fraction/Size*

SrcT0	*
_output_shapes
: *

DstT0
�
Ulinear/zero_fraction/total_zero/zero_count_1/zero_fraction/counts_to_fraction/truedivRealDivRlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/counts_to_fraction/CastTlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/counts_to_fraction/Cast_1*
T0*
_output_shapes
: 
�
Clinear/zero_fraction/total_zero/zero_count_1/zero_fraction/fractionIdentityUlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/counts_to_fraction/truediv*
T0*
_output_shapes
: 
�
4linear/zero_fraction/total_zero/zero_count_1/ToFloatCast;linear/zero_fraction/total_zero/zero_count_1/ToFloat/Switch*

SrcT0	*
_output_shapes
: *

DstT0
�
;linear/zero_fraction/total_zero/zero_count_1/ToFloat/SwitchSwitch&linear/zero_fraction/total_size/Size_14linear/zero_fraction/total_zero/zero_count_1/pred_id*
T0	*9
_class/
-+loc:@linear/zero_fraction/total_size/Size_1*
_output_shapes
: : 
�
0linear/zero_fraction/total_zero/zero_count_1/mulMulClinear/zero_fraction/total_zero/zero_count_1/zero_fraction/fraction4linear/zero_fraction/total_zero/zero_count_1/ToFloat*
T0*
_output_shapes
: 
�
2linear/zero_fraction/total_zero/zero_count_1/MergeMerge0linear/zero_fraction/total_zero/zero_count_1/mul2linear/zero_fraction/total_zero/zero_count_1/Const*
T0*
N*
_output_shapes
: : 
i
'linear/zero_fraction/total_zero/Const_2Const*
value	B	 R *
dtype0	*
_output_shapes
: 
�
'linear/zero_fraction/total_zero/Equal_2Equal&linear/zero_fraction/total_size/Size_2'linear/zero_fraction/total_zero/Const_2*
T0	*
_output_shapes
: 
�
3linear/zero_fraction/total_zero/zero_count_2/SwitchSwitch'linear/zero_fraction/total_zero/Equal_2'linear/zero_fraction/total_zero/Equal_2*
T0
*
_output_shapes
: : 
�
5linear/zero_fraction/total_zero/zero_count_2/switch_tIdentity5linear/zero_fraction/total_zero/zero_count_2/Switch:1*
T0
*
_output_shapes
: 
�
5linear/zero_fraction/total_zero/zero_count_2/switch_fIdentity3linear/zero_fraction/total_zero/zero_count_2/Switch*
T0
*
_output_shapes
: 
�
4linear/zero_fraction/total_zero/zero_count_2/pred_idIdentity'linear/zero_fraction/total_zero/Equal_2*
T0
*
_output_shapes
: 
�
2linear/zero_fraction/total_zero/zero_count_2/ConstConst6^linear/zero_fraction/total_zero/zero_count_2/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Ilinear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOpReadVariableOpPlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOp/Switch*
dtype0*
_output_shapes

:
�
Plinear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOp/SwitchSwitch.linear/linear_model/direction.z/weights/part_04linear/zero_fraction/total_zero/zero_count_2/pred_id*
T0*A
_class7
53loc:@linear/linear_model/direction.z/weights/part_0*
_output_shapes
: : 
�
?linear/zero_fraction/total_zero/zero_count_2/zero_fraction/SizeConst6^linear/zero_fraction/total_zero/zero_count_2/switch_f*
value	B	 R*
dtype0	*
_output_shapes
: 
�
Flinear/zero_fraction/total_zero/zero_count_2/zero_fraction/LessEqual/yConst6^linear/zero_fraction/total_zero/zero_count_2/switch_f*
valueB	 R����*
dtype0	*
_output_shapes
: 
�
Dlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/LessEqual	LessEqual?linear/zero_fraction/total_zero/zero_count_2/zero_fraction/SizeFlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/LessEqual/y*
T0	*
_output_shapes
: 
�
Flinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/SwitchSwitchDlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/LessEqualDlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/LessEqual*
T0
*
_output_shapes
: : 
�
Hlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/switch_tIdentityHlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/Switch:1*
T0
*
_output_shapes
: 
�
Hlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/switch_fIdentityFlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/Switch*
T0
*
_output_shapes
: 
�
Glinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/pred_idIdentityDlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/LessEqual*
T0
*
_output_shapes
: 
�
Slinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/zerosConstI^linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Vlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/NotEqualNotEqual_linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/NotEqual/Switch:1Slinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/zeros*
T0*
_output_shapes

:
�
]linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/NotEqual/SwitchSwitchIlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOpGlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/pred_id*
T0*\
_classR
PNloc:@linear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOp*(
_output_shapes
::
�
Rlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/CastCastVlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/NotEqual*

SrcT0
*
_output_shapes

:*

DstT0
�
Slinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/ConstConstI^linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/switch_t*
valueB"       *
dtype0*
_output_shapes
:
�
[linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/nonzero_countSumRlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/CastSlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/Const*
T0*
_output_shapes
: 
�
Dlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/CastCast[linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/nonzero_count*

SrcT0*
_output_shapes
: *

DstT0	
�
Ulinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/zerosConstI^linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/switch_f*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Xlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/NotEqualNotEqual_linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/NotEqual/SwitchUlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/zeros*
T0*
_output_shapes

:
�
_linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/NotEqual/SwitchSwitchIlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOpGlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/pred_id*
T0*\
_classR
PNloc:@linear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOp*(
_output_shapes
::
�
Tlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/CastCastXlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/NotEqual*

SrcT0
*
_output_shapes

:*

DstT0	
�
Ulinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/ConstConstI^linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/switch_f*
valueB"       *
dtype0*
_output_shapes
:
�
]linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/nonzero_countSumTlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/CastUlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/Const*
T0	*
_output_shapes
: 
�
Elinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/MergeMerge]linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/nonzero_countDlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/Cast*
T0	*
N*
_output_shapes
: : 
�
Qlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/counts_to_fraction/subSub?linear/zero_fraction/total_zero/zero_count_2/zero_fraction/SizeElinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/Merge*
T0	*
_output_shapes
: 
�
Rlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/counts_to_fraction/CastCastQlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/counts_to_fraction/sub*

SrcT0	*
_output_shapes
: *

DstT0
�
Tlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/counts_to_fraction/Cast_1Cast?linear/zero_fraction/total_zero/zero_count_2/zero_fraction/Size*

SrcT0	*
_output_shapes
: *

DstT0
�
Ulinear/zero_fraction/total_zero/zero_count_2/zero_fraction/counts_to_fraction/truedivRealDivRlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/counts_to_fraction/CastTlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/counts_to_fraction/Cast_1*
T0*
_output_shapes
: 
�
Clinear/zero_fraction/total_zero/zero_count_2/zero_fraction/fractionIdentityUlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/counts_to_fraction/truediv*
T0*
_output_shapes
: 
�
4linear/zero_fraction/total_zero/zero_count_2/ToFloatCast;linear/zero_fraction/total_zero/zero_count_2/ToFloat/Switch*

SrcT0	*
_output_shapes
: *

DstT0
�
;linear/zero_fraction/total_zero/zero_count_2/ToFloat/SwitchSwitch&linear/zero_fraction/total_size/Size_24linear/zero_fraction/total_zero/zero_count_2/pred_id*
T0	*9
_class/
-+loc:@linear/zero_fraction/total_size/Size_2*
_output_shapes
: : 
�
0linear/zero_fraction/total_zero/zero_count_2/mulMulClinear/zero_fraction/total_zero/zero_count_2/zero_fraction/fraction4linear/zero_fraction/total_zero/zero_count_2/ToFloat*
T0*
_output_shapes
: 
�
2linear/zero_fraction/total_zero/zero_count_2/MergeMerge0linear/zero_fraction/total_zero/zero_count_2/mul2linear/zero_fraction/total_zero/zero_count_2/Const*
T0*
N*
_output_shapes
: : 
i
'linear/zero_fraction/total_zero/Const_3Const*
value	B	 R *
dtype0	*
_output_shapes
: 
�
'linear/zero_fraction/total_zero/Equal_3Equal&linear/zero_fraction/total_size/Size_3'linear/zero_fraction/total_zero/Const_3*
T0	*
_output_shapes
: 
�
3linear/zero_fraction/total_zero/zero_count_3/SwitchSwitch'linear/zero_fraction/total_zero/Equal_3'linear/zero_fraction/total_zero/Equal_3*
T0
*
_output_shapes
: : 
�
5linear/zero_fraction/total_zero/zero_count_3/switch_tIdentity5linear/zero_fraction/total_zero/zero_count_3/Switch:1*
T0
*
_output_shapes
: 
�
5linear/zero_fraction/total_zero/zero_count_3/switch_fIdentity3linear/zero_fraction/total_zero/zero_count_3/Switch*
T0
*
_output_shapes
: 
�
4linear/zero_fraction/total_zero/zero_count_3/pred_idIdentity'linear/zero_fraction/total_zero/Equal_3*
T0
*
_output_shapes
: 
�
2linear/zero_fraction/total_zero/zero_count_3/ConstConst6^linear/zero_fraction/total_zero/zero_count_3/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Ilinear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOpReadVariableOpPlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOp/Switch*
dtype0*
_output_shapes

:
�
Plinear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOp/SwitchSwitch+linear/linear_model/origin.x/weights/part_04linear/zero_fraction/total_zero/zero_count_3/pred_id*
T0*>
_class4
20loc:@linear/linear_model/origin.x/weights/part_0*
_output_shapes
: : 
�
?linear/zero_fraction/total_zero/zero_count_3/zero_fraction/SizeConst6^linear/zero_fraction/total_zero/zero_count_3/switch_f*
value	B	 R*
dtype0	*
_output_shapes
: 
�
Flinear/zero_fraction/total_zero/zero_count_3/zero_fraction/LessEqual/yConst6^linear/zero_fraction/total_zero/zero_count_3/switch_f*
valueB	 R����*
dtype0	*
_output_shapes
: 
�
Dlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/LessEqual	LessEqual?linear/zero_fraction/total_zero/zero_count_3/zero_fraction/SizeFlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/LessEqual/y*
T0	*
_output_shapes
: 
�
Flinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/SwitchSwitchDlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/LessEqualDlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/LessEqual*
T0
*
_output_shapes
: : 
�
Hlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/switch_tIdentityHlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/Switch:1*
T0
*
_output_shapes
: 
�
Hlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/switch_fIdentityFlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/Switch*
T0
*
_output_shapes
: 
�
Glinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/pred_idIdentityDlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/LessEqual*
T0
*
_output_shapes
: 
�
Slinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/zerosConstI^linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Vlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/NotEqualNotEqual_linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/NotEqual/Switch:1Slinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/zeros*
T0*
_output_shapes

:
�
]linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/NotEqual/SwitchSwitchIlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOpGlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/pred_id*
T0*\
_classR
PNloc:@linear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOp*(
_output_shapes
::
�
Rlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/CastCastVlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/NotEqual*

SrcT0
*
_output_shapes

:*

DstT0
�
Slinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/ConstConstI^linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/switch_t*
valueB"       *
dtype0*
_output_shapes
:
�
[linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/nonzero_countSumRlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/CastSlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/Const*
T0*
_output_shapes
: 
�
Dlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/CastCast[linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/nonzero_count*

SrcT0*
_output_shapes
: *

DstT0	
�
Ulinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/zerosConstI^linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/switch_f*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Xlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/NotEqualNotEqual_linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/NotEqual/SwitchUlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/zeros*
T0*
_output_shapes

:
�
_linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/NotEqual/SwitchSwitchIlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOpGlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/pred_id*
T0*\
_classR
PNloc:@linear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOp*(
_output_shapes
::
�
Tlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/CastCastXlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/NotEqual*

SrcT0
*
_output_shapes

:*

DstT0	
�
Ulinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/ConstConstI^linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/switch_f*
valueB"       *
dtype0*
_output_shapes
:
�
]linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/nonzero_countSumTlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/CastUlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/Const*
T0	*
_output_shapes
: 
�
Elinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/MergeMerge]linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/nonzero_countDlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/Cast*
T0	*
N*
_output_shapes
: : 
�
Qlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/counts_to_fraction/subSub?linear/zero_fraction/total_zero/zero_count_3/zero_fraction/SizeElinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/Merge*
T0	*
_output_shapes
: 
�
Rlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/counts_to_fraction/CastCastQlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/counts_to_fraction/sub*

SrcT0	*
_output_shapes
: *

DstT0
�
Tlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/counts_to_fraction/Cast_1Cast?linear/zero_fraction/total_zero/zero_count_3/zero_fraction/Size*

SrcT0	*
_output_shapes
: *

DstT0
�
Ulinear/zero_fraction/total_zero/zero_count_3/zero_fraction/counts_to_fraction/truedivRealDivRlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/counts_to_fraction/CastTlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/counts_to_fraction/Cast_1*
T0*
_output_shapes
: 
�
Clinear/zero_fraction/total_zero/zero_count_3/zero_fraction/fractionIdentityUlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/counts_to_fraction/truediv*
T0*
_output_shapes
: 
�
4linear/zero_fraction/total_zero/zero_count_3/ToFloatCast;linear/zero_fraction/total_zero/zero_count_3/ToFloat/Switch*

SrcT0	*
_output_shapes
: *

DstT0
�
;linear/zero_fraction/total_zero/zero_count_3/ToFloat/SwitchSwitch&linear/zero_fraction/total_size/Size_34linear/zero_fraction/total_zero/zero_count_3/pred_id*
T0	*9
_class/
-+loc:@linear/zero_fraction/total_size/Size_3*
_output_shapes
: : 
�
0linear/zero_fraction/total_zero/zero_count_3/mulMulClinear/zero_fraction/total_zero/zero_count_3/zero_fraction/fraction4linear/zero_fraction/total_zero/zero_count_3/ToFloat*
T0*
_output_shapes
: 
�
2linear/zero_fraction/total_zero/zero_count_3/MergeMerge0linear/zero_fraction/total_zero/zero_count_3/mul2linear/zero_fraction/total_zero/zero_count_3/Const*
T0*
N*
_output_shapes
: : 
i
'linear/zero_fraction/total_zero/Const_4Const*
value	B	 R *
dtype0	*
_output_shapes
: 
�
'linear/zero_fraction/total_zero/Equal_4Equal&linear/zero_fraction/total_size/Size_4'linear/zero_fraction/total_zero/Const_4*
T0	*
_output_shapes
: 
�
3linear/zero_fraction/total_zero/zero_count_4/SwitchSwitch'linear/zero_fraction/total_zero/Equal_4'linear/zero_fraction/total_zero/Equal_4*
T0
*
_output_shapes
: : 
�
5linear/zero_fraction/total_zero/zero_count_4/switch_tIdentity5linear/zero_fraction/total_zero/zero_count_4/Switch:1*
T0
*
_output_shapes
: 
�
5linear/zero_fraction/total_zero/zero_count_4/switch_fIdentity3linear/zero_fraction/total_zero/zero_count_4/Switch*
T0
*
_output_shapes
: 
�
4linear/zero_fraction/total_zero/zero_count_4/pred_idIdentity'linear/zero_fraction/total_zero/Equal_4*
T0
*
_output_shapes
: 
�
2linear/zero_fraction/total_zero/zero_count_4/ConstConst6^linear/zero_fraction/total_zero/zero_count_4/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Ilinear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOpReadVariableOpPlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOp/Switch*
dtype0*
_output_shapes

:
�
Plinear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOp/SwitchSwitch+linear/linear_model/origin.y/weights/part_04linear/zero_fraction/total_zero/zero_count_4/pred_id*
T0*>
_class4
20loc:@linear/linear_model/origin.y/weights/part_0*
_output_shapes
: : 
�
?linear/zero_fraction/total_zero/zero_count_4/zero_fraction/SizeConst6^linear/zero_fraction/total_zero/zero_count_4/switch_f*
value	B	 R*
dtype0	*
_output_shapes
: 
�
Flinear/zero_fraction/total_zero/zero_count_4/zero_fraction/LessEqual/yConst6^linear/zero_fraction/total_zero/zero_count_4/switch_f*
valueB	 R����*
dtype0	*
_output_shapes
: 
�
Dlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/LessEqual	LessEqual?linear/zero_fraction/total_zero/zero_count_4/zero_fraction/SizeFlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/LessEqual/y*
T0	*
_output_shapes
: 
�
Flinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/SwitchSwitchDlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/LessEqualDlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/LessEqual*
T0
*
_output_shapes
: : 
�
Hlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/switch_tIdentityHlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/Switch:1*
T0
*
_output_shapes
: 
�
Hlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/switch_fIdentityFlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/Switch*
T0
*
_output_shapes
: 
�
Glinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/pred_idIdentityDlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/LessEqual*
T0
*
_output_shapes
: 
�
Slinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/zerosConstI^linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Vlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/NotEqualNotEqual_linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/NotEqual/Switch:1Slinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/zeros*
T0*
_output_shapes

:
�
]linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/NotEqual/SwitchSwitchIlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOpGlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/pred_id*
T0*\
_classR
PNloc:@linear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOp*(
_output_shapes
::
�
Rlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/CastCastVlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/NotEqual*

SrcT0
*
_output_shapes

:*

DstT0
�
Slinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/ConstConstI^linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/switch_t*
valueB"       *
dtype0*
_output_shapes
:
�
[linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/nonzero_countSumRlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/CastSlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/Const*
T0*
_output_shapes
: 
�
Dlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/CastCast[linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/nonzero_count*

SrcT0*
_output_shapes
: *

DstT0	
�
Ulinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/zerosConstI^linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/switch_f*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Xlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/NotEqualNotEqual_linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/NotEqual/SwitchUlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/zeros*
T0*
_output_shapes

:
�
_linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/NotEqual/SwitchSwitchIlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOpGlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/pred_id*
T0*\
_classR
PNloc:@linear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOp*(
_output_shapes
::
�
Tlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/CastCastXlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/NotEqual*

SrcT0
*
_output_shapes

:*

DstT0	
�
Ulinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/ConstConstI^linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/switch_f*
valueB"       *
dtype0*
_output_shapes
:
�
]linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/nonzero_countSumTlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/CastUlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/Const*
T0	*
_output_shapes
: 
�
Elinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/MergeMerge]linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/nonzero_countDlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/Cast*
T0	*
N*
_output_shapes
: : 
�
Qlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/counts_to_fraction/subSub?linear/zero_fraction/total_zero/zero_count_4/zero_fraction/SizeElinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/Merge*
T0	*
_output_shapes
: 
�
Rlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/counts_to_fraction/CastCastQlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/counts_to_fraction/sub*

SrcT0	*
_output_shapes
: *

DstT0
�
Tlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/counts_to_fraction/Cast_1Cast?linear/zero_fraction/total_zero/zero_count_4/zero_fraction/Size*

SrcT0	*
_output_shapes
: *

DstT0
�
Ulinear/zero_fraction/total_zero/zero_count_4/zero_fraction/counts_to_fraction/truedivRealDivRlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/counts_to_fraction/CastTlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/counts_to_fraction/Cast_1*
T0*
_output_shapes
: 
�
Clinear/zero_fraction/total_zero/zero_count_4/zero_fraction/fractionIdentityUlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/counts_to_fraction/truediv*
T0*
_output_shapes
: 
�
4linear/zero_fraction/total_zero/zero_count_4/ToFloatCast;linear/zero_fraction/total_zero/zero_count_4/ToFloat/Switch*

SrcT0	*
_output_shapes
: *

DstT0
�
;linear/zero_fraction/total_zero/zero_count_4/ToFloat/SwitchSwitch&linear/zero_fraction/total_size/Size_44linear/zero_fraction/total_zero/zero_count_4/pred_id*
T0	*9
_class/
-+loc:@linear/zero_fraction/total_size/Size_4*
_output_shapes
: : 
�
0linear/zero_fraction/total_zero/zero_count_4/mulMulClinear/zero_fraction/total_zero/zero_count_4/zero_fraction/fraction4linear/zero_fraction/total_zero/zero_count_4/ToFloat*
T0*
_output_shapes
: 
�
2linear/zero_fraction/total_zero/zero_count_4/MergeMerge0linear/zero_fraction/total_zero/zero_count_4/mul2linear/zero_fraction/total_zero/zero_count_4/Const*
T0*
N*
_output_shapes
: : 
i
'linear/zero_fraction/total_zero/Const_5Const*
value	B	 R *
dtype0	*
_output_shapes
: 
�
'linear/zero_fraction/total_zero/Equal_5Equal&linear/zero_fraction/total_size/Size_5'linear/zero_fraction/total_zero/Const_5*
T0	*
_output_shapes
: 
�
3linear/zero_fraction/total_zero/zero_count_5/SwitchSwitch'linear/zero_fraction/total_zero/Equal_5'linear/zero_fraction/total_zero/Equal_5*
T0
*
_output_shapes
: : 
�
5linear/zero_fraction/total_zero/zero_count_5/switch_tIdentity5linear/zero_fraction/total_zero/zero_count_5/Switch:1*
T0
*
_output_shapes
: 
�
5linear/zero_fraction/total_zero/zero_count_5/switch_fIdentity3linear/zero_fraction/total_zero/zero_count_5/Switch*
T0
*
_output_shapes
: 
�
4linear/zero_fraction/total_zero/zero_count_5/pred_idIdentity'linear/zero_fraction/total_zero/Equal_5*
T0
*
_output_shapes
: 
�
2linear/zero_fraction/total_zero/zero_count_5/ConstConst6^linear/zero_fraction/total_zero/zero_count_5/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Ilinear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOpReadVariableOpPlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOp/Switch*
dtype0*
_output_shapes

:
�
Plinear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOp/SwitchSwitch+linear/linear_model/origin.z/weights/part_04linear/zero_fraction/total_zero/zero_count_5/pred_id*
T0*>
_class4
20loc:@linear/linear_model/origin.z/weights/part_0*
_output_shapes
: : 
�
?linear/zero_fraction/total_zero/zero_count_5/zero_fraction/SizeConst6^linear/zero_fraction/total_zero/zero_count_5/switch_f*
value	B	 R*
dtype0	*
_output_shapes
: 
�
Flinear/zero_fraction/total_zero/zero_count_5/zero_fraction/LessEqual/yConst6^linear/zero_fraction/total_zero/zero_count_5/switch_f*
valueB	 R����*
dtype0	*
_output_shapes
: 
�
Dlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/LessEqual	LessEqual?linear/zero_fraction/total_zero/zero_count_5/zero_fraction/SizeFlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/LessEqual/y*
T0	*
_output_shapes
: 
�
Flinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/SwitchSwitchDlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/LessEqualDlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/LessEqual*
T0
*
_output_shapes
: : 
�
Hlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/switch_tIdentityHlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/Switch:1*
T0
*
_output_shapes
: 
�
Hlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/switch_fIdentityFlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/Switch*
T0
*
_output_shapes
: 
�
Glinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/pred_idIdentityDlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/LessEqual*
T0
*
_output_shapes
: 
�
Slinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/zerosConstI^linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Vlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/NotEqualNotEqual_linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/NotEqual/Switch:1Slinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/zeros*
T0*
_output_shapes

:
�
]linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/NotEqual/SwitchSwitchIlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOpGlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/pred_id*
T0*\
_classR
PNloc:@linear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOp*(
_output_shapes
::
�
Rlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/CastCastVlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/NotEqual*

SrcT0
*
_output_shapes

:*

DstT0
�
Slinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/ConstConstI^linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/switch_t*
valueB"       *
dtype0*
_output_shapes
:
�
[linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/nonzero_countSumRlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/CastSlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/Const*
T0*
_output_shapes
: 
�
Dlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/CastCast[linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/nonzero_count*

SrcT0*
_output_shapes
: *

DstT0	
�
Ulinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/zerosConstI^linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/switch_f*
valueB
 *    *
dtype0*
_output_shapes
: 
�
Xlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/NotEqualNotEqual_linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/NotEqual/SwitchUlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/zeros*
T0*
_output_shapes

:
�
_linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/NotEqual/SwitchSwitchIlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOpGlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/pred_id*
T0*\
_classR
PNloc:@linear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOp*(
_output_shapes
::
�
Tlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/CastCastXlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/NotEqual*

SrcT0
*
_output_shapes

:*

DstT0	
�
Ulinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/ConstConstI^linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/switch_f*
valueB"       *
dtype0*
_output_shapes
:
�
]linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/nonzero_countSumTlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/CastUlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/Const*
T0	*
_output_shapes
: 
�
Elinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/MergeMerge]linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/nonzero_countDlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/Cast*
T0	*
N*
_output_shapes
: : 
�
Qlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/counts_to_fraction/subSub?linear/zero_fraction/total_zero/zero_count_5/zero_fraction/SizeElinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/Merge*
T0	*
_output_shapes
: 
�
Rlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/counts_to_fraction/CastCastQlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/counts_to_fraction/sub*

SrcT0	*
_output_shapes
: *

DstT0
�
Tlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/counts_to_fraction/Cast_1Cast?linear/zero_fraction/total_zero/zero_count_5/zero_fraction/Size*

SrcT0	*
_output_shapes
: *

DstT0
�
Ulinear/zero_fraction/total_zero/zero_count_5/zero_fraction/counts_to_fraction/truedivRealDivRlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/counts_to_fraction/CastTlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/counts_to_fraction/Cast_1*
T0*
_output_shapes
: 
�
Clinear/zero_fraction/total_zero/zero_count_5/zero_fraction/fractionIdentityUlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/counts_to_fraction/truediv*
T0*
_output_shapes
: 
�
4linear/zero_fraction/total_zero/zero_count_5/ToFloatCast;linear/zero_fraction/total_zero/zero_count_5/ToFloat/Switch*

SrcT0	*
_output_shapes
: *

DstT0
�
;linear/zero_fraction/total_zero/zero_count_5/ToFloat/SwitchSwitch&linear/zero_fraction/total_size/Size_54linear/zero_fraction/total_zero/zero_count_5/pred_id*
T0	*9
_class/
-+loc:@linear/zero_fraction/total_size/Size_5*
_output_shapes
: : 
�
0linear/zero_fraction/total_zero/zero_count_5/mulMulClinear/zero_fraction/total_zero/zero_count_5/zero_fraction/fraction4linear/zero_fraction/total_zero/zero_count_5/ToFloat*
T0*
_output_shapes
: 
�
2linear/zero_fraction/total_zero/zero_count_5/MergeMerge0linear/zero_fraction/total_zero/zero_count_5/mul2linear/zero_fraction/total_zero/zero_count_5/Const*
T0*
N*
_output_shapes
: : 
�
$linear/zero_fraction/total_zero/AddNAddN0linear/zero_fraction/total_zero/zero_count/Merge2linear/zero_fraction/total_zero/zero_count_1/Merge2linear/zero_fraction/total_zero/zero_count_2/Merge2linear/zero_fraction/total_zero/zero_count_3/Merge2linear/zero_fraction/total_zero/zero_count_4/Merge2linear/zero_fraction/total_zero/zero_count_5/Merge*
T0*
N*
_output_shapes
: 
�
)linear/zero_fraction/compute/float32_sizeCast$linear/zero_fraction/total_size/AddN*

SrcT0	*
_output_shapes
: *

DstT0
�
$linear/zero_fraction/compute/truedivRealDiv$linear/zero_fraction/total_zero/AddN)linear/zero_fraction/compute/float32_size*
T0*
_output_shapes
: 
|
)linear/zero_fraction/zero_fraction_or_nanIdentity$linear/zero_fraction/compute/truediv*
T0*
_output_shapes
: 
�
$linear/fraction_of_zero_weights/tagsConst*0
value'B% Blinear/fraction_of_zero_weights*
dtype0*
_output_shapes
: 
�
linear/fraction_of_zero_weightsScalarSummary$linear/fraction_of_zero_weights/tags)linear/zero_fraction/zero_fraction_or_nan*
T0*
_output_shapes
: 
�
$linear/head/predictions/logits/ShapeShape:linear/linear_model/linear_model/linear_model/weighted_sum*
T0*
_output_shapes
:
z
8linear/head/predictions/logits/assert_rank_at_least/rankConst*
value	B :*
dtype0*
_output_shapes
: 
j
blinear/head/predictions/logits/assert_rank_at_least/assert_type/statically_determined_correct_typeNoOp
[
Slinear/head/predictions/logits/assert_rank_at_least/static_checks_determined_all_okNoOp
�
 linear/head/predictions/logisticSigmoid:linear/linear_model/linear_model/linear_model/weighted_sum*
T0*'
_output_shapes
:���������
�
"linear/head/predictions/zeros_like	ZerosLike:linear/linear_model/linear_model/linear_model/weighted_sum*
T0*'
_output_shapes
:���������
x
-linear/head/predictions/two_class_logits/axisConst*
valueB :
���������*
dtype0*
_output_shapes
: 
�
(linear/head/predictions/two_class_logitsConcatV2"linear/head/predictions/zeros_like:linear/linear_model/linear_model/linear_model/weighted_sum-linear/head/predictions/two_class_logits/axis*
T0*
N*'
_output_shapes
:���������
�
%linear/head/predictions/probabilitiesSoftmax(linear/head/predictions/two_class_logits*
T0*'
_output_shapes
:���������
v
+linear/head/predictions/class_ids/dimensionConst*
valueB :
���������*
dtype0*
_output_shapes
: 
�
!linear/head/predictions/class_idsArgMax(linear/head/predictions/two_class_logits+linear/head/predictions/class_ids/dimension*
T0*#
_output_shapes
:���������
q
&linear/head/predictions/ExpandDims/dimConst*
valueB :
���������*
dtype0*
_output_shapes
: 
�
"linear/head/predictions/ExpandDims
ExpandDims!linear/head/predictions/class_ids&linear/head/predictions/ExpandDims/dim*
T0	*'
_output_shapes
:���������
�
linear/head/predictions/ShapeShape:linear/linear_model/linear_model/linear_model/weighted_sum*
T0*
_output_shapes
:
u
+linear/head/predictions/strided_slice/stackConst*
valueB: *
dtype0*
_output_shapes
:
w
-linear/head/predictions/strided_slice/stack_1Const*
valueB:*
dtype0*
_output_shapes
:
w
-linear/head/predictions/strided_slice/stack_2Const*
valueB:*
dtype0*
_output_shapes
:
�
%linear/head/predictions/strided_sliceStridedSlicelinear/head/predictions/Shape+linear/head/predictions/strided_slice/stack-linear/head/predictions/strided_slice/stack_1-linear/head/predictions/strided_slice/stack_2*
shrink_axis_mask*
Index0*
T0*
_output_shapes
: 
e
#linear/head/predictions/range/startConst*
value	B : *
dtype0*
_output_shapes
: 
e
#linear/head/predictions/range/limitConst*
value	B :*
dtype0*
_output_shapes
: 
e
#linear/head/predictions/range/deltaConst*
value	B :*
dtype0*
_output_shapes
: 
�
linear/head/predictions/rangeRange#linear/head/predictions/range/start#linear/head/predictions/range/limit#linear/head/predictions/range/delta*
_output_shapes
:
j
(linear/head/predictions/ExpandDims_1/dimConst*
value	B : *
dtype0*
_output_shapes
: 
�
$linear/head/predictions/ExpandDims_1
ExpandDimslinear/head/predictions/range(linear/head/predictions/ExpandDims_1/dim*
T0*
_output_shapes

:
j
(linear/head/predictions/Tile/multiples/1Const*
value	B :*
dtype0*
_output_shapes
: 
�
&linear/head/predictions/Tile/multiplesPack%linear/head/predictions/strided_slice(linear/head/predictions/Tile/multiples/1*
T0*
N*
_output_shapes
:
�
linear/head/predictions/TileTile$linear/head/predictions/ExpandDims_1&linear/head/predictions/Tile/multiples*
T0*'
_output_shapes
:���������
�
linear/head/predictions/Shape_1Shape:linear/linear_model/linear_model/linear_model/weighted_sum*
T0*
_output_shapes
:
w
-linear/head/predictions/strided_slice_1/stackConst*
valueB: *
dtype0*
_output_shapes
:
y
/linear/head/predictions/strided_slice_1/stack_1Const*
valueB:*
dtype0*
_output_shapes
:
y
/linear/head/predictions/strided_slice_1/stack_2Const*
valueB:*
dtype0*
_output_shapes
:
�
'linear/head/predictions/strided_slice_1StridedSlicelinear/head/predictions/Shape_1-linear/head/predictions/strided_slice_1/stack/linear/head/predictions/strided_slice_1/stack_1/linear/head/predictions/strided_slice_1/stack_2*
shrink_axis_mask*
Index0*
T0*
_output_shapes
: 
g
%linear/head/predictions/range_1/startConst*
value	B : *
dtype0*
_output_shapes
: 
g
%linear/head/predictions/range_1/limitConst*
value	B :*
dtype0*
_output_shapes
: 
g
%linear/head/predictions/range_1/deltaConst*
value	B :*
dtype0*
_output_shapes
: 
�
linear/head/predictions/range_1Range%linear/head/predictions/range_1/start%linear/head/predictions/range_1/limit%linear/head/predictions/range_1/delta*
_output_shapes
:
r
 linear/head/predictions/AsStringAsStringlinear/head/predictions/range_1*
T0*
_output_shapes
:
j
(linear/head/predictions/ExpandDims_2/dimConst*
value	B : *
dtype0*
_output_shapes
: 
�
$linear/head/predictions/ExpandDims_2
ExpandDims linear/head/predictions/AsString(linear/head/predictions/ExpandDims_2/dim*
T0*
_output_shapes

:
l
*linear/head/predictions/Tile_1/multiples/1Const*
value	B :*
dtype0*
_output_shapes
: 
�
(linear/head/predictions/Tile_1/multiplesPack'linear/head/predictions/strided_slice_1*linear/head/predictions/Tile_1/multiples/1*
T0*
N*
_output_shapes
:
�
linear/head/predictions/Tile_1Tile$linear/head/predictions/ExpandDims_2(linear/head/predictions/Tile_1/multiples*
T0*'
_output_shapes
:���������
�
#linear/head/predictions/str_classesAsString"linear/head/predictions/ExpandDims*
T0	*'
_output_shapes
:���������
f
linear/head/ShapeShape%linear/head/predictions/probabilities*
T0*
_output_shapes
:
i
linear/head/strided_slice/stackConst*
valueB: *
dtype0*
_output_shapes
:
k
!linear/head/strided_slice/stack_1Const*
valueB:*
dtype0*
_output_shapes
:
k
!linear/head/strided_slice/stack_2Const*
valueB:*
dtype0*
_output_shapes
:
�
linear/head/strided_sliceStridedSlicelinear/head/Shapelinear/head/strided_slice/stack!linear/head/strided_slice/stack_1!linear/head/strided_slice/stack_2*
shrink_axis_mask*
Index0*
T0*
_output_shapes
: 
Y
linear/head/range/startConst*
value	B : *
dtype0*
_output_shapes
: 
Y
linear/head/range/limitConst*
value	B :*
dtype0*
_output_shapes
: 
Y
linear/head/range/deltaConst*
value	B :*
dtype0*
_output_shapes
: 
�
linear/head/rangeRangelinear/head/range/startlinear/head/range/limitlinear/head/range/delta*
_output_shapes
:
X
linear/head/AsStringAsStringlinear/head/range*
T0*
_output_shapes
:
\
linear/head/ExpandDims/dimConst*
value	B : *
dtype0*
_output_shapes
: 

linear/head/ExpandDims
ExpandDimslinear/head/AsStringlinear/head/ExpandDims/dim*
T0*
_output_shapes

:
^
linear/head/Tile/multiples/1Const*
value	B :*
dtype0*
_output_shapes
: 
�
linear/head/Tile/multiplesPacklinear/head/strided_slicelinear/head/Tile/multiples/1*
T0*
N*
_output_shapes
:
~
linear/head/TileTilelinear/head/ExpandDimslinear/head/Tile/multiples*
T0*'
_output_shapes
:���������

initNoOp

init_all_tablesNoOp

init_1NoOp
4

group_depsNoOp^init^init_1^init_all_tables
Y
save/filename/inputConst*
valueB Bmodel*
dtype0*
_output_shapes
: 
n
save/filenamePlaceholderWithDefaultsave/filename/input*
shape: *
dtype0*
_output_shapes
: 
e

save/ConstPlaceholderWithDefaultsave/filename*
shape: *
dtype0*
_output_shapes
: 
|
save/Read/ReadVariableOpReadVariableOp'linear/linear_model/bias_weights/part_0*
dtype0*
_output_shapes
:
X
save/IdentityIdentitysave/Read/ReadVariableOp*
T0*
_output_shapes
:
^
save/Identity_1Identitysave/Identity"/device:CPU:0*
T0*
_output_shapes
:
�
save/Read_1/ReadVariableOpReadVariableOp.linear/linear_model/direction.x/weights/part_0*
dtype0*
_output_shapes

:
`
save/Identity_2Identitysave/Read_1/ReadVariableOp*
T0*
_output_shapes

:
d
save/Identity_3Identitysave/Identity_2"/device:CPU:0*
T0*
_output_shapes

:
�
save/Read_2/ReadVariableOpReadVariableOp.linear/linear_model/direction.y/weights/part_0*
dtype0*
_output_shapes

:
`
save/Identity_4Identitysave/Read_2/ReadVariableOp*
T0*
_output_shapes

:
d
save/Identity_5Identitysave/Identity_4"/device:CPU:0*
T0*
_output_shapes

:
�
save/Read_3/ReadVariableOpReadVariableOp.linear/linear_model/direction.z/weights/part_0*
dtype0*
_output_shapes

:
`
save/Identity_6Identitysave/Read_3/ReadVariableOp*
T0*
_output_shapes

:
d
save/Identity_7Identitysave/Identity_6"/device:CPU:0*
T0*
_output_shapes

:
�
save/Read_4/ReadVariableOpReadVariableOp+linear/linear_model/origin.x/weights/part_0*
dtype0*
_output_shapes

:
`
save/Identity_8Identitysave/Read_4/ReadVariableOp*
T0*
_output_shapes

:
d
save/Identity_9Identitysave/Identity_8"/device:CPU:0*
T0*
_output_shapes

:
�
save/Read_5/ReadVariableOpReadVariableOp+linear/linear_model/origin.y/weights/part_0*
dtype0*
_output_shapes

:
a
save/Identity_10Identitysave/Read_5/ReadVariableOp*
T0*
_output_shapes

:
f
save/Identity_11Identitysave/Identity_10"/device:CPU:0*
T0*
_output_shapes

:
�
save/Read_6/ReadVariableOpReadVariableOp+linear/linear_model/origin.z/weights/part_0*
dtype0*
_output_shapes

:
a
save/Identity_12Identitysave/Read_6/ReadVariableOp*
T0*
_output_shapes

:
f
save/Identity_13Identitysave/Identity_12"/device:CPU:0*
T0*
_output_shapes

:
�
save/StringJoin/inputs_1Const*<
value3B1 B+_temp_87854219b5924ebd838ac2c8e09c5fa7/part*
dtype0*
_output_shapes
: 
d
save/StringJoin
StringJoin
save/Constsave/StringJoin/inputs_1*
N*
_output_shapes
: 
Q
save/num_shardsConst*
value	B :*
dtype0*
_output_shapes
: 
k
save/ShardedFilename/shardConst"/device:CPU:0*
value	B : *
dtype0*
_output_shapes
: 
�
save/ShardedFilenameShardedFilenamesave/StringJoinsave/ShardedFilename/shardsave/num_shards"/device:CPU:0*
_output_shapes
: 
{
save/SaveV2/tensor_namesConst"/device:CPU:0* 
valueBBglobal_step*
dtype0*
_output_shapes
:
t
save/SaveV2/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:
�
save/SaveV2SaveV2save/ShardedFilenamesave/SaveV2/tensor_namessave/SaveV2/shape_and_slicesglobal_step"/device:CPU:0*
dtypes
2	
�
save/control_dependencyIdentitysave/ShardedFilename^save/SaveV2"/device:CPU:0*
T0*'
_class
loc:@save/ShardedFilename*
_output_shapes
: 
m
save/ShardedFilename_1/shardConst"/device:CPU:0*
value	B :*
dtype0*
_output_shapes
: 
�
save/ShardedFilename_1ShardedFilenamesave/StringJoinsave/ShardedFilename_1/shardsave/num_shards"/device:CPU:0*
_output_shapes
: 
�
save/Read_7/ReadVariableOpReadVariableOp'linear/linear_model/bias_weights/part_0"/device:CPU:0*
dtype0*
_output_shapes
:
l
save/Identity_14Identitysave/Read_7/ReadVariableOp"/device:CPU:0*
T0*
_output_shapes
:
b
save/Identity_15Identitysave/Identity_14"/device:CPU:0*
T0*
_output_shapes
:
�
save/Read_8/ReadVariableOpReadVariableOp.linear/linear_model/direction.x/weights/part_0"/device:CPU:0*
dtype0*
_output_shapes

:
p
save/Identity_16Identitysave/Read_8/ReadVariableOp"/device:CPU:0*
T0*
_output_shapes

:
f
save/Identity_17Identitysave/Identity_16"/device:CPU:0*
T0*
_output_shapes

:
�
save/Read_9/ReadVariableOpReadVariableOp.linear/linear_model/direction.y/weights/part_0"/device:CPU:0*
dtype0*
_output_shapes

:
p
save/Identity_18Identitysave/Read_9/ReadVariableOp"/device:CPU:0*
T0*
_output_shapes

:
f
save/Identity_19Identitysave/Identity_18"/device:CPU:0*
T0*
_output_shapes

:
�
save/Read_10/ReadVariableOpReadVariableOp.linear/linear_model/direction.z/weights/part_0"/device:CPU:0*
dtype0*
_output_shapes

:
q
save/Identity_20Identitysave/Read_10/ReadVariableOp"/device:CPU:0*
T0*
_output_shapes

:
f
save/Identity_21Identitysave/Identity_20"/device:CPU:0*
T0*
_output_shapes

:
�
save/Read_11/ReadVariableOpReadVariableOp+linear/linear_model/origin.x/weights/part_0"/device:CPU:0*
dtype0*
_output_shapes

:
q
save/Identity_22Identitysave/Read_11/ReadVariableOp"/device:CPU:0*
T0*
_output_shapes

:
f
save/Identity_23Identitysave/Identity_22"/device:CPU:0*
T0*
_output_shapes

:
�
save/Read_12/ReadVariableOpReadVariableOp+linear/linear_model/origin.y/weights/part_0"/device:CPU:0*
dtype0*
_output_shapes

:
q
save/Identity_24Identitysave/Read_12/ReadVariableOp"/device:CPU:0*
T0*
_output_shapes

:
f
save/Identity_25Identitysave/Identity_24"/device:CPU:0*
T0*
_output_shapes

:
�
save/Read_13/ReadVariableOpReadVariableOp+linear/linear_model/origin.z/weights/part_0"/device:CPU:0*
dtype0*
_output_shapes

:
q
save/Identity_26Identitysave/Read_13/ReadVariableOp"/device:CPU:0*
T0*
_output_shapes

:
f
save/Identity_27Identitysave/Identity_26"/device:CPU:0*
T0*
_output_shapes

:
�
save/SaveV2_1/tensor_namesConst"/device:CPU:0*�
value�B�B linear/linear_model/bias_weightsB'linear/linear_model/direction.x/weightsB'linear/linear_model/direction.y/weightsB'linear/linear_model/direction.z/weightsB$linear/linear_model/origin.x/weightsB$linear/linear_model/origin.y/weightsB$linear/linear_model/origin.z/weights*
dtype0*
_output_shapes
:
�
save/SaveV2_1/shape_and_slicesConst"/device:CPU:0*h
value_B]B1 0,1B1 1 0,1:0,1B1 1 0,1:0,1B1 1 0,1:0,1B1 1 0,1:0,1B1 1 0,1:0,1B1 1 0,1:0,1*
dtype0*
_output_shapes
:
�
save/SaveV2_1SaveV2save/ShardedFilename_1save/SaveV2_1/tensor_namessave/SaveV2_1/shape_and_slicessave/Identity_15save/Identity_17save/Identity_19save/Identity_21save/Identity_23save/Identity_25save/Identity_27"/device:CPU:0*
dtypes
	2
�
save/control_dependency_1Identitysave/ShardedFilename_1^save/SaveV2_1"/device:CPU:0*
T0*)
_class
loc:@save/ShardedFilename_1*
_output_shapes
: 
�
+save/MergeV2Checkpoints/checkpoint_prefixesPacksave/ShardedFilenamesave/ShardedFilename_1^save/control_dependency^save/control_dependency_1"/device:CPU:0*
T0*
N*
_output_shapes
:
u
save/MergeV2CheckpointsMergeV2Checkpoints+save/MergeV2Checkpoints/checkpoint_prefixes
save/Const"/device:CPU:0
�
save/Identity_28Identity
save/Const^save/MergeV2Checkpoints^save/control_dependency^save/control_dependency_1"/device:CPU:0*
T0*
_output_shapes
: 
~
save/RestoreV2/tensor_namesConst"/device:CPU:0* 
valueBBglobal_step*
dtype0*
_output_shapes
:
w
save/RestoreV2/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:
�
save/RestoreV2	RestoreV2
save/Constsave/RestoreV2/tensor_namessave/RestoreV2/shape_and_slices"/device:CPU:0*
dtypes
2	*
_output_shapes
:
s
save/AssignAssignglobal_stepsave/RestoreV2*
T0	*
_class
loc:@global_step*
_output_shapes
: 
(
save/restore_shardNoOp^save/Assign
�
save/RestoreV2_1/tensor_namesConst"/device:CPU:0*�
value�B�B linear/linear_model/bias_weightsB'linear/linear_model/direction.x/weightsB'linear/linear_model/direction.y/weightsB'linear/linear_model/direction.z/weightsB$linear/linear_model/origin.x/weightsB$linear/linear_model/origin.y/weightsB$linear/linear_model/origin.z/weights*
dtype0*
_output_shapes
:
�
!save/RestoreV2_1/shape_and_slicesConst"/device:CPU:0*h
value_B]B1 0,1B1 1 0,1:0,1B1 1 0,1:0,1B1 1 0,1:0,1B1 1 0,1:0,1B1 1 0,1:0,1B1 1 0,1:0,1*
dtype0*
_output_shapes
:
�
save/RestoreV2_1	RestoreV2
save/Constsave/RestoreV2_1/tensor_names!save/RestoreV2_1/shape_and_slices"/device:CPU:0*
dtypes
	2*V
_output_shapesD
B:::::::
b
save/Identity_29Identitysave/RestoreV2_1"/device:CPU:0*
T0*
_output_shapes
:
�
save/AssignVariableOpAssignVariableOp'linear/linear_model/bias_weights/part_0save/Identity_29"/device:CPU:0*
dtype0
h
save/Identity_30Identitysave/RestoreV2_1:1"/device:CPU:0*
T0*
_output_shapes

:
�
save/AssignVariableOp_1AssignVariableOp.linear/linear_model/direction.x/weights/part_0save/Identity_30"/device:CPU:0*
dtype0
h
save/Identity_31Identitysave/RestoreV2_1:2"/device:CPU:0*
T0*
_output_shapes

:
�
save/AssignVariableOp_2AssignVariableOp.linear/linear_model/direction.y/weights/part_0save/Identity_31"/device:CPU:0*
dtype0
h
save/Identity_32Identitysave/RestoreV2_1:3"/device:CPU:0*
T0*
_output_shapes

:
�
save/AssignVariableOp_3AssignVariableOp.linear/linear_model/direction.z/weights/part_0save/Identity_32"/device:CPU:0*
dtype0
h
save/Identity_33Identitysave/RestoreV2_1:4"/device:CPU:0*
T0*
_output_shapes

:
�
save/AssignVariableOp_4AssignVariableOp+linear/linear_model/origin.x/weights/part_0save/Identity_33"/device:CPU:0*
dtype0
h
save/Identity_34Identitysave/RestoreV2_1:5"/device:CPU:0*
T0*
_output_shapes

:
�
save/AssignVariableOp_5AssignVariableOp+linear/linear_model/origin.y/weights/part_0save/Identity_34"/device:CPU:0*
dtype0
h
save/Identity_35Identitysave/RestoreV2_1:6"/device:CPU:0*
T0*
_output_shapes

:
�
save/AssignVariableOp_6AssignVariableOp+linear/linear_model/origin.z/weights/part_0save/Identity_35"/device:CPU:0*
dtype0
�
save/restore_shard_1NoOp^save/AssignVariableOp^save/AssignVariableOp_1^save/AssignVariableOp_2^save/AssignVariableOp_3^save/AssignVariableOp_4^save/AssignVariableOp_5^save/AssignVariableOp_6"/device:CPU:0
2
save/restore_all/NoOpNoOp^save/restore_shard
E
save/restore_all/NoOp_1NoOp^save/restore_shard_1"/device:CPU:0
J
save/restore_allNoOp^save/restore_all/NoOp^save/restore_all/NoOp_1"&?
save/Const:0save/Identity_28:0save/restore_all (5 @F8"�
trainable_variables��
�
0linear/linear_model/direction.x/weights/part_0:05linear/linear_model/direction.x/weights/part_0/AssignDlinear/linear_model/direction.x/weights/part_0/Read/ReadVariableOp:0"5
'linear/linear_model/direction.x/weights  "(2Blinear/linear_model/direction.x/weights/part_0/Initializer/zeros:08
�
0linear/linear_model/direction.y/weights/part_0:05linear/linear_model/direction.y/weights/part_0/AssignDlinear/linear_model/direction.y/weights/part_0/Read/ReadVariableOp:0"5
'linear/linear_model/direction.y/weights  "(2Blinear/linear_model/direction.y/weights/part_0/Initializer/zeros:08
�
0linear/linear_model/direction.z/weights/part_0:05linear/linear_model/direction.z/weights/part_0/AssignDlinear/linear_model/direction.z/weights/part_0/Read/ReadVariableOp:0"5
'linear/linear_model/direction.z/weights  "(2Blinear/linear_model/direction.z/weights/part_0/Initializer/zeros:08
�
-linear/linear_model/origin.x/weights/part_0:02linear/linear_model/origin.x/weights/part_0/AssignAlinear/linear_model/origin.x/weights/part_0/Read/ReadVariableOp:0"2
$linear/linear_model/origin.x/weights  "(2?linear/linear_model/origin.x/weights/part_0/Initializer/zeros:08
�
-linear/linear_model/origin.y/weights/part_0:02linear/linear_model/origin.y/weights/part_0/AssignAlinear/linear_model/origin.y/weights/part_0/Read/ReadVariableOp:0"2
$linear/linear_model/origin.y/weights  "(2?linear/linear_model/origin.y/weights/part_0/Initializer/zeros:08
�
-linear/linear_model/origin.z/weights/part_0:02linear/linear_model/origin.z/weights/part_0/AssignAlinear/linear_model/origin.z/weights/part_0/Read/ReadVariableOp:0"2
$linear/linear_model/origin.z/weights  "(2?linear/linear_model/origin.z/weights/part_0/Initializer/zeros:08
�
)linear/linear_model/bias_weights/part_0:0.linear/linear_model/bias_weights/part_0/Assign=linear/linear_model/bias_weights/part_0/Read/ReadVariableOp:0"+
 linear/linear_model/bias_weights "(2;linear/linear_model/bias_weights/part_0/Initializer/zeros:08"A
	summaries4
2
linear/bias:0
!linear/fraction_of_zero_weights:0"�
	variables��
Z
global_step:0global_step/Assignglobal_step/read:02global_step/Initializer/zeros:0H
�
0linear/linear_model/direction.x/weights/part_0:05linear/linear_model/direction.x/weights/part_0/AssignDlinear/linear_model/direction.x/weights/part_0/Read/ReadVariableOp:0"5
'linear/linear_model/direction.x/weights  "(2Blinear/linear_model/direction.x/weights/part_0/Initializer/zeros:08
�
0linear/linear_model/direction.y/weights/part_0:05linear/linear_model/direction.y/weights/part_0/AssignDlinear/linear_model/direction.y/weights/part_0/Read/ReadVariableOp:0"5
'linear/linear_model/direction.y/weights  "(2Blinear/linear_model/direction.y/weights/part_0/Initializer/zeros:08
�
0linear/linear_model/direction.z/weights/part_0:05linear/linear_model/direction.z/weights/part_0/AssignDlinear/linear_model/direction.z/weights/part_0/Read/ReadVariableOp:0"5
'linear/linear_model/direction.z/weights  "(2Blinear/linear_model/direction.z/weights/part_0/Initializer/zeros:08
�
-linear/linear_model/origin.x/weights/part_0:02linear/linear_model/origin.x/weights/part_0/AssignAlinear/linear_model/origin.x/weights/part_0/Read/ReadVariableOp:0"2
$linear/linear_model/origin.x/weights  "(2?linear/linear_model/origin.x/weights/part_0/Initializer/zeros:08
�
-linear/linear_model/origin.y/weights/part_0:02linear/linear_model/origin.y/weights/part_0/AssignAlinear/linear_model/origin.y/weights/part_0/Read/ReadVariableOp:0"2
$linear/linear_model/origin.y/weights  "(2?linear/linear_model/origin.y/weights/part_0/Initializer/zeros:08
�
-linear/linear_model/origin.z/weights/part_0:02linear/linear_model/origin.z/weights/part_0/AssignAlinear/linear_model/origin.z/weights/part_0/Read/ReadVariableOp:0"2
$linear/linear_model/origin.z/weights  "(2?linear/linear_model/origin.z/weights/part_0/Initializer/zeros:08
�
)linear/linear_model/bias_weights/part_0:0.linear/linear_model/bias_weights/part_0/Assign=linear/linear_model/bias_weights/part_0/Read/ReadVariableOp:0"+
 linear/linear_model/bias_weights "(2;linear/linear_model/bias_weights/part_0/Initializer/zeros:08"m
global_step^\
Z
global_step:0global_step/Assignglobal_step/read:02global_step/Initializer/zeros:0H"̴
cond_context����
�
4linear/zero_fraction/total_zero/zero_count/cond_text4linear/zero_fraction/total_zero/zero_count/pred_id:05linear/zero_fraction/total_zero/zero_count/switch_t:0 *�
2linear/zero_fraction/total_zero/zero_count/Const:0
4linear/zero_fraction/total_zero/zero_count/pred_id:0
5linear/zero_fraction/total_zero/zero_count/switch_t:0l
4linear/zero_fraction/total_zero/zero_count/pred_id:04linear/zero_fraction/total_zero/zero_count/pred_id:0
�.
6linear/zero_fraction/total_zero/zero_count/cond_text_14linear/zero_fraction/total_zero/zero_count/pred_id:05linear/zero_fraction/total_zero/zero_count/switch_f:0*�
0linear/linear_model/direction.x/weights/part_0:0
&linear/zero_fraction/total_size/Size:0
;linear/zero_fraction/total_zero/zero_count/ToFloat/Switch:0
4linear/zero_fraction/total_zero/zero_count/ToFloat:0
0linear/zero_fraction/total_zero/zero_count/mul:0
4linear/zero_fraction/total_zero/zero_count/pred_id:0
5linear/zero_fraction/total_zero/zero_count/switch_f:0
Flinear/zero_fraction/total_zero/zero_count/zero_fraction/LessEqual/y:0
Dlinear/zero_fraction/total_zero/zero_count/zero_fraction/LessEqual:0
Plinear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOp/Switch:0
Ilinear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOp:0
?linear/zero_fraction/total_zero/zero_count/zero_fraction/Size:0
Dlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/Cast:0
Elinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/Merge:0
Elinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/Merge:1
Flinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/Switch:0
Flinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/Switch:1
Rlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/Cast:0
Slinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/Const:0
]linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/NotEqual/Switch:1
Vlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/NotEqual:0
[linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/nonzero_count:0
Slinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/zeros:0
Tlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/Cast:0
Ulinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/Const:0
_linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
Xlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/NotEqual:0
]linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/nonzero_count:0
Ulinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/zeros:0
Glinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/pred_id:0
Hlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/switch_f:0
Hlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/switch_t:0
Rlinear/zero_fraction/total_zero/zero_count/zero_fraction/counts_to_fraction/Cast:0
Tlinear/zero_fraction/total_zero/zero_count/zero_fraction/counts_to_fraction/Cast_1:0
Qlinear/zero_fraction/total_zero/zero_count/zero_fraction/counts_to_fraction/sub:0
Ulinear/zero_fraction/total_zero/zero_count/zero_fraction/counts_to_fraction/truediv:0
Clinear/zero_fraction/total_zero/zero_count/zero_fraction/fraction:0l
4linear/zero_fraction/total_zero/zero_count/pred_id:04linear/zero_fraction/total_zero/zero_count/pred_id:0�
0linear/linear_model/direction.x/weights/part_0:0Plinear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOp/Switch:0e
&linear/zero_fraction/total_size/Size:0;linear/zero_fraction/total_zero/zero_count/ToFloat/Switch:02�

�

Glinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/cond_textGlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/pred_id:0Hlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/switch_t:0 *�
Ilinear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOp:0
Dlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/Cast:0
Rlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/Cast:0
Slinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/Const:0
]linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/NotEqual/Switch:1
Vlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/NotEqual:0
[linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/nonzero_count:0
Slinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/zeros:0
Glinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/pred_id:0
Hlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/switch_t:0�
Ilinear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOp:0]linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero/NotEqual/Switch:1�
Glinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/pred_id:0Glinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/pred_id:02�

�

Ilinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/cond_text_1Glinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/pred_id:0Hlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/switch_f:0*�
Ilinear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOp:0
Tlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/Cast:0
Ulinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/Const:0
_linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
Xlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/NotEqual:0
]linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/nonzero_count:0
Ulinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/zeros:0
Glinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/pred_id:0
Hlinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/switch_f:0�
Ilinear/zero_fraction/total_zero/zero_count/zero_fraction/ReadVariableOp:0_linear/zero_fraction/total_zero/zero_count/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0�
Glinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/pred_id:0Glinear/zero_fraction/total_zero/zero_count/zero_fraction/cond/pred_id:0
�
6linear/zero_fraction/total_zero/zero_count_1/cond_text6linear/zero_fraction/total_zero/zero_count_1/pred_id:07linear/zero_fraction/total_zero/zero_count_1/switch_t:0 *�
4linear/zero_fraction/total_zero/zero_count_1/Const:0
6linear/zero_fraction/total_zero/zero_count_1/pred_id:0
7linear/zero_fraction/total_zero/zero_count_1/switch_t:0p
6linear/zero_fraction/total_zero/zero_count_1/pred_id:06linear/zero_fraction/total_zero/zero_count_1/pred_id:0
�0
8linear/zero_fraction/total_zero/zero_count_1/cond_text_16linear/zero_fraction/total_zero/zero_count_1/pred_id:07linear/zero_fraction/total_zero/zero_count_1/switch_f:0*�
0linear/linear_model/direction.y/weights/part_0:0
(linear/zero_fraction/total_size/Size_1:0
=linear/zero_fraction/total_zero/zero_count_1/ToFloat/Switch:0
6linear/zero_fraction/total_zero/zero_count_1/ToFloat:0
2linear/zero_fraction/total_zero/zero_count_1/mul:0
6linear/zero_fraction/total_zero/zero_count_1/pred_id:0
7linear/zero_fraction/total_zero/zero_count_1/switch_f:0
Hlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/LessEqual/y:0
Flinear/zero_fraction/total_zero/zero_count_1/zero_fraction/LessEqual:0
Rlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOp/Switch:0
Klinear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOp:0
Alinear/zero_fraction/total_zero/zero_count_1/zero_fraction/Size:0
Flinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/Cast:0
Glinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/Merge:0
Glinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/Merge:1
Hlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/Switch:0
Hlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/Switch:1
Tlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/Cast:0
Ulinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/Const:0
_linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/NotEqual/Switch:1
Xlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/NotEqual:0
]linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/nonzero_count:0
Ulinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/zeros:0
Vlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/Cast:0
Wlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/Const:0
alinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
Zlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/NotEqual:0
_linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/nonzero_count:0
Wlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/switch_f:0
Jlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/switch_t:0
Tlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/counts_to_fraction/Cast:0
Vlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/counts_to_fraction/Cast_1:0
Slinear/zero_fraction/total_zero/zero_count_1/zero_fraction/counts_to_fraction/sub:0
Wlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/counts_to_fraction/truediv:0
Elinear/zero_fraction/total_zero/zero_count_1/zero_fraction/fraction:0i
(linear/zero_fraction/total_size/Size_1:0=linear/zero_fraction/total_zero/zero_count_1/ToFloat/Switch:0p
6linear/zero_fraction/total_zero/zero_count_1/pred_id:06linear/zero_fraction/total_zero/zero_count_1/pred_id:0�
0linear/linear_model/direction.y/weights/part_0:0Rlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOp/Switch:02�

�

Ilinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/cond_textIlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/pred_id:0Jlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/switch_t:0 *�	
Klinear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOp:0
Flinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/Cast:0
Tlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/Cast:0
Ulinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/Const:0
_linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/NotEqual/Switch:1
Xlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/NotEqual:0
]linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/nonzero_count:0
Ulinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/switch_t:0�
Klinear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOp:0_linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero/NotEqual/Switch:1�
Ilinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/pred_id:0Ilinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/pred_id:02�

�

Klinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/cond_text_1Ilinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/pred_id:0Jlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/switch_f:0*�
Klinear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOp:0
Vlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/Cast:0
Wlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/Const:0
alinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
Zlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/NotEqual:0
_linear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/nonzero_count:0
Wlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/switch_f:0�
Klinear/zero_fraction/total_zero/zero_count_1/zero_fraction/ReadVariableOp:0alinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0�
Ilinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/pred_id:0Ilinear/zero_fraction/total_zero/zero_count_1/zero_fraction/cond/pred_id:0
�
6linear/zero_fraction/total_zero/zero_count_2/cond_text6linear/zero_fraction/total_zero/zero_count_2/pred_id:07linear/zero_fraction/total_zero/zero_count_2/switch_t:0 *�
4linear/zero_fraction/total_zero/zero_count_2/Const:0
6linear/zero_fraction/total_zero/zero_count_2/pred_id:0
7linear/zero_fraction/total_zero/zero_count_2/switch_t:0p
6linear/zero_fraction/total_zero/zero_count_2/pred_id:06linear/zero_fraction/total_zero/zero_count_2/pred_id:0
�0
8linear/zero_fraction/total_zero/zero_count_2/cond_text_16linear/zero_fraction/total_zero/zero_count_2/pred_id:07linear/zero_fraction/total_zero/zero_count_2/switch_f:0*�
0linear/linear_model/direction.z/weights/part_0:0
(linear/zero_fraction/total_size/Size_2:0
=linear/zero_fraction/total_zero/zero_count_2/ToFloat/Switch:0
6linear/zero_fraction/total_zero/zero_count_2/ToFloat:0
2linear/zero_fraction/total_zero/zero_count_2/mul:0
6linear/zero_fraction/total_zero/zero_count_2/pred_id:0
7linear/zero_fraction/total_zero/zero_count_2/switch_f:0
Hlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/LessEqual/y:0
Flinear/zero_fraction/total_zero/zero_count_2/zero_fraction/LessEqual:0
Rlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOp/Switch:0
Klinear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOp:0
Alinear/zero_fraction/total_zero/zero_count_2/zero_fraction/Size:0
Flinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/Cast:0
Glinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/Merge:0
Glinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/Merge:1
Hlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/Switch:0
Hlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/Switch:1
Tlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/Cast:0
Ulinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/Const:0
_linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/NotEqual/Switch:1
Xlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/NotEqual:0
]linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/nonzero_count:0
Ulinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/zeros:0
Vlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/Cast:0
Wlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/Const:0
alinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
Zlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/NotEqual:0
_linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/nonzero_count:0
Wlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/switch_f:0
Jlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/switch_t:0
Tlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/counts_to_fraction/Cast:0
Vlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/counts_to_fraction/Cast_1:0
Slinear/zero_fraction/total_zero/zero_count_2/zero_fraction/counts_to_fraction/sub:0
Wlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/counts_to_fraction/truediv:0
Elinear/zero_fraction/total_zero/zero_count_2/zero_fraction/fraction:0i
(linear/zero_fraction/total_size/Size_2:0=linear/zero_fraction/total_zero/zero_count_2/ToFloat/Switch:0�
0linear/linear_model/direction.z/weights/part_0:0Rlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOp/Switch:0p
6linear/zero_fraction/total_zero/zero_count_2/pred_id:06linear/zero_fraction/total_zero/zero_count_2/pred_id:02�

�

Ilinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/cond_textIlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/pred_id:0Jlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/switch_t:0 *�	
Klinear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOp:0
Flinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/Cast:0
Tlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/Cast:0
Ulinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/Const:0
_linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/NotEqual/Switch:1
Xlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/NotEqual:0
]linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/nonzero_count:0
Ulinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/switch_t:0�
Ilinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/pred_id:0Ilinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/pred_id:0�
Klinear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOp:0_linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero/NotEqual/Switch:12�

�

Klinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/cond_text_1Ilinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/pred_id:0Jlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/switch_f:0*�
Klinear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOp:0
Vlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/Cast:0
Wlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/Const:0
alinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
Zlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/NotEqual:0
_linear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/nonzero_count:0
Wlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/switch_f:0�
Ilinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/pred_id:0Ilinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/pred_id:0�
Klinear/zero_fraction/total_zero/zero_count_2/zero_fraction/ReadVariableOp:0alinear/zero_fraction/total_zero/zero_count_2/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
�
6linear/zero_fraction/total_zero/zero_count_3/cond_text6linear/zero_fraction/total_zero/zero_count_3/pred_id:07linear/zero_fraction/total_zero/zero_count_3/switch_t:0 *�
4linear/zero_fraction/total_zero/zero_count_3/Const:0
6linear/zero_fraction/total_zero/zero_count_3/pred_id:0
7linear/zero_fraction/total_zero/zero_count_3/switch_t:0p
6linear/zero_fraction/total_zero/zero_count_3/pred_id:06linear/zero_fraction/total_zero/zero_count_3/pred_id:0
�0
8linear/zero_fraction/total_zero/zero_count_3/cond_text_16linear/zero_fraction/total_zero/zero_count_3/pred_id:07linear/zero_fraction/total_zero/zero_count_3/switch_f:0*�
-linear/linear_model/origin.x/weights/part_0:0
(linear/zero_fraction/total_size/Size_3:0
=linear/zero_fraction/total_zero/zero_count_3/ToFloat/Switch:0
6linear/zero_fraction/total_zero/zero_count_3/ToFloat:0
2linear/zero_fraction/total_zero/zero_count_3/mul:0
6linear/zero_fraction/total_zero/zero_count_3/pred_id:0
7linear/zero_fraction/total_zero/zero_count_3/switch_f:0
Hlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/LessEqual/y:0
Flinear/zero_fraction/total_zero/zero_count_3/zero_fraction/LessEqual:0
Rlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOp/Switch:0
Klinear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOp:0
Alinear/zero_fraction/total_zero/zero_count_3/zero_fraction/Size:0
Flinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/Cast:0
Glinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/Merge:0
Glinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/Merge:1
Hlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/Switch:0
Hlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/Switch:1
Tlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/Cast:0
Ulinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/Const:0
_linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/NotEqual/Switch:1
Xlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/NotEqual:0
]linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/nonzero_count:0
Ulinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/zeros:0
Vlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/Cast:0
Wlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/Const:0
alinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
Zlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/NotEqual:0
_linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/nonzero_count:0
Wlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/switch_f:0
Jlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/switch_t:0
Tlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/counts_to_fraction/Cast:0
Vlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/counts_to_fraction/Cast_1:0
Slinear/zero_fraction/total_zero/zero_count_3/zero_fraction/counts_to_fraction/sub:0
Wlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/counts_to_fraction/truediv:0
Elinear/zero_fraction/total_zero/zero_count_3/zero_fraction/fraction:0i
(linear/zero_fraction/total_size/Size_3:0=linear/zero_fraction/total_zero/zero_count_3/ToFloat/Switch:0�
-linear/linear_model/origin.x/weights/part_0:0Rlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOp/Switch:0p
6linear/zero_fraction/total_zero/zero_count_3/pred_id:06linear/zero_fraction/total_zero/zero_count_3/pred_id:02�

�

Ilinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/cond_textIlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/pred_id:0Jlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/switch_t:0 *�	
Klinear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOp:0
Flinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/Cast:0
Tlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/Cast:0
Ulinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/Const:0
_linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/NotEqual/Switch:1
Xlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/NotEqual:0
]linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/nonzero_count:0
Ulinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/switch_t:0�
Klinear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOp:0_linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero/NotEqual/Switch:1�
Ilinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/pred_id:0Ilinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/pred_id:02�

�

Klinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/cond_text_1Ilinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/pred_id:0Jlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/switch_f:0*�
Klinear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOp:0
Vlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/Cast:0
Wlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/Const:0
alinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
Zlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/NotEqual:0
_linear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/nonzero_count:0
Wlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/switch_f:0�
Klinear/zero_fraction/total_zero/zero_count_3/zero_fraction/ReadVariableOp:0alinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0�
Ilinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/pred_id:0Ilinear/zero_fraction/total_zero/zero_count_3/zero_fraction/cond/pred_id:0
�
6linear/zero_fraction/total_zero/zero_count_4/cond_text6linear/zero_fraction/total_zero/zero_count_4/pred_id:07linear/zero_fraction/total_zero/zero_count_4/switch_t:0 *�
4linear/zero_fraction/total_zero/zero_count_4/Const:0
6linear/zero_fraction/total_zero/zero_count_4/pred_id:0
7linear/zero_fraction/total_zero/zero_count_4/switch_t:0p
6linear/zero_fraction/total_zero/zero_count_4/pred_id:06linear/zero_fraction/total_zero/zero_count_4/pred_id:0
�0
8linear/zero_fraction/total_zero/zero_count_4/cond_text_16linear/zero_fraction/total_zero/zero_count_4/pred_id:07linear/zero_fraction/total_zero/zero_count_4/switch_f:0*�
-linear/linear_model/origin.y/weights/part_0:0
(linear/zero_fraction/total_size/Size_4:0
=linear/zero_fraction/total_zero/zero_count_4/ToFloat/Switch:0
6linear/zero_fraction/total_zero/zero_count_4/ToFloat:0
2linear/zero_fraction/total_zero/zero_count_4/mul:0
6linear/zero_fraction/total_zero/zero_count_4/pred_id:0
7linear/zero_fraction/total_zero/zero_count_4/switch_f:0
Hlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/LessEqual/y:0
Flinear/zero_fraction/total_zero/zero_count_4/zero_fraction/LessEqual:0
Rlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOp/Switch:0
Klinear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOp:0
Alinear/zero_fraction/total_zero/zero_count_4/zero_fraction/Size:0
Flinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/Cast:0
Glinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/Merge:0
Glinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/Merge:1
Hlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/Switch:0
Hlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/Switch:1
Tlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/Cast:0
Ulinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/Const:0
_linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/NotEqual/Switch:1
Xlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/NotEqual:0
]linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/nonzero_count:0
Ulinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/zeros:0
Vlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/Cast:0
Wlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/Const:0
alinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
Zlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/NotEqual:0
_linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/nonzero_count:0
Wlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/switch_f:0
Jlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/switch_t:0
Tlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/counts_to_fraction/Cast:0
Vlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/counts_to_fraction/Cast_1:0
Slinear/zero_fraction/total_zero/zero_count_4/zero_fraction/counts_to_fraction/sub:0
Wlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/counts_to_fraction/truediv:0
Elinear/zero_fraction/total_zero/zero_count_4/zero_fraction/fraction:0�
-linear/linear_model/origin.y/weights/part_0:0Rlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOp/Switch:0i
(linear/zero_fraction/total_size/Size_4:0=linear/zero_fraction/total_zero/zero_count_4/ToFloat/Switch:0p
6linear/zero_fraction/total_zero/zero_count_4/pred_id:06linear/zero_fraction/total_zero/zero_count_4/pred_id:02�

�

Ilinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/cond_textIlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/pred_id:0Jlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/switch_t:0 *�	
Klinear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOp:0
Flinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/Cast:0
Tlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/Cast:0
Ulinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/Const:0
_linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/NotEqual/Switch:1
Xlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/NotEqual:0
]linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/nonzero_count:0
Ulinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/switch_t:0�
Klinear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOp:0_linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero/NotEqual/Switch:1�
Ilinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/pred_id:0Ilinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/pred_id:02�

�

Klinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/cond_text_1Ilinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/pred_id:0Jlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/switch_f:0*�
Klinear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOp:0
Vlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/Cast:0
Wlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/Const:0
alinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
Zlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/NotEqual:0
_linear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/nonzero_count:0
Wlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/switch_f:0�
Klinear/zero_fraction/total_zero/zero_count_4/zero_fraction/ReadVariableOp:0alinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0�
Ilinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/pred_id:0Ilinear/zero_fraction/total_zero/zero_count_4/zero_fraction/cond/pred_id:0
�
6linear/zero_fraction/total_zero/zero_count_5/cond_text6linear/zero_fraction/total_zero/zero_count_5/pred_id:07linear/zero_fraction/total_zero/zero_count_5/switch_t:0 *�
4linear/zero_fraction/total_zero/zero_count_5/Const:0
6linear/zero_fraction/total_zero/zero_count_5/pred_id:0
7linear/zero_fraction/total_zero/zero_count_5/switch_t:0p
6linear/zero_fraction/total_zero/zero_count_5/pred_id:06linear/zero_fraction/total_zero/zero_count_5/pred_id:0
�0
8linear/zero_fraction/total_zero/zero_count_5/cond_text_16linear/zero_fraction/total_zero/zero_count_5/pred_id:07linear/zero_fraction/total_zero/zero_count_5/switch_f:0*�
-linear/linear_model/origin.z/weights/part_0:0
(linear/zero_fraction/total_size/Size_5:0
=linear/zero_fraction/total_zero/zero_count_5/ToFloat/Switch:0
6linear/zero_fraction/total_zero/zero_count_5/ToFloat:0
2linear/zero_fraction/total_zero/zero_count_5/mul:0
6linear/zero_fraction/total_zero/zero_count_5/pred_id:0
7linear/zero_fraction/total_zero/zero_count_5/switch_f:0
Hlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/LessEqual/y:0
Flinear/zero_fraction/total_zero/zero_count_5/zero_fraction/LessEqual:0
Rlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOp/Switch:0
Klinear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOp:0
Alinear/zero_fraction/total_zero/zero_count_5/zero_fraction/Size:0
Flinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/Cast:0
Glinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/Merge:0
Glinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/Merge:1
Hlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/Switch:0
Hlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/Switch:1
Tlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/Cast:0
Ulinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/Const:0
_linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/NotEqual/Switch:1
Xlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/NotEqual:0
]linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/nonzero_count:0
Ulinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/zeros:0
Vlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/Cast:0
Wlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/Const:0
alinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
Zlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/NotEqual:0
_linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/nonzero_count:0
Wlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/switch_f:0
Jlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/switch_t:0
Tlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/counts_to_fraction/Cast:0
Vlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/counts_to_fraction/Cast_1:0
Slinear/zero_fraction/total_zero/zero_count_5/zero_fraction/counts_to_fraction/sub:0
Wlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/counts_to_fraction/truediv:0
Elinear/zero_fraction/total_zero/zero_count_5/zero_fraction/fraction:0p
6linear/zero_fraction/total_zero/zero_count_5/pred_id:06linear/zero_fraction/total_zero/zero_count_5/pred_id:0i
(linear/zero_fraction/total_size/Size_5:0=linear/zero_fraction/total_zero/zero_count_5/ToFloat/Switch:0�
-linear/linear_model/origin.z/weights/part_0:0Rlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOp/Switch:02�

�

Ilinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/cond_textIlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/pred_id:0Jlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/switch_t:0 *�	
Klinear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOp:0
Flinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/Cast:0
Tlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/Cast:0
Ulinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/Const:0
_linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/NotEqual/Switch:1
Xlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/NotEqual:0
]linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/nonzero_count:0
Ulinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/switch_t:0�
Klinear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOp:0_linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero/NotEqual/Switch:1�
Ilinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/pred_id:0Ilinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/pred_id:02�

�

Klinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/cond_text_1Ilinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/pred_id:0Jlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/switch_f:0*�
Klinear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOp:0
Vlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/Cast:0
Wlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/Const:0
alinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0
Zlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/NotEqual:0
_linear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/nonzero_count:0
Wlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/zeros:0
Ilinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/pred_id:0
Jlinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/switch_f:0�
Klinear/zero_fraction/total_zero/zero_count_5/zero_fraction/ReadVariableOp:0alinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/count_nonzero_1/NotEqual/Switch:0�
Ilinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/pred_id:0Ilinear/zero_fraction/total_zero/zero_count_5/zero_fraction/cond/pred_id:0"%
saved_model_main_op


group_deps*�
classification�
3
inputs)
input_example_tensor:0���������H
scores>
'linear/head/predictions/probabilities:0���������4
classes)
linear/head/Tile:0���������tensorflow/serving/classify*�

regression�
3
inputs)
input_example_tensor:0���������D
outputs9
"linear/head/predictions/logistic:0���������tensorflow/serving/regress*�
serving_default�
3
inputs)
input_example_tensor:0���������H
scores>
'linear/head/predictions/probabilities:0���������4
classes)
linear/head/Tile:0���������tensorflow/serving/classify*�
predict�
5
examples)
input_example_tensor:0���������H
	class_ids;
$linear/head/predictions/ExpandDims:0	���������G
classes<
%linear/head/predictions/str_classes:0���������F
all_class_ids5
linear/head/predictions/Tile:0���������F
all_classes7
 linear/head/predictions/Tile_1:0���������E
logistic9
"linear/head/predictions/logistic:0���������O
probabilities>
'linear/head/predictions/probabilities:0���������]
logitsS
<linear/linear_model/linear_model/linear_model/weighted_sum:0���������tensorflow/serving/predict