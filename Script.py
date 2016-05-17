import mathutils
import bpy
import numpy
import scipy
import scipy.spatial
import math



context = bpy.context;
e = (context.active_object.dimensions.x + context.active_object.dimensions.y + context.active_object.dimensions.z) / 100
vertices = context.active_object.data.vertices
points = list()

for v in vertices:
    points.append(v.co)
    

normals = list()

for n in context.active_object['vertex_normal_list']:
    normal = mathutils.Vector((n[0], n[1], n[2]))
    normals.append(normal)

pointsPlusN = list()
pointsPlus2N = list()
i = 0
while i < len(points):
    point = points[i] + e * normals[i]
    pointsPlusN.append(point) #normals maal e toevoegen
    distances = scipy.spatial.distance.cdist(points,[point]) 
    for distance in distances:
        if distance < e and distance > 0:
            pointsPlusN = list()
            pointsPlus2N = list()
            i = -1
            e = e / 2
    if i > -1:
        point = points[i] - e * normals[i] #normals maal e aftrekken
        pointsPlus2N.append(point) 
        distances = scipy.spatial.distance.cdist(points,[point]) 
        for distance in distances:
            if distance < e and distance > 0:
                pointsPlusN = list()
                pointsPlus2N = list()
                i = -1
                e = e / 2
    i = i+1

dvector = list()

for d in points:
    dvector.append(0)

for d in pointsPlusN:
    dvector.append(e)

for d in pointsPlus2N:
    dvector.append(-e)

cvalues = list()

i = 0
while i < len(points): 
    cvalues.append(points[i])
    cvalues.append(pointsPlusN[i])
    cvalues.append(pointsPlus2N[i])
    i = i + 1

constraintMatrix = list(list())
i = 0
while i < len(cvalues): 
    j = 0
    constraintRow = list()
    while j < len(cvalues):
        constraintRow.append(math.pow(scipy.spatial.distance.euclidean(cvalues[i], cvalues[j]),3))
        j = j+1
    constraintMatrix.append(constraintRow)
    i = i+1

realConstraintMatrix = numpy.matrix(constraintMatrix)
weightValues = numpy.linalg.lstsq(realConstraintMatrix,dvector)


# calculate bounding area 
minX = 100000
minY = 100000
minZ = 100000
maxX = -100000
maxY = -100000
maxZ = -100000

for p in points:
    if (p[x] > maxX):
        maxX = p[x]
    if (p[x] < minX):
        minX = p[x]
    if (p[y] > maxY):
        maxY = p[y]
    if (p[y] < minY):
        minY = p[y]
    if (p[z] > maxZ):
        maxZ = p[z]
    if (p[z] < minZ):
        minZ = p[z]

minX = minX * 1.01
minY = minY * 1.01
minZ = minZ * 1.01
maxX = maxX * 1.01
maxY = maxY * 1.01
maxZ = maxZ * 1.01

#divide de bounding area in cubes
x = 0
y = 0
z = 0