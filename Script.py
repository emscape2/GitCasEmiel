import mathutils
import bpy
import numpy
import scipy
import scipy.spatial
import math
import mcubes

class Spatial:   

    # initialize spatial index
    def __init__(self, res):
        self.resolution = res
        self.data = {}

    # converts the three-tuple spatial index to a string that can be used as a key
    def indexConvert(self, indexTuple):
        return str(indexTuple[0] + (indexTuple[1] * 128) + (indexTuple[2] * (128**2)))

    # returns the spatial index of a given point
    def getIndex(self, point):
        indX = math.floor((point.x - boundingArea['minX']) / context.active_object.dimensions.x * self.resolution)
        indY = math.floor((point.y - boundingArea['minY']) / context.active_object.dimensions.y * self.resolution)
        indZ = math.floor((point.z - boundingArea['minZ']) / context.active_object.dimensions.z * self.resolution)
        return (indX, indY, indZ)

    # searches the spatial index for a certain point, returning all points that lie within the same cube, EXCLUDING the point itself
    def search(self, point, provideIndices = False):
        ind = self.indexConvert(self.getIndex(point))


        if ind in self.data:
            if provideIndices:
                return self.data[ind]
            else: 
                rList = list()
                for item in self.data[ind]:
                    rList.append(item)
                    return rList

        return list()

    # creates a spatial index from point list
    def create(self, points):
        for i in range(len(points)):
            self.add(i, points[i])

    # adds a point to the spatial index
    def add(self, index, point):
        ind = self.indexConvert(self.getIndex(point))
        if not ind in self.data:
            self.data[ind] = list()
        self.data[ind].append((index, point))
 
    # returns the neighbouring cubes of index
    def neighbours(self, index, outerRange):
        neighbours = list()

        # loop through neighbours
        for x in range(max(0, index[0] - outerRange), min(self.resolution, index[0] + outerRange) + 1):
            for y in range(max(0, index[1] - outerRange), min(self.resolution, index[1] + outerRange) + 1):
                for z in range(max(0, index[2] - outerRange), min(self.resolution, index[2] + outerRange) + 1):
                    if not (x == index[0] and y == index[1] and z == index[2]):
                        neighbours.append((x,y,z))

        return neighbours

    # returns the neighbouring points of a point, EXCLUDING the point itself
    def neighbouringPoints(self, point, provideIndices = False, minResults = 0):
        pList = list()

        # retrieve points local cube
        for p in self.search(point):
            if not (p == point):
                pList.append(p)
        
        range = 0
        while len(pList) < minResults and range < 3:
            pList = list()
            range = range + 1
            # retrieve points from neighbouring cubes
            for ns in self.neighbours(self.getIndex(point), range):
                ind = self.indexConvert(ns)
                if ind in self.data:
                    for p in self.data[ind]:
                        pList.append(p)
        
        if provideIndices:
            return pList

        rList = list()
        for p in pList:
            rList.append(p[1])

        return rList

# calculates the bounding area of the object
def getBoundingArea(inputData):
    minX = 100000
    minY = 100000
    minZ = 100000
    maxX = -100000
    maxY = -100000
    maxZ = -100000
    
    for vertex in inputData:
        minX = min(vertex.x, minX)
        maxX = max(vertex.x, maxX)

        minY = min(vertex.y, minY)
        maxY = max(vertex.y, maxY)

        minZ = min(vertex.z, minZ)
        maxZ = max(vertex.z, maxZ)
    
    return {"minX": minX * 1.01, "minY": minY * 1.01, "minZ": minZ * 1.01, "maxX": maxX * 1.01, "maxY": maxY * 1.01, "maxZ": maxZ * 1.01}

# implicit function to be used by the marching cubes algorithm
def implicit(x, y, z):
    vector = mathutils.Vector((x, y, z))
    degree = Degree
    neighbours = spatialIndex.neighbouringPoints(vector, True, 9)

    if len(neighbours) == 0:
        return 10000

    indices = list()
    for n in neighbours:
        indices.append(n[0])


    a = MlsFunction(vector, indices, degree)
    b = BasePolynomial(vector,0,True,degree)

    return numpy.matrix.sum( a.transpose() * b )

#calculates the MlsFunction
def MlsFunction(point, controlindices, degree = 1):
    
    wendlandValues = numpy.zeros(len(controlindices))
    i = 0
    for index in controlindices:
        wendlandValue =  Wendland(scipy.spatial.distance.euclidean(point,cpValues[index]))
        wendlandValues[i] = wendlandValue
        i = i + 1

    part1 = buildIdealPart1(point,wendlandValues, controlindices,degree)
    part2 = buildIdealPart2(point,wendlandValues, controlindices,degree)
    aVector = numpy.linalg.lstsq(part1,part2)
    return numpy.matrix( aVector[0])

    
    

    #Leest de waarde van polynomial bij punt af, position maakt niet uit bij needEntirePolygon
def BasePolynomial(point, position, needEntirePolygon = False, degree = 1):#nu voor degree 0 en 1, TODO: 2 implementeren en entire polygon achterhalen.
    if needEntirePolygon == False:
        if degree == 0 or position == 0:
            return 1
        if degree == 1:
            switcher = {
                1: point.x,
                2: point.y,
                3: point.z
                }
            return switcher.get(position, 1)
        if degree == 2:
            switcher = {
                1: point.x,
                2: point.y,
                3: point.z,
                4: point.x * point.y,
                5: point.x * point.z,
                6: point.z * point.y,
                7: point.x * point.x,
                8: point.y * point.y,
                9: point.z * point.z
                }
            return switcher.get(position, 1)
    if degree == 0:
        return [1]
    if degree == 1:
        vector = list()
        vector.append(1)
        vector.append(point.x)
        vector.append(point.y)
        vector.append(point.z)
        return vector
    if degree == 2:
        vector = list()
        vector.append(1)
        vector.append(point.x)
        vector.append(point.y)
        vector.append(point.z)
        vector.append(point.y * point.x)
        vector.append(point.x * point.z)
        vector.append(point.z * point.y)
        vector.append(point.x * point.x)
        vector.append(point.y * point.y)
        vector.append(point.z * point.z)
        return vector


    #returns the amount of elements present in the base polynomial
def BasePolynomialLength(degree = 1):
    switcher = {
        0: 1,
        1: 4,
        2: 10
        }
    return switcher.get(degree, 1)

def findAllNeighbours(x,y,z):
    min = mathutils.Vector((boundingArea['minX'],boundingArea['minY'],boundingArea['minZ']))
    Range = mathutils.Vector((boundingArea['maxX'], boundingArea['maxY'], boundingArea['maxZ'])) - min
    steps = mathutils.Vector((1/x,1/y,1/z))
    for X in range(0,x):
        for Y in range(0,y):
            for Z in range(0,z):
                pos = steps * mathutils.Vector((X,Y,Z))
                coords = pos * Range + min
                neighbours = spatialIndex.neighbouringPoints(coords,False,9)


def findAllValues(x,y,z):
    min = mathutils.Vector((boundingArea['minX'],boundingArea['minY'],boundingArea['minZ']))
    Range = mathutils.Vector((boundingArea['maxX'], boundingArea['maxY'], boundingArea['maxZ'])) - min
    steps = mathutils.Vector((1/x,1/y,1/z))
    for X in range(0,x):
        for Y in range(0,y):
            for Z in range(0,z):
                pos = steps * mathutils.Vector((X,Y,Z))
                coords = pos * Range + min
                implicit(coords.x,coords.y,coords.z)

#Creates the optimal A matrix
def buildIdealAMatrix(point, controlindices = list(), degree = 1):
    idealAMatrix = list()
    for index in controlindices:
        i = 0 
        matrixRow = list()
        while i < BasePolynomialLength(degree):
            matrixRow.append(math.sqrt(Wendland(scipy.spatial.distance.euclidean(point,cpValues[index]))) * basePolynomialList[index][i])#basepolynomial cp, i cp = controlpoint als niet werkt
            i = i+1
        idealAMatrix.append(matrixRow)
    return numpy.matrix(idealAMatrix)

#Creates matrix A * A^T, much more efficient
def buildIdealPart1(point,wendlandValues , controlindices = list(), degree = 1):
    idealMatrix = numpy.zeros((len(wendlandValues),BasePolynomialLength(degree)))
    for index in range(0,len(wendlandValues)):
        i = 0 
        wendlandValue = wendlandValues[index]
        while i < BasePolynomialLength(degree):
            bp =  basePolynomialList[controlindices[index]][i] * sum(basePolynomialList[controlindices[index]])
            idealMatrix[index,i] = wendlandValue * bp#
            i = i+1
    return numpy.matrix(idealMatrix)

#Creates matrix A^T * r, much more efficient
def buildIdealPart2(point,wendlandValues, controlindices = list(), degree = 1):
    idealMatrix = numpy.zeros(len(wendlandValues))
    for index in range(0,len(wendlandValues)):
        wendlandValue = wendlandValues[index]
        bp =  sum(basePolynomialList[controlindices[index]])
        idealMatrix[index] = wendlandValue * bp * dvector[controlindices[index]]
        
    return idealMatrix

#creates the optimal r Vector
def buildIdealrVector(point, controlIndices = list(), degree = 1):
    rVector = list()
    for i in controlIndices:
        rVector.append(math.sqrt(Wendland(scipy.spatial.distance.euclidean(point,cpValues[i]))) * dvector[i])
    return rVector

#creates the optimal a Vector
def buildIdealaVector(aMatrix, rVector):
    part1 = aMatrix.transpose() * numpy.matrix(rVector).transpose()
    part2 = aMatrix.transpose() * aMatrix
    #iets met singular verwerking doen
    a = numpy.linalg.lstsq(part2,part1)
    return numpy.matrix(a[0])



     
#Weendland function, smoothness still needs defining
def Wendland(inputValue, smoothness = 0.5):
    if inputValue > wendlandRadius:
        return 0
    else:
        return (((1-(inputValue/wendlandRadius)) ** 4) * (4*inputValue/wendlandRadius))


## zie boven de defs voor de grote TODO lijst bounding box die half af is staat helemaal onderaan in bounding area
context = bpy.context;
e = (context.active_object.dimensions.x + context.active_object.dimensions.y + context.active_object.dimensions.z) / 100
wendlandRadius = e * 6
vertices = context.active_object.data.vertices
points = list()
Degree = 1 
for v in vertices:
    points.append(v.co)

# calculate object bounding area
boundingArea = getBoundingArea(points)

# create spatial index
tmpSpatial = Spatial(16)
tmpSpatial.create(points)

normals = list()

# add normals to normal list
for n in context.active_object['vertex_normal_list']:
    normal = mathutils.Vector((n[0], n[1], n[2]))
    normals.append(normal)

pointsPlusN = list()
pointsPlus2N = list()

i = 0
while i < len(points): #Deze nog optimaliseren
    point = points[i] + e * normals[i]
    pointsPlusN.append(point) #normals maal e toevoegen
    neighbours = tmpSpatial.neighbouringPoints(point)

    if len(neighbours) > 0:
        distances = scipy.spatial.distance.cdist(neighbours, [point]) 
        for distance in distances:
            if distance < e and distance > 0:
                pointsPlusN = list()
                pointsPlus2N = list()
                i = -1
                e = e / 2
    if i > -1:
        point = points[i] - e * normals[i] #normals maal e aftrekken
        pointsPlus2N.append(point) 
        neighbours = tmpSpatial.neighbouringPoints(point)

        if len(neighbours) > 0:
            distances = scipy.spatial.distance.cdist(neighbours, [point]) 
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

cpValues = list()

i = 0
while i < len(points): 
    cpValues.append(points[i])
    i = i + 1

i = 0
while i < len(points): 
    cpValues.append(pointsPlusN[i])
    i = i + 1

i = 0
while i < len(points): 
    cpValues.append(pointsPlus2N[i])
    i = i + 1

spatialIndex = Spatial(16)
spatialIndex.create(cpValues)

def newBasePolynomialList(degree):
    Degree = degree
    basePolynomialList = list()
    for p in cpValues:
        basePolynomialList.append(BasePolynomial(p,0,True,Degree))
    return basePolynomialList

basePolynomialList = newBasePolynomialList(Degree)

def evaluate(dimX,dimY,dimZ):
    vertices, triangles = mcubes.marching_cubes_func(
            (boundingArea['minX'],boundingArea['minY'],boundingArea['minZ']), 
            (boundingArea['maxX'], boundingArea['maxY'], boundingArea['maxZ']),  # Bounds
            dimX, dimY, dimZ,                                                    # Number of samples in each dimension
            implicit,                                                            # Implicit function
            0)                                                                   # Isosurface value                                                                

    # Export the result to sphere2.dae
    mcubes.export_mesh(vertices, triangles, "/Users/Emscape/Documents/Blender Projects/result.dae", "MLS result")
    print("Done. Result saved in 'result.dae'.")




