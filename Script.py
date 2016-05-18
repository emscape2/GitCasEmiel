import mathutils
import bpy
import numpy
import scipy
import scipy.spatial
import math

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

    def search(self, point):
        ind = self.indexConvert(self.getIndex(point))

        if ind in self.data:
            return self.data[ind]

        return list()

    # creates a spatial index from point list
    def create(self, points):
        for p in points:
            ind = self.indexConvert(self.getIndex(p))
            if not ind in self.data:
                self.data[ind] = list()
            self.data[ind].append(p)

    # returns the neighbouring cubes of index
    def neighbours(self, index):
        neighbours = list()

        # loop through neighbours
        for x in range(max(0, index[0] - 1), min(self.resolution, index[0] + 1) + 1):
            for y in range(max(0, index[1] - 1), min(self.resolution, index[1] + 1) + 1):
                for z in range(max(0, index[2] - 1), min(self.resolution, index[2] + 1) + 1):
                    if not (x == index[0] and y == index[1] and z == index[2]):
                        neighbours.append((x,y,z))

        return neighbours

    # returns the neighbouring points of a point, EXCLUDING the point itself
    def neighbouringPoints(self, point):
        pList = list()

        # retrieve points local cube
        for p in self.search(point):
            if not (p == point):
                pList.append(p)

        # retrieve points from neighbouring cubes
        for ns in self.neighbours(self.getIndex(point)):
            ind = self.indexConvert(ns)
            if ind in self.data:
                for p in self.data[ind]:
                    pList.append(p)

        return pList

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


## zie boven de defs voor de grote TODO lijst bounding box die half af is staat helemaal onderaan in bounding area
context = bpy.context;
e = (context.active_object.dimensions.x + context.active_object.dimensions.y + context.active_object.dimensions.z) / 100
vertices = context.active_object.data.vertices
points = list()

for v in vertices:
    points.append(v.co)

# calculate object bounding area
boundingArea = getBoundingArea(points)

# create spatial index
spatialIndex = Spatial(100)
spatialIndex.create(points)

normals = list()

# add normals to normal list
for n in context.active_object['vertex_normal_list']:
    normal = mathutils.Vector((n[0], n[1], n[2]))
    normals.append(normal)

pointsPlusN = list()
pointsPlus2N = list()

i = 0
while i < len(points): 
    point = points[i] + e * normals[i]
    pointsPlusN.append(point) #normals maal e toevoegen
    neighbours = spatialIndex.neighbouringPoints(point)

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
        neighbours = spatialIndex.neighbouringPoints(point)

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
    cpValues.append(pointsPlusN[i])
    cpValues.append(pointsPlus2N[i])
    i = i + 1



def MlsFunction(point, controlindices, smoothness = 4, degree = 1):
    controlPoints = list()
    dValues = list()
    basePoly = numpy.matrix( BasePolynomial(point, 0, True, degree))
    basePoly.transpose();
    #pick wich points and values to use for interpolation
    for i in controlindices:
        controlPoints.append(cpValues[i])
        dValues.append(dvector[i])
    
    #Calculates the B matrix
    BMatrix = buildIdealBMatrix(point, controlPoints, degree)
    
    dValues = numpy.matrix(dValues).transpose()
    
    #Gets the inverse of the A matrix, hier kijken of de 3 dimensionale matrix wel goed geaccepteerd wordt 
    AMatrixInverse = buildIdealAMatrix(point, controlPoints, degree)
    for matrix in AMatrixInverse:
        print(matrix)
        matrix.I
        #berekent voor iedere matrix in A de juiste waarden en telt deze op
        matrix = matrix * BMatrix * dValues * basePoly
        matrix = numpy.sum(matrix)

    result = (AMatrixInverse )
    return numpy.sum(result)
    

    #Leest de waarde van polynomial bij punt af, position maakt niet uit bij needEntirePolygon
def BasePolynomial(point, position, needEntirePolygon = False, degree = 1):#nu voor degree 0 en 1, TODO: 2 implementeren en entire polygon achterhalen.
    if needEntirePolygon == False:
        if degree == 0 or position == 0:
            return 1
        if degree == 1:
            switcher = {
                1: point[x],
                2: point[y],
                3: point[z]
                }
            return switcher.get(position, 1)
    if degree == 0:
        return [1]
    if degree == 1:
        vector = list()
        vector.append(1)
        vector.append(point[x])
        vector.append(point[y])
        vector.append(point[z])
        return vector


    #returns the amount of elements present in the base polynomial
def BasePolynomialLength(degree = 1):
    switcher = {
        1: 4,
        2: 10
        }
    return switcher.get(degree, 1)

#creates the optimal B matrix
def buildIdealBMatrix(point, controlpoints = list(), degree = 1):
    idealBMatrix = list()
    i = 0
    while i < BasePolynomialLength(degree):
        matrixRow = list()
        for cp in controlpoints:
            matrixRow.append(Wendland(scipy.spatial.distance.euclidean(point,cp)) * BasePolynomial(cp,i))
        idealBMatrix.append(matrixRow)
        i = i +1
    return numpy.matrix(idealBMatrix)

#Creates the optimal A matrix 
def buildIdealAMatrix(point, controlpoints = list(), degree = 1):
    idealAMatrix = list()
    
    for cp in controlpoints:
        i = 0
        matrixRow = list()
        while i < BasePolynomialLength(degree):
            j = 0
            innerVector = list()
            while j < BasePolynomialLength(degree):
                innerVector.append(Wendland(scipy.spatial.distance.euclidean(point,cp)) * BasePolynomial(cp,i) * BasePolynomial(cp,j))
                j = j + 1
            matrixRow.append(innerVector)
            i = i + 1
        idealAMatrix.append(numpy.matrix(matrixRow))
    return idealAMatrix
     
#Weendland function, smoothness still needs defining
def Wendland(inputValue, smoothness = 1):
    if inputValue > smoothness:
        return 0
    else:
        return (math.pow(1-(inputValue/smoothness),4) * (4*inputValue/smoothness))
