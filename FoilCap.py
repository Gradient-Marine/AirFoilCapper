import bpy
import numpy as np
import mathutils

# def curveFit(vertices, resolution):
#     x = vertices[:,0]
#     y = vertices[:,1]
#     z = vertices[:,2]
#     f = scipy.interpolate.interp2d(x, y, z, kind='cubic')


def main():

    topFoil = bpy.data.objects["topFoil"].data
    bottomFoil = bpy.data.objects["bottomFoil"].data
    centerLine = bpy.data.objects["centerLine"].data
    trailingEdge = bpy.data.objects["trailingEdge"].data
    outputCollection = bpy.data.collections["outputs"]

    # Create a new mesh object
    mesh = bpy.data.meshes.new(name="RandomMesh")
    obj = bpy.data.objects.new("RandomObject", mesh)
    outputCollection.objects.link(obj)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)

    # Generate vertices
    curves = [
        topFoil,
        bottomFoil,
        trailingEdge
    ]

    # parameters:
    interpRes = 18  # interpolation along foil curves (will replace with just the csv point data later)
    verticalIters = 8 # number of loops going upwards
    verticalHeight = 1 # distance to go upwards

    # create mesh for top, bottom, and trailing edge
    interpPoints = []
    bezier_points = 0
    for curve in curves:
        bezier = curve.splines.active.bezier_points
        for pts in range(1, len(bezier)):
            firstPt = bezier[pts-1]
            secondPt = bezier[pts]
            # leaving off last vertex of interp, to loop together curves
            interpPoints = interpPoints + mathutils.geometry.interpolate_bezier(firstPt.co, firstPt.handle_right, secondPt.handle_left, secondPt.co, interpRes+1)[0:-1]
            bezier_points += 1
    outlinePoints = np.array(interpPoints)

    # create ceneterline mesh
    centerBez = centerLine.splines.active.bezier_points
    centerPts = []
    for pts in range(1, len(centerBez)):
        firstPt = centerBez[pts-1]
        secondPt = centerBez[pts]
        if (pts < len(centerBez)-1): # need both endpoints since centerline isn't a loop like other curves
            centerPts = centerPts + mathutils.geometry.interpolate_bezier(firstPt.co, firstPt.handle_right, secondPt.handle_left, secondPt.co, interpRes+1)[0:-1]
        else:
            centerPts = centerPts + mathutils.geometry.interpolate_bezier(firstPt.co, firstPt.handle_right, secondPt.handle_left, secondPt.co, interpRes+1)
    
    # trim centerline (this is a terrible way to set the size of the top edge, dependent on resolution)
    centerPts = np.array(centerPts[10:-10])

    # find nearest points for base outline and centerline, the mapping betwix the two
    outlineToCenter = [] # index of centerline that is nearest to outline index
    for i in range(outlinePoints.shape[0]):
        pt = outlinePoints[i]
        dists = ((centerPts - pt)**2).sum(1)
        outlineToCenter.append(np.argmin(dists))
    outlineToCenter = np.array(outlineToCenter)

    # make this part a loop and move subsequent outlines closer to centerline
    # create outline and lift
    numVertices = bezier_points*interpRes
    vertexArray = outlinePoints
    squishDirections = ((centerPts[outlineToCenter] - outlinePoints)*np.array([1,1,0]))
    edgeArray = np.array([[i-1, i] for i in range(1,numVertices)] + [[(numVertices - 1), 0]])
    edgeLoop = edgeArray.copy()
    faces = []

    for v in range(verticalIters-1):
        # MAIN CODE: lifts the outline and squishes it in towards the center line and raises it up
        factor = (v+1)/(verticalIters-1) # travel parameter
        newVertices = outlinePoints + np.array([0,0,factor*verticalHeight]) 
        newVertices = newVertices + np.power(factor,2.0)*squishDirections
        vertexArray = np.append(vertexArray, newVertices, 0)
        # top loop
        edgeArray = np.append(edgeArray, edgeLoop + [numVertices*(v+1),numVertices*(v+1)],0)
        # vertical edges
        edgeArray = np.append(edgeArray, [[i + numVertices*v, i+(numVertices*(v+1))] for i in range(numVertices)],0)
        # square faces
        faces = faces + ([[i-1+numVertices*v, i+numVertices*v, i+numVertices*(v+1), i-1+numVertices*(v+1)] for i in range(1, numVertices)] + [[numVertices*(v+1)-1, numVertices*v, numVertices*(v+1), numVertices*(v+2)-1]])

    numCenterPts = len(centerPts)
    # add edges and faces that connect to the centerline of the foil
    vertexArray = np.append(vertexArray, centerPts + np.array([0,0,vertexArray[-1,2]]),0)
    # edges along spine
    edgeArray = np.append(edgeArray, np.array([[i-1, i] for i in range(1,numCenterPts)]) + [numVertices*verticalIters, numVertices*verticalIters],0)
    # edges from topmost outline to nearest spine vertex
    edgeArray = np.append(edgeArray, [[i+(verticalIters-1)*numVertices, outlineToCenter[i]+(verticalIters)*numVertices] for i in range(numVertices)], 0)

    # faces from topmost outline to spine, need for loop and if statement
    vertexOffset = (verticalIters-1)*numVertices
    for i in range(1, numVertices):
        if outlineToCenter[i] == outlineToCenter[i-1]:
            # face has only 3 vertices
            faces = faces + [[i - 1 + vertexOffset, i + vertexOffset, outlineToCenter[i] + vertexOffset + numVertices]]
        else:
            # face has 4 vertices and one edge along spine
            faces = faces + [[i - 1 + vertexOffset, i + vertexOffset, outlineToCenter[i] + vertexOffset + numVertices, outlineToCenter[i-1] + vertexOffset + numVertices]]
        pass

    # add last missing face from wrapping around profile, from 0 to numVertices-1
    if outlineToCenter[0] == outlineToCenter[numVertices-1]:
        faces = faces + [[vertexOffset, numVertices - 1 + vertexOffset, outlineToCenter[numVertices-1] + vertexOffset + numVertices]]
    else:
        faces = faces + [[vertexOffset, numVertices - 1 + vertexOffset, outlineToCenter[numVertices-1] + vertexOffset + numVertices, outlineToCenter[0] + vertexOffset + numVertices]]
    pass

    afPoints = getPoints()
    numAFPoints = afPoints.shape[0]
    halfPointsU = int(np.ceil(numAFPoints/2.0))
    halfPointsL = int(np.floor(numAFPoints/2.0))
    meanLine = (afPoints[:halfPointsU] + afPoints[halfPointsL:])/2.0
    # vertexArray = np.append(afPoints, meanLine, 0)
    # vertexArray = afPoints

    # mesh.from_pydata(vertexArray, [], [])
    mesh.from_pydata(vertexArray, edgeArray, faces)

    # Update and display the mesh
    mesh.update()

    print("yas queen")
    pass

def squishFactorH(i, iMax):
    return 0.1*(i/iMax)


def squishFactorV(i, iMax):
    return 0.1*(i/iMax)
    # return 0.1*np.sin((i/iMax)*(np.pi/2.0))

def getPoints():
    return np.array([
            [1.000084,  0.001257, 0.0],
            [0.999405,  0.001399, 0.0],
            [0.997370,  0.001822, 0.0],
            [0.993984,  0.002524, 0.0],
            [0.989256,  0.003499, 0.0],
            [0.983197,  0.004738, 0.0],
            [0.975825,  0.006231, 0.0],
            [0.967156,  0.007967, 0.0],
            [0.957215,  0.009932, 0.0],
            [0.946027,  0.012110, 0.0],
            [0.933621,  0.014485, 0.0],
            [0.920029,  0.017038, 0.0],
            [0.905287,  0.019752, 0.0],
            [0.889434,  0.022606, 0.0],
            [0.872512,  0.025580, 0.0],
            [0.854565,  0.028653, 0.0],
            [0.835642,  0.031805, 0.0],
            [0.815792,  0.035015, 0.0],
            [0.795069,  0.038260, 0.0],
            [0.773528,  0.041520, 0.0],
            [0.751228,  0.044774, 0.0],
            [0.728228,  0.048000, 0.0],
            [0.704592,  0.051177, 0.0],
            [0.680382,  0.054285, 0.0],
            [0.655665,  0.057302, 0.0],
            [0.630509,  0.060209, 0.0],
            [0.604982,  0.062985, 0.0],
            [0.579155,  0.065609, 0.0],
            [0.553099,  0.068063, 0.0],
            [0.526886,  0.070326, 0.0],
            [0.500588,  0.072381, 0.0],
            [0.474279,  0.074210, 0.0],
            [0.448032,  0.075795, 0.0],
            [0.421921,  0.077122, 0.0],
            [0.395987,  0.078173, 0.0],
            [0.370157,  0.078879, 0.0],
            [0.344680,  0.079198, 0.0],
            [0.319630,  0.079125, 0.0],
            [0.295081,  0.078659, 0.0],
            [0.271106,  0.077802, 0.0],
            [0.247774,  0.076558, 0.0],
            [0.225154,  0.074936, 0.0],
            [0.203313,  0.072947, 0.0],
            [0.182315,  0.070606, 0.0],
            [0.162221,  0.067930, 0.0],
            [0.143088,  0.064941, 0.0],
            [0.124973,  0.061660, 0.0],
            [0.107927,  0.058112, 0.0],
            [0.091996,  0.054325, 0.0],
            [0.077226,  0.050327, 0.0],
            [0.063657,  0.046145, 0.0],
            [0.051324,  0.041808, 0.0],
            [0.040261,  0.037346, 0.0],
            [0.030495,  0.032785, 0.0],
            [0.022051,  0.028152, 0.0],
            [0.014950,  0.023471, 0.0],
            [0.009206,  0.018764, 0.0],
            [0.004833,  0.014049, 0.0],
            [0.001838,  0.009343, 0.0],
            [0.000227,  0.004657, 0.0],
            [0.000000,  0.000000, 0.0],
            [0.001143, -0.004520, 0.0],
            [0.003640, -0.008797, 0.0],
            [0.007479, -0.012828, 0.0],
            [0.012647, -0.016609, 0.0],
            [0.019125, -0.020136, 0.0],
            [0.026892, -0.023408, 0.0],
            [0.035924, -0.026419, 0.0],
            [0.046194, -0.029168, 0.0],
            [0.057669, -0.031651, 0.0],
            [0.070318, -0.033869, 0.0],
            [0.084103, -0.035820, 0.0],
            [0.098987, -0.037507, 0.0],
            [0.114928, -0.038931, 0.0],
            [0.131882, -0.040098, 0.0],
            [0.149805, -0.041013, 0.0],
            [0.168649, -0.041686, 0.0],
            [0.188365, -0.042126, 0.0],
            [0.208902, -0.042346, 0.0],
            [0.230207, -0.042360, 0.0],
            [0.252226, -0.042183, 0.0],
            [0.274904, -0.041834, 0.0],
            [0.298182, -0.041331, 0.0],
            [0.322002, -0.040693, 0.0],
            [0.346303, -0.039941, 0.0],
            [0.371024, -0.039095, 0.0],
            [0.396102, -0.038177, 0.0],
            [0.421644, -0.037174, 0.0],
            [0.447439, -0.036049, 0.0],
            [0.473385, -0.034816, 0.0],
            [0.499412, -0.033493, 0.0],
            [0.525450, -0.032095, 0.0],
            [0.551429, -0.030639, 0.0],
            [0.577279, -0.029138, 0.0],
            [0.602929, -0.027607, 0.0],
            [0.628310, -0.026057, 0.0],
            [0.653352, -0.024500, 0.0],
            [0.677986, -0.022945, 0.0],
            [0.702145, -0.021403, 0.0],
            [0.725762, -0.019880, 0.0],
            [0.748772, -0.018385, 0.0],
            [0.771111, -0.016922, 0.0],
            [0.792716, -0.015499, 0.0],
            [0.813528, -0.014119, 0.0],
            [0.833489, -0.012788, 0.0],
            [0.852541, -0.011510, 0.0],
            [0.870633, -0.010289, 0.0],
            [0.887712, -0.009129, 0.0],
            [0.903730, -0.008033, 0.0],
            [0.918642, -0.007006, 0.0],
            [0.932405, -0.006052, 0.0],
            [0.944979, -0.005174, 0.0],
            [0.956330, -0.004376, 0.0],
            [0.966424, -0.003662, 0.0],
            [0.975232, -0.003035, 0.0],
            [0.982729, -0.002498, 0.0],
            [0.988892, -0.002055, 0.0],
            [0.993705, -0.001708, 0.0],
            [0.997152, -0.001458, 0.0],
            [0.999225, -0.001307, 0.0],
            [0.999916, -0.001257, 0.0]
    ])