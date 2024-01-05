import arcpy
import numpy
import queue
import math
import copy

#Parameters   you just need to adjust these parameters
path_inputdem = r"E:\Doctor\Paper\surrounding_dir\data\testdem.tif"      #the path of DEMs. when you text this code,you can use the DEMs we provide in our repository
path_result = r"E:\Doctor\Paper\surrounding_dir\data\D8.tif"                                    #the path for saveing flow direction result
numpy_dem = arcpy.RasterToNumPyArray(path_inputdem,nodata_to_value = -9999)
startacc = 117500                                                            #the threshold for extracting backbone flow path in the first round. the maximum acc value in flats boundary
tolerance = 500                                                               #50 or 500;when you test this code,maybe you can use an higher tolerance value such as 500. the compute time will be decrease.
#Parameters
round = int((startacc-tolerance)/tolerance)
#environments
arcpy.env.overwriteOutput = True
arcpy.env.snapRaster = path_inputdem
arcpy.env.extent = path_inputdem
temp_cellsizex = arcpy.GetRasterProperties_management(path_inputdem, "CELLSIZEX")
cellsizex = float(temp_cellsizex.getOutput(0))
cellsizey = cellsizex
inraster = arcpy.sa.Raster(path_inputdem, True)
lowerleft_point = arcpy.Point(inraster.extent.XMin, inraster.extent.YMin)
#objects
class cell:
    def __init__(self,row,col,ele):
        self.row = row
        self.col = col
        self.ele = ele
#functions
def D8(numpy_dem):
    row, col = numpy.shape(numpy_dem)                                       #record the shape of DEMs
    numpy_D8result = numpy.zeros((row,col)).astype("uint8")                 #create a numpy to record the flow direction result
    cellnodirinD8 = queue.Queue()                                           #create a queue to record those cells without unique downstream slope(in this queue all cells have downstream slope)
    cellinflat = queue.Queue()                                              #create a queue to record those cells without unique downstream slope(in this queue all cells have no downstream slope)
    for i in range(0,row):
        for j in range(0,col):
            #step1 all the raster cell on the edge of dem would flow out
            if i == 0 and j == 0:
                numpy_D8result[i][j] = 32
            elif i == 0 and j != col - 1:
                numpy_D8result[i][j] = 64
            elif i == 0 and j == col - 1:
                numpy_D8result[i][j] = 128
            elif j == col - 1 and i != row - 1:
                numpy_D8result[i][j] = 1
            elif j == col - 1 and i == row - 1:
                numpy_D8result[i][j] = 2
            elif j != 0 and i == row - 1:
                numpy_D8result[i][j] = 4
            elif j == 0 and i == row - 1:
                numpy_D8result[i][j] = 8
            elif j == 0:
                numpy_D8result[i][j] = 16
            else:
                #step 2 traditional D8 algorithm
                surrounding = [(numpy_dem[i][j]-numpy_dem[i-1][j-1])/1.414,(numpy_dem[i][j]-numpy_dem[i-1][j]),(numpy_dem[i][j]-numpy_dem[i-1][j+1])/1.414,(numpy_dem[i][j]-numpy_dem[i][j-1]),(numpy_dem[i][j]-numpy_dem[i][j+1]),(numpy_dem[i][j]-numpy_dem[i+1][j-1])/1.414,(numpy_dem[i][j]-numpy_dem[i+1][j]),(numpy_dem[i][j]-numpy_dem[i+1][j+1])/1.414]
                numofds,location,maxvalue,_ = findnumofdownstream(surrounding,1)        # 1 calculate the number of maximum in surrounding
                if numofds == 1 and maxvalue > 0:
                    if location == 0:
                        numpy_D8result[i][j] = 32
                    elif location == 1:
                        numpy_D8result[i][j] = 64
                    elif location == 2:
                        numpy_D8result[i][j] = 128
                    elif location == 3:
                        numpy_D8result[i][j] = 16
                    elif location == 4:
                        numpy_D8result[i][j] = 1
                    elif location == 5:
                        numpy_D8result[i][j] = 8
                    elif location == 6:
                        numpy_D8result[i][j] = 4
                    elif location == 7:
                        numpy_D8result[i][j] = 2
                elif maxvalue > 0:                                                       #raster cells with equal staus in traditional D8 algorithm
                    cellnodirinD8.put(cell(i,j,numpy_dem[i][j]))
                elif maxvalue == 0:                                                      #cells in flats
                    cellinflat.put(cell(i,j,numpy_dem[i,j]))
    return numpy_D8result,cellnodirinD8,cellinflat

def findnumofdownstream(array_surele,flag):
    row = numpy.shape(array_surele)[0]
    numofds = 0
    location = -1
    sur = [0, 0, 0, 0, 0, 0, 0, 0]                                                       #Record the position of the extreme value
    if flag == 1:                                                                        #maximum
        maxvalue = numpy.max(array_surele)
        for i in range(0,row):
            if array_surele[i] == maxvalue:
                numofds += 1
                location = i
                sur[i] = 1
        return numofds,location,maxvalue,sur
    if flag == 2:                                                                       #minimum
        minvalue = numpy.min(array_surele)
        for i in range(0,row):
            if array_surele[i] == minvalue:
                numofds += 1
                location = i
                sur[i] = 1
        return numofds, location, minvalue, sur


def improveD8(numpy_dem,result_D8,cellnodirinD8):                                       #calculate the flow direction value of those cells without unique downstream slope(in this queue all cells have downstream slope)
    while(not queue.Queue.empty(cellnodirinD8)):
        cell = cellnodirinD8.get()
        row = cell.row
        col = cell.col
        dir = dirwithsurele(cell,numpy_dem,1,None)
        result_D8[row][col] = dir
    return result_D8

def improvePF(numpy_dem,result_D8,cellinflat,numpy_riverinfo):
    while (not queue.Queue.empty(cellinflat)):
        cell = cellinflat.get()
        row = cell.row
        col = cell.col
        ele = cell.ele
        if result_D8[row][col] > 0:
            continue
        surroundingele = [numpy_dem[row-1][col-1], numpy_dem[row-1][col], numpy_dem[row-1][col+1],
                          numpy_dem[row][col-1], numpy_dem[row][col+1], numpy_dem[row+1][col-1],
                          numpy_dem[row+1][col], numpy_dem[row+1][col+1]]
        surroundingdir = [result_D8[row - 1][col - 1], result_D8[row - 1][col], result_D8[row - 1][col + 1],
                          result_D8[row][col - 1], result_D8[row][col + 1], result_D8[row + 1][col - 1],
                          result_D8[row + 1][col], result_D8[row + 1][col + 1]]
        surroundingriverinfo = [numpy_riverinfo[row - 1][col - 1], numpy_riverinfo[row - 1][col], numpy_riverinfo[row - 1][col + 1],
                          numpy_riverinfo[row][col - 1], numpy_riverinfo[row][col + 1], numpy_riverinfo[row + 1][col - 1],
                          numpy_riverinfo[row + 1][col], numpy_riverinfo[row + 1][col + 1]]

        flowright = 0
        sur = [0,0,0,0,0,0,0,0]
        for i in range(0,8):
            if surroundingele[i] <= ele and surroundingdir[i] > 0:
                flowright = 1 #This raster cell has right to flow into downstream raster cell.
            else:
                surroundingriverinfo[i] = 0
        maxriverinfo = numpy.max(surroundingriverinfo)
        for i in range(0,8):
            if surroundingriverinfo[i] == maxriverinfo and surroundingele[i] <= ele and surroundingdir[i] > 0:
                sur[i] = 1
        if flowright == 0:
            cellinflat.put(cell)
        else:
            dir = dirwithsurele(cell, numpy_dem,2,sur)
            result_D8[row][col] = dir
    return result_D8

def dirwithsurele(cell,numpy_dem,flag,sur2):                                           #calculate flow direction for these grid cells without unique downstream slope
    window_size = 1
    (row,col) = numpy.shape(numpy_dem)
    i = cell.row
    j = cell.col
    if flag == 1:
        surrounding = [(numpy_dem[i][j] - numpy_dem[i - 1][j - 1]) / 1.414, (numpy_dem[i][j] - numpy_dem[i - 1][j]),
                       (numpy_dem[i][j] - numpy_dem[i - 1][j + 1]) / 1.414, (numpy_dem[i][j] - numpy_dem[i][j - 1]),
                       (numpy_dem[i][j] - numpy_dem[i][j + 1]), (numpy_dem[i][j] - numpy_dem[i + 1][j - 1]) / 1.414,
                       (numpy_dem[i][j] - numpy_dem[i + 1][j]), (numpy_dem[i][j] - numpy_dem[i + 1][j + 1]) / 1.414]
    else:
        surrounding = [0,0,0,0,0,0,0,0]
        for m in range(0,8):
            if sur2[m] == 0:
                surrounding[m] = -9999
            else:
                if m == 0:
                    surrounding[m] = (numpy_dem[i][j] - numpy_dem[i - 1][j - 1]) / 1.414
                elif m == 1:
                    surrounding[m] = (numpy_dem[i][j] - numpy_dem[i - 1][j])
                elif m == 2:
                    surrounding[m] = (numpy_dem[i][j] - numpy_dem[i - 1][j + 1]) / 1.414
                elif m == 3:
                    surrounding[m] = (numpy_dem[i][j] - numpy_dem[i][j - 1])
                elif m == 4:
                    surrounding[m] = (numpy_dem[i][j] - numpy_dem[i][j + 1])
                elif m == 5:
                    surrounding[m] = (numpy_dem[i][j] - numpy_dem[i + 1][j - 1]) / 1.414
                elif m == 6:
                    surrounding[m] = (numpy_dem[i][j] - numpy_dem[i + 1][j])
                elif m == 7:
                    surrounding[m] = (numpy_dem[i][j] - numpy_dem[i + 1][j + 1]) / 1.414

    numofds, location, maxvalue, sur = findnumofdownstream(surrounding, 1)                   #calculate the number of grid cells with maximum downstream slope

    while numofds > 1:
        for i2 in range(0, numpy.shape(sur)[0]):
            if sur[i2] == 0:
                surrounding[i2] = 9999
            else:
                sumele = 0
                if i2 == 0:
                    x = i - 1
                    y = j -1
                elif i2 == 1:
                    x = i - 1
                    y = j
                elif i2 == 2:
                    x = i - 1
                    y = j + 1
                elif i2 == 3:
                    x = i
                    y = j - 1
                elif i2 == 4:
                    x = i
                    y = j + 1
                elif i2 == 5:
                    x = i + 1
                    y = j - 1
                elif i2 == 6:
                    x = i + 1
                    y = j
                elif i2 == 7:
                    x = i + 1
                    y = j + 1
                for m in range(-1*window_size,window_size+1):
                    for n in range(-1*window_size,window_size+1):
                        if (x + m >= row) or (y + n >= col) or (x + m < 0) or (y + n < 0):  #if windows over the edge of window,ele of center cell would be input.
                            sumele += numpy_dem[x][y]
                        else:
                            sumele += numpy_dem[x+m][y+n]
                avgele = sumele/((2*window_size+1)*(2*window_size+1))
                surrounding[i2] = avgele
        numofds, location, _, _ = findnumofdownstream(surrounding, 2)                             #calculate the number of candidate grid cells with lowest elevation
        window_size += 1

    if location == 0:
        return 32
    elif location == 1:
        return 64
    elif location == 2:
        return 128
    elif location == 3:
        return 16
    elif location == 4:
        return 1
    elif location == 5:
        return 8
    elif location == 6:
        return 4
    elif location == 7:
        return 2


def riverinflat(numpy_D8withnoflat,numpy_D8withflat,numpy_river,numpy_riverorder):
    (row,col) = numpy.shape(numpy_D8withnoflat)
    for i in range(0,row):
        for j in range(0,col):
            if numpy_river[i][j] > 0:
                numpy_D8withnoflat[i][j] = numpy_D8withflat[i][j]
                numpy_riverorder[i][j] +=1
    return numpy_D8withnoflat,numpy_riverorder

def copyqueue(queue1):
    queue2 = queue.Queue()
    queue3 = queue.Queue()
    while (not queue.Queue.empty(queue1)):
        cell = queue1.get()
        queue2.put(cell)
        queue3.put(cell)
    return queue2,queue3

def main():
    numpy_D8result,cellnodirinD8,cellinflat1 = D8(numpy_dem)                #calculate the flow direction value with unique downstream slope
    numpy_D8result2 = improveD8(numpy_dem, numpy_D8result, cellnodirinD8)
    numpy_riverorder = numpy.zeros(numpy.shape(numpy_D8result))
    numpy_D8result3 =  copy.deepcopy(numpy_D8result2)
    numpy_river = numpy.zeros(numpy.shape(numpy_D8result2))

    for i in range(0,round):
        print(i)
        cellinflat2, cellinflat1 = copyqueue(cellinflat1)
        numpy_D8result2copy = copy.deepcopy(numpy_D8result2)
        if i > 0:
            numpy_D8result2copy, numpy_riverorder = riverinflat(numpy_D8result2copy, numpy_D8result3, numpy_river,
                                                            numpy_riverorder)

        numpy_D8result3 = improvePF(numpy_dem, numpy_D8result2copy, cellinflat2,numpy_riverorder)

        result_D81 = arcpy.NumPyArrayToRaster(numpy_D8result3, lowerleft_point, cellsizex, cellsizey, -9999)
        result_acc1 = arcpy.sa.FlowAccumulation(result_D81)
        threshold = startacc-i*tolerance

        river = arcpy.sa.SetNull(result_acc1,1,"VALUE < {0}".format(threshold))
        river.save(r"in_memory\river"+str(threshold)+".tif")
        numpy_river = arcpy.RasterToNumPyArray(river,nodata_to_value=-9999)
        order = arcpy.NumPyArrayToRaster(numpy_riverorder, lowerleft_point, cellsizex, cellsizey, -9999)
        order.save(r"in_memory\order" + str(threshold) + ".tif")

    result_D8 = arcpy.NumPyArrayToRaster(numpy_D8result3, lowerleft_point, cellsizex, cellsizey, -9999)
    result_D8.save(path_result)
main()



