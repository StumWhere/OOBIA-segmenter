"""
Created on Tue Apr 17 19:29:56 2012

@author: Alex Stum
"""

# -*- coding: utf-8 -*-

class obia:
    def __init__(self,obID,value,thresh,relative=False):
        self.obID=obID
        self.N=1
        self.var=0
        self.mean=value
        self.vals=[value]
        self.sum=float(value)
        self.thresh=thresh 
    
    def compare(self,value):
        '''evaluates individual pixels to be incorporated into self'''
        self.vals.append(value)
        newVar=np.var(self.vals)
        if newVar<=(self.thresh*self.mean)**2:
            self.N+=1
            self.var=newVar
            self.sum+=value
            self.mean=self.sum/self.N
            return self.obID
        else:
            self.vals.pop(-1)
            return 0


    def merge(self,right):
        '''compares to obia objects, and merges their values if comprable'''
        newVals=self.vals+right.vals
        newVar=np.var(newVals)
        newN=self.N+right.N
        newSum=self.sum+right.sum
        newMean=newSum/newN
        if newVar<=(self.thresh*newMean)**2:
#                print str(self.N)+':'+ str(right.N)+':'+str(newN)
            self.N=newN
            self.var=newVar
            self.sum=newSum
            self.mean=newMean
            self.vals=newVals
            return self.obID
        else: return 0


import numpy as np
import arcpy
from arcpy import env
import arcpy.sa as sa

arcpy.CheckOutExtension("Spatial")
env.workspace = r"D:\Dissertation\testing"
env.overwriteOutput = 1

#create raster objects from input grids

slp100 = sa.Raster("obia_degree")
rows = slp100.height
cols = slp100.width

cellSize = slp100.meanCellHeight
#get lower left point of extent to write output after analysis
llpnt = slp100.extent.lowerLeft
#get spatial reference object from the raster object
spref = slp100.spatialReference

#convert raster object to numpy array

slp = arcpy.RasterToNumPyArray(slp100)
print "starting segmentation"

slp_diff = .1
slp_var=.1
groupCount = 2

#make a zeros grid with same dimentions as the DEM
obID = np.zeros_like(slp,dtype=np.int)
imgObjects=[] #this will hold the list of obia objects

#first object upper left hand corner
imgObjects.append(obia(1,slp[0,0],slp_var))
obID[0,0]=1

#evaluate first row
for x in range(1,cols):
    if abs(slp[0,x]-slp[0,x-1])<( slp_diff*slp[0,x]):
        if imgObjects[obID[0,x-1]-1].compare(slp[0,x]):
            obID[0,x]=obID[0,x-1]
        else: #not accepted into object, new object
            obID[0,x]=groupCount
            imgObjects.append(obia(groupCount,slp[0,x],slp_var))
            groupCount+=1
    else: #new object
        obID[0,x]=groupCount
        imgObjects.append(obia(groupCount,slp[0,x],slp_var))
        groupCount+=1

for y in range(1,rows):
    #evaluate the first column of each row
    if abs(slp[y,0]-slp[y-1,0])<( slp_diff*slp[y,0]):
        if imgObjects[obID[y-1,0]-1].compare(slp[y,0]):
            obID[y,0]=obID[y-1,0]
        else: #not accepted into object, new object
            obID[y,0]=groupCount
            imgObjects.append(obia(groupCount,slp[y,0],slp_var))
            groupCount+=1
    else:
        obID[y,0]=groupCount
        imgObjects.append(obia(groupCount,slp[y,0],slp_var))
        groupCount+=1
    for x in range(1,cols):
        limit=slp_diff*slp[y,x] 
        if abs(slp[y,x]-slp[y-1,x])<limit and abs(slp[y,x]-slp[y,x-1])<limit:  #what if the same...
            up=abs(slp[y,x]-imgObjects[obID[y-1,x]-1].mean)
            left=abs(slp[y,x]-imgObjects[obID[y,x-1]-1].mean)
            if up==left: #both objects have the same mean find compare object of most similar pixel
                up=abs(slp[y,x]-slp[y-1,x])
                left=abs(slp[y,x]-slp[y,x-1])
            if up<=left:
                if imgObjects[obID[y-1,x]-1].compare(slp[y,x]):
                    obID[y,x]=obID[y-1,x]
                    if obID[y,x]!=obID[y,x-1]: #if both objects aren't the same try merging
                        if imgObjects[obID[y,x]-1].merge(imgObjects[obID[y,x-1]-1]):         
                            old=obID[y,x-1]
                            obID=np.where(obID==old,obID[y,x],obID)
                            imgObjects[old-1]=[]
                elif imgObjects[obID[y,x-1]-1].compare(slp[y,x]): #compare to the other object then
                    obID[y,x]=obID[y,x-1]
                else: #not accepted into object, new object
                    obID[y,x]=groupCount
                    imgObjects.append(obia(groupCount,slp[y,x],slp_var))
                    groupCount+=1
            else: #left object more similar
                if imgObjects[obID[y,x-1]-1].compare(slp[y,x]):
                    obID[y,x]=obID[y,x-1]
                    if obID[y,x]!=obID[y-1,x]:
                        if imgObjects[obID[y,x]-1].merge(imgObjects[obID[y-1,x]-1]):         
                            old=obID[y-1,x]
                            obID=np.where(obID==old,obID[y,x],obID)
                            imgObjects[old-1]=[]
                elif imgObjects[obID[y-1,x]-1].compare(slp[y,x]):
                    obID[y,x]=obID[y-1,x]
                else: #not accepted into object, new object
                    obID[y,x]=groupCount
                    imgObjects.append(obia(groupCount,slp[y,x],slp_var))
                    groupCount+=1
        elif abs(slp[y,x]-slp[y-1,x])<limit: #is the pixel above similar
            if imgObjects[obID[y-1,x]-1].compare(slp[y,x]):
                obID[y,x]=obID[y-1,x]
            else: #not accepted into object, new object
                obID[y,x]=groupCount
                imgObjects.append(obia(groupCount,slp[y,x],slp_var))
                groupCount+=1
        elif abs(slp[y,x]-slp[y,x-1])<limit: #how about the one to left?
            if imgObjects[obID[y,x-1]-1].compare(slp[y,x]):
                    obID[y,x]=obID[y,x-1]
            else: #not accepted into object, new object
                obID[y,x]=groupCount
                imgObjects.append(obia(groupCount,slp[y,x],slp_var))
                groupCount+=1
        else: #both pixels are different, new object
            obID[y,x]=groupCount
            imgObjects.append(obia(groupCount,slp[y,x],slp_var))
            groupCount+=1

resultRast = arcpy.NumPyArrayToRaster(obID,llpnt,cellSize,cellSize)

#define the projection for the output raster object

arcpy.DefineProjection_management(resultRast,spref)

resultRast.save('new_slp_obi2') #save the output raster as a permanent file
    
#u=np.unique(obID)
#pop_grid=np.zeros_like(obID)
#
#for v in u: pop_grid[obID==v]=imgObjects[v-1].N
#
#population = arcpy.NumPyArrayToRaster(pop_grid,llpnt,cellSize,cellSize)
#arcpy.DefineProjection_management(population,spref)
#population.save('pop_grid') #save the output raster as a permanent file
#    

