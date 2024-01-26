#!/usr/bin/env python

import cv2
import os
import sys
import makeCompositeImg
import concurrent.futures


def formatAdd0Number(number,nDigits):
    number=str(number)
    while len(number)<nDigits:
        number="0"+number
    return number

num_threads=4

inputFolder=sys.argv[1] + "/"
outputFolder=sys.argv[2]
num_threads=int(sys.argv[3])

outputComposite = outputFolder + "/outComposite/"
outputWellPics = outputFolder + "/outWellPics/"
outPutWellSmall = outputFolder + "/outWellSmall/"

if not os.path.exists(outputWellPics): os.mkdir(outputWellPics)
if not os.path.exists(outputComposite): os.mkdir(outputComposite)
if not os.path.exists(outPutWellSmall): os.mkdir(outPutWellSmall)

files=os.listdir(inputFolder)
filePerField_Well = {}
for file in files:
    splitted=file.split("_")
    well=splitted[1][1]+formatAdd0Number(splitted[1][2:],2)
    field="fld"+formatAdd0Number(splitted[3][1:],2)
    if not well in filePerField_Well: filePerField_Well[well] = {}
    if not field in filePerField_Well[well]: filePerField_Well[well][field] = []
    filePerField_Well[well][field].append(file)

fieldList=[*list(range(2,9)),*[14,13,12,1,11,10,9],*list(range(15,22)),*list(range(28,21,-1))]
nFields=max(fieldList)
imagePerField=[0]*nFields
for i in range(len(fieldList)): imagePerField[fieldList[i]-1]=i

#call makeCompositeImg.py
filePerWell = {}
argsForMakeCompositeImg = {}
for well in filePerField_Well:
    for field in filePerField_Well[well]:
        paths = [inputFolder+file for file in filePerField_Well[well][field]]
        compositeImgPath = outputComposite+well+field+".png"
        argsForMakeCompositeImg[well + field] = (paths, compositeImgPath)
        fieldInt = int(field[3:])
        if not well in filePerWell: filePerWell[well]=[0]*nFields
        filePerWell[well][imagePerField[fieldInt-1]]=compositeImgPath

num_threads = 4
with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:

    # Submit tasks to the thread pool
    futures = []
    for wellField in argsForMakeCompositeImg:
        future = executor.submit(makeCompositeImg.makeCompositeImg, *argsForMakeCompositeImg[wellField])
        futures.append(future)

    # Wait for all tasks to complete
    concurrent.futures.wait(futures)


H=2160-1
W=2160-1
overlap=105

def concat_tile(im_list_2d):
    return cv2.vconcat([cv2.hconcat(im_list_h) for im_list_h in im_list_2d])

def mergeField(well):
    imgList = list((cv2.imread(filePerWell[well][i]) for i in range(nFields)))
    for i in range(nFields):
        maxH= H-overlap if i<21 else H
        maxW= W-overlap if (i+1)%7 > 0 else W
        imgList[i]=imgList[i][0:maxH,0:maxW]
    
    output=concat_tile([imgList[0:6],imgList[7:13],imgList[14:20],imgList[21:27]])
    makeCompositeImg.overwriteImage(outputWellPics+"/"+well+".png", output)
    outputSmall=cv2.resize(output, (0,0), fx=0.25, fy=0.25)
    makeCompositeImg.overwriteImage(outPutWellSmall+"/"+well+".png", outputSmall)

with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:

    # Submit tasks to the thread pool
    futures = []
    for well in filePerWell.keys():
        future = executor.submit(mergeField, well)
        futures.append(future)

    # Wait for all tasks to complete
    concurrent.futures.wait(futures)

