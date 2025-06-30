import cv2
import numpy as np
import os

def overwriteImage(name, image):
    if os.path.exists(name):
        os.remove(name)
    cv2.imwrite(name, image)

def convertHexToBGR(hex):
    return np.array(tuple(int(hex[i:i+2], 16) for i in (5, 3, 1)))

def makeCompositeImg(inputPaths, outputComposite, outputRescaled, img_basename):
    
    channels=[path.split("/")[-1].split("_")[4] for path in inputPaths]
    if not "c1" in channels: raise Exception("No cyto channel")
    if not "c2" in channels: raise Exception("No mito channel")
    if not "c3" in channels: raise Exception("No nucl channel")
    if not "c4" in channels: raise Exception("No endo channel")
    
    cytoChannel = cv2.imread(inputPaths[channels.index("c1")], cv2.IMREAD_GRAYSCALE)
    mitoChannel = cv2.imread(inputPaths[channels.index("c2")], cv2.IMREAD_GRAYSCALE)
    nuclChannel = cv2.imread(inputPaths[channels.index("c3")], cv2.IMREAD_GRAYSCALE)
    endoChannel = cv2.imread(inputPaths[channels.index("c4")], cv2.IMREAD_GRAYSCALE)
    
    # Define the new minimum and maximum values for each channel
    # first value : cut low intensity, second: higher=less bright
    min_cyto, max_cyto = 1, 15
    min_mito, max_mito = 2, 18
    min_nucl, max_nucl = 0, 13
    min_endo, max_endo = 5, 40
    
    # Rescale the intensities of each channel
    cytoRescaled = np.interp(cytoChannel, (min_cyto, max_cyto), (0, 255)).astype(np.uint8)
    mitoRescaled = np.interp(mitoChannel, (min_mito, max_mito), (0, 255)).astype(np.uint8)
    nuclRescaled = np.interp(nuclChannel, (min_nucl, max_nucl), (0, 255)).astype(np.uint8)
    endoRescaled = np.interp(endoChannel, (min_endo, max_endo), (0, 255)).astype(np.uint8)
    
    color_cyto_hex = '#0000FF'
    colored_mito_hex = '#FF0000'
    colored_nucl_hex = '#FFFFFF'
    colored_endo_hex = '#00FF00'
    
    
    color_cyto = convertHexToBGR(color_cyto_hex) / 255.0
    color_mito = convertHexToBGR(colored_mito_hex) / 255.0
    color_nucl = convertHexToBGR(colored_nucl_hex) / 255.0
    color_endo = convertHexToBGR(colored_endo_hex) / 255.0
    
    # multiply greyscale image by color to get colored image
    colored_cyto = np.zeros((cytoRescaled.shape[0], cytoRescaled.shape[1], 3), dtype=np.uint8)
    colored_cyto = color_cyto * cytoRescaled[:, :, np.newaxis]
    colored_mito = np.zeros((mitoRescaled.shape[0], mitoRescaled.shape[1], 3), dtype=np.uint8)
    colored_mito = color_mito * mitoRescaled[:, :, np.newaxis]
    colored_nucl = np.zeros((nuclRescaled.shape[0], nuclRescaled.shape[1], 3), dtype=np.uint8)
    colored_nucl = color_nucl * nuclRescaled[:, :, np.newaxis]
    colored_endo = np.zeros((endoRescaled.shape[0], endoRescaled.shape[1], 3), dtype=np.uint8)
    colored_endo = color_endo * endoRescaled[:, :, np.newaxis]
    
    overwriteImage(outputRescaled + img_basename + "_cyto.png", colored_cyto)
    overwriteImage(outputRescaled + img_basename + "_mito.png", colored_mito)
    overwriteImage(outputRescaled + img_basename + "_nucl.png", colored_nucl)
    overwriteImage(outputRescaled + img_basename + "_endo.png", colored_endo)
    
    # Create a composite image with each channel represented by a color
    composite_image = colored_cyto + colored_mito + colored_nucl + colored_endo
    composite_image[composite_image > 255] = 255
    overwriteImage(outputComposite + img_basename + ".png", composite_image)

