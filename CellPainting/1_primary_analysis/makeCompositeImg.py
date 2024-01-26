import cv2
import numpy as np
import os

def overwriteImage(name, image):
    if os.path.exists(name):
        os.remove(name)
    cv2.imwrite(name, image)

def convertHexToBGR(hex):
    return np.array(tuple(int(hex[i:i+2], 16) for i in (5, 3, 1)))

# inputPaths=["test/p1_wB5_t1_m12_c1_z0_l1_o0.tiff", "test/p1_wB5_t1_m12_c2_z0_l1_o0.tiff", "test/p1_wB5_t1_m12_c3_z0_l1_o0.tiff", "test/p1_wB5_t1_m12_c4_z0_l1_o0.tiff"]

def makeCompositeImg(inputPaths, outputPath):
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
    min_cyto, max_cyto = 0, 30
    min_mito, max_mito = 2, 20
    min_nucl, max_nucl = 0, 10
    min_endo, max_endo = 5, 50

    # Rescale the intensities of each channel
    cytoRescaled = np.interp(cytoChannel, (min_cyto, max_cyto), (0, 255)).astype(np.uint8)
    mitoRescaled = np.interp(mitoChannel, (min_mito, max_mito), (0, 255)).astype(np.uint8)
    nuclRescaled = np.interp(nuclChannel, (min_nucl, max_nucl), (0, 255)).astype(np.uint8)
    endoRescaled = np.interp(endoChannel, (min_endo, max_endo), (0, 255)).astype(np.uint8)

    #overwriteImage(outputFolder+'cytoB5_m12.png', cytoRescaled)
    #overwriteImage(outputFolder+'mitoB5_m12.png', mitoRescaled)
    #overwriteImage(outputFolder+'nuclB5_m12.png', nuclRescaled)
    #overwriteImage(outputFolder+'endoB5_m12.png', endoRescaled)

    color_cyto_hex = '#E0BB00'
    colored_mito_hex = '#FF0000'
    colored_nucl_hex = '#0000FF'
    colored_endo_hex = '#00C354'

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

    #overwriteImage(outputFolder+'cytoB5_m12.png', colored_cyto)
    #overwriteImage(outputFolder+'mitoB5_m12.png', colored_mito)
    #overwriteImage(outputFolder+'nuclB5_m12.png', colored_nucl)
    #overwriteImage(outputFolder+'endoB5_m12.png', colored_endo)

    # Create a composite image with each channel represented by a color
    composite_image = colored_cyto + colored_mito + colored_nucl + colored_endo
    composite_image[composite_image > 255] = 255
    overwriteImage(outputPath, composite_image)

