# -*- coding: utf-8 -*-
"""
Chessboard clibration for CopterCurrents

Example from:
 https://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_calib3d/py_calibration/py_calibration.html
"""

import numpy as np
import cv2
import glob
import scipy.io as sio


# # define grid crorners
n_corners_vert = 6 # number of corners in vertical direction
n_corners_horz = 10 # number of corners in Horizontal direction 
square_size = 0.105 # square size in meters
subpixel_search_winSize = 21 # subpixels resoluion
path2images = '/media/d1/Drone_current_fit/githup/CopterCurrents/Calibration_examples/video/'

# images = glob.glob('*.jpg')
images = glob.glob( path2images + '*.tif')

# termination criteria
criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)

# prepare object points, like (0,0,0), (1,0,0), (2,0,0) ....,(6,5,0)
objp = np.zeros((n_corners_horz*n_corners_vert,3), np.float32)
objp[:,:2] = np.mgrid[0:n_corners_vert,0:n_corners_horz].T.reshape(-1,2) *square_size

# Arrays to store object points and image points from all the images.
objpoints = [] # 3d point in real world space
imgpoints = [] # 2d points in image plane.



n = 1

for fname in images:
    
    print fname
    

    img = cv2.imread(fname)
    gray = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)

    # Find the chess board corners
    ret, corners = cv2.findChessboardCorners(gray, (n_corners_vert,n_corners_horz),None)

    # If found, add object points, image points (after refining them)
    if ret == True:
        
        objpoints.append(objp)

        corners2 = cv2.cornerSubPix(gray,corners,(subpixel_search_winSize,subpixel_search_winSize),(-1,-1),criteria)
        imgpoints.append(corners2)

        # Draw and display the corners
        img2 = cv2.drawChessboardCorners(img, (n_corners_vert,n_corners_horz), corners2,ret)
 
        # save picture
        cv2.imwrite('corners_'+ str(n) + '.png',img2)
        
        print 'corners_'+ str(n) + '.png'
        n = n + 1

cv2.destroyAllWindows()


# Calibration ( mtx = camera matrix; dist =distortion coeff)
ret, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(objpoints, imgpoints, gray.shape[::-1],None,None)

# save Calibration to CopterCurrents_CamCalib format
shape_img = img.shape

fc = np.array([[mtx[0][0]],
               [mtx[1][1]]])      
cc = np.array([[mtx[0][2]],
               [mtx[1][2]]])
kc = dist
alpha_c = mtx[0][1] / mtx[0][0]
nx = np.float64(shape_img[1])
ny = np.float64(shape_img[0])
camera_offset_Z = np.float64(0)
source = 'OpenCV_script'

CopterCurrents_CamCalib = {'fc' : fc, 'cc' : cc, 'kc': kc, 'alpha_c': alpha_c, 'nx' : nx, 'ny' : ny, 'camera_offset_Z': camera_offset_Z, 'source' : source}

filename = 'CopterCurrents_CamCalib_' + str(shape_img[1]) + 'x' + str(shape_img[0]) + '.mat'

sio.savemat(filename, {'CopterCurrents_CamCalib': CopterCurrents_CamCalib})


# get mean error (reprojecting)
tot_error = 0
for i in xrange(len(objpoints)):
    imgpoints2, _ = cv2.projectPoints(objpoints[i], rvecs[i], tvecs[i], mtx, dist)
    error = cv2.norm(imgpoints[i],imgpoints2, cv2.NORM_L2)/len(imgpoints2)
    tot_error += error

print "total error: ", tot_error/len(objpoints)


 
 
