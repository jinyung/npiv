#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
# from scipy.ndimage import gaussian_filter1d

# ------------------

def readtps(input):
    """
    Function to read a .TPS file

    Args:
        input (str): path to the .TPS file

    Returns:
        lm (str list): info extracted from 'LM=' field
        im (str list): info extracted from 'IMAGE=' field
        id (str list): info extracted from 'ID=' filed
        coords: returns a 3D numpy array if all the individuals have same
                number of landmarks, otherwise returns a list containing 2d
                matrices of landmarks
    """

    # open the file
    tps_file = open(input, 'r')  # 'r' = read
    tps = tps_file.read().splitlines()  # read as lines and split by new lines
    tps_file.close()

    # initiate lists to take fields of "LM=","IMAGE=", "ID=" and the coords
    lm, im, ID, coords_array = [], [], [], []

    # looping thru the lines
    for i, ln in enumerate(tps):

        # Each individual starts with "LM="
        if ln.startswith("LM"):
            # number of landmarks of this ind
            lm_num = int(ln.split('=')[1])
            # fill the info to the list for all inds
            lm.append(lm_num)
            # initiate a list to take 2d coordinates
            coords_mat = []

            # fill the coords list by reading next lm_num of lines
            for j in range(i + 1, i + 1 + lm_num):
                coords_mat.append(tps[j].split(' '))  # split lines into values

            # change the list into a numpy matrix storing float vals
            coords_mat = np.array(coords_mat, dtype=float)
            # fill the ind 2d matrix into the 3D coords array of all inds
            coords_array.append(coords_mat)
            # coords_array.append(coords_mat)

        # Get info of IMAGE= and ID= fields
        if ln.startswith("IMAGE"):
            im.append(ln.split('=')[1])

        if ln.startswith("ID"):
            ID.append(ln.split('=')[1])

    # check if all inds contains same number of landmarks
    all_lm_same = all(x == lm[0] for x in lm)
    # if all same change the list into a 3d numpy array
    if all_lm_same:
        coords_array = np.dstack(coords_array)

    # return results in dictionary form
    return {'lm': lm, 'im': im, 'id': ID, 'coords': coords_array}

# ------------------

def cent(contour):
    """
    Function to calculate centroids

    Args:
        contours (numpy ndarray): xy(z) coordinates

    Returns:
        centroid (tuple)
    """
    return np.apply_along_axis(np.mean, 0, contour)

# ------------------

def vangle(v1, v2, degree = True):
    """
    Returns angle between two vectors

    Args:
        v1 (tuple): vector 1
        v2 (tuple): vector 2
        degree (boolean): return degree is True, else radian

    Returns:
        angle (float)
    """
    # modified from https://stackoverflow.com/questions/2827393/
    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)
    angle = np.arccos(np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0))
    if degree:
        angle = angle * 180. / np.pi
    return angle

# ------------------

def langle(body_land, limb_land, degree = True):
    """
    Function to calculate angles of nauplius limbs

    Args:
        body_land (2d numpy ndarray): body landmarks or contours to calculate
            centroid
        limb_land (2d numpy ndarray): limb landmarks

    Returns:
        limb_angle (list): 1st value is ant1 angle, 2nd is ant2, 3rd is mand
    """

    # get the centroid
    cen = cent(body_land)

    # define vectors for calculation of angles, using centroid as origin
    ref_vect = body_land[1, ] - cen
    limb_vect = limb_land - cen

    # calculate angle based on cosine law
    angle = np.apply_along_axis(vangle, 1, limb_vect, ref_vect, degree = degree)

    return np.array(angle, dtype = float)

# ------------------

def tlangle(body_land, limb_land, degree = True):
    """
    Wrapper for langle for time series

    Args:
        body_land (3d numpy ndarray): body landmarks or contours to calculate
            centroid, 3rd dimension is the frames
        limb_land (3d numpy ndarray): limb landmarks,
            3rd dimension is the frames

    Returns:
        limb_angle (2d numpy ndarray): 1st column is ant1 angle,
            2nd is ant2, 3rd is mand
    """
    nframe = body_land.shape[2]
    result = np.empty((nframe, 3))
    for i in range(nframe):
        result[i,:] = langle(body_land[:, :, i], limb_land[:, :, i], 
                              degree = degree)
    return result

# ------------------

def sangle(body_land, limb_land, degree = True):
  """
  Calculate angular separation among the limbs
  
  Returns:
    list of 1.ant1-ant2 angle, 2.ant1-mand angle and 3.ant2-mand angle
  """
  # centroid
  cen = cent(body_land)
  
  # limb vectors
  limb_vect = limb_land - cen
  
  # calculate for all limb combinations
  result = []
  for pair in [(0, 1), (0, 2), (1, 2)]:
    result.append(vangle(limb_vect[pair[0], ], limb_vect[pair[1], ], 
                  degree = degree))
  # return
  return np.array(result)
  
# ------------------

def tsangle(body_land, limb_land, degree = True):
    """
    Wrapper for sangle for time series

    Args:
        body_land (3d numpy ndarray): body landmarks or contours to calculate
            centroid, 3rd dimension is the frames
        limb_land (3d numpy ndarray): limb landmarks,
            3rd dimension is the frames

    Returns:
        limb_angle (2d numpy ndarray): 1st column is ant1 angle,
            2nd is ant2, 3rd is mand
    """
    nframe = body_land.shape[2]
    result = np.empty((nframe, 3))
    for i in range(nframe):
        result[i,:] = sangle(body_land[:, :, i], limb_land[:, :, i], 
                              degree = degree)
    return result

# ------------------

# https://stackoverflow.com/a/26337730
def smooth(y, win_size = 5):
  """
  Calculate simple moving average smoothing
  
  Args:
    y (list/ 1d ndarray): series to be smoothed
    win_size (int): moving window size
  """
  box = np.ones(win_size) / win_size
  y_smooth = np.convolve(y, box, mode = 'same')
  return y_smooth
  
# ------------------
def plot_angle(body_tps_file, limb_tps_file, t = 'id', scale = 0.5, cut = None, 
               win_size = 5, col = [u'#1f77b4', u'#ff7f0e', u'#2ca02c'], 
               degree = True, sub = None, axlab = True, **kwargs):
  """
  Wrapper for to plot tlangle
  
  Args:
      body_tps_file(char): tps file path containing body landmarks
      limb_tps_file (char): tps file path containing limb landmarks
      t (list/1d ndarray): time label in ms/ 'id'
      scale (float): used for scaling 't' if 'id' is used
      win_size (int): to be passed to smooth() for smoothing
      col: control color of markers/ lines
      cut (int): index to cut first n-th frame to plot, useful for animation
      degree (bool): to express in degree if true, else in radian
      sub (matplotlib object): for easier subplots integration 
      axes (bool): draw axes labels if True
  
  Returns:
    A plot
  """
  # reading files
  body_tps = readtps(body_tps_file)
  limb_tps = readtps(limb_tps_file)
  body_land = body_tps['coords']
  limb_land = limb_tps['coords']
  
  # sanity check
  if len(body_tps['im']) is not len(limb_tps['im']):
    raise ValueError('length of frames in body and limb are different. ' + 
                     'There are %s frames in body and %s frames in limb' 
                     % (len(body_tps['im']), len(limb_tps['im'])))
  if not body_tps['id']:
    if body_tps['id'] is not limb_tps['id']:
      print('Warnings: body_tps and limb_tps has different ID! '+ 
            'Did you use different time series?')
    
  # calculate
  result = tlangle(body_land, limb_land, degree = degree)
  result_smooth = np.apply_along_axis(smooth, 0, result, win_size)
  result_smooth[-(win_size/2):, :] = None  # NA the head and tail
  result_smooth[:(win_size/2), :] = None

  # cut frames out
  if cut is not None:
    result = result[:cut, :]
    result_smooth = result_smooth[:cut, :]

  # add in time axis
  if t is 'id':
    t = np.array(body_tps['id']).astype(int)
    t = (t - t[0]) * scale
  
  # set plot limit  
  xlim = [min(t), max(t)]
  ylim = [0, 190 if degree == True else 190./180.*np.pi] 
  
  if cut is not None:
    t = t[:cut]
      
  # plot
  if sub is None:
    sub = plt
  for i in range(3):
    sub.scatter(t, result[:, i], facecolors = 'none', edgecolors = col[i], **kwargs)
    sub.plot(t, result_smooth[:, i], '-', color = col[i])
  sub.set_xlim(xlim)
  sub.set_ylim(ylim)
  
  # label axes
  if axlab is True:
    sub.set_ylabel('Angle (%s)' % ('radian' if degree == False else 'degree'))
    sub.set_xlabel('time (ms)')

# ------------------
def kinemate(body_tps_file, limb_tps_file, img_dir, r = '5', degree = True):
  
  # read the tps files
  body_tps = readtps(body_tps_file)
  limb_tps = readtps(limb_tps_file)
  body_land = body_tps['coords']
  limb_land = limb_tps['coords']
  
  # create result folder
  out_dir = os.path.join(img_dir, "kinemate_output")
  if not os.path.exists(out_dir):
    os.makedirs(out_dir)
  
  # loop thru images
  for i, img_list in enumerate(body_tps['im']):
    img = mpimg.imread(os.path.join(img_dir, img_list))
    fig = plt.figure(figsize = (14, 5))
    
    # subplot-1: the video
    sub1 = fig.add_subplot(1, 2, 1)
    sub1.axis('off')
    sub1.imshow(img, cmap = 'gray')
    plt.scatter(limb_land[:, 0, i], img.shape[0] - limb_land[:, 1, i],  # invert y-axis 
                facecolors = 'none', 
                edgecolors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c'])
    
    # subplot-2: the angle plot
    sub2 = fig.add_subplot(1, 2, 2)
    plot_angle(body_tps_file, limb_tps_file, degree = degree, 
               cut = i + 1, sub = sub2)
    tmp_img_name = os.path.join(out_dir, 'tmp_%04d.jpg' % (i))
    fig.savefig(os.path.join(out_dir, tmp_img_name), dpi = 167, 
                bbox_inches = 'tight')
    plt.close(fig)

  # create video
  vid_name = os.path.basename(os.path.normpath(img_dir))
  ffmpeg_input = os.path.join(out_dir, 'tmp_%04d.jpg')
  ffmpeg_output = os.path.join(out_dir,
                              'kinemate_%s.mp4' % (vid_name))
  ffmpeg_option = '-c:v libx264 -vf "fps=25,format=yuv420p"'
  ffmpeg_cmd = 'ffmpeg -r %s -i %s %s -r 25 %s' % (
               r, ffmpeg_input, ffmpeg_option, ffmpeg_output)
  os.system(ffmpeg_cmd)

  # remove tmp imgs
  tmp_imgs = [x for x in os.listdir(out_dir) if x.startswith('tmp_')]
  for i in tmp_imgs:
    os.remove(os.path.join(out_dir, i))

# #----------------
# def find_maxima (x):
# 	return (np.diff(np.sign(np.diff(x))) < 0).nonzero()[0] + 1 
# 
# def find_minima (x):
# 	return (np.diff(np.sign(np.diff(x))) > 0).nonzero()[0] + 1 
