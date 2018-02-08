#!/usr/bin/env python

'''
---------------------------------------------
Title: .png to .mp4
Author: IJ Rodriguez
Date: 2018.02.08
Available at https://www.github.com/space-isa
---------------------------------------------

This code pulls a set of images stored in a directory, sorts them numerically, 
and compiles them into an mp4 file. Other supported output files include .mov and .gif.

(Images pulled in this case were of the form: 'Frame_XXX.png')

Feel free to use and modify this code!

'''

import imageio
import os, sys
import re

dir_path = '/Users/Isabel/Documents/Git/PIC/Fulltest4/'

images = []

for names in os.listdir(folder_path):
    if names.startswith('Frame'):
        path = dir_path + names
        images.append(path)
        images.sort(key = lambda f: int(filter(str.isdigit, f)))

writer = imageio.get_writer('Fulltest4.1.mp4', fps = 10)

for im in images:
    writer.append_data(imageio.imread(im))
writer.close()

#print(images)
