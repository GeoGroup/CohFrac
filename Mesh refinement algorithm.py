# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import gmsh
import dxfgrabber as df
import math
import sys
import string
from matplotlib import path

# To avoid the prblem of the coding format
# Change from ANSI_1252 to ANSI_936


file_name = 'smooth_t_06_gai.dxf'
# gmsh.finalize()
try:
    f = open(file_name, 'r+')
    flist = f.readlines()
except:
    f = open(file_name, 'r+', encoding='UTF-8')
    flist = f.readlines()

f.close()
flist[15] = 'ANSI_936\n'
f = open(file_name, 'w+')
f.writelines(flist)
f.close()

dxf = df.readfile(file_name)
# Get the vertex on the closed line
polypos = {}  # Store the poistion
points = []
# lc=[]
group = []
lc_point = {}
cable_group = []
cable_pos = []
# max_cable_group=0
I = 0  # store the number of rock and fault
for e in dxf.entities:
    # print(e.dxftype,e.layer)
    if e.dxftype == 'POLYLINE':
        # I=I+1
        t = e.points
        pt = [[i[0]-5, i[1]-5] for i in t]
        polypos[I] = pt  # np.vstack((ptx,pty))
        I = I + 1


def point_in_list(p, l, err):
    if l:
        if (len(l[0]) == 2):
            ld = [((i[0] - p[0]) ** 2 + (i[1] - p[1]) ** 2) ** 0.5 for i in l]
            if min(ld) < err:
                return [True, ld.index(min(ld))]
            else:
                return [False, 0]
    else:
        return [False, 0]

err = 0.08
for i in range(len(polypos)):
    pt = polypos[i]
    for ipt in pt:
        tpos = point_in_list(ipt, points, err)
        if not tpos[0]:
            if abs(abs(ipt[0])-5) < err:
                ipt[0] = ipt[0]/abs(ipt[0])*5
            if abs(abs(ipt[1])-5) < err:
                ipt[1] = ipt[1]/abs(ipt[1])*5
            points.append(ipt)
#
# ptx = [i[0] for i in points]
# pty = [i[1] for i in points]
# xmin = min(ptx)
# xmax = max(ptx)
# ymin = min(pty)
# ymax = max(pty)

# plt.plot(ptx, pty,'*')
# plt.show()

filename='dfn_model_gen_t_06.py'
f=open(filename,'w')
f.write('# encoding: utf-8\n')
f.write('from abaqus import *\n')
f.write('from abaqusConstants import *\n')
f.write('from caeModules import *\n')
f.write('s = mdb.models["Model-1"].ConstrainedSketch(name="__profile__", sheetSize=10.0)\n')
f.write('s.rectangle(point1=(-5*2, -5*2), point2=(5*2, 5*2))\n')
f.write('p = mdb.models["Model-1"].Part(name="Part-1",dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)\n')
f.write('p.BaseShell(sketch=s)\n')
f.write('s0 = mdb.models["Model-1"].ConstrainedSketch(name="__profile__",sheetSize=200.0)\n')
f.write('g, v, d, c = s0.geometry, s0.vertices, s0.dimensions, s0.constraints\n')
# Deterimine the point in the range of refinement
for i in range(I):
    pt = polypos[i]
    for j in range(len(pt)-1):
        tmp0 = pt[j]
        idx = point_in_list(tmp0, points, err)
        pt1 = points[idx[1]]
        tmp1 = pt[j+1]
        idx = point_in_list(tmp1, points, err)
        pt2 = points[idx[1]]
        dis = math.sqrt((pt1[0]-pt2[0])**2+(pt1[1]-pt2[1])**2)
        if (pt1 != pt2) & (dis > 0.01):
            line = 's0.Line(point1=(' + str(pt1[0]*2) + ',' + str(pt1[1]*2) + '), point2=(' + str(pt2[0]*2) + ',' + str(pt2[1]*2) + '))'
            f.write('%s\n' % line)

f.write("p1 = mdb.models['Model-1'].parts['Part-1']\n")
f.write("pickedFaces = p1.faces[0:1]\n")
f.write("p1.PartitionFaceBySketch(faces=pickedFaces, sketch=s0)")

f.close()
