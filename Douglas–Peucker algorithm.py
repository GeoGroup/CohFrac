import dxfgrabber as df
import math
import ezdxf


def get_point_line_distance(point,line):
    point_x = point[0]
    point_y = point[1]
    line_s_x = line[0][0]
    line_s_y = line[0][1]
    line_e_x = line[1][0]
    line_e_y = line[1][1]
    # 若直线与y轴平行，则距离为点的x坐标与直线上任意一点的x坐标差值的绝对值
    if line_e_x - line_s_x == 0:
        return math.fabs(point_x - line_s_x)
    # 若直线与x轴平行，则距离为点的y坐标与直线上任意一点的y坐标差值的绝对值
    if line_e_y - line_s_y == 0:
        return math.fabs(point_y - line_s_y)
    # 斜率
    k = (line_e_y - line_s_y) / (line_e_x - line_s_x)
    # 截距
    b = line_s_y - k * line_s_x
    # 带入公式得到距离dis
    dis = math.fabs(k * point_x - point_y + b) / math.pow(k * k + 1, 0.5)
    return dis


def get_line_points(lp1, rp1):
    list_distance = []
    for j in range(rp1, lp1+1):
        point1 = i.points[j]
        line1 = [i.points[rp1], i.points[lp1]]
        dist = get_point_line_distance(point1, line1)
        list_distance.append(dist)
    mp = list_distance.index(max(list_distance))+rp1
    if max(list_distance) <= h:
        return list_mp.append(mp)
    else:
        list_mp.append(mp)
        line_r_points = [mp, rp1]
        line_l_points = [lp1, mp]
        return get_line_points(line_r_points[0], line_r_points[1]),get_line_points(line_l_points[0], line_l_points[1])


f = open('set10_poly.txt','w')
dxf = df.readfile("four_fractal.dxf")

dwg = ezdxf.new('R2010')
modespace = dwg.modelspace()

list_dis=[]
# list_mp = []

list_polyline=[]
h = 0.008
for i in dxf.entities:
    if i.dxftype == 'POLYLINE' and len(i.points)>=4:
        list_mp = []
        list_pl = []
        lp = len(i.points) - 1
        rp = 0
        get_line_points(lp, rp)

        list_mp.append(0)
        list_mp.append(len(i.points) - 1)

        list_mp = list(set(list_mp))
        list_mp.sort()
        for k in list_mp:
            lst = i.points[k]
            list_pl.append(lst)
            f.writelines(str(lst[0]) + ',' + str(lst[1]) + ',' + str(lst[2]) + '\n')
        list_polyline.append(list_pl)
        modespace.add_lwpolyline(list_pl)
    else:
        list_pl = []
        for p in range(0,len(i.points)):
            lst = i.points[p]
            list_pl.append(lst)
            f.writelines(str(lst[0]) + ',' + str(lst[1]) + ',' + str(lst[2]) + '\n')
        list_polyline.append(list_pl)
        modespace.add_lwpolyline(list_pl)
print(list_polyline,len(list_polyline))
dwg.saveas("smooth_four_fractal.dxf")
#
# filename='D:\\abaqus2021\\temp\\dfn_model_gen.py'
# f=open(filename,'w')
# f.write('# encoding: utf-8\n')
# f.write('from abaqus import *\n')
# f.write('from abaqusConstants import *\n')
# f.write('from caeModules import *\n')
# f.write('s = mdb.models["Model-1"].ConstrainedSketch(name="__profile__", sheetSize=10.0)\n')
# f.write('s.rectangle(point1=(-0.5, -0.5), point2=(0.5, 0.5))\n')
# f.write('p = mdb.models["Model-1"].Part(name="Part-1",dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)\n')
# f.write('p.BaseShell(sketch=s)\n')
# f.write('s0 = mdb.models["Model-1"].ConstrainedSketch(name="__profile__",sheetSize=200.0)\n')
# f.write('g, v, d, c = s0.geometry, s0.vertices, s0.dimensions, s0.constraints\n')
# error =1e-3
# for i in list_polyline:
#     if len(i) > 2:
#         for j in range(len(i)-1):
#             pt1 = [i[j][0]-0.5,i[j][1]-0.5]
#             pt2 = [i[j+1][0]-0.5,i[j+1][1]-0.5]
#             if abs(abs(pt1[0])-0.5) < error:
#                 pt1[0] = pt1[0]/abs(pt1[0])*0.5
#             if abs(abs(pt1[1])-0.5) < error:
#                 pt1[1] = pt1[1]/abs(pt1[1])*0.5
#             if abs(abs(pt2[0])-0.5) < error:
#                 pt2[0] = pt2[0]/abs(pt2[0])*0.5
#             if abs(abs(pt2[1])-0.5) < error:
#                 pt2[1] = pt2[1]/abs(pt2[1])*0.5
#             dist_line = math.sqrt((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2 )
#             if dist_line >5e-4:
#                 line= 's0.Line(point1=(' + str(pt1[0])+','+ str(pt1[1]) + '), point2=('+ str(pt2[0])+','+ str(pt2[1]) + '))'
#                 f.write('%s\n' % line)
#
# f.write('p1 = mdb.models["Model-1"].parts["Part-1"]\n')
# f.write('pickedFaces = p1.faces[0:1]\n')
# f.write('p1.PartitionFaceBySketch(faces=pickedFaces, sketch=s0)')
# f.close()
