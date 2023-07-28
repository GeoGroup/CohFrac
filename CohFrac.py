for var in dir():
    if not var.startswith('_'):
        del globals()[var]
# encoding: utf-8
from abaqus import *
from odbAccess import *
from abaqusConstants import *
from caeModules import *
from abaqus import mdb
from driverUtils import executeOnCaeStartup
import displayGroupOdbToolset as dgo
import os
import numpy as np
import re
import string

#Set parameters
current_directory = 'D://abaqus2021//temp//CZM'
ini_crack_p1 = [0.0, 0.1]
ini_crack_p2 = [0.0, -0.1]
MeshSizeSeed = 1
GroupInp = 'Model-Group'
CZMInp = "Model-CZM"
jobname ='Job-CZM3'
Rock_Elastic_Young_Modulus = 15e9
Rock_Elastic_Poission = 0.21
Rock_Permeability_Weight = 9800.0
Rock_Permeability_k = 1e-07
Rock_Permeability_Void_Ratio = 0.1
Cohe1_Elastic_E = 15e9
Cohe1_Elastic_G1 = 15e9
Cohe1_Elastic_G2 = 15e9
Cohe1_MaxsDamage_Nominal_Stress_Normal_only = 2e6
Cohe1_MaxsDamage_Nominal_Stress_First_Direction = 10e6
Cohe1_MaxsDamage_Nominal_Stress_Second_Direction = 10e6
Cohe1_MaxsDamage_DamageEvolution_Displacement = 0.001
Cohe1_FuildLeakoff = 1e-14
Cohe1_GapFlow = 0.001
Cohe2_MaxsDamage_Nominal_Stress_Normal_only = 6e6
Cohe2_MaxsDamage_Nominal_Stress_First_Direction = 20e6
Cohe2_MaxsDamage_Nominal_Stress_Second_Direction = 20e6
Step1_time = 0.1
Step2_time = 20
Model_xmax = 10
Model_xmin = -10
Model_ymax = 10
Model_ymin = -10
sxx = 0.8e6
syy = 3.2e6
szz = 7.8e6

#Generate initial cracks
os.chdir(current_directory)
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['Part-1']
s = p.features['Partition face-1'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s)
s1 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
g_num=0
for i in g.keys():
    g_num = g_num+1
s1.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s1,
    upToFeature=p.features['Partition face-1'], filter=COPLANAR_EDGES)
s1.Line(point1=((ini_crack_p1[0]+ini_crack_p2[0])/2, (ini_crack_p1[1]+ini_crack_p2[1])/2), point2=(ini_crack_p1[0], ini_crack_p1[1]))
s1.VerticalConstraint(entity=g[g_num+2], addUndoState=False)
s1.Line(point1=((ini_crack_p1[0]+ini_crack_p2[0])/2, (ini_crack_p1[1]+ini_crack_p2[1])/2), point2=(ini_crack_p2[0], ini_crack_p2[1]))
s1.VerticalConstraint(entity=g[g_num+3], addUndoState=False)
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
p.features['Partition face-1'].setValues(sketch=s1)
del mdb.models['Model-1'].sketches['__edit__']
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
p = mdb.models['Model-1'].parts['Part-1']
p.features['Partition face-1'].restore()
p = mdb.models['Model-1'].parts['Part-1']
print('Initial crack has generated')

# Generate mesh
p.seedPart(size=MeshSizeSeed, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Part-1']
f1 = p.faces
p.setMeshControls(regions=f1, elemShape=QUAD, allowMapped=False)
p = mdb.models['Model-1'].parts['Part-1']
p.generateMesh()
elemType1 = mesh.ElemType(elemCode=CPE4P, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=UNKNOWN_TRI, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Part-1']
f2 = p.faces
pickedRegions =(f2, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
print('Mesh has generated')

#Grouping
modelName = session.sessionState[session.currentViewportName]['modelName']
partName = "Part-1"
model = mdb.models[modelName]
p = model.parts[partName]
elemPair = []
totalNodel = []
if len(p.elements)==0:
    getWarningReply(message='Part without mesh...',buttons=(CANCEL,))
oldEdges = {}
for edge in p.edges:
    edgeNodes = edge.getNodes()
    edgeElems = edge.getElements()
    if len(edgeNodes)>len(edgeElems):
        continue
    else:
        edgeNodeLabels = [n.label for n in edgeNodes]
        totalNodel.append([edgeNodeLabels])
        for n in edgeNodes:
            for elemEdge in n.getElemEdges():
                n1,n2 = elemEdge.getNodes()
                if n.label==n1.label:
                    otherLabel = n2.label
                else:
                    otherLabel = n1.label
                if otherLabel in edgeNodeLabels:
                    edgeKey = "%d-%d"%(max([n1.label,n2.label]),min([n1.label,n2.label]))
                    if oldEdges.has_key(edgeKey):continue
                    oldEdges[edgeKey] = ""
                    e1,e2 = elemEdge.getElements()
                    elemPair.append([e1.label,e2.label])

mergeNodel = [num for sublist in totalNodel for inner_list in sublist for num in inner_list]
newNodel = list(set(mergeNodel))
node_group = p.nodes.sequenceFromLabels(newNodel)
p.Set(nodes=node_group, name='dfn_node')
newEle = list(set([num for sublist in elemPair for num in sublist]))
ele_group = p.elements.sequenceFromLabels(newEle)
p.Set(elements=ele_group, name='dfn_elem')
print('Model has been grouped')

#Out Inp
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=ON)
mdb.Job(name=GroupInp, model='Model-1', description='', type=ANALYSIS,
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
    scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN,
    numDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT,
    numCpus=1, numGPUs=0)
mdb.jobs[GroupInp].writeInput(consistencyChecking=OFF)


def judge_in(dat, vert):
    vert = np.sort(vert)
    dat = np.sort(dat, axis=1)
    if len(vert) == 2:
        T1 = [np.abs(dat[:, 0] - vert[0]),np.abs(dat[:, 1] - vert[1])]
        T1 = np.sum(T1, axis=0)
        pos = np.where(T1 < 0.1)[0]
    elif len(vert) == 3:
        T1 = [np.abs(dat[:, 0] - vert[0]) + np.abs(dat[:, 1] - vert[1]) + np.abs(dat[:, 2] - vert[2])]
        T1 = np.sum(T1, axis=0)
        pos = np.where(T1 < 0.1)[0]
    elif len(vert) == 4:
        T1 = [np.abs(dat[:, 0] - vert[0, 0]) + np.abs(dat[:, 1] - vert[0, 1]) + np.abs(dat[:, 2] - vert[0, 2]) + np.abs(dat[:, 3] - vert[0, 3])]
        T1 = np.sum(T1, axis=0)
        pos = np.where(T1 < 0.1)[0]
    return pos

filename = GroupInp+'.inp'
finp = open(filename, 'rt')
dim = 2
flag = ''
k = 0
ieletype = 0
igroup = 0
ingroup = 0
nnode = 0
nelement = 0
model={'node':[],'element':{'type': [], 'data': []},'group':{'name': [], 'data': []},'nodegroup':{'name': [], 'data': []}}

#read inp
while True:
    s = finp.readline()
    if not s:
        break
    k += 1
    if k % 100000 == 0:
        print(k)
    if '*' in s:
        if '*Node' in s:
            flag = 'node'
            empty_node=np.zeros((10000000, dim+1))
            model['node']=empty_node
        elif '*Element, type=' in s:
            ieletype += 1
            ti = s.find('=')
            flag = 'element'
            name = s[ti + 1:].strip()
            t = re.findall(r'\d*\.?\d*', name)
            no_empty= [x for x in t if x != '']
            nn = int(no_empty[0])
            model['element']['type'].append(name)
            model['element']['data'].append(np.zeros((1000000,nn+1)))
            ilength = 0
        elif '*Elset, elset=' in s:
            tmp = s
            if 'generate' not in s:
                ti = s.find('=')
                name = s[ti+1:].strip()
                ig = 0
            else:
                ti1 = s.find('=')
                ti = s.find(',')
                # name = s[ti1+1:ti[1]-1].strip()
                name = s[ti1 + 1:ti - 1].strip()
                ig = 1
            igroup += 1
            flag = 'group'
            model['group']['name'].append(name)
            model['group']['data'].append([])
        elif '*Nset, nset' in s:
            tmp = s
            if 'generate' not in s:
                ti = s.find('=')
                name = s[ti+1:].strip()
                ig = 0
            else:
                ti1 = s.find('=')
                ti = s.find(',')
                # name = s[ti1+1:ti[1]-1].strip()
                name = s[ti1 + 1:ti - 1].strip()
                ig = 1
            ingroup += 1
            flag = 'nodegroup'
            model['nodegroup']['name'].append(name)
            model['nodegroup']['data'].append([])
        else:
            flag = 'comment'
    # add inp date into model
    if flag == 'node' and '*' not in s:
        tmp = [float(num) if '.' in num else int(num) for num in s.split(',')]
        model['node'][nnode, :] = tmp
        nnode = nnode + 1
    elif flag == 'element' and '*' not in s:
        nelement += 1
        strtmp = '%d, ' * nn + '%d '
        tmp = [int(x) for x in re.findall(r'\d+', s)]
        ilength += 1
        model['element']['data'][ieletype-1][ilength-1, :] = tmp
    elif flag == 'group' and '*' not in s:
        if ig == 1:
            tmp = [int(x) for x in re.findall(r'\d+', s)]
            t = list(range(tmp[0], tmp[1] + 1, tmp[2]))
            model['group']['data'][igroup-1] = t
        else:
            tmp = [float(x) for x in re.findall(r'\d*\.?\d+', s)]
            tmp = [int(x) for x in tmp]
            model['group']['data'][igroup-1] += tmp
    elif flag == 'nodegroup' and '*' not in s:
        if ig == 1:
            tmp = [int(x) for x in re.findall(r'\d+', s)]
            t = list(range(tmp[0], tmp[1] + 1, tmp[2]))
            model['nodegroup']['data'][igroup-1] = t
        else:
            tmp = [float(x) for x in re.findall(r'\d*\.?\d+', s)]
            tmp = [int(x) for x in tmp]
            model['nodegroup']['data'][igroup-1] += tmp
finp.close()

# #remove the empty value
model['node'] = model['node'][:nnode, :]  # node
for i in range(len(model['element']['type'])):  # element
    remove_empty_element=model['element']['data'][i]
    model['element']['data'][i]=remove_empty_element[remove_empty_element[:,0] != 0]

nno = model['node'].shape[0]  # node number
nel = 0
for i in range(len(model['element']['type'])):
    nel += model['element']['data'][i].shape[0]  # element number

# allocation information
modu = {}
modu['node'] = model['node']
modu['element'] = {}
t = len(model['group']['name'])
nn = 0
set = {'pro':[],'setid':[]}
set['pro']=np.empty((nel,1),dtype=object)
set['pro']=np.array(set['pro']).reshape(-1, 1)
set['setid']=np.zeros((nel,1))

for i in range(t):
    tmp = model['group']['data'][i]
    for j in range(len(tmp)):
        nn += 1
        set['pro'][tmp[j]-1,0] = model['group']['name'][i]
        set['setid'][tmp[j]-1,0]= i+1
tdfn = model['nodegroup']['data'][0]

t = model['element']['type']
nn = 0
for i in range(len(t)):
    tl = model['element']['data'][i].shape[0]

t = model['element']['type']
nn = 0
for i in range(len(t)):
    tl = model['element']['data'][i].shape[0]
    for j in range(tl):
        nn += 1
        vert = model['element']['data'][i][j,:]
        id = int(vert[0])
        vert = vert[1:len(vert)]
        ele = {}
        ele['id'] = np.ones(len(vert), dtype=int) * id
        ele['setid'] = np.ones(len(vert), dtype=int) * set['setid'][id-1,0]
        ele['vert'] = vert
        ele['line'] = np.zeros((len(vert), 2), dtype=int)
        for ii in range(len(vert)-1):
            ele['line'][ii,:] = [vert[ii], vert[ii+1]]
        ele['line'][ii+1,:] = [vert[ii+1], vert[0]]
        ele['prop'] = {}
        for ii in range(len(vert)):
            tvert = ele['line'][ii,:]
            if len(np.setdiff1d(tvert, tdfn)) == 0:
                ele['prop'][ii] = 'DFN'
            else:
                ele['prop'][ii] = 'NDFN'
        pvert = np.arange(nno+(nn-1)*4+1, nno+(nn-1)*4+1+len(vert))
        ele['pvert'] = pvert
        ele['pline'] = np.zeros((len(pvert), 2), dtype=int)
        for ii in range(len(pvert)-1):
            ele['pline'][ii,:] = [pvert[ii], pvert[ii+1]]
        ele['pline'][ii+1,:] = [pvert[ii+1], pvert[0]]
        modu['element'][nn] = ele
print('Allocate information finish')

pool = {
    "id": np.zeros((1,1),dtype=int),
    "line": np.zeros((1,2),dtype=int),
    "pline": np.zeros((1,2),dtype=int),
    "prop": [],
    "setid": {0: '0'}
}
pool['prop'].append('0')

newnum = 0
update = {
    "line": {},
    "pnode": [],
    "linetype": {},
    "node": [np.array([]) for i in range(nno)],
}
index = np.zeros((nno,1),dtype=int)
delete_pos=[]

for i in range(nel):
    tid = modu['element'][i+1]['id']
    tline = modu['element'][i+1]["line"]
    tpline = modu['element'][i+1]["pline"]
    tprop = modu['element'][i+1]["prop"]
    tsetid = modu['element'][i+1]["setid"]
    for j in range(tline.shape[0]):
        tid2 = tid[j]
        tline2 = tline[j]
        tpline2 = tpline[j]
        tprop2 = tprop[j]
        tsetid2 = tsetid[j]
        pos = judge_in(pool["line"], tline2)
        if len(pos) != 0:
            if tprop2 == 'DFN':
                newnum = newnum + 1
                update["pnode"]= np.concatenate((update["pnode"], tline2), axis=0).astype(int)
                update["line"][newnum-1] = [tline2, tpline2, pool["pline"][pos,:][0]]
                update["linetype"][newnum-1] = 1
                for k in range(tline2.shape[0]):
                    index[tline2[k]-1] = 1
                    if len(update['node'][tline2[k] - 1]) == 0:
                        update['node'][tline2[k] - 1] = np.array([[1, tline2[k], tpline2[k], pool['pline'][pos[0], 1 - k]]])
                    else:
                        update['node'][tline2[k] - 1] = np.vstack((update['node'][tline2[k] - 1],
                                                                   [1, tline2[k], tpline2[k], pool['pline'][pos[0], 1 - k]]))
            else:
                newnum = newnum + 1
                update["pnode"]= np.concatenate((update["pnode"], tline2), axis=0)
                update["line"][newnum-1] = [tline2, tpline2, pool["pline"][pos,:][0]]
                update["linetype"][newnum-1] = 2
                for k in range(tline2.shape[0]):
                    index[tline2[k]-1] = 1
                    if len(update['node'][tline2[k] - 1]) == 0:
                        update['node'][tline2[k] - 1] = np.array([[1, tline2[k], tpline2[k], pool['pline'][pos[0], 1 - k]]])
                    else:
                        update['node'][tline2[k] - 1] = np.vstack((update['node'][tline2[k] - 1],
                                                                   [1, tline2[k], tpline2[k], pool['pline'][pos[0], 1 - k]]))
            pool['line'] = np.delete(pool['line'], pos[0], axis=0)
            pool['pline'] = np.delete(pool['pline'], pos[0], axis=0)
            del pool['prop'][pos[0]]
            pool['id'] = np.delete(pool['id'], pos[0], axis=0)
        else:
            pool["line"] = np.vstack((pool["line"], tline2))
            pool["pline"] = np.vstack((pool["pline"], tpline2))
            pool["prop"].append(tprop2)
            pool["id"] = np.vstack((pool["id"], tid2))

# Create new model
newnode = model["node"][np.squeeze(index) == 0, :]
tl = newnode.shape[0]
k = 0
nm = 0
for i in range(len(update["node"])):
    if update["node"][i].size != 0:
        nm += 1
        l = update["node"][i].shape[0]
        td = update["node"][i][:, 2:4]
        t = np.reshape(td, (l * 2, 1))
        t = np.unique(t)
        tt = update["node"][i][0, 1]
        for j in range(len(t)):
            k += 1
            tp = [t[j], model["node"][tt-1, 1], model["node"][tt-1, 2]]
            newnode = np.vstack((newnode, tp))

update["pnode"] = np.unique(update["pnode"])
for i in range(len(update["pnode"])):
    k += 1
    tp = [update["pnode"][i], model["node"][update["pnode"][i]-1, 1], model["node"][update["pnode"][i]-1, 2]]
    newnode = np.vstack((newnode, tp))

newele = model['element']
t = len(model["element"]["type"])
nn = 0
for i in range(t):
    tmp = newele["data"][i]
    for j in range(tmp.shape[0]):
        nn += 1
        tmp1 = tmp[j, :]
        for k in range(1, len(tmp1)):
            if index[int(tmp1[k])-1,0] == 1:
                newele["data"][i][j, k] = modu["element"][nn]["pvert"][k - 1]
print('Update element finish')

newele["type"].append('COH2D4P')
newele["data"].append(np.zeros((len(update["line"]), 7)))
kne = nel
ig = 0
for i in range(len(update["linetype"])):
    if update["linetype"][i] == 1:
        kne += 1
        ig += 1
        newele["data"][len(newele['type'])-1][ig-1, :] = [kne, update["line"][i][1][1], update["line"][i][1][0], update["line"][i][2][1], update["line"][i][2][0], update["line"][i][0][1], update["line"][i][0][0]]
gn1 = kne
for i in range(len(update["linetype"])):
    if update["linetype"][i] == 2:
        kne += 1
        ig += 1
        newele["data"][len(newele['type'])-1][ig-1, :] = [kne, update["line"][i][1][1], update["line"][i][1][0], update["line"][i][2][1], update["line"][i][2][0], update["line"][i][0][1], update["line"][i][0][0]]
gn2 = kne

# Node reorder
newmodel={}
newmodel['node'] = newnode
newmodel['element'] = newele
newmodel['group'] = model['group']
newmodel['group']['name'].append('rock')
newmodel['group']['data'].append(list(range(1, nel+1)))
newmodel['group']['name'].append('pjoint_inter_set')
newmodel['group']['data'].append(list(range(nel+1, gn1+1)))
newmodel['group']['name'].append('pjoint_adj_set')
newmodel['group']['data'].append(list(range(gn1+1, gn2+1)))
newmodel['group']['name'].append('pdfn')
newmodel['group']['data'].append(list(range(nel+1, gn2+1)))

# Add node
newmodel['nodegroup'] = model['nodegroup']
newmodel['nodegroup']['name'].append('pnode')
newmodel['nodegroup']['data'].append(list(range(tt+1, len(newnode)+1)))

indexnode = newnode[:, 0].astype(int)
for i in range(newmodel['node'].shape[0]):
    newmodel['node'][i, :] = [i+1, newnode[i, 1], newnode[i, 2]]


indexnode_dict = {v: i for i, v in enumerate(indexnode)}

for i in range(len(newmodel['element']['type'])):
    for j in range(newmodel['element']['data'][i].shape[0]):
        for k in range(1, newmodel['element']['data'][i].shape[1]):
            tmpnode = newmodel['element']['data'][i][j, k]
            newmodel['element']['data'][i][j, k] = indexnode_dict[int(tmpnode)] + 1

#Output inp
filename = CZMInp+'.inp'
with open(filename, 'w') as fid:
    fid.write('*Heading\n')
    fid.write('** Cohesive Test model for Jointed Rock\n')
    fid.write('** PARTS\n')
    fid.write('**\n')
    fid.write('*Part, name=Joint_rock\n')
    fid.write('*Node\n')

    for i in range(len(newmodel['node'])):
        if len(newmodel['node'][i]) == 3:
            fid.write('      %d,  %f,  %f\n' % (int(newmodel["node"][i][0]), newmodel["node"][i][1], newmodel["node"][i][2]))
        else:
            fid.write('      %d,  %f,  %f,  %f\n' % (int(newmodel["node"][i][0]), newmodel["node"][i][1], newmodel["node"][i][2], newmodel["node"][i][3]))

    for i in range(len(newmodel['element']['type'])):
        fid.write('*Element, type=%s\n' % newmodel["element"]["type"][i])
        for j in range(len(newmodel["element"]["data"][i])):
            data = newmodel["element"]["data"][i][j]
            fid.write('      ')
            for k in range(len(data)-1):
                fid.write('%d, ' % int(data[k]))
            fid.write('%d\n' % int(data[-1]))

    for i in range(len(newmodel['group']['name'])):
        grname = newmodel["group"]["name"][i]
        grdata = newmodel["group"]["data"][i]
        test = list(range(min(grdata), max(grdata)+1))
        if test == grdata:
            fid.write('*Elset, elset=%s, generate\n' % grname)
            fid.write('      %d,  %d,  1\n' % (int(min(grdata)), int(max(grdata))))
        else:
            fid.write('*Elset, elset=%s\n' % grname)
            if len(grdata) <= 10:
                j = 0
            else:
                for j in range(len(grdata)//10):
                    tmp = grdata[10*j:10*(j+1)]
                    fid.write('      %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, \n' % (int(tmp[0]), int(tmp[1]), int(tmp[2]), int(tmp[3]), int(tmp[4]), int(tmp[5]), int(tmp[6]),int(tmp[7]), int(tmp[8]), int(tmp[9])))

            if len(grdata) % 10 != 0:
                tmp = grdata[10*j:]
                str_ = 'fid.write("      %d' % int(tmp[0])
                for k in range(1, len(tmp)):
                    str_ += ', %d' % int(tmp[k])
                str_ += '\\n")'
                exec(str_)

    for i in range(len(newmodel["nodegroup"]["name"])):
        grname = newmodel["nodegroup"]["name"][i]
        grdata = newmodel["nodegroup"]["data"][i]
        test = list(range(int(min(grdata)), int(max(grdata))+1))
        if test == grdata:
            fid.write('*Nset, nset=%s, generate\n' % grname)
            fid.write('      %d, %d,  1\n' % (int(min(grdata)), int(max(grdata))))
        else:
            fid.write('*Nset, nset=%s\n' % grname)
            if len(grdata) < 10:
                j = 0
            else:
                for j in range(1, len(grdata)//10+1):
                    tmp = grdata[10*(j-1):10*j]
                    fid.write('      %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, \n' % (int(tmp[0]), int(tmp[1]), int(tmp[2]), int(tmp[3]), int(tmp[4]), int(tmp[5]), int(tmp[6]),int(tmp[7]), int(tmp[8]), int(tmp[9])))

            if len(grdata)%10 != 0:
                tmp = grdata[10*j:]
                str_ = 'fid.write("      %d' % int(tmp[0])
                for k in range(1, len(tmp)):
                    str_ += ', %d' % int(tmp[k])
                str_ += '\\n")'
                exec(str_)
    fid.write('*End Part\n' % ())

#Read inp and create new model and part
newmodelName = modelName+ '_czm'
newpartName = partName + '_czm'
model_name = mdb.ModelFromInputFile(inputFileName=filename,name=newmodelName)
mdb.models['Model-1_czm'].parts.changeKey(fromName='JOINT_ROCK', toName=newpartName)
p = mdb.models[newmodelName].parts[newpartName]
session.viewports['Viewport: 1'].setValues(displayedObject=p)
print('New model with CZM has generated')

p = mdb.models[newmodelName].parts[newpartName]
a=mdb.models[newmodelName].parts[newpartName]
a.Set(name='inject-node', objectToCopy=a.sets['PNODE'])
p = mdb.models[newmodelName].parts[newpartName]
nodes = p.nodes
max_node_label = -1
for node in nodes:
    x, y, z = node.coordinates
    if x == 0 and y == 0:
        if node.label > max_node_label:
            max_node_label = node.label
            max_node = node

n=p.nodes
node = max_node
idt=node.label
lab=tuple([idt])
nl=n.sequenceFromLabels(lab)
p.Set(nodes=nl, name='inject-node')
el=node.getElements()
xlist=[]
for ie in el:
    nodes=ie.getNodes()
    x=nodes[0].coordinates[0] + nodes[1].coordinates[0]
    xlist.append(x)

xlist2=xlist[:]
xlist2=[abs(i) for i in xlist2]
xlist3=[abs(i) for i in xlist2]
xlist_abs=xlist3
xlist_abs.sort()
id=0
idlist=[]
for i in xlist2:
    if i<xlist_abs[2]:
        idlist.append(id)
    id=id+1

labs=[]
for i in idlist:
    labs.append(el[i].label)

labs=tuple(labs)
e=p.elements
e2=e.sequenceFromLabels(labs)
p.Set(elements=e2, name='ini-elem')
elemType1 = mesh.ElemType(elemCode=COH2D4P, elemLibrary=STANDARD,
    viscosity=0.05)
p = mdb.models[newmodelName].parts[newpartName]
z1 = p.sets['PDFN'].elements
pickedRegions =(z1, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))

mdb.models[newmodelName].Material(name='Material-rock')
mdb.models[newmodelName].materials['Material-rock'].Elastic(table=((Rock_Elastic_Young_Modulus,
    Rock_Elastic_Poission), ))
mdb.models[newmodelName].materials['Material-rock'].Permeability(
    specificWeight=Rock_Permeability_Weight, inertialDragCoefficient=0.142887, table=((Rock_Permeability_k,
    Rock_Permeability_Void_Ratio), ))
mdb.models[newmodelName].Material(name='Material-cohe1')
mdb.models[newmodelName].materials['Material-cohe1'].Elastic(type=TRACTION,
    table=((Cohe1_Elastic_E, Cohe1_Elastic_G1, Cohe1_Elastic_G2), ))
mdb.models[newmodelName].materials['Material-cohe1'].MaxsDamageInitiation(table=((
    Cohe1_MaxsDamage_Nominal_Stress_Normal_only, Cohe1_MaxsDamage_Nominal_Stress_First_Direction, Cohe1_MaxsDamage_Nominal_Stress_Second_Direction), ))
mdb.models[newmodelName].materials['Material-cohe1'].maxsDamageInitiation.DamageEvolution(
    type=DISPLACEMENT, table=((Cohe1_MaxsDamage_DamageEvolution_Displacement, ), ))
mdb.models[newmodelName].materials['Material-cohe1'].FluidLeakoff(table=((Cohe1_FuildLeakoff,
    Cohe1_FuildLeakoff), ))
mdb.models[newmodelName].materials['Material-cohe1'].GapFlow(table=((Cohe1_GapFlow, ), ))
mdb.models[newmodelName].Material(name='Material-cohe2',
    objectToCopy=mdb.models[newmodelName].materials['Material-cohe1'])
mdb.models[newmodelName].materials['Material-cohe2'].maxsDamageInitiation.setValues(
    table=((Cohe2_MaxsDamage_Nominal_Stress_Normal_only, Cohe2_MaxsDamage_Nominal_Stress_First_Direction, Cohe2_MaxsDamage_Nominal_Stress_Second_Direction), ))

mdb.models[newmodelName].HomogeneousSolidSection(name='Section-rock',
    material='Material-rock', thickness=None)
mdb.models[newmodelName].CohesiveSection(name='Section-cohe1',
    material='Material-cohe1', response=TRACTION_SEPARATION,
    initialThicknessType=SPECIFY, initialThickness=0.001,
    outOfPlaneThickness=None)
mdb.models[newmodelName].CohesiveSection(name='Section-cohe2',
    material='Material-cohe2', response=TRACTION_SEPARATION,
    initialThicknessType=SPECIFY, initialThickness=0.001,
    outOfPlaneThickness=None)

p = mdb.models[newmodelName].parts[newpartName]
region = p.sets['ROCK']
p = mdb.models[newmodelName].parts[newpartName]
p.SectionAssignment(region=region, sectionName='Section-rock', offset=0.0,
    offsetType=MIDDLE_SURFACE, offsetField='',
    thicknessAssignment=FROM_SECTION)
p = mdb.models[newmodelName].parts[newpartName]
region = p.sets['PJOINT_INTER_SET']
p = mdb.models[newmodelName].parts[newpartName]
p.SectionAssignment(region=region, sectionName='Section-cohe1', offset=0.0,
    offsetType=MIDDLE_SURFACE, offsetField='',
    thicknessAssignment=FROM_SECTION)
p = mdb.models[newmodelName].parts[newpartName]
region = p.sets['PJOINT_ADJ_SET']
p = mdb.models[newmodelName].parts[newpartName]
p.SectionAssignment(region=region, sectionName='Section-cohe2', offset=0.0,
    offsetType=MIDDLE_SURFACE, offsetField='',
    thicknessAssignment=FROM_SECTION)

a = mdb.models[newmodelName].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models[newmodelName].parts[newpartName]
a.Instance(name=newpartName, part=p, dependent=ON)
mdb.models[newmodelName].GeostaticStep(name='Step-1', previous='Initial',
    maxNumInc=10000, timeIncrementationMethod=AUTOMATIC, minInc=1e-05,
    timePeriod=Step1_time, maxInc=0.1, utol=1e-05, matrixSolver=DIRECT, matrixStorage=UNSYMMETRIC,
    nlgeom=ON)
mdb.models[newmodelName].SoilsStep(name='Step-2', previous='Step-1',
    timePeriod=Step2_time, maxNumInc=10000, stabilizationMagnitude=0.0002,
    stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
    continueDampingFactors=False, adaptiveDampingRatio=0.05, initialInc=0.01,
    minInc=1e-08, maxInc=1.0, utol=100000000.0, cetol=None,
    matrixSolver=DIRECT, matrixStorage=UNSYMMETRIC)
mdb.models[newmodelName].steps['Step-1'].control.setValues(allowPropagation=OFF,
    resetDefaultValues=OFF, timeIncrementation=(4.0, 8.0, 9.0, 16.0, 10.0, 4.0,
    12.0, 10.0, 6.0, 3.0, 50.0))
mdb.models[newmodelName].FieldOutputRequest(name='F-Output-1',
    createStepName='Step-1', variables=(
    'S', 'LE', 'U', 'SDEG', 'FLVEL', 'POR'))
regionDef=mdb.models[newmodelName].rootAssembly.allInstances[newpartName].sets['inject-node']
mdb.models[newmodelName].HistoryOutputRequest(name='H-Output-1',
    createStepName='Step-1', variables=('POR', ), frequency=1,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)

mdb.models[newmodelName].TabularAmplitude(name='Amp-1', timeSpan=STEP,
    smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))
a = mdb.models[newmodelName].rootAssembly
region = a.instances[newpartName].sets['inject-node']
mdb.models[newmodelName].ConcPoreFluid(name='Load-1', createStepName='Step-2',
    region=region, magnitude=-0.001)

a = mdb.models[newmodelName].rootAssembly
n1 = a.instances[newpartName].nodes
fx1 = n1.getByBoundingBox(xMin=Model_xmin, xMax=Model_xmin, yMin=Model_ymin, yMax=Model_ymax, zMin=0, zMax=0)
fx2 = n1.getByBoundingBox(xMin=Model_xmax, xMax=Model_xmax, yMin=Model_ymin, yMax=Model_ymax, zMin=0, zMax=0)
region = a.Set(nodes=fx1+fx2, name='Set-x')
mdb.models[newmodelName].DisplacementBC(name='BC-x', createStepName='Step-1',
    region=region, u1=0.0, u2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF,
    distributionType=UNIFORM, fieldName='', localCsys=None)
a = mdb.models[newmodelName].rootAssembly
n1 = a.instances[newpartName].nodes
fy1 = n1.getByBoundingBox(xMin=Model_xmin, xMax=Model_xmax, yMin=Model_ymax, yMax=Model_ymax, zMin=0, zMax=0)
fy2 = n1.getByBoundingBox(xMin=Model_xmin, xMax=Model_xmax, yMin=Model_ymin, yMax=Model_ymin, zMin=0, zMax=0)
region = a.Set(nodes=fy1+fy2, name='Set-y')
mdb.models[newmodelName].DisplacementBC(name='BC-y', createStepName='Step-1',
    region=region, u1=UNSET, u2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF,
    distributionType=UNIFORM, fieldName='', localCsys=None)
a = mdb.models[newmodelName].rootAssembly
n1 = a.instances[newpartName].nodes
region = a.Set(nodes=fx1+fx2+fy1+fy2, name='Set-pp')
mdb.models[newmodelName].PorePressureBC(name='BC-pp', createStepName='Step-1',
    region=region, fixed=OFF, distributionType=UNIFORM, fieldName='',
    magnitude=0.0, amplitude=UNSET)
a = mdb.models[newmodelName].rootAssembly
e1 = a.instances[newpartName].elements
region = a.Set(elements=e1, name='Set-all')
mdb.models[newmodelName].Stress(name='Predefined Field-1', region=region,
    distributionType=UNIFORM, sigma11=-sxx, sigma22=-syy,
    sigma33=-szz, sigma12=0.0, sigma13=None, sigma23=None)

a = mdb.models[newmodelName].rootAssembly
n1 = a.instances[newpartName].nodes
region = a.Set(nodes=n1, name='Set-all-node')
mdb.models[newmodelName].VoidsRatio(name='Predefined Field-2', region=region,
    voidsRatio1=0.1, distributionType=UNIFORM, variation=CONSTANT_RATIO)
mdb.models[newmodelName].steps['Step-1'].control.setValues(discontinuous=ON)
mdb.models[newmodelName].keywordBlock.synchVersions()

def GetBlockPosition(modelName, blockPrefix):
    if blockPrefix == '':
        return len(mdb.models[modelName].keywordBlock.sieBlocks)-1
    pos = 0
    for block in mdb.models[modelName].keywordBlock.sieBlocks:
        if string.lower(block[0:len(blockPrefix)])==string.lower(blockPrefix):
            return pos
        pos=pos+1
    return -1

pos=GetBlockPosition(newmodelName, '*Step')-1
mdb.models[newmodelName].keywordBlock.insert(pos, """
    *initial conditions,type=initial gap
    %s.ini-elem""" % newpartName)

pos = 0
for block in mdb.models[newmodelName].keywordBlock.sieBlocks:
  if string.lower(block[0:len('*Element Output')])==string.lower('*Element Output'):
    #npos.append(pos)
    mdb.models[newmodelName].keywordBlock.replace(pos, """
      *Element Output, directions=YES
      FLVEL, LE, S, SDEG, PFOPEN""")
  pos=pos+1

print('Parameter setting completed')

del mdb.jobs[GroupInp]
mdb.models[newmodelName].keywordBlock.synchVersions()
mdb.Job(name=jobname, model=newmodelName, description='', type=ANALYSIS,
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=6,
    numDomains=6, numGPUs=0)

mdb.jobs[jobname].submit(consistencyChecking=OFF)
print('Job has subimitted')

