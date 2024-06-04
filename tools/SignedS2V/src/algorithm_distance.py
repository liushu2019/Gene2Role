# -*- coding: utf-8 -*-
from time import time
from collections import deque
import numpy as np
import math,logging
from fastdtw import fastdtw
from concurrent.futures import ProcessPoolExecutor, as_completed
from utils import *
import os
# import cmath

limiteDist = 20
maxA = 2.0 # for calculating the distance

def getDegreeListsVertices_complex(gp,gn,vertices,calcUntilLayer):
    degreeList = {}

    for v in vertices:
        degreeList[v] = getDegreeLists_complex(gp,gn,v,calcUntilLayer)

    return degreeList

def getCompactDegreeListsVertices_complex(Gp,Gn,vertices,maxDegree,calcUntilLayer):
    degreeList = {}

    logging.info('TEST getCompactDegreeListsVertices_complex in') # DEBUG
    for v in vertices:
        degreeList[v] = getCompactDegreeLists_complex(Gp,Gn,v,maxDegree,calcUntilLayer)

    return degreeList

def getCompactDegreeLists_complex(Gp,Gn, root, maxDegree,calcUntilLayer):
    t0 = time()

    listas = {}
    vetor_marcacao = [0] * (max(list(set(list(Gp.keys()) + list(Gn.keys())))) + 1)
    # maxNodeId = max(list(set(list(Gp.keys()) + list(Gn.keys()))))
    # vetor_marcacao = np.zeros((maxNodeId+1,maxNodeId+1))

    # Marcar s e inserir s na fila Q
    queue = deque()
    queue.append(root)
    vetor_marcacao[root] = 1
    l = {}

    ## Variáveis de controle de distância
    depth = 0
    pendingDepthIncrease = 0
    timeToDepthIncrease = 1

    logging.info('BFS vertex {}. in !'.format(root))
    while queue:
        vertex = queue.popleft()
        timeToDepthIncrease -= 1

        d = complex(len(Gp[vertex]), len(Gn[vertex]))
        if(d not in l):
            l[d] = 0
        l[d] += 1

        for v in Gp[vertex]:
            if(vetor_marcacao[v] == 0):
                vetor_marcacao[v] = 1
                queue.append(v)
                pendingDepthIncrease += 1   
        for v in Gn[vertex]:
            if(vetor_marcacao[v] == 0):
                vetor_marcacao[v] = 1
                queue.append(v)
                pendingDepthIncrease += 1  

        if(timeToDepthIncrease == 0):

            list_d = []
            for degree,freq in l.items():
                list_d.append((degree,freq))
            list_d.sort(key=lambda x: x[0].real+x[0].imag)
            listas[depth] = np.array(list_d)#,dtype=np.int32)

            l = {}

            if(calcUntilLayer == depth):
                break

            depth += 1
            timeToDepthIncrease = pendingDepthIncrease
            pendingDepthIncrease = 0


    t1 = time()
    logging.info('BFS vertex {}. Time: {}s'.format(root,(t1-t0)))

    return listas

def getDegreeLists_complex(gp, gn, root, calcUntilLayer):
    # print (gp, gn, max(gn))
    t0 = time()

    listas = {}
    
    list_gp = list(gp.keys()) if len(gp) > 0 else []
    list_gn = list(gn.keys()) if len(gn) > 0 else []
    vertices = list(set(list_gp + list_gn))
    vetor_marcacao = [0] * (max(vertices) + 1)
    # Marcar s e inserir s na fila Q
    queue = deque()
    queue.append(root)
    vetor_marcacao[root] = 1


    l = deque()

    ## Variáveis de controle de distância
    depth = 0
    pendingDepthIncrease = 0
    timeToDepthIncrease = 1

    while queue:
        vertex = queue.popleft()
        timeToDepthIncrease -= 1

        # l.append(len(g[vertex]))
        l.append(complex(len(gp[vertex]), len(gn[vertex])))

        for v in gp[vertex]:
            if(vetor_marcacao[v] == 0):
                vetor_marcacao[v] = 1
                queue.append(v)
                pendingDepthIncrease += 1    
        for v in gn[vertex]:
            if(vetor_marcacao[v] == 0):
                vetor_marcacao[v] = 1
                queue.append(v)
                pendingDepthIncrease += 1   

        if(timeToDepthIncrease == 0):

            # lp = np.array(l,dtype='float')
            # lp = np.sort(lp)
            lp = np.array(sorted(l, key= lambda x: (x.real + x.imag, x.real)))
            listas[depth] = lp
            l = deque()

            if(calcUntilLayer == depth):
                break

            depth += 1
            timeToDepthIncrease = pendingDepthIncrease
            pendingDepthIncrease = 0


    t1 = time()
    logging.info('BFS vertex {}. Time: {}s'.format(root,(t1-t0)))
    
    return listas

def cost(a,b):
    ep = 0.5
    m = max(a,b) + ep
    mi = min(a,b) + ep
    return ((m/mi) - 1)

def cost_max(a,b):
    ep = 0.5
    m = max(a[0],b[0]) + ep
    mi = min(a[0],b[0]) + ep
    return ((m/mi) - 1) * max(a[1],b[1])

def cost_complex(a,b):
    return abs(a*b.conjugate()-a.conjugate()*b)/4

def cost_complex_vector(a,b):
    return abs(a[0]*b[1]-a[1]*b[0])/2

def cost_complex_vector_2area_logscale(a,b):
    a_ = np.log(np.array(a)+1)
    b_ = np.log(np.array(b)+1)    
    return abs(a_[0]*b_[1]-a_[1]*b_[0]) + (np.linalg.norm(a_) - np.linalg.norm(b_))**2

def cost_complex_vector_Euclidean_logscale(a,b):
    a_ = np.log(np.array(a)+1)
    b_ = np.log(np.array(b)+1)    
    return np.linalg.norm(a_-b_)

def cost_complex_vector_Euclidean_logscale_max(a,b):
    a_ = np.log(np.array(a[:-1])+1)
    b_ = np.log(np.array(b[:-1])+1)    
    return np.linalg.norm(a_-b_)*max(a[-1],b[-1])

def cost_complex_vector_Euclidean(a,b):
    # a_ = np.log(np.array(a)+1)
    # b_ = np.log(np.array(b)+1)    
    return np.linalg.norm(a-b)

def cost_complex_vector_Euclidean_max(a,b):
    # a_ = np.log(np.array(a[0])+1)
    # b_ = np.log(np.array(b[0])+1)    
    # return np.linalg.norm(a_-b_)*max(a[1],b[1])
    # return np.linalg.norm(a[0]-b[0])*max(a[1][0],b[1][0])
    return np.linalg.norm(a[:-1]-b[:-1])*max(a[-1],b[-1])

def cost_complex_sinusoidalWave_logscale_vector(a,b): # complex numbers NG for fastDTW
    maxA = get_maxA()
    # print (maxA)
    if (a==b).all():
        return 0
    if maxA == 1:
        return 0
    a1 = np.log(np.sqrt(a[0]**2+a[1]**2)+1)
    a1 = (a1-1)/(maxA-1) + 1
    a2 = np.log(np.sqrt(b[0]**2+b[1]**2)+1)
    a2 = (a2-1)/(maxA-1) + 1
    # a1 = np.sqrt(np.log(a[0]+1)**2+np.log(a[1]+1)**2)
    # a2 = np.sqrt(np.log(b[0]+1)**2+np.log(b[1]+1)**2)
    # print ('000 111')
    alpha1 = math.atan2(a[1],a[0])
    alpha2 = math.atan2(b[1],b[0])
    # print ('000 112')
    p1 = np.arctan((a2*np.sin(alpha2)-a1*np.sin(alpha1))/(a1*np.cos(alpha1)-a2*np.cos(alpha2)))
    # print ('000 113')
    return abs(4*a2*np.cos(p1+alpha2)-4*a1*np.cos(p1+alpha1))

def cost_min(a,b):
    ep = 0.5
    m = max(a[0],b[0]) + ep
    mi = min(a[0],b[0]) + ep
    return ((m/mi) - 1) * min(a[1],b[1])


def verifyDegrees(degrees,degree_v_root,degree_a,degree_b):

    if(degree_b == -1):
        degree_now = degree_a
    elif(degree_a == -1):
        degree_now = degree_b
    elif(abs(degree_b - degree_v_root) < abs(degree_a - degree_v_root)):
        degree_now = degree_b
    else:
        degree_now = degree_a

    return degree_now 

def get_vertices(v,degree_v,degrees,a_vertices):
    a_vertices_selected = 2 * math.log(a_vertices,2)
    #logging.info("Selecionando {} próximos ao vértice {} ...".format(int(a_vertices_selected),v))
    vertices = deque()

    try:
        c_v = 0  

        for v2 in degrees[degree_v]['vertices']:
            if(v != v2):
                vertices.append(v2)
                c_v += 1
                if(c_v > a_vertices_selected):
                    raise StopIteration

        if('before' not in degrees[degree_v]):
            degree_b = -1
        else:
            degree_b = degrees[degree_v]['before']
        if('after' not in degrees[degree_v]):
            degree_a = -1
        else:
            degree_a = degrees[degree_v]['after']
        if(degree_b == -1 and degree_a == -1):
            raise StopIteration
        degree_now = verifyDegrees(degrees,degree_v,degree_a,degree_b)

        while True:
            for v2 in degrees[degree_now]['vertices']:
                if(v != v2):
                    vertices.append(v2)
                    c_v += 1
                    if(c_v > a_vertices_selected):
                        raise StopIteration

            if(degree_now == degree_b):
                if('before' not in degrees[degree_b]):
                    degree_b = -1
                else:
                    degree_b = degrees[degree_b]['before']
            else:
                if('after' not in degrees[degree_a]):
                    degree_a = -1
                else:
                    degree_a = degrees[degree_a]['after']

            if(degree_b == -1 and degree_a == -1):
                raise StopIteration

            degree_now = verifyDegrees(degrees,degree_v,degree_a,degree_b)

    except StopIteration:
        #logging.info("Vértice {} - próximos selecionados.".format(v))
        return list(vertices)

    return list(vertices)

def verifyDegrees_matrix(degrees,degree_v_root, degree_list_p, degree_list_n,degree_left_bottom, degree_right_top ):#,degree_a_p,degree_b_p,degree_a_n,degree_b_n):
    # print (degree_v_root, degree_left_bottom, degree_right_top)
    pass_flag_p = False
    pass_flag_n = False
    if (degree_list_p.index(degree_left_bottom.real) == 0) and (degree_list_p.index(degree_right_top.real)+1 == len(degree_list_p)) and (degree_list_n.index(degree_left_bottom.imag) == 0) and (degree_list_n.index(degree_right_top.imag)+1 == len(degree_list_n)):
        raise StopIteration
    if (degree_list_p.index(degree_left_bottom.real) == 0) and (degree_list_p.index(degree_right_top.real)+1 == len(degree_list_p)):
        pass_flag_p = True
        pass
    elif (degree_list_p.index(degree_left_bottom.real) == 0):
        degree_now_p = degree_list_p[degree_list_p.index(degree_right_top.real)+1]
    elif (degree_list_p.index(degree_right_top.real)+1 == len(degree_list_p)):
        degree_now_p = degree_list_p[degree_list_p.index(degree_left_bottom.real)-1]
    elif (abs(degree_list_p[degree_list_p.index(degree_right_top.real)+1] - degree_v_root.real) < abs(degree_list_p[degree_list_p.index(degree_left_bottom.real)-1] - degree_v_root.real)):
        degree_now_p = degree_list_p[degree_list_p.index(degree_right_top.real)+1]
    else:
        degree_now_p = degree_list_p[degree_list_p.index(degree_left_bottom.real)-1]
    
    if (degree_list_n.index(degree_left_bottom.imag) == 0) and (degree_list_n.index(degree_right_top.imag)+1 == len(degree_list_n)):
        pass_flag_n = True
        pass
    elif (degree_list_n.index(degree_left_bottom.imag) == 0):
        degree_now_n = degree_list_n[degree_list_n.index(degree_right_top.imag)+1]
    elif (degree_list_n.index(degree_right_top.imag)+1 == len(degree_list_n)):
        degree_now_n = degree_list_n[degree_list_n.index(degree_left_bottom.imag)-1]
    elif (abs(degree_list_n[degree_list_n.index(degree_right_top.imag)+1] - degree_v_root.imag) < abs(degree_list_n[degree_list_n.index(degree_left_bottom.imag)-1] - degree_v_root.imag)):
        degree_now_n = degree_list_n[degree_list_n.index(degree_right_top.imag)+1]
    else:
        degree_now_n = degree_list_n[degree_list_n.index(degree_left_bottom.imag)-1]
    # print (degree_now_n, degree_now_p)
    # if(degree_b_p.real == -1):
    #     degree_now_p = degree_a_p.real
    # elif(degree_a_p.real == -1):
    #     degree_now_p = degree_b_p.real
    # elif(abs(degree_b_p.real - degree_v_root.real) < abs(degree_a_p.real - degree_v_root.real)):
    #     degree_now_p = degree_b_p.real
    # else:
    #     degree_now_p = degree_a_p.real
    

    # if(degree_b_n.imag == -1):
    #     degree_now_n = degree_a_n.imag
    # elif(degree_a_n.imag == -1):
    #     degree_now_n = degree_b_n.imag
    # elif(abs(degree_b_n.imag - degree_v_root.imag) < abs(degree_a_n.imag - degree_v_root.imag)):
    #     degree_now_n = degree_b_n.imag
    # else:
    #     degree_now_n = degree_a_n.imag
    
    if ((not pass_flag_n) and (not pass_flag_p) and (abs(degree_now_n - degree_v_root.imag) < abs(degree_now_p - degree_v_root.real))) or (pass_flag_p):
        # print ('N')
        assert( (degree_now_n >= degree_right_top.imag) or (degree_now_n <= degree_left_bottom.imag) ), 'Search ERROR in verifyDegrees_matrix'
        degree_now = [complex(x, degree_now_n) for x in degree_list_p[degree_list_p.index(degree_left_bottom.real): degree_list_p.index(degree_right_top.real)+1 ]]
        if (degree_now_n > degree_right_top.imag):
            degree_right_top = complex(degree_right_top.real, degree_now_n)
        else:
            degree_left_bottom = complex(degree_left_bottom.real, degree_now_n)
    else:
        # print ('P',degree_right_top)
        assert( (degree_now_p >= degree_right_top.real) or (degree_now_p <= degree_left_bottom.real) ), 'Search ERROR in verifyDegrees_matrix'
        degree_now = [complex(degree_now_p, x) for x in degree_list_n[degree_list_n.index(degree_left_bottom.imag): degree_list_n.index(degree_right_top.imag)+1 ]]
        if (degree_now_p > degree_right_top.real):
            degree_right_top = complex(degree_now_p, degree_right_top.imag)
            # print ('P',degree_right_top, )
        else:
            degree_left_bottom = complex(degree_now_p, degree_left_bottom.imag)
    degree_now.sort(key=lambda x:abs(x - degree_v_root))

    # print (degree_now, degree_left_bottom, degree_right_top)
    return degree_now, degree_left_bottom, degree_right_top

def get_vertices_matrix(v,degree_v,degrees,a_vertices,degrees_sorted_n,degrees_sorted_p):
    '''
    degree_v: v's degree+ and v's degree- complex style
    degrees: matrix style
    a_vertices: # nodes
    '''
    a_vertices_selected = 2 * math.log(a_vertices,2)
    #logging.info("Selecionando {} próximos ao vértice {} ...".format(int(a_vertices_selected),v))
    vertices = deque()
    try:
        c_v = 0  

        for v2 in degrees[degree_v]['vertices']:
            if(v != v2):
                vertices.append(v2)
                c_v += 1
                if(c_v > a_vertices_selected):
                    raise StopIteration

        # if('before_p' not in degrees[degree_v]):
        #     degree_b_p = -1
        # else:
        #     degree_b_p = degrees[degree_v]['before_p']
        # if('after_p' not in degrees[degree_v]):
        #     degree_a_p = -1
        # else:
        #     degree_a_p = degrees[degree_v]['after_p']

        # if('before_n' not in degrees[degree_v]):
        #     degree_b_n = -1
        # else:
        #     degree_b_n = degrees[degree_v]['before_n']
        # if('after_n' not in degrees[degree_v]):
        #     degree_a_n = -1
        # else:
        #     degree_a_n = degrees[degree_v]['after_n']

        # if(degree_b_p == -1 and degree_a_p == -1 and degree_b_n == -1 and degree_a_n == -1):
        #     raise StopIteration
        degree_left_bottom = degree_v
        degree_right_top = degree_v
        degree_now_list, degree_left_bottom, degree_right_top = verifyDegrees_matrix(degrees,degree_v, degrees_sorted_p, degrees_sorted_n,degree_left_bottom, degree_right_top )#,degree_a_p,degree_b_p,degree_a_n,degree_b_n)

        while True:
            # print (v, c_v, a_vertices_selected, degree_now_list, degree_left_bottom, degree_right_top)
            for degree_now in degree_now_list:
                for v2 in degrees[degree_now]['vertices']:
                    if(v != v2):
                        vertices.append(v2)
                        c_v += 1
                        if(c_v > a_vertices_selected):
                            raise StopIteration
            degree_now_list, degree_left_bottom, degree_right_top = verifyDegrees_matrix(degrees,degree_v, degrees_sorted_p, degrees_sorted_n,degree_left_bottom, degree_right_top )#,degree_a_p,degree_b_p,degree_a_n,degree_b_n)

            # if(degree_now == degree_b):
            #     if('before' not in degrees[degree_b]):
            #         degree_b = -1
            #     else:
            #         degree_b = degrees[degree_b]['before']
            # else:
            #     if('after' not in degrees[degree_a]):
            #         degree_a = -1
            #     else:
            #         degree_a = degrees[degree_a]['after']

            # if(degree_b == -1 and degree_a == -1):
            #     raise StopIteration

            # degree_now = verifyDegrees(degrees,degree_v,degree_a,degree_b)

    except StopIteration:
        #logging.info("Vértice {} - próximos selecionados.".format(v))
        return list(vertices)

    return list(vertices)

def splitDegreeList(part,c,G,compactDegree):
    if(compactDegree):
        logging.info("Recovering compactDegreeList from disk...")
        degreeList = restoreVariableFromDisk('compactDegreeList')
    else:
        logging.info("Recovering degreeList from disk...")
        degreeList = restoreVariableFromDisk('degreeList')

    logging.info("Recovering degree vector from disk...")
    degrees = restoreVariableFromDisk('degrees_vector')

    degreeListsSelected = {}
    vertices = {}
    a_vertices = len(G)

    for v in c:
        nbs = get_vertices(v,len(G[v]),degrees,a_vertices)
        vertices[v] = nbs
        degreeListsSelected[v] = degreeList[v]
        for n in nbs:
            degreeListsSelected[n] = degreeList[n]

    saveVariableOnDisk(vertices,'split-vertices-'+str(part))
    saveVariableOnDisk(degreeListsSelected,'split-degreeList-'+str(part))


def splitDegreeList_complex(part,c,Gp,Gn,compactDegree):
    if(compactDegree):
        logging.info("Recovering compactDegreeList from disk...")
        degreeList = restoreVariableFromDisk('compactDegreeList')
    else:
        logging.info("Recovering degreeList from disk...")
        degreeList = restoreVariableFromDisk('degreeList')

    logging.info("Recovering degree vector from disk...")
    degrees = restoreVariableFromDisk('degrees_vector')
    degrees_sorted_n = restoreVariableFromDisk('degrees_vector_negativeList')
    degrees_sorted_p = restoreVariableFromDisk('degrees_vector_positiveList')

    degreeListsSelected = {}
    vertices = {}
    a_vertices = len(list(set(list(Gp.keys()) + list(Gn.keys()))))
    for v in c:
        nbs = get_vertices_matrix(v,complex(len(Gp[v]), len(Gn[v])),degrees,a_vertices,degrees_sorted_n,degrees_sorted_p) # TODO
        # print (v, nbs)
        vertices[v] = nbs
        degreeListsSelected[v] = degreeList[v]
        for n in nbs:
            degreeListsSelected[n] = degreeList[n]

    saveVariableOnDisk(vertices,'split-vertices-'+str(part))
    saveVariableOnDisk(degreeListsSelected,'split-degreeList-'+str(part))

def splitDegreeList_signed(part,c,G,compactDegree,type):
    if(compactDegree):
        logging.info("Recovering compactDegreeList from disk...")
        degreeList = restoreVariableFromDisk(type+'compactDegreeList')
    else:
        logging.info("Recovering degreeList from disk...")
        degreeList = restoreVariableFromDisk(type+'degreeList')

    logging.info("Recovering degree vector from disk...")
    degrees = restoreVariableFromDisk(type+'degrees_vector')

    degreeListsSelected = {}
    vertices = {}
    a_vertices = len(G)

    for v in c:
        nbs = get_vertices(v,len(G[v]),degrees,a_vertices)
        vertices[v] = nbs
        degreeListsSelected[v] = degreeList[v]
        for n in nbs:
            degreeListsSelected[n] = degreeList[n]

    saveVariableOnDisk(vertices,type+'split-vertices-'+str(part))
    saveVariableOnDisk(degreeListsSelected,type+'split-degreeList-'+str(part))


def calc_distances(part, compactDegree = False):

    vertices = restoreVariableFromDisk('split-vertices-'+str(part))
    degreeList = restoreVariableFromDisk('split-degreeList-'+str(part))

    distances = {}

    if compactDegree:
        dist_func = cost_max
    else:
        dist_func = cost

    for v1,nbs in vertices.items():
        lists_v1 = degreeList[v1]

        for v2 in nbs:
            t00 = time()
            lists_v2 = degreeList[v2]

            max_layer = min(len(lists_v1),len(lists_v2))
            distances[v1,v2] = {}

            for layer in range(0,max_layer):
                dist, path = fastdtw(lists_v1[layer],lists_v2[layer],radius=1,dist=dist_func)
# start dist_minus calculation

#  end  dist_minus calculation
                distances[v1,v2][layer] = dist

            t11 = time()
            logging.info('fastDTW between vertices ({}, {}). Time: {}s'.format(v1,v2,(t11-t00)))


    preprocess_consolides_distances(distances)
    saveVariableOnDisk(distances,'distances-'+str(part))
    return

def calc_distances_complex(part, compactDegree = False, scale_free = False):
    # scale_free = False
    # scale_free = True
    vertices = restoreVariableFromDisk('split-vertices-'+str(part))
    degreeList = restoreVariableFromDisk('split-degreeList-'+str(part))

    distances = {}

    if compactDegree:
        dist_func = cost_complex_vector_Euclidean_logscale_max if scale_free else cost_complex_vector_Euclidean_max
        for v1,nbs in vertices.items():
            lists_v1 = degreeList[v1]
            lists_v1_new = {}
            for key_,value_ in lists_v1.items():
                lists_v1_new.update({key_:[(x[0].real, x[0].imag, x[1].real) for x in value_]}) # TODO
            for v2 in nbs:
                lists_v2 = degreeList[v2]
                lists_v2_new = {}
                for key_,value_ in lists_v2.items():
                    lists_v2_new.update({key_:[(x[0].real, x[0].imag, x[1].real) for x in value_]}) # TODO

                max_layer = min(len(lists_v1),len(lists_v2))
                distances[v1,v2] = {}
                for layer in range(0,max_layer):
                    #t0 = time()
                    
                    # print ('101 111')
                    # print (layer,lists_v1_new[layer],lists_v2_new[layer])
                    dist, path = fastdtw(lists_v1_new[layer],lists_v2_new[layer],radius=1,dist=dist_func)
                    # print ('100 111')
                    # dist, path = fastdtw(lists_v1[layer],lists_v2[layer],radius=1,dist=dist_func)
    # start dist_minus calculation

    #  end  dist_minus calculation
                    #t1 = time()
                    #logging.info('D ({} , {}), Tempo fastDTW da camada {} : {}s . Distância: {}'.format(v1,v2,layer,(t1-t0),dist))    
                    distances[v1,v2][layer] = np.exp(dist) # Link to # edges might be better!!! TODO
                    # distances[v1,v2][layer] = dist 
                    # distances[v1,v2][layer] = np.exp(dist)**np.e
    else:
        dist_func = cost_complex_vector_Euclidean_logscale if scale_free else cost_complex_vector_Euclidean
        # dist_func = cost_complex_vector_2area_logscale
        # dist_func = cost_complex_sinusoidalWave_logscale_vector
        for v1,nbs in vertices.items():
            lists_v1 = degreeList[v1]
            lists_v1_new = {}
            for key_,value_ in lists_v1.items():
                lists_v1_new.update({key_:[(x.real, x.imag) for x in value_]})
            for v2 in nbs:
                lists_v2 = degreeList[v2]
                lists_v2_new = {}
                for key_,value_ in lists_v2.items():
                    lists_v2_new.update({key_:[(x.real, x.imag) for x in value_]})

                max_layer = min(len(lists_v1),len(lists_v2))
                distances[v1,v2] = {}
                for layer in range(0,max_layer):
                    #t0 = time()
                    
                    # print ('101 111')
                    # print (layer,lists_v1_new[layer],lists_v2_new[layer])
                    dist, path = fastdtw(lists_v1_new[layer],lists_v2_new[layer],radius=1,dist=dist_func)
                    # print ('100 111')
                    # dist, path = fastdtw(lists_v1[layer],lists_v2[layer],radius=1,dist=dist_func)
    # start dist_minus calculation

    #  end  dist_minus calculation
                    #t1 = time()
                    #logging.info('D ({} , {}), Tempo fastDTW da camada {} : {}s . Distância: {}'.format(v1,v2,layer,(t1-t0),dist))    
                    # distances[v1,v2][layer] = np.exp(dist) # Link to # edges might be better!!! TODO
                    distances[v1,v2][layer] = dist 
                    # distances[v1,v2][layer] = np.exp(dist)**np.e

#     if compactDegree:
#         dist_func = cost_max
#     else:
#         dist_func = cost

#     for v1,nbs in vertices.items():
#         lists_v1 = degreeList[v1]

#         for v2 in nbs:
#             t00 = time()
#             lists_v2 = degreeList[v2]

#             max_layer = min(len(lists_v1),len(lists_v2))
#             distances[v1,v2] = {}

#             for layer in range(0,max_layer):
#                 dist, path = fastdtw(lists_v1[layer],lists_v2[layer],radius=1,dist=dist_func)
# # start dist_minus calculation

# #  end  dist_minus calculation
#                 distances[v1,v2][layer] = dist

            # t11 = time()
            # logging.info('fastDTW between vertices ({}, {}). Time: {}s'.format(v1,v2,(t11-t00)))


    preprocess_consolides_distances(distances)
    saveVariableOnDisk(distances,'distances-'+str(part))
    return

def calc_distances_signed(part,type, compactDegree = False):

    vertices = restoreVariableFromDisk(type+'split-vertices-'+str(part))
    degreeList = restoreVariableFromDisk(type+'split-degreeList-'+str(part))

    distances = {}

    if compactDegree:
        dist_func = cost_max
    else:
        dist_func = cost

    for v1,nbs in vertices.items():
        lists_v1 = degreeList[v1]

        for v2 in nbs:
            t00 = time()
            lists_v2 = degreeList[v2]

            max_layer = min(len(lists_v1),len(lists_v2))
            distances[v1,v2] = {}

            for layer in range(0,max_layer):
                dist, path = fastdtw(lists_v1[layer],lists_v2[layer],radius=1,dist=dist_func)
# start dist_minus calculation

#  end  dist_minus calculation
                distances[v1,v2][layer] = dist

            t11 = time()
            logging.info('fastDTW between vertices ({}, {}). Time: {}s'.format(v1,v2,(t11-t00)))


    preprocess_consolides_distances(distances)
    saveVariableOnDisk(distances,type+'distances-'+str(part))
    return

def calc_distances_all(vertices,list_vertices,degreeList,part, compactDegree = False):

    distances = {}
    cont = 0

    if compactDegree:
        dist_func = cost_max
    else:
        dist_func = cost

    for v1 in vertices:
        lists_v1 = degreeList[v1]

        for v2 in list_vertices[cont]:
            lists_v2 = degreeList[v2]

            max_layer = min(len(lists_v1),len(lists_v2))
            distances[v1,v2] = {}

            for layer in range(0,max_layer):
                #t0 = time()
                dist, path = fastdtw(lists_v1[layer],lists_v2[layer],radius=1,dist=dist_func)
# start dist_minus calculation

#  end  dist_minus calculation
                #t1 = time()
                #logging.info('D ({} , {}), Tempo fastDTW da camada {} : {}s . Distância: {}'.format(v1,v2,layer,(t1-t0),dist))    
                distances[v1,v2][layer] = dist


        cont += 1

    preprocess_consolides_distances(distances)
    saveVariableOnDisk(distances,'distances-'+str(part))
    return

def get_maxA():
    global maxA
    return maxA

def set_maxA(value):
    global maxA
    maxA = value
    return

def calc_distances_all_complex(vertices,list_vertices,degreeList,part, compactDegree = False, scale_free = False):
    # scale_free = False
    # scale_free = True
    distances = {}
    cont = 0

    # TODO
    if compactDegree:
        dist_func = cost_complex_vector_Euclidean_logscale_max if scale_free else cost_complex_vector_Euclidean_max
        for v1 in vertices:
            lists_v1 = degreeList[v1]
            lists_v1_new = {}
            for key_,value_ in lists_v1.items():
                lists_v1_new.update({key_:[(x[0].real, x[0].imag, x[1].real) for x in value_]})
            for v2 in list_vertices[cont]:
                lists_v2 = degreeList[v2]
                lists_v2_new = {}
                for key_,value_ in lists_v2.items():
                    lists_v2_new.update({key_:[(x[0].real, x[0].imag, x[1].real) for x in value_]})

                max_layer = min(len(lists_v1),len(lists_v2))
                distances[v1,v2] = {}
                for layer in range(0,max_layer):
                    dist, path = fastdtw(lists_v1_new[layer],lists_v2_new[layer],radius=1,dist=dist_func)
                    distances[v1,v2][layer] = np.exp(dist)

            cont += 1

    else:
        dist_func = cost_complex_vector_Euclidean_logscale if scale_free else cost_complex_vector_Euclidean
        # dist_func = cost_complex_vector_2area_logscale
        # dist_func = cost_complex_sinusoidalWave_logscale_vector
        for v1 in vertices:
            lists_v1 = degreeList[v1]
            lists_v1_new = {}
            for key_,value_ in lists_v1.items():
                lists_v1_new.update({key_:[(x.real, x.imag) for x in value_]})
            for v2 in list_vertices[cont]:
                lists_v2 = degreeList[v2]
                lists_v2_new = {}
                for key_,value_ in lists_v2.items():
                    lists_v2_new.update({key_:[(x.real, x.imag) for x in value_]})

                max_layer = min(len(lists_v1),len(lists_v2))
                distances[v1,v2] = {}
                for layer in range(0,max_layer):
                    #t0 = time()
                    
                    # print ('101 111')
                    # print (layer,lists_v1_new[layer],lists_v2_new[layer])
                    dist, path = fastdtw(lists_v1_new[layer],lists_v2_new[layer],radius=1,dist=dist_func)
                    # print ('100 111')
                    # dist, path = fastdtw(lists_v1[layer],lists_v2[layer],radius=1,dist=dist_func)
    # start dist_minus calculation

    #  end  dist_minus calculation
                    #t1 = time()
                    #logging.info('D ({} , {}), Tempo fastDTW da camada {} : {}s . Distância: {}'.format(v1,v2,layer,(t1-t0),dist))    
                    # distances[v1,v2][layer] = np.exp(dist) # Link to # edges might be better!!! TODO
                    distances[v1,v2][layer] = np.exp(dist)
                    # distances[v1,v2][layer] = np.exp(dist)**np.e


            cont += 1

    preprocess_consolides_distances(distances)
    saveVariableOnDisk(distances,'distances-'+str(part))
    return

def selectVertices(layer,fractionCalcDists):
    previousLayer = layer - 1

    logging.info("Recovering distances from disk...")
    distances = restoreVariableFromDisk('distances')

    threshold = calcThresholdDistance(previousLayer,distances,fractionCalcDists)

    logging.info('Selecting vertices...')

    vertices_selected = deque()

    for vertices,layers in distances.items():
        if(previousLayer not in layers):
            continue
        if(layers[previousLayer] <= threshold):
            vertices_selected.append(vertices)

    distances = {}

    logging.info('Vertices selected.')

    return vertices_selected


def preprocess_consolides_distances(distances, startLayer = 1):

    logging.info('Consolidating distances...')

    for vertices,layers in distances.items():
        keys_layers = sorted(list(layers.keys()))
        startLayer = min(len(keys_layers),startLayer)
        for layer in range(0,startLayer):
            keys_layers.pop(0)


        for layer in keys_layers:
            layers[layer] += layers[layer - 1]

    logging.info('Distances consolidated.')

def exec_bfs_compact_complex(Gp,Gn,workers,calcUntilLayer):

    futures = {}
    degreeList = {}

    t0 = time()
    vertices = list(set(list(Gp.keys()) + list(Gn.keys())))
    parts = workers
    chunks = partition(vertices,parts)

    logging.info('Capturing larger degree...')
    maxDegree = 0
    for v in vertices:
        if(len(Gp[v])+len(Gn[v]) > maxDegree):
            maxDegree = len(Gp[v])+len(Gn[v])
    logging.info('Larger degree captured')

    # # print (chunks) #DEBUG
    # part = 1
    # for c in chunks:
    #     dl = getCompactDegreeListsVertices_complex(Gp,Gn,c,maxDegree,calcUntilLayer)
    #     degreeList.update(dl)
    # # print (degreeList) #DEBUG

    with ProcessPoolExecutor(max_workers=workers) as executor:

        part = 1
        for c in chunks:
            job = executor.submit(getCompactDegreeListsVertices_complex,Gp,Gn,c,maxDegree,calcUntilLayer)
            futures[job] = part
            part += 1

        for job in as_completed(futures):
            dl = job.result()
            v = futures[job]
            degreeList.update(dl)

    logging.info("Saving degreeList on disk...")
    saveVariableOnDisk(degreeList,'compactDegreeList')
    t1 = time()
    logging.info('Execution time - BFS: {}m'.format((t1-t0)/60))


    return

def exec_bfs_complex(Gp,Gn,workers,calcUntilLayer):
    futures = {}
    degreeList = {}

    t0 = time()
    list_gp = list(Gp.keys()) if len(Gp) > 0 else []
    list_gn = list(Gn.keys()) if len(Gn) > 0 else []
    vertices = list(set(list_gp + list_gn))
    parts = workers
    chunks = partition(vertices,parts)
    # for c in chunks:
    #     dl = getDegreeListsVertices_complex(Gp,Gn,c,calcUntilLayer)
    #     degreeList.update(dl)
    with ProcessPoolExecutor(max_workers=workers) as executor:

        part = 1
        for c in chunks:
            job = executor.submit(getDegreeListsVertices_complex,Gp,Gn,c,calcUntilLayer)
            futures[job] = part
            part += 1

        for job in as_completed(futures):
            dl = job.result()
            v = futures[job]
            degreeList.update(dl)

    logging.info("Saving degreeList on disk...")
    saveVariableOnDisk(degreeList,'degreeList')
    t1 = time()
    logging.info('Execution time - BFS: {}m'.format((t1-t0)/60))


    return

def generate_distances_network_part1(workers):
    parts = workers
    weights_distances = {}
    for part in range(1,parts + 1):    

        logging.info('Executing part {}...'.format(part))
        distances = restoreVariableFromDisk('distances-'+str(part))
        for vertices,layers in distances.items():
            for layer,distance in layers.items():
                vx = vertices[0]
                vy = vertices[1]
                if(layer not in weights_distances):
                    weights_distances[layer] = {}
                weights_distances[layer][vx,vy] = distance

        logging.info('Part {} executed.'.format(part))

    for layer,values in weights_distances.items():
        saveVariableOnDisk(values,'weights_distances-layer-'+str(layer))
    return

def generate_distances_network_part2(workers):
    parts = workers
    graphs = {}
    for part in range(1,parts + 1):

        logging.info('Executing part {}...'.format(part))
        distances = restoreVariableFromDisk('distances-'+str(part))

        for vertices,layers in distances.items():
            for layer,distance in layers.items():
                vx = vertices[0]
                vy = vertices[1]
                if(layer not in graphs):
                    graphs[layer] = {}
                if(vx not in graphs[layer]):
                   graphs[layer][vx] = [] 
                if(vy not in graphs[layer]):
                   graphs[layer][vy] = [] 
                graphs[layer][vx].append(vy)
                graphs[layer][vy].append(vx)
        logging.info('Part {} executed.'.format(part))

    for layer,values in graphs.items():
        saveVariableOnDisk(values,'graphs-layer-'+str(layer))

    return

def generate_distances_network_part3():

    layer = 0
    while(isPickle('graphs-layer-'+str(layer))):
        graphs = restoreVariableFromDisk('graphs-layer-'+str(layer))
        weights_distances = restoreVariableFromDisk('weights_distances-layer-'+str(layer))

        logging.info('Executing layer {}...'.format(layer))
        alias_method_j = {}
        alias_method_q = {}
        weights = {}

        for v,neighbors in graphs.items():
            e_list = deque()
            sum_w = 0.0


            for n in neighbors:
                if (v,n) in weights_distances:
                    wd = weights_distances[v,n]
                else:
                    wd = weights_distances[n,v]
                w = np.exp(-float(wd))
                e_list.append(w)
                sum_w += w

            e_list = [x / sum_w for x in e_list]
            weights[v] = e_list
            J, q = alias_setup(e_list)
            alias_method_j[v] = J
            alias_method_q[v] = q

        saveVariableOnDisk(weights,'distances_nets_weights-layer-'+str(layer))
        saveVariableOnDisk(alias_method_j,'alias_method_j-layer-'+str(layer))
        saveVariableOnDisk(alias_method_q,'alias_method_q-layer-'+str(layer))
        logging.info('Layer {} executed.'.format(layer))
        layer += 1

    logging.info('Weights created.')

    return


def generate_distances_network_part4():
    logging.info('Consolidating graphs...')
    graphs_c = {}
    layer = 0
    while(isPickle('graphs-layer-'+str(layer))):
        logging.info('Executing layer {}...'.format(layer))
        graphs = restoreVariableFromDisk('graphs-layer-'+str(layer))
        graphs_c[layer] = graphs
        logging.info('Layer {} executed.'.format(layer))
        layer += 1


    logging.info("Saving distancesNets on disk...")
    saveVariableOnDisk(graphs_c,'distances_nets_graphs')
    logging.info('Graphs consolidated.')
    return

def generate_distances_network_part5():
    alias_method_j_c = {}
    layer = 0
    while(isPickle('alias_method_j-layer-'+str(layer))):
        logging.info('Executing layer {}...'.format(layer))          
        alias_method_j = restoreVariableFromDisk('alias_method_j-layer-'+str(layer))
        alias_method_j_c[layer] = alias_method_j
        logging.info('Layer {} executed.'.format(layer))
        layer += 1

    logging.info("Saving nets_weights_alias_method_j on disk...")
    saveVariableOnDisk(alias_method_j_c,'nets_weights_alias_method_j')

    return

def generate_distances_network_part6():
    alias_method_q_c = {}
    layer = 0
    while(isPickle('alias_method_q-layer-'+str(layer))):
        logging.info('Executing layer {}...'.format(layer))          
        alias_method_q = restoreVariableFromDisk('alias_method_q-layer-'+str(layer))
        alias_method_q_c[layer] = alias_method_q
        logging.info('Layer {} executed.'.format(layer))
        layer += 1

    logging.info("Saving nets_weights_alias_method_q on disk...")
    saveVariableOnDisk(alias_method_q_c,'nets_weights_alias_method_q')

    return

def generate_distances_network(workers):
    t0 = time()
    logging.info('Creating distance network...')

    os.system("rm "+returnPathSignedS2V()+"/../pickles/weights_distances-layer-*.pickle")
    with ProcessPoolExecutor(max_workers=1) as executor:
        job = executor.submit(generate_distances_network_part1,workers)
        job.result()
    t1 = time()
    t = t1-t0
    logging.info('- Time - part 1: {}s'.format(t))

    t0 = time()
    os.system("rm "+returnPathSignedS2V()+"/../pickles/graphs-layer-*.pickle")
    with ProcessPoolExecutor(max_workers=1) as executor:
        job = executor.submit(generate_distances_network_part2,workers)
        job.result()
    t1 = time()
    t = t1-t0
    logging.info('- Time - part 2: {}s'.format(t))
    logging.info('distance network created.')

    logging.info('Transforming distances into weights...')

    t0 = time()
    os.system("rm "+returnPathSignedS2V()+"/../pickles/distances_nets_weights-layer-*.pickle")
    os.system("rm "+returnPathSignedS2V()+"/../pickles/alias_method_j-layer-*.pickle")
    os.system("rm "+returnPathSignedS2V()+"/../pickles/alias_method_q-layer-*.pickle")
    with ProcessPoolExecutor(max_workers=1) as executor:
        job = executor.submit(generate_distances_network_part3)
        job.result()
    t1 = time()
    t = t1-t0
    logging.info('- Time - part 3: {}s'.format(t))

    t0 = time()
    with ProcessPoolExecutor(max_workers=1) as executor:
        job = executor.submit(generate_distances_network_part4)
        job.result()
    t1 = time()
    t = t1-t0
    logging.info('- Time - part 4: {}s'.format(t))

    t0 = time()
    with ProcessPoolExecutor(max_workers=1) as executor:
        job = executor.submit(generate_distances_network_part5)
        job.result()
    t1 = time()
    t = t1-t0
    logging.info('- Time - part 5: {}s'.format(t))

    t0 = time()
    with ProcessPoolExecutor(max_workers=1) as executor:
        job = executor.submit(generate_distances_network_part6)
        job.result()
    t1 = time()
    t = t1-t0
    logging.info('- Time - part 6: {}s'.format(t))

    return

def alias_setup(probs):
    '''
    Compute utility lists for non-uniform sampling from discrete distributions.
    Refer to https://hips.seas.harvard.edu/blog/2013/03/03/the-alias-method-efficient-sampling-with-many-discrete-outcomes/
    for details
    '''
    K = len(probs)
    q = np.zeros(K)
    J = np.zeros(K, dtype=int)

    smaller = []
    larger = []
    for kk, prob in enumerate(probs):
        q[kk] = K*prob
        if q[kk] < 1.0:
            smaller.append(kk)
        else:
            larger.append(kk)

    while len(smaller) > 0 and len(larger) > 0:
        small = smaller.pop()
        large = larger.pop()

        J[small] = large
        q[large] = q[large] + q[small] - 1.0
        if q[large] < 1.0:
            smaller.append(large)
        else:
            larger.append(large)

    return J, q

def merge_distance_by_sech(height):
    '''
    $$f(x,y)=(x+y)^2+ h*\frac{cosh^2(\sqrt{2}|\frac{x-y}{2}|) - sinh^2(\sqrt{2}|\frac{x-y}{2}|))}{cosh^2(\sqrt{2}|\frac{x-y}{2}|)}$$
    '''
    return lambda x,y: (x+y+1)**2 + height*(np.cosh(np.sqrt(2)*abs(x-y)/2)**2 - np.sinh(np.sqrt(2)*abs(x-y)/2)**2)/np.cosh(np.sqrt(2)*abs(x-y)/2)**2

def mylambda():
    return []