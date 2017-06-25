# Subject to BSD-3-Clause
# Copyright (c) 2016-2017 Filippo Ranalli, Forest Flager, Pratyush Havelia
# See LICENSE.txt and NOTICE.txt

import os
import sys
import shutil
import numpy as np
import pandas as pd
import comtypes.client
from SAP_wrapper_class import SAPAnalysis
from sizing_class import size
from detailing_class import detail
from continuity_class import prioritize
from coster_class import cost
from optimizer_class import *
import time

def initializeData(model_name,overwrite_files, overwrite_model, load_checkpoint, checkpoint_iteration):
    if not load_checkpoint:
        print ('\n No checkpoint loaded')
        path = './starter_files/' + model_name + '/'
        files_src = (path)
        files_dst = ("./")
        model_path = path + '/model/'
        model_src = (model_path)
        model_dst = ("./SAPMODEL/")
    else:
        print ('\n loading checkpoint # ', str(checkpoint_iteration))
        path = './checkpoints/' + str(checkpoint_iteration) + '/'
        files_src = (path)
        files_dst = ("./")
        model_path = path
        model_src = (model_path)
        model_dst = ("./SAPMODEL/")

    if overwrite_files:
        for filename in os.listdir(files_src):
            if filename.endswith('.txt'):
                shutil.copy( files_src + filename, files_dst)
                print('initializing file: ',filename)

    if overwrite_model:
        for filename in os.listdir(model_src):
            if filename.endswith('.sdb'):
                shutil.copy( model_src + filename, model_dst)
                print('initializing file: ',filename)

    sapImember = pd.read_csv("SAP_I_member.txt")
    geomData = pd.read_csv("member_geometry.txt")
    sapInode = pd.read_csv("SAP_I_Node.txt")
    nodeData = pd.read_csv("node_geometry.txt")
    parameters = pd.read_csv("inputs.txt")
    sec_info_W = pd.read_csv('W.txt',sep=',',encoding='utf-8')
    sec_info_Round = pd.read_csv('Round.txt',sep=',',encoding='utf-8')
    sec_info_Square = pd.read_csv('Square.txt',sep=',',encoding='utf-8')
    all_sections = pd.read_csv('AllFeaturedSections.txt',sep=',')
    return sapImember,geomData,sapInode,nodeData,parameters,sec_info_W,sec_info_Round,sec_info_Square,all_sections



(overwrite_model, overwrite_files, load_checkpoint, checkpoint_iteration) = (True, True, False, 1)
sapImember,geomData,sapInode,nodeData,parameters,sec_info_W,sec_info_Round,sec_info_Square,all_sections = initializeData('3x6 truss',overwrite_files,overwrite_model,load_checkpoint,checkpoint_iteration)
lastRemoved,no_remove_list,cost_history,total_reactions = initializeOptimizerVars()

analysis = SAPAnalysis('./SAPMODEL/optimization.sdb',logging=False)
analysis.initialize()

mdo_converged = False
while not mdo_converged:
    frequency,virtual_forces,mem_force,node_forces,node_displacements = analysis.run(sapImember,sapInode,sizingLoop=False)
    nodeToElementMap = detail(sapImember,sapInode)
    continuity = prioritize(nodeToElementMap,parameters)
    sizing_converged = False
    while not sizing_converged:
        sapImember,continuity,member_cost,proposed_design,memberStiffnessContribution,sizing_converged,parameters = size(continuity,virtual_forces,sapImember,sapInode,geomData,node_displacements,mem_force,frequency,parameters,sec_info_W,sec_info_Round,sec_info_Square, logging=False)
        frequency,virtual_forces,mem_force,node_forces,node_displacements = analysis.run(sapImember,sapInode,sizingLoop=True)
    if parameters['Values'][13]:
        totalCosts = cost(nodeToElementMap,continuity,sapImember,member_cost,connection_ratio=0.0,logging=False)
    else:
        totalCosts = -1
    lastRemoved,no_remove_list,cost_history,total_reactions,parameters,sapImember,sapINode,mdo_converged = optimize(lastRemoved,no_remove_list,cost_history,total_reactions,parameters,sapImember,sapInode,proposed_design,memberStiffnessContribution,totalCosts,mem_force,node_forces,nodeToElementMap,frequency,all_sections)