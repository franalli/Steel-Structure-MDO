# Subject to BSD-3-Clause
# Copyright (c) 2016-2017 Filippo Ranalli, Forest Flager, Pratyush Havelia
# See LICENSE.txt and NOTICE.txt

import numpy as np
import pandas as pd
import subprocess
import sys
import os
import shutil
from continuity_class import prioritize

def initializeOptimizerVars():
    lastRemoved = pd.read_csv("last_removed.txt")
    no_remove_list = pd.read_csv("no_remove_list.txt")
    cost_history = pd.read_csv("cost_history.txt")
    total_reactions = pd.read_csv("Total_Reactions.txt")
    return lastRemoved,no_remove_list,cost_history,total_reactions


def optimize(lastRemoved,no_remove_list,cost_history,total_reactions,parameters,sapIMember,sapINode,memberSizes,memberStiffnessContribution,memberCosts,memberForces,nodeReactions,nodeToMemberMap,modeFrequencies,sec_info, logging=False, checkpoint=True):
    print('OPTIMIZING')

    def findLoadedNodes():
        loadedNodes = []
        nodeIds = sapINode["node_ID"]
        loadBearing = sapINode["load_bearing?"]
        for i in range(0, len(nodeIds)):
            if bool(loadBearing[i]) == True:
                loadedNodes.append(nodeIds[i])
        return loadedNodes

    def rankAndRemoveMembers(sapIMember,memberCosts,memberSizes,memberStiffnessContribution,memberForces,no_remove_list,removalMultiplier,parameters):
        Fcs = []
        # find Fcs and order by descending values
        memberIDs = memberStiffnessContribution['member_ID']
        memberStiffnessContribution = memberStiffnessContribution.sort_values(by = "member_ID", ascending = True).reset_index(drop = True)
        for i in range(0, len(memberIDs)):
            Fcs.append( memberCosts['cost'][i]/ memberStiffnessContribution.get_value(i, "VW") )  # cost based optimization
        memberForces = memberForces.sort_values(by = "member_ID", ascending = True)     # because all the outputs from sizing are ordered alphabetically
        memberForces = memberForces.reset_index(drop = True)

        memberFcs = memberStiffnessContribution
        memberFcs["cost"] = list(memberCosts["cost"])
        memberFcs["member_Fcs"] = Fcs
        memberFcs["cross_section"] = list(memberSizes["cross_section"])
        memberFcs["start_node"] = memberForces["Start_Node"]
        memberFcs["end_node"] = memberForces["End_Node"]

        memberFcsSorted = memberFcs.sort_values(by = ['member_Fcs','member_ID'], ascending = False)   # sorted by values of Fc
        memberFcsSorted = memberFcsSorted.reset_index(drop = True)
        membersRemaining = memberFcsSorted

        # do not try to remove members that are already on the do-not-remove list
        for no_remove_member in no_remove_list['member_ID']:
            membersRemaining = membersRemaining[membersRemaining['member_ID'] != no_remove_member]
        membersRemaining = membersRemaining.reset_index(drop = True)
        membersRemoved = membersRemaining


        # find x% members to remove and update input file
        threshold = parameters.get_value(1, "Values")
        percent_to_remove = int(parameters.get_value(2, "Values"))
        iteration = int(parameters.get_value(3, "Values"))
        decrease_rate = 0

        parameters['Values'][2] = int(percent_to_remove)
        numberToRemove = int(np.ceil(percent_to_remove * len(membersRemaining["member_ID"]) / 100.0))
        desiredNumberToRemove = int(numberToRemove)
        numberToRemove = min(int(np.ceil(percent_to_remove * len(membersRemaining["member_ID"]) / 100.0) * removalMultiplier),len(membersRemaining["member_ID"]))

        if numberToRemove % 2 != 0:    # ensure the number of elements removed is even
            numberToRemove += 1

        # Find the zero-force members
        zero_member_count = 0
        zero_member_list = []

        ''' P threshold'''
        checklist = no_remove_list['member_ID'].tolist()

        VW_threshold = 0.0001
        for i in range(0,len(memberStiffnessContribution['member_ID'])):
            member = memberStiffnessContribution.get_value(i,'member_ID')
            if (memberStiffnessContribution[memberStiffnessContribution['member_ID']==member].iloc[0]['VW'] <= VW_threshold or memberForces[memberForces['member_ID']==member].iloc[0]['Governing_Load_Case'] == 'ZERO_FORCE') and member not in checklist:
                zero_member_count = zero_member_count + 1
                zero_member_list.append(member)

        # remove all the zero force members first
        for zero_VW_member in zero_member_list:
            membersRemaining = membersRemaining[membersRemaining['member_ID'] != zero_VW_member]
            membersRemaining = membersRemaining.reset_index(drop = True)


        for i in range(0,len(membersRemaining['member_ID'])):
            member_remaining = membersRemaining.get_value(i,'member_ID')
            membersRemoved = membersRemoved[membersRemoved['member_ID'] != member_remaining]
        membersRemoved = membersRemoved.reset_index(drop = True)

        # remove some additional members
        if numberToRemove > zero_member_count:
            difference = numberToRemove - zero_member_count
            if difference % 2 != 0:
                difference += 1
            newMembersRemaining = membersRemaining[difference:]
            newMembersRemaining = newMembersRemaining.reset_index(drop = True)
            additionalRemovals = membersRemaining[:difference]
            additionalRemovals = additionalRemovals.reset_index(drop = True)
            membersRemaining = newMembersRemaining
            membersRemoved = membersRemoved.append(additionalRemovals)
            membersRemoved = membersRemoved.reset_index(drop = True)

        # add back the members on the do-not-remove list
        for no_remove_member in no_remove_list['member_ID']:
            d = memberFcsSorted[memberFcsSorted['member_ID'] == no_remove_member]
            membersRemaining = membersRemaining.append(d)
        membersRemaining = membersRemaining.reset_index(drop = True)

        '''# update SAP_I_Member in order to enforce the continuity of the new member lists
        membersToRemove = membersRemoved['member_ID'].tolist()
        for member in membersToRemove:
            sapIMember.loc[sapIMember['member_ID'] == member,'size'] == '0'''
        return membersRemaining,membersRemoved,memberFcsSorted,sapIMember,parameters,desiredNumberToRemove

    def resetFixity(sapIMember):
        for i in range(0,len(sapIMember['member_ID'])):
            member = sapIMember.get_value(i,'member_ID')
            if str(sapIMember.loc[sapIMember['member_ID'] == member,'size']) != '0':
                if member[0:2] == 'BE' or member[0:2] == 'CO':
                    sapIMember.loc[sapIMember['member_ID'] == member,'M2I'] = True
                    sapIMember.loc[sapIMember['member_ID'] == member,'M2J'] = True
                    sapIMember.loc[sapIMember['member_ID'] == member,'M3I'] = True
                    sapIMember.loc[sapIMember['member_ID'] == member,'M3J'] = True
                else:
                    sapIMember.loc[sapIMember['member_ID'] == member,'M2I'] = True
                    sapIMember.loc[sapIMember['member_ID'] == member,'M2J'] = True
                    sapIMember.loc[sapIMember['member_ID'] == member,'M3I'] = True
                    sapIMember.loc[sapIMember['member_ID'] == member,'M3J'] = True
        return sapIMember
        
    def enforceMemberContinuity(sapIMember,nodeToMemberMap,parameters):
        sapIMember = resetFixity(sapIMember)
        
        sizing_continuity = eval(parameters['Values'][9])
        if sizing_continuity:
            element_continuity = prioritize(nodeToMemberMap,parameters)
        else:
            element_continuity = pd.DataFrame({'Member_A':[],'Member_B':[]})
         
        if eval(parameters['Values'][10]) == False:
            continuous_type = ['NA']
        else:
            continuous_type = ['CO','BE','BR']

        ''' Find collinear elements and change their fixity based on continuity rules'''
        continuous_elements = {}

        for i in range(0,len(element_continuity['Member_A'])):
            member_a = element_continuity.get_value(i,'Member_A')
            member_b = element_continuity.get_value(i,'Member_B')
            continuous_elements[member_a] = member_b

        for element in continuous_elements:
            element_a = element
            element_b = continuous_elements[element_a]
            a_start = sapIMember[sapIMember['member_ID'] == element_a].iloc[0]['start_node']
            a_end = sapIMember[sapIMember['member_ID'] == element_a].iloc[0]['end_node']
            b_start = sapIMember[sapIMember['member_ID'] == element_b].iloc[0]['start_node']
            b_end = sapIMember[sapIMember['member_ID'] == element_b].iloc[0]['end_node']
            if a_start == b_start and element[0:2] in continuous_type:
                sapIMember.loc[sapIMember['member_ID'] == element_a,'M2I'] = False
                sapIMember.loc[sapIMember['member_ID'] == element_a,'M3I'] = False
                sapIMember.loc[sapIMember['member_ID'] == element_b,'M2I'] = False
                sapIMember.loc[sapIMember['member_ID'] == element_b,'M3I'] = False
            elif a_end == b_end and element[0:2] in continuous_type:
                sapIMember.loc[sapIMember['member_ID'] == element_a,'M2J'] = False
                sapIMember.loc[sapIMember['member_ID'] == element_a,'M3J'] = False
                sapIMember.loc[sapIMember['member_ID'] == element_b,'M2J'] = False
                sapIMember.loc[sapIMember['member_ID'] == element_b,'M3J'] = False
            elif a_start == b_end and element[0:2] in continuous_type:
                sapIMember.loc[sapIMember['member_ID'] == element_a,'M2I'] = False
                sapIMember.loc[sapIMember['member_ID'] == element_a,'M3I'] = False
                sapIMember.loc[sapIMember['member_ID'] == element_b,'M2J'] = False
                sapIMember.loc[sapIMember['member_ID'] == element_b,'M3J'] = False
            elif a_end == b_start and element[0:2] in continuous_type:
                sapIMember.loc[sapIMember['member_ID'] == element_a,'M2J'] = False
                sapIMember.loc[sapIMember['member_ID'] == element_a,'M3J'] = False
                sapIMember.loc[sapIMember['member_ID'] == element_b,'M2I'] = False
                sapIMember.loc[sapIMember['member_ID'] == element_b,'M3I'] = False
        return sapIMember,element_continuity

    def removeVerticalAndHorizontalMechanisms(membersRemaining,membersRemoved):
        '''MECHANISM CHECK: avoiding vertical and horizontal mechanisms within membersRemaining by removing unstable elements'''
        verticalElements = []
        for i in range(0,len(membersRemaining['member_ID'])):   # finding all vertical elements
            start_node = membersRemaining.get_value(i,'start_node')
            end_node = membersRemaining.get_value(i,'end_node')
            start_node_Z = sapINode[sapINode['node_ID'] == start_node].iloc[0]['Z']
            end_node_Z = sapINode[sapINode['node_ID'] == end_node].iloc[0]['Z']
            if start_node_Z != end_node_Z:
                verticalElements.append(membersRemaining.get_value(i,'member_ID'))

        mechanism_elements = []
        for i in range(0,len(sapINode['node_ID'])):
            currentNode = sapINode.get_value(i,'node_ID')
            membersConnected = []
            for j in range(0,len(membersRemaining['member_ID'])):
                if membersRemaining.get_value(j,'start_node')==currentNode or membersRemaining.get_value(j,'end_node')==currentNode:
                    membersConnected.append(membersRemaining.get_value(j,'member_ID'))
            flag = True     # no vertical elements connected to the node
            for element in membersConnected:
                if element in verticalElements:
                    flag = False
            if (flag == True and len(membersConnected) != 0) or ( (currentNode not in loadedNodes) and (len(membersConnected) <= 2) and (len(membersConnected) != 0) ):
                for member in membersConnected:
                    mechanism_elements.append(member)

                    mechanism_element = membersRemaining[membersRemaining['member_ID'] == member]
                    mechanism_element = mechanism_element.reset_index(drop = True)
                    membersRemoved = membersRemoved.append(mechanism_element)
                    membersRemoved = membersRemoved.reset_index(drop = True)

        for member in mechanism_elements:
                membersRemaining = membersRemaining[membersRemaining['member_ID'] != member]
                membersRemaining = membersRemaining.reset_index(drop = True)
        return membersRemaining,membersRemoved

    def addBackMembersAtLoadedNodes(membersRemaining,membersRemoved,loadedNodes,nodeToMemberMap):
        ''' find lone members connected to loaded nodes and add back'''
        memberIDs = sapIMember["member_ID"]
        startNodes = sapIMember["start_node"]
        endNodes = sapIMember["end_node"]
        memSizes = sapIMember["size"]   # the difference between memSizes and memberSizes is that the array also has the 0 (removed) members
        memberIDsRemaining = membersRemaining["member_ID"].tolist()
        membersConsidered = membersRemoved["member_ID"]

        memsAtNodeMap = {}   # all the members from SAP_I_Member connected to the loaded nodes
        minNumber = 1

        for i in range(0, len(loadedNodes)):
            thisNode = loadedNodes[i]
            membersLoadedAtNode = []
            for element in nodeToMemberMap[thisNode]:
                if element in memberIDsRemaining:
                    membersLoadedAtNode.append(element)

            if len(membersLoadedAtNode) < minNumber:
                numberToAddBack = minNumber - len(membersLoadedAtNode)
                for element in membersConsidered:
                    if element in nodeToMemberMap[thisNode] and numberToAddBack > 0:
                        d = membersRemoved[membersRemoved["member_ID"] == element]
                        d = d.reset_index(drop = True)
                        membersRemaining = membersRemaining.append(d)
                        membersRemaining = membersRemaining.reset_index(drop = True)
                        membersRemoved = membersRemoved[membersRemoved["member_ID"] != element]
                        membersRemoved = membersRemoved.reset_index(drop = True)
                        numberToAddBack = numberToAddBack - 1

                    if numberToAddBack == 0:
                        break
        return membersRemaining,membersRemoved

    def ensureContinuityInBraceRemoval(membersRemaining,membersRemoved,nodeToMemberMap,element_continuity):
        ''' Ensure brace continuity in member removal (will only work with a consistent labeling convention)'''

        member_continuity = {}
        for i in range(0,len(element_continuity['Member_A'])):
            member_a = element_continuity.get_value(i,'Member_A')
            member_b = element_continuity.get_value(i,'Member_B')
            member_continuity[member_a] = member_b
            member_continuity[member_b] = member_a

        for i in range(0,len(membersRemoved)):
            member1 = membersRemoved.get_value(i,'member_ID')
            node1 = membersRemoved.get_value(i,'start_node')
            node2 = membersRemoved.get_value(i,'end_node')
            if member1[0:2] == 'BR':
                other_members = nodeToMemberMap[node1]
                other_members += nodeToMemberMap[node2]
                # print other_members
                for member2 in other_members:
                    if member2[0:2] == 'BR' and member2 in list(membersRemaining['member_ID']) and member2 in member_continuity.keys() and member1 in member_continuity.keys() and (member_continuity[member1] == member2 or member_continuity[member2] == member1):
                        secondMemberData = membersRemaining[membersRemaining['member_ID'] == member2]
                        secondMemberData = secondMemberData.reset_index(drop = True)
                        membersRemoved = membersRemoved.append(secondMemberData)
                        membersRemoved = membersRemoved.reset_index(drop = True)
                        membersRemaining = membersRemaining[membersRemaining['member_ID'] != member2]
                        membersRemaining = membersRemaining.reset_index(drop = True)
        return membersRemaining,membersRemoved

    def ensureContinuityInBraceAddBack(membersRemaining,membersRemoved,nodeToMemberMap):
        ''' Ensure brace continuity in member removal (will only work with a consistent labeling convention)'''
        for i in range(0,len(membersRemoved)):
            member1 = membersRemaining.get_value(i,'member_ID')
            node1 = membersRemaining.get_value(i,'start_node')
            node2 = membersRemaining.get_value(i,'end_node')
            if member1[0:2] == 'BR':
                other_members = nodeToMemberMap[node1]
                other_members += (nodeToMemberMap[node2])
                for member2 in other_members:
                    if member2[0:2] == 'BR' and member2 in list(membersRemoved['member_ID']) and ((int(member1[2:len(member1)]) - int(member2[2:len(member1)])) == 1 or (int(member1[2:len(member1)]) - int(member2[2:len(member1)])) == -1):
                        secondMemberData = membersRemoved[membersRemoved['member_ID'] == member2]
                        secondMemberData = secondMemberData.reset_index(drop = True)
                        membersRemaining = membersRemaining.append(secondMemberData)
                        membersRemaining = membersRemaining.reset_index(drop = True)
                        membersRemoved = membersRemoved[membersRemoved['member_ID'] != member2]
                        membersRemoved = membersRemoved.reset_index(drop = True)

        return membersRemaining,membersRemoved

    def removeFloatingAndSingleElements(membersRemaining,membersRemoved,no_remove_list,dataNone):
        ''' Ensure floating elements and single-element-mechanisms are removed'''
        checklist = no_remove_list['member_ID'].tolist()
        for k in range(0,1):
            floating_members = []
            for i in range(0,len(membersRemaining['member_ID'])):
                member = membersRemaining.get_value(i,"member_ID")
                if member not in checklist:
                    member_start = membersRemaining.get_value(i,'start_node')
                    member_end = membersRemaining.get_value(i,'end_node')

                    membersConnectedtoStart = 0
                    membersConnectedtoEnd = 0

                    for j in range(0,len(membersRemaining['member_ID'])):
                        if member != membersRemaining.get_value(j,'member_ID'):
                            if membersRemaining.get_value(j,'start_node') == member_start or membersRemaining.get_value(j,'end_node') == member_start:
                                membersConnectedtoStart = membersConnectedtoStart + 1
                            if membersRemaining.get_value(j,'start_node')== member_end or membersRemaining.get_value(j,'end_node')== member_end:
                                membersConnectedtoEnd = membersConnectedtoEnd + 1

                    start_support = False
                    end_support = False

                    if dataNode.isSupport[member_start] == 1:
                        start_support = True
                    elif dataNode.isSupport[member_end] == 1:
                        end_support = True

                    if (membersConnectedtoStart == 0 and start_support == False) or (membersConnectedtoEnd == 0 and end_support == False):
                        floating_members.append(member)

                        floating_member = membersRemaining[membersRemaining['member_ID'] == member]
                        floating_member = floating_member.reset_index(drop = True)
                        membersRemoved = membersRemoved.append(floating_member)
                        membersRemoved = membersRemoved.reset_index(drop = True)


            for element in floating_members:
                membersRemaining = membersRemaining[membersRemaining['member_ID'] != element]
                membersRemaining = membersRemaining.reset_index(drop = True)
        return membersRemaining,membersRemoved

    def finalizeNormalRunOutputs(nodeToMemberMap,membersRemoved,membersRemaining,totalCost,parameters,sapIMember,no_remove_list):
        '''emptyD = {"member_ID": [], "checked?": [], "cross_section": []}
        lastRemovedNew = pd.DataFrame(data  = emptyD)'''
        lastRemovedNew = membersRemoved

        d = {"member_ID": membersRemaining["member_ID"], "cross_section": membersRemaining["cross_section"]}
        newMemberList = pd.DataFrame(data = d)
        d = {"member_ID": membersRemoved["member_ID"], "cross_section": ["0"] * len(membersRemoved["member_ID"])}
        newMemberList = newMemberList.append(pd.DataFrame(data = d))
        newMemberList = newMemberList.reset_index(drop = True)
        lastRemoved = lastRemovedNew

        memberIDs = newMemberList["member_ID"]
        sizes = newMemberList["cross_section"]
        for i in range(0, len(memberIDs)):
            sapIMember.loc[sapIMember["member_ID"] == memberIDs[i], "size"] = sizes[i]  # get rid of removed (0) members from the list

        sapIMember,element_continuity = enforceMemberContinuity(sapIMember,nodeToMemberMap,parameters)
        iteration = int(parameters.get_value(3, "Values"))
        iteration = iteration + 1
        parameters['Values'][3] = iteration

        return parameters,sapIMember,lastRemoved,no_remove_list

    def totalBaseReactionsCheck(total_reactions,unstable):
        current_iteration = int(parameters['Values'][3])
        totalReaction = 0
        for i in range(0,len(nodeReactions['vertical_reaction_at_Column_Base'])):
            totalReaction = totalReaction + np.ceil( nodeReactions.get_value(i,'horizontal_reaction_x_at_Column_Base'))
            #totalReaction = totalReaction + np.ceil( nodeReactions.get_value(i,'vertical_reaction_at_Column_Base') + nodeReactions.get_value(i,'horizontal_reaction_x_at_Column_Base') + nodeReactions.get_value(i,'horizontal_reaction_y_at_Column_Base'))
        d = pd.DataFrame({"iteration": current_iteration,"Total_Reactions": totalReaction},index=[current_iteration]).set_index('iteration')
        total_reactions = total_reactions.append(d,ignore_index = True)
        total_reactions.reset_index(drop = True)

        eliminated_load_path = False
        if unstable == False and current_iteration >= 1 and abs(total_reactions.iloc[-1]["Total_Reactions"] - total_reactions.iloc[0]["Total_Reactions"]) >= abs(total_reactions.iloc[-1]["Total_Reactions"])*0.5:
            eliminated_load_path = True
        eliminated_load_path = False #temporary
        return total_reactions, eliminated_load_path

    def costIncreaseCheck(cost_history,parameters,unstable):
        totalCost = 0
        current_iteration = int(parameters['Values'][3])

        if not unstable:
            for cost in memberCosts['cost']:
                totalCost += cost
        else:
            totalCost = -1

        costIncrease = False
        if eval(parameters['Values'][14]) == True:
            tolerance_index = 1.5
        else:
            tolerance_index = 1.05

        if current_iteration >= 10 and totalCost > tolerance_index*cost_history.iloc[-1]["cost"] and unstable == False:
            costIncrease = True

        if unstable == False and costIncrease == False:
            # only record feasible solutions
            d = pd.DataFrame({"iteration": current_iteration,"cost": totalCost},index=[current_iteration]).set_index('iteration')
            cost_history = cost_history.append(d)
            cost_history = cost_history.reset_index(drop=True)

        return cost_history, totalCost, costIncrease

    def instabilityCheck(parameters,modeFrequencies):
        current_iteration = int(parameters['Values'][3])
        minimum_iteration = 1
        unstable = False
        if current_iteration >= minimum_iteration and np.abs(float(modeFrequencies['Frequency'][0])) <= 0.01:
            unstable = True
        if parameters['Values'][13] == False:
            unstable = True
        '''elif current_iteration < minimum_iteration and np.abs(modeFrequencies['Frequency'][0]) <= 0.001:
            unstable = True'''
        return unstable

    def coplanarCheck(node,members,dataMemb,dataNode):
        "returns True if all members are coplanar"
        nMem=len(members)
        if nMem<=2:
            return True
            #break
        Vector=np.zeros((nMem,3))
        for i in range(0,nMem):
            if dataMemb.start_node[members[i]]==node:
                Vector[i,0]=dataNode.x_coord[dataMemb.end_node[members[i]]]-dataNode.x_coord[node]
                Vector[i,1]=dataNode.y_coord[dataMemb.end_node[members[i]]]-dataNode.y_coord[node]
                Vector[i,2]=dataNode.z_coord[dataMemb.end_node[members[i]]]-dataNode.z_coord[node]
            elif dataMemb.end_node[members[i]]==node:
                Vector[i,0]=dataNode.x_coord[dataMemb.start_node[members[i]]]-dataNode.x_coord[node]
                Vector[i,1]=dataNode.y_coord[dataMemb.start_node[members[i]]]-dataNode.y_coord[node]
                Vector[i,2]=dataNode.z_coord[dataMemb.start_node[members[i]]]-dataNode.z_coord[node]
            else:
                return 'members no connectied to node'
        Base=np.cross(Vector[0,:],Vector[1,:])
        i = 0
        #while np.linalg.norm(Base)==0 and i<=nMem:
        while abs(np.linalg.norm(Base))<400 and i<=nMem:
            i +=1 
            Vector2=Vector[1:nMem,:]
            Vector2=np.vstack([Vector2, Vector[0,:]])
            Vector=Vector2
            Base=np.cross(Vector[0,:],Vector[1,:])

        Planar=True
        for i in range(2,nMem):
            #if abs(np.dot(Base,Vector[i,:]))>1:
            if abs(np.dot(Base,Vector[i,:]))>150:
                Planar=False
        return Planar

    def colinearCheck(node,members,sapIMember,sapINode):
        dataMemb = sapIMember.set_index('member_ID')
        dataNode = sapINode.set_index('node_ID')

        nodes = []
        if dataMemb.start_node[members[0]] != node:
            nodes.append(dataMemb.start_node[members[0]])
        if dataMemb.end_node[members[0]] != node:
            nodes.append(dataMemb.end_node[members[0]])
        if dataMemb.start_node[members[1]] != node:
            nodes.append(dataMemb.start_node[members[1]])
        if dataMemb.end_node[members[1]] != node:
           nodes.append(dataMemb.end_node[members[1]])
        node1 = nodes[0]
        node2 = node
        node3 = nodes[1]
        # ASSUMED 2D STRUCTURE IN X-Z PLANE

        if dataNode.X[node2] == dataNode.X[node1] and dataNode.X[node1] != dataNode.X[node3]:
            return False
        elif dataNode.X[node2] == dataNode.X[node1] and dataNode.X[node1] == dataNode.X[node3]:
            return True


        m = (dataNode.Z[node2] - dataNode.Z[node1]) / (dataNode.Z[node2] - dataNode.Z[node1])
        if np.floor(m*(dataNode.X[node3] - dataNode.X[node1]) + dataNode.Z[node1]) == np.floor(dataNode.Z[node3]):
            return True
        else:
            return False

    def isNodeStable(parameters,node,nodeToElementRemainingMapList,sapIMember,sapINode):
        dimensionality = int(parameters['Values'][11])
        # The node is stable if all attached elements are removed, or if it is a support node
        if len(nodeToElementRemainingMapList) == 0 or sapINode[sapINode['node_ID'] == node].iloc[0]['ground?'] == True:
            return True

        # stability due to fixity at adjacent nodes


        # for element in nodeToElementRemainingMapList:
        #     #print element
        #     element_start = sapIMember[sapIMember['member_ID'] == element].iloc[0]['start_node']
        #     element_end = sapIMember[sapIMember['member_ID'] == element].iloc[0]['end_node']
        #     if (sapIMember[sapIMember['member_ID'] == element].iloc[0]['M3I'] == False) and (sapIMember[sapIMember['member_ID'] == element].iloc[0]['M3J'] == False):
        #         return True
        #     else:
        #         print ' - no through continuity'

        # if stability through member fixity isn't available for this node, we check it's out-of-plane stability:
        # less than 3 elements cannot be stable out of plane
        if len(nodeToElementRemainingMapList) < dimensionality:
            print (' - less than',dimensionality, 'elements')
            return False
        # if all the members are coplanar, the node is not stable
        if dimensionality == 3 and coplanarCheck(node,nodeToElementRemainingMapList,sapIMember,sapINode) == True:
            print (' - coplanar')
            return False
        elif dimensionality == 2 and len(nodeToElementRemainingMapList) == 2 and colinearCheck(node,nodeToElementRemainingMapList,sapIMember,sapINode) == True:
            if eval(parameters['Values'][10]) == True:
                return True
            else:
                return False

        return True

    def addBestElementBack(parameters,node,nodeToElementRemovedMap,nodeToElementRemainingMap,membersRemoved,membersRemaining,nodesToCheck,sapIMember,sapINode):
        candidates = nodeToElementRemovedMap[node]
        stableCandidateList = list(candidates)
        
        for candidate in candidates:
            stability_check_group = nodeToElementRemainingMap[node] + [candidate]
            if isNodeStable(parameters,node,stability_check_group,sapIMember,sapINode) == False:
                stableCandidateList.remove(candidate)

        if len(stableCandidateList) == 0:
            best = candidates[0]
            for element in candidates:
                if membersRemoved[membersRemoved['member_ID'] == element].iloc[0]["member_Fcs"] < membersRemoved[membersRemoved['member_ID'] == best].iloc[0]["member_Fcs"]:
                    best=element
        else:
            best = stableCandidateList[0]
            for element in stableCandidateList:
                if membersRemoved[membersRemoved['member_ID'] == element].iloc[0]["member_Fcs"] < membersRemoved[membersRemoved['member_ID'] == best].iloc[0]["member_Fcs"]:
                    best=element


        # update member lists accordingly
        add_back = membersRemoved[membersRemoved['member_ID'] == best]
        add_back = add_back.reset_index(drop = True)
        membersRemoved = membersRemoved[membersRemoved['member_ID'] != best]
        membersRemoved = membersRemoved.reset_index(drop = True)
        membersRemaining = membersRemaining.append(add_back)
        membersRemaining = membersRemaining.reset_index(drop = True)

        # add member back in SAP_I_Member
        sapIMember.loc[sapIMember['member_ID'] == best,'size'] = membersRemaining[membersRemaining['member_ID'] == best].iloc[0]['cross_section']


        # update the loop variables and data structures to reflect the member being added back
        new_start = membersRemaining[membersRemaining['member_ID'] == best].iloc[0]["start_node"]
        new_end = membersRemaining[membersRemaining['member_ID'] == best].iloc[0]["end_node"]

        if node == new_start:
            new_node_to_check = new_end
        elif node == new_end:
            new_node_to_check = new_start

        nodesToCheck.append(new_node_to_check)

        nodeToElementRemainingMap[node].append(best)
        nodeToElementRemovedMap[node].remove(best)

        # We now need to check if the new element is shared with any other node to check. If so, add that element to the nodetoelementremainingmap as well
        for otherNode in nodesToCheck:
            if otherNode == new_node_to_check and best not in nodeToElementRemainingMap[otherNode]:
                nodeToElementRemainingMap[otherNode].append(best)
            if otherNode == new_node_to_check and best in nodeToElementRemovedMap[otherNode]:
                nodeToElementRemovedMap[otherNode].remove(best)

        return membersRemoved,membersRemaining,nodesToCheck,nodeToElementRemovedMap,nodeToElementRemainingMap,sapIMember

    def preventMultiMemberMechanisms(loadedNodes,parameters,membersRemaining,membersRemoved,nodeToMemberMap,sapIMember,sapINode,no_remove_list,desiredNumberToRemove,removalMultiplier):
        # at every iteration, try to remove all the members in the temporary no remove list, plus the ones removed by the ranking function

        desired_remove_list = membersRemoved['member_ID'].tolist()
        print ('desired remove list',desired_remove_list)

        current_iteration = int(parameters['Values'][3])
        nodesToCheck = []

        for i in range(0,len(membersRemoved['member_ID'])):
            if membersRemoved.get_value(i,'start_node') not in nodesToCheck:
                nodesToCheck.append(membersRemoved.get_value(i,'start_node'))
            if membersRemoved.get_value(i,'end_node') not in nodesToCheck:
                nodesToCheck.append(membersRemoved.get_value(i,'end_node'))
        print ('intial nodes to check', len(nodesToCheck),nodesToCheck)

        nodeToElementRemainingMap = {}
        nodeToElementRemovedMap = {}
        membersRemainingList = membersRemaining['member_ID'].tolist()
        for node in nodesToCheck:
            nodeToElementRemainingMap[node] = []
            nodeToElementRemovedMap[node] = []
            for member in nodeToMemberMap[node]:
                if member in membersRemainingList:
                    nodeToElementRemainingMap[node].append(member)
                else:
                    nodeToElementRemovedMap[node].append(member)

        # If the node is removed completely (has no remaining elements) or if it is unstable to begin with, do not attempt to stabilize it
        nodeList = list(nodesToCheck)

        initially_unstable = []
        for node in nodeList:
            if isNodeStable(parameters,node,nodeToMemberMap[node],sapIMember,sapINode) == False:
                initially_unstable.append(node)
                nodesToCheck.remove(node)
        
        print ('final nodes to check',len(nodesToCheck),nodesToCheck)
        # assume all nodes are unstable. Check each one of them and stabilize by adding a single member back iteratively when necessary
        while len(nodesToCheck) != 0:
            nodeList = list(nodesToCheck)
            #print nodeList
            for node in nodeList:
                #print node
                if isNodeStable(parameters,node,nodeToElementRemainingMap[node],sapIMember,sapINode) == True or node in initially_unstable:
                    nodesToCheck.remove(node)
                    print ('stable',node)
                else:
                    print ('not stable', node)
                    membersRemoved,membersRemaining,nodesToCheck,nodeToElementRemovedMap,nodeToElementRemainingMap,sapIMember = addBestElementBack(parameters,node,nodeToElementRemovedMap,nodeToElementRemainingMap,membersRemoved,membersRemaining,nodesToCheck,sapIMember,sapINode)


        #membersRemaining,membersRemoved = removeFloatingAndSingleElements(membersRemaining,membersRemoved,no_remove_list,dataNode)
        membersRemaining,membersRemoved = addBackMembersAtLoadedNodes(membersRemaining,membersRemoved,loadedNodes,nodeToMemberMap)

        print ('attempted number to remove',len(membersRemoved['member_ID']))
        print ('desired to remove', desiredNumberToRemove)
        print ('removed members', membersRemoved['member_ID'])
        print ('length of members remaining', len(membersRemaining))
        print ('should be backfilling?', len(membersRemoved) < desiredNumberToRemove)

        if len(membersRemoved) == 0 and len(desired_remove_list) == len(membersRemaining) - len(no_remove_list):
            for element in desired_remove_list:
                d = pd.DataFrame([str(element)],columns=['member_ID'])
                no_remove_list = no_remove_list.append(d,ignore_index = True)
                no_remove_list = no_remove_list.reset_index(drop = True)
            backFill = False
        elif len(membersRemoved) != 0 and len(desired_remove_list) == len(membersRemaining) - len(no_remove_list):
            backFill = False
        elif len(membersRemoved) >= desiredNumberToRemove:
            backFill = False
        elif len(membersRemoved) < desiredNumberToRemove:
            backFill = True
            removalMultiplier = removalMultiplier + 1

        #print 'NODE STABLE CHECK', isNodeStable(dataNode,dataMemb,'P020306',nodeToElementRemainingMap['P020306'],sapIMember)

        return membersRemaining,membersRemoved,sapIMember,no_remove_list,backFill,removalMultiplier

    def initializeLocalSearch(nodeToMemberMap,lastRemoved,sapIMember,no_remove_list,parameters):
        current_iteration = parameters['Values'][3]
        localSearch = 'False'

        # Add all the members back
        for removedMember in lastRemoved['member_ID']:
            sapIMember.loc[sapIMember['member_ID'] == removedMember,'size'] = lastRemoved[lastRemoved['member_ID'] == removedMember].iloc[0]['cross_section']

        if len(lastRemoved['member_ID']) == 1:
            last_removed_member = lastRemoved['member_ID'][0]
            d = pd.DataFrame({'member_ID':last_removed_member},index = [int(current_iteration)])
            no_remove_list = no_remove_list.append(d)
            no_remove_list = no_remove_list.reset_index(drop = True)
            # no_remove_list.to_csv("no_remove_list.txt", index = False)     
        else:
            localSearch = 'True'
            # remove only the first one on the list, and update the list
            toRemove = lastRemoved['member_ID'][0]
            sapIMember.loc[sapIMember['member_ID'] == toRemove,'size'] = '0'

        parameters['Values'][6] = localSearch
        sapIMember,element_continuity = enforceMemberContinuity(sapIMember,nodeToMemberMap,parameters)
        return sapIMember,lastRemoved,parameters,no_remove_list

    def manageLocalSearchIteration(nodeToMemberMap,cost_history,total_reactions,no_remove_list,lastRemoved,sapIMember,parameters,eliminated_load_path,unstable,costIncrease):

        current_iteration = int(parameters['Values'][3])

        #no_remove_list = pd.DataFrame({'member_ID':[]})
        #no_remove_list.to_csv('no_remove_list.txt',index = False)
        # The member removed is necessary: add it back and try the next one
        proceed = False
        if (eliminated_load_path == True or unstable == True or costIncrease == True):
            # if current_iteration > 20:
            last_removed_member = lastRemoved['member_ID'][0]
            d = pd.DataFrame({'member_ID':last_removed_member},index = [current_iteration])
            no_remove_list = no_remove_list.append(d)
            no_remove_list = no_remove_list.reset_index(drop = True)

            # add the member back and remove it from the list
            sapIMember.loc[sapIMember['member_ID'] == last_removed_member,'size'] = lastRemoved[lastRemoved['member_ID'] == last_removed_member].iloc[0]['cross_section']
            lastRemoved = lastRemoved[lastRemoved['member_ID'] != last_removed_member]
            lastRemoved = lastRemoved.reset_index(drop = True)
            # If there are no more members to check, local search is over
            if len(lastRemoved['member_ID']) == 0:
                parameters['Values'][6] = 'False'
            else:
                toRemove = lastRemoved['member_ID'][0]
                sapIMember.loc[sapIMember['member_ID'] == toRemove,'size'] = '0'

        # If a single member that can be removed is found, exit local search
        else:
            proceed = True
            parameters['Values'][6] = 'False'

        sapIMember,element_continuity = enforceMemberContinuity(sapIMember,nodeToMemberMap,parameters)
        return sapIMember,parameters,lastRemoved,no_remove_list,proceed,parameters,no_remove_list

    def checkForConvergence(no_remove_list,memberSizes,parameters):
        localSearch = eval(parameters['Values'][6])
        convergedBool = False
        if (len(no_remove_list['member_ID']) == len(list(memberSizes.index)) ) and localSearch == False:
            convergedBool = True
        return convergedBool

    def reset_no_remove_list(parameters,convergedOutput,no_remove_list):
        first_time_convergence = eval(parameters['Values'][14])
        mdo_converged = False
        if first_time_convergence == True:
            print ('FIRST TIME CONVERGENCE - RESETTING NO REMOVE LIST')
            parameters['Values'][14] = 'False'
            no_remove_list = pd.DataFrame({'member_ID':[]})
        else:
            print ('\t CONVERGED!')
            mdo_converged = True
        return mdo_converged,parameters,no_remove_list

    def saveCheckpoint(lastRemoved,no_remove_list,cost_history,total_reactions,parameters,sapIMember,sapINode,mdo_converged):

        current_iteration = int(parameters['Values'][3])
        path = './checkpoints/' + str(current_iteration) + '/'
        directory = os.path.dirname(path)
        if not os.path.exists(directory):
            os.makedirs(directory)

        lastRemoved.to_csv(path + 'last_removed.txt',index=False)
        no_remove_list.to_csv(path+ 'no_remove_list.txt',index=False)
        cost_history.to_csv(path + 'cost_history.txt',index=False)
        total_reactions.to_csv(path + 'Total_Reactions.txt',index=False)
        parameters.to_csv(path + 'inputs.txt',index=False)
        sapIMember.to_csv(path + 'SAP_I_Member.txt',index=False)
        sapINode.to_csv(path + 'SAP_I_Node.txt',index=False)

        dir_src = ("./SAPMODEL/")
        dir_dst = (path)
        for filename in os.listdir(dir_src):
            if filename.endswith('.sdb'):
                shutil.copy( dir_src + filename, dir_dst)
        print ('Checkpoint' + str(current_iteration) + 'saved!')


    ### MAIN EXECUTABLE PROGRAM ###

    mdo_converged = checkForConvergence(no_remove_list,memberSizes,parameters)
    if not mdo_converged:
        mdo_converged = False
        unstable = instabilityCheck(parameters,modeFrequencies)
        cost_history,totalCost,costIncrease = costIncreaseCheck(cost_history,parameters,unstable)
        total_reactions,eliminated_load_path = totalBaseReactionsCheck(total_reactions,unstable)

        print ('costIncrease: ',costIncrease)
        print ('unstable: ',unstable)
        print ('eliminated load path: ',eliminated_load_path)

        localSearch = eval(parameters['Values'][6])
        proceed = True


        if localSearch == True:
            print ('MANAGING LOCAL SEARCH ITERATION')
            sapIMember,parameters,lastRemoved,no_remove_list,proceed,parameters,no_remove_list = manageLocalSearchIteration(nodeToMemberMap,cost_history,total_reactions,no_remove_list,lastRemoved,sapIMember,parameters,eliminated_load_path,unstable,costIncrease)
        elif localSearch == False and (eliminated_load_path == True or costIncrease == True or unstable == True):
            print ('INITIALIZING LOCAL SEARCH')
            proceed = False
            sapIMember,lastRemoved,parameters,no_remove_list = initializeLocalSearch(nodeToMemberMap,lastRemoved,sapIMember,no_remove_list,parameters)

        if proceed == True:
            loadedNodes = findLoadedNodes()
            removalMultiplier = 1
            backFill = True
            while backFill == True:
                membersRemaining,membersRemoved,memberFcsSorted,sapIMember,parameters,desiredNumberToRemove = rankAndRemoveMembers(sapIMember,memberCosts,memberSizes,memberStiffnessContribution,memberForces,no_remove_list,removalMultiplier,parameters)
                membersRemaining,membersRemoved,sapIMember,no_remove_list,backFill,removalMultiplier = preventMultiMemberMechanisms(loadedNodes,parameters,membersRemaining,membersRemoved,nodeToMemberMap,sapIMember,sapINode,no_remove_list,desiredNumberToRemove,removalMultiplier)
                print ('---- BACKFILL: ---- ',backFill)
            #membersRemaining,membersRemoved = removeFloatingAndSingleElements(membersRemaining,membersRemoved,no_remove_list)
            #membersRemaining,membersRemoved = addBackMembersAtLoadedNodes(membersRemaining,membersRemoved,loadedNodes,nodeToMemberMap)
            sapIMember,element_continuity = enforceMemberContinuity(sapIMember,nodeToMemberMap,parameters)
            #membersRemaining,membersRemoved = ensureContinuityInBraceRemoval(membersRemaining,membersRemoved,nodeToMemberMap,element_continuity)
            parameters,sapIMember,lastRemoved,no_remove_list = finalizeNormalRunOutputs(nodeToMemberMap,membersRemoved,membersRemaining,totalCost,parameters,sapIMember,no_remove_list)
    else:
        cost_history,totalCost,costIncrease = costIncreaseCheck(cost_history,parameters,False) # cost the structure at convergence
        mdo_converged,parameters,no_remove_list = reset_no_remove_list(parameters,mdo_converged,no_remove_list)

    if checkpoint: saveCheckpoint(lastRemoved,no_remove_list,cost_history,total_reactions,parameters,sapIMember,sapINode,mdo_converged)
    return lastRemoved,no_remove_list,cost_history,total_reactions,parameters,sapIMember,sapINode,mdo_converged

# optimize()