# Subject to BSD-3-Clause
# Copyright (c) 2016-2017 Filippo Ranalli, Forest Flager, Pratyush Havelia
# See LICENSE.txt and NOTICE.txt

def cost(nodeToElementMap,continuity,sapImember,weightCosts, connection_ratio = 0.0,logging= False):
    print('COSTING')
    import pandas as pd
    import numpy as np

    def assembleFullContinuityList(continuity):
        continuity_lists = []
        for i in range(0,len(continuity['Member_A'])):
            continuity_lists.append([continuity.get_value(i,'Member_A'),continuity.get_value(i,'Member_B')])

        for cont_list in continuity_lists:
            for member in cont_list:
                for other_cont_list in continuity_lists:
                    if member in other_cont_list and other_cont_list != cont_list:
                        other_cont_list.remove(member)
                        cont_list += other_cont_list
                        continuity_lists.remove(other_cont_list)

        return continuity_lists

    def continuityPairs(continuity):
        member_continuity = {}
        for i in range(0,len(continuity['Member_A'])):
            member_a = continuity.get_value(i,'Member_A')
            member_b = continuity.get_value(i,'Member_B')
            if member_a not in member_continuity:
                member_continuity[member_a] = [member_b]
            else:
                member_continuity[member_a].append(member_b)
            if member_b not in member_continuity:
                member_continuity[member_b] = [member_a]
            else:
                member_continuity[member_b].append(member_a)
        return member_continuity
    def jointCost(members,continuity_pairs,weightCosts,connection_ratio):
        connection_cost = 0.0
        if (len(members) == 2 and (members[0] in continuity_pairs.keys() and members[1] in continuity_pairs[members[0]])) or len(members) <2:
            return connection_cost
        else:
            for element in members:
                connection_cost += connection_ratio * weightCosts['cost'][element]
            return connection_cost
    def distributeConnectionCostLocally(joint_cost,members,continuity_pairs):
        memberJointCost = {}
        continuous_members = []
        for member in members:
            if member in continuity_pairs.keys() and continuity_pairs[member] in members:
                continuous_members = [member,continuity_pairs[member]]
                break
        for member in continuous_members:
            members.remove(member)
        try:
            cost_fraction = joint_cost/ (len(members)+max(0,len(continuous_members)-1))
        except:
            embed()
        for member in members:
            memberJointCost[member] = cost_fraction
        for member in continuous_members:
            memberJointCost[member] = cost_fraction/2.0
        return memberJointCost
    def distributeConnectionCostGlobally(connectionCosts,memberJointCost,continuity_lists):
        for member in memberJointCost.keys():
            member_continuity = False
            for cont_list in continuity_lists:
                if member in cont_list:
                    for element in cont_list:
                        connectionCosts['cost'][element] += memberJointCost[member]/len(cont_list)
                        member_continuity = True
                    break
            if member_continuity == False:
                connectionCosts['cost'][member] += memberJointCost[member]
        return connectionCosts
    def totalCosts(weightCosts,connectionCosts):
        weightCosts['cost'] += connectionCosts['cost']
        return weightCosts

    continuity_lists = assembleFullContinuityList(continuity)
    continuity_pairs = continuityPairs(continuity)
    # print(weightCosts)
    # connectionCosts = pd.DataFrame({'member_ID':list(weightCosts.index),'cost':np.zeros(len(weightCosts['member_ID']))}).set_index('member_ID')
    # weightCosts = weightCosts.set_index('member_ID')

    connectionCosts = pd.DataFrame({'member_ID':list(weightCosts.index),'cost':np.zeros(len(list(weightCosts.index)))}).set_index('member_ID')

    for joint in nodeToElementMap.keys():
        if len(nodeToElementMap[joint]) != 0:
            members = nodeToElementMap[joint]
            joint_cost = jointCost(members,continuity_pairs,weightCosts,connection_ratio)
            memberJointCost = distributeConnectionCostLocally(joint_cost,members,continuity_pairs)
            connectionCosts = distributeConnectionCostGlobally(connectionCosts,memberJointCost,continuity_lists)

    totalCosts = totalCosts(weightCosts,connectionCosts)
    if logging: 
        totalCosts.to_csv("./logs/member_cost.txt")
    return totalCosts






