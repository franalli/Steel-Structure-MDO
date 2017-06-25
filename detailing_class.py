# Subject to BSD-3-Clause
# Copyright (c) 2016-2017 Filippo Ranalli, Forest Flager, Pratyush Havelia
# See LICENSE.txt and NOTICE.txt

import pandas as pd 
import numpy as np 

def detail(sapImember,sapInode):
	print('\t DETAILING')
	nodeToElementMap = {}
	for i in range(len(sapInode['node_ID'])):
		node = sapInode.get_value(i,'node_ID')
		nodeToElementMap[node] = []
		for j in range(len(sapImember['member_ID'])):
			member = sapImember.get_value(j,'member_ID')
			start_node = sapImember.get_value(j,'start_node')
			end_node = sapImember.get_value(j,'end_node')
			if (start_node == node or end_node == node) and str(sapImember.get_value(j,'size')) != '0':
				nodeToElementMap[node].append(member)
	return nodeToElementMap
