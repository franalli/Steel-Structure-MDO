# Subject to BSD-3-Clause
# Copyright (c) 2016-2017 Filippo Ranalli, Forest Flager, Pratyush Havelia
# See LICENSE.txt and NOTICE.txt

import scipy.integrate as integ
import pandas as pd
def calculate_virtual_work(member_info,sec_info,forces,real,virtual,verbose=False):

	if verbose: print ('------------- Calculating Virtual Work, COMBINATIONS:',real,virtual, '-------------')

	'####### ONLY WORKS IN 2D AND FOR W SECTIONS FOR NOW ########'
	# print(member_info)
	member_length_map = {}
	for i in range(0,len(member_info['member_length'])):
		if member_info['member_ID'][i] not in member_length_map:
			member_length_map[member_info['member_ID'][i]] = []
		member_length_map[member_info['member_ID'][i]] = member_info['member_length'][i]


	force_map = {}
	for i in range(0,len(forces['member_ID'])):
		if forces['Load_case'][i] not in force_map:
			force_map[forces['Load_case'][i]] = []
		member_forces = {}
		member_forces[forces['member_ID'][i]] = []
		member_forces[forces['member_ID'][i]].append(forces['P'][i])
		member_forces[forces['member_ID'][i]].append(forces['V'][i])
		member_forces[forces['member_ID'][i]].append(forces['Mx_start'][i])
		member_forces[forces['member_ID'][i]].append(forces['Mx_end'][i])
		force_map[forces['Load_case'][i]].append(member_forces)

	real_areas = {}
	for member in force_map[real]:
		label = list(member.keys())[0]
		L = member_length_map[label]
		#L = 9144
		real_areas[label] = []
		real_areas[label].append(member[label][0] * L) # P
		real_areas[label].append(member[label][1] * L) # V
		real_areas[label].append(member[label][2])
		real_areas[label].append(member[label][3])

	virtual_areas = {}
	for member in force_map[virtual]:
		label = list(member.keys())[0]
		L = member_length_map[label]
		#L = 9144
		virtual_areas[label] = []
		virtual_areas[label].append(member[label][0] * L) # P
		virtual_areas[label].append(member[label][1]* L) # V
		virtual_areas[label].append(member[label][2])
		virtual_areas[label].append(member[label][3])


	# Vtot = 0
	# for member in member_info['member_ID']:
	# 	section = sapImember['size'][member]
	# 	A = W_info[W_info["AISC_Manual_Label"] == section].iloc[0]['A']
	# 	Vtot += member_length_map[member] * A


	virtual_work = {}
	for member in member_info['member_ID']:

		L = member_length_map[member]
		#L = 9144
		#section = sap_I_member['size'][member]
		#section = member_info['size'][member]
		section = member_info[member_info['member_ID']==member].iloc[0]['size']
		member_type = member_info[member_info['member_ID']==member].iloc[0]['member_type']
		info = sec_info[member_type]
		
		E = 200000
		A = info['A'][section]
		Av = A - 2 * info['bf'][section] * info['tf'][section]
		I = info['Ix'][section]
		V = L * A
		virtual_work[member] = 0
		virtual_work[member] = max(0,real_areas[member][0]*virtual_areas[member][0]/(E*A)) + max(0,real_areas[member][1]*virtual_areas[member][1]/(E*Av)) # P + V
		B = real_areas[member][2]
		C = real_areas[member][3]
		D = virtual_areas[member][2]
		F = virtual_areas[member][3]
		m1 = (C-B)/L
		m2 = (F-D)/L
		virtual_work[member] += integ.quad(lambda x: (m1*x+B) * (m2*x+D),0,L)[0] / (E*I) # M

	normalized_VW = pd.DataFrame()
	for member in virtual_work:
		normalized_VW_line = pd.DataFrame({'VW':[virtual_work[member]],'member_ID':[member]})
		normalized_VW = normalized_VW.append(normalized_VW_line)
		normalized_VW = normalized_VW.reset_index(drop=True)

	normalized_VW = normalized_VW.sort_values(by = ['VW'], ascending = [False]).reset_index(drop=True)
	max_VW =  max(normalized_VW['VW'])
	normalized_VW['VW'] = normalized_VW['VW']/max_VW
	normalized_VW = normalized_VW

	return normalized_VW