# Subject to BSD-3-Clause
# Copyright (c) 2016-2017 Filippo Ranalli, Forest Flager, Pratyush Havelia
# See LICENSE.txt and NOTICE.txt

import os
import sys
import comtypes.client
import numpy as np
import pandas as pd
import time

class SAPAnalysis(object):
    def __init__(self, path, logging=True):
        self.path = path
        self.logging = logging

    def initialize(self):
        #set the following flag to True to attach to an existing instance of the program
        #otherwise a new instance of the program will be started
        writeToFile = True

        #set the following flag to True to manually specify the path to SAP2000.exe
        #this allows for a connection to a version of SAP2000 other than the latest installation
        #otherwise the latest installed version of SAP2000 will be launched
        SpecifyPath = True
        #if the above flag is set to True, specify the path to SAP2000 below
        ProgramPath = 'C:/Program Files/Computers and Structures/SAP2000 18/SAP2000.exe'


        #full path to the model
        #set it to the desired path of your model

        try:
            #get the active SapObject
            attachToInstance = True
            mySapObject = comtypes.client.GetActiveObject("CSI.SAP2000.API.SapObject")
        except (OSError, comtypes.COMError):
            attachToInstance = False
            print("No running instance of the program found - starting new instance")
            pass
        if not attachToInstance:
            #create API helper object
            helper = comtypes.client.CreateObject('SAP2000v18.Helper')
            helper = helper.QueryInterface(comtypes.gen.SAP2000v18.cHelper)
            mySapObject = helper.CreateObject(ProgramPath)
            mySapObject.ApplicationStart()

        self.SM = mySapObject.SapModel
        self.SM.File.OpenFile(self.path)
        self.SM.SetPresentUnits(9) #unit: N,mm,C
        self.SM.SetModelIsLocked(False)

    def run(self,sapImember,sapInode,sizingLoop=False):
      if sizingLoop:
        print('\t ANALYSIS')
      else: 
        print('ANALYSIS')

      SM = self.SM
      self.SM.SetPresentUnits(9)
      SM.SetModelIsLocked(False)

      # 'LOAD IN THE NODE AND MEMBER INFO'
      # node_info = np.loadtxt('SAP_I_Node.txt', dtype=str, comments='#', delimiter=',', converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0)
      # member_info = np.loadtxt('SAP_I_Member.txt', dtype=str, comments='#', delimiter=',', converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0)

      'GENERATE THE NODES'
      nodeList = []
      groundNodeList = []
      for i in range(len(sapInode['node_ID'])):
        x = float(sapInode.get_value(i, 'X'))
        y = float(sapInode.get_value(i, 'Y'))
        z = float(sapInode.get_value(i, 'Z'))
        node_ID = str(sapInode.get_value(i, 'node_ID'))
        U1 = sapInode.get_value(i, 'U1')
        U2 = sapInode.get_value(i, 'U2')
        U3 = sapInode.get_value(i, 'U3')
        R1 = sapInode.get_value(i, 'R1')
        R2 = sapInode.get_value(i, 'R2')
        R3 = sapInode.get_value(i, 'R3')

        SM.PointObj.AddCartesian(x,y,z,node_ID,node_ID)
        restraints = [U1, U2, U3, R1, R2, R3]
        SM.PointObj.SetRestraint(node_ID, restraints)
        nodeList.append(node_ID)
        SM.PointObj.SetSpecialPoint(node_ID, True)
        if sapInode.get_value(i, 'ground?'):
          groundNodeList.append(node_ID)
      

      'GENERATE THE MEMBERS'
      memberList = []
      memberStart = []
      memberEnd = []
      for i in range(len(sapImember['member_ID'])):
        size = sapImember.get_value(i,'size')
        member_ID = sapImember.get_value(i,'member_ID')
        if size == '0':
          SM.FrameObj.Delete(member_ID)
        else:
          start = sapImember.get_value(i,'start_node')
          end = sapImember.get_value(i,'end_node')
          orientation = float(sapImember.get_value(i,'orientation'))
          SM.FrameObj.AddByPoint(start, end, '', size, member_ID)
          SM.FrameObj.SetLocalAxes(member_ID, orientation)
          SM.FrameObj.SetSection(member_ID, size)

          startValue = [0,0,0,0,0,0]
          endValue = [0,0,0,0,0,0]
          U1_I = sapImember.get_value(i,'PI')
          U2_I = sapImember.get_value(i,'V2I')
          U3_I = sapImember.get_value(i,'V3I')
          R1_I = sapImember.get_value(i,'TI')
          R2_I = sapImember.get_value(i,'M2I')
          R3_I = sapImember.get_value(i,'M3I')
          U1_J = sapImember.get_value(i,'PJ')
          U2_J = sapImember.get_value(i,'V2J')
          U3_J = sapImember.get_value(i,'V3J')
          R1_J = sapImember.get_value(i,'TJ')
          R2_J = sapImember.get_value(i,'M2J')
          R3_J = sapImember.get_value(i,'M3J')
          ii = [U1_I, U2_I, U3_I, R1_I, R2_I, R3_I]
          jj = [U1_J, U2_J, U3_J, R1_J, R2_J, R3_J]
          SM.FrameObj.SetReleases(member_ID, ii, jj, startValue, endValue)
          memberList.append(member_ID)
          memberStart.append(start)
          memberEnd.append(end)

      'RUN THE ANALYSIS'
      SM.File.Save(self.path)
      SM.Analyze.RunAnalysis()
      SM.DesignSteel.StartDesign()

      'RETREIVE THE MODAL FREQUENCY'
      SM.Results.Setup.DeselectAllCasesAndCombosForOutput()
      SM.Results.Setup.SetCaseSelectedForOutput("MODAL")

      results= SM.Results.ModalPeriod()
      numberResults = results[0]

      while numberResults == 0: #Temporary fix to avoid bug of empty result list
        print ('FAILURE TO RETREIVE RESULTS - RERUNNING')
        SM.Analyze.RunAnalysis()
        SM.DesignSteel.StartDesign()
        SM.Results.Setup.DeselectAllCasesAndCombosForOutput()
        SM.Results.Setup.SetCaseSelectedForOutput("MODAL")
        results= SM.Results.ModalPeriod()
        numberResults = results[0]
      numberResults = results[0]
      frequency = results[5][0]
    
      data = {'Mode':numberResults,'Frequency':frequency}
      modal_info = pd.DataFrame(data,index=[0]).reindex(columns=['Mode','Frequency'])

      'RETREIVE ALL COMBINATIONS'
      names = SM.RespCombo.GetNameList()
      numberCombos = names[0]
      combos = names[1]

      'RETRIEVE VIRTUAL WORK INFORMATION'
      virtualLoadCases = ["P02_z", "P03_z", "P05_z", "P02_x_pos", "P02_x_neg", "P03_x_pos", "P03_x_neg", "P05_x_pos", "P05_x_neg"]
      member_ID = []
      load_case = []
      P = []
      V2 = []
      M3_start = []
      M3_end = []
      for i in range(len(memberList)):
        for j in range(len(virtualLoadCases)):
          SM.Results.Setup.DeselectAllCasesAndCombosForOutput()
          SM.Results.Setup.SetCaseSelectedForOutput(virtualLoadCases[j])
          results = SM.Results.FrameForce(memberList[i], 0)
          member_ID.append(memberList[i])
          elmSta = results[4]
          load_case.append(results[5][0])
          P.append(results[8][0])
          V2.append(results[9][0])
          # V3 = results[10]
          # T.append(results[11][0])
          # M2 = results[12]
          M3_start.append(results[13][0])
          M3_end.append(results[13][len(elmSta)-1])

        for k in range(numberCombos):
          SM.Results.Setup.DeselectAllCasesAndCombosForOutput()
          SM.Results.Setup.SetComboSelectedForOutput(combos[k])
          results = SM.Results.FrameForce(memberList[i], 0)
          member_ID.append(memberList[i])
          elmSta = results[4]
          load_case.append(results[5][0])
          P.append(results[8][0])
          V2.append(results[9][0])
          # V3 = results[10]
          # T.append(results[11][0])
          # M2 = results[12]
          M3_start.append(results[13][0])
          M3_end.append(results[13][len(elmSta)-1])

      data = {'member_ID':member_ID,'Load_case':load_case,'P':P,'V':V2,'Mx_start':M3_start,'Mx_end':M3_end}
      virtual_info = pd.DataFrame(data).reindex(columns=['member_ID','Load_case','P','V','Mx_start','Mx_end'])

      'RETRIEVE MEMBER FORCES FOR ALL LOAD CASES'
      member_ID = []
      maxPs = []
      maxM2s = []
      maxM3s = []
      combonames = []
      ratios = []

      for i in range(len(memberList)):
        member_ID.append(memberList[i])
        results = SM.DesignSteel.GetSummaryResults(memberList[i])
        comboname = results[5][0]
        ratio = results[2][0]
        if comboname != None:
          SM.Results.Setup.DeselectAllCasesAndCombosForOutput()
          SM.Results.Setup.SetComboSelectedForOutput(comboname)
          results = SM.Results.FrameForce(memberList[i], 0)
          maxP = np.max(results[8])
          maxM2 = np.max(results[12])
          maxM3 = np.max(results[13])
        else:
          maxP = 0
          maxM2 = 0
          maxM3 = 0
          comboname = "ZERO_FORCE"
          ratio = 0
        maxPs.append(maxP)
        maxM2s.append(maxM2)
        maxM3s.append(maxM3)
        combonames.append(comboname)
        ratios.append(ratio)

      data = {'member_ID':member_ID,'Start_Node':memberStart,'End_Node':memberEnd,'Governing_Load_Case':combonames,'P':maxPs,'Mx':maxM2s,'My':maxM3s,'DC_Ratio':ratios}
      member_forces = pd.DataFrame(data).reindex(columns=['member_ID','Start_Node','End_Node','Governing_Load_Case','P','Mx','My','DC_Ratio'])
        

      'RETRIEVE THE NODE REACTIONS'
      F1 = []
      F2 = []
      F3 = []
      node_ID = []
      combonames = []
      for i in range(len(groundNodeList)):
        for j in range(numberCombos):
          SM.Results.Setup.DeselectAllCasesAndCombosForOutput()
          SM.Results.Setup.SetComboSelectedForOutput(combos[j])
          results = SM.Results.JointReact(groundNodeList[i], 0)
          if results[0] != 0:
            node_ID.append(groundNodeList[i])
            combonames.append(combos[j])
            F1.append(results[6][0])
            F2.append(results[7][0])
            F3.append(results[8][0])
            # M1 = results[9]
            # M2 = results[10]
            # M3 = results[11]
      data = {'Node_ID':node_ID,'Load_Combination':combonames,'vertical_reaction_at_Column_Base':F3,'horizontal_reaction_x_at_Column_Base':F1,'horizontal_reaction_y_at_Column_Base':F2}
      node_forces = pd.DataFrame(data).reindex(columns=['Node_ID','Load_Combination','vertical_reaction_at_Column_Base','horizontal_reaction_x_at_Column_Base','horizontal_reaction_y_at_Column_Base'])

      'RETRIEVE THE NODE DISPLACEMENTS'
      node_ID = []
      combonames = []
      U1 = []
      U2 = []
      U3 = []
      R1 = []
      R2 = []
      R3 = []
      for i in range(len(nodeList)):
        for j in range(numberCombos):
          SM.Results.Setup.DeselectAllCasesAndCombosForOutput()
          SM.Results.Setup.SetComboSelectedForOutput(combos[j])
          results = SM.Results.JointDispl(nodeList[i], 0)
          if results[0] != 0:
            node_ID.append(nodeList[i])
            combonames.append(combos[j])
            U1.append(results[6][0])
            U2.append(results[7][0])
            U3.append(results[8][0])
            R1.append(results[9][0])
            R2.append(results[10][0])
            R3.append(results[11][0])

      data = {'Node_ID':node_ID,'Load_Combination':combonames,'u1':U1,'u2':U2,'u3':U3,'r1':R1,'r2':R2,'r3':R3}
      node_displacements = pd.DataFrame(data).reindex(columns=['Node_ID','Load_Combination','u1','u2','u3','r1','r2','r3'])

      if self.logging:
        modal_info.to_csv('SAP_O_ModeFrequency.txt',index=False)
        virtual_info.to_csv('SAP_O_VirtualForce.txt',index=False)
        member_forces.to_csv('SAP_O_MemberForce.txt',index=False)
        node_forces.to_csv('SAP_O_NodeForce.txt',index=False)
        node_displacements.to_csv('SAP_O_NodeDisplacement.txt',index=False)
      return modal_info, virtual_info, member_forces, node_forces, node_displacements