# Subject to BSD-3-Clause
# Copyright (c) 2016-2017 Filippo Ranalli, Forest Flager, Pratyush Havelia
# See LICENSE.txt and NOTICE.txt

def size(continuity,virtual_forces,sapImember,sapInode,geomData,node_displacements,mem_force,frequency,parameters,sec_info_W,sec_info_Round,sec_info_Square, logging=False):
    print('\t SIZING')
    import numpy as np
    import pandas as pd
    from detailing_class import detail
    from virtual_work import calculate_virtual_work
    import sys

    def read_member_data(geomData,sapImember,mem_force):
        remaining_members = list(mem_force['member_ID'])
        sizes = pd.DataFrame({'member_ID':sapImember['member_ID'], 'size': sapImember['size']})
        for member in geomData['member_ID']:
            if member not in remaining_members:
                geomData = geomData[geomData['member_ID'] != member]
                geomData = geomData.reset_index(drop = True)
                sizes = sizes[sizes['member_ID'] != member]
                size = sizes.reset_index(drop=True)
        # print mem_force
        # print geomData
        mem_info=pd.merge(mem_force,geomData).drop_duplicates()
        mem_info = mem_info.drop('End_Node',1)
        mem_info = mem_info.drop('end_node',1)
        mem_info = mem_info.drop('Start_Node',1)
        mem_info = mem_info.drop('start_node',1)
        mem_info = mem_info.drop('radial_tab',1)
        size_column = sizes['size'].reset_index(drop=True)
        mem_info= mem_info.join(size_column,how='outer')
        return mem_info
    def check_feasibility(mem_id,sec_info,sec,section_type,mem_info):
        E = 200000.0
        K = 1.2
        l = mem_info[mem_info['member_ID'] == mem_id].iloc[0]['member_length']

        pr = float(mem_info[mem_info['member_ID'] == mem_id].iloc[0]['P']) # take care of sign matching in pc equation
        mrx = abs(float(mem_info[mem_info['member_ID'] == mem_id].iloc[0]['My']))
        mry = abs(float(mem_info[mem_info['member_ID'] == mem_id].iloc[0]['Mx']))

        # if section_type == 'Square':
        #   Fy = 315.0
     #        A=sec_info["Square"].A[sec]
     #        I=sec_info["Square"].I[sec]
     #        Z=sec_info["Square"].Z[sec]
     #        B=sec_info["Square"].B[sec]
     #        S=sec_info["Square"].S[sec]
     #        b_t=sec_info["Square"].b_t[sec]
     #        r = sec_info["Square"].r[sec]
     #        t = sec_info["Square"].t[sec]
     #        h_t = sec_info["Square"].h_t[sec]
     #        b = sec_info["Square"].b[sec]
     #        # Pc
     #        Q = 1.0
     #        if sec_info["Square"].C_s[sec] == "S":
     #            be=min(b,1.92*t*np.sqrt(E/Fy)*(1 - 0.38/b_t * np.sqrt(E/Fy)))
     #            Be=be+2*t
     #            Ae=Be*t*2 + (Be-(2*t))*t*2
     #            Q= min(1,A/Ae)
                
     #        fe=np.pi**2*E/(K*l/r)**2
     #        if (Q*Fy/fe)<=2.25:
     #            fcr=Q*Fy*[0.658]**(Q*Fy/fe)
     #        else:
     #            fcr=0.877*fe
     #        pc=(fcr*A)*0.9

     #        # Mcx & Mcy
     #        if sec_info["Square"].FF_s[sec] == "C":     # Compact Flange
     #            mcx = mcy = 0.9*Z*Fy
     #        elif sec_info["Square"].FF_s[sec] == "NC":  # Non-Compact Flange
     #            mp = Z*Fy
     #            mcx = mcy = 0.9* min(mp, mp-(mp-Fy*S)*(3.57*b_t*np.sqrt(Fy/E)-4.0))
     #        elif sec_info["Square"].FF_s[sec] == "S":   # Slender Flange
     #            be = min(b,1.92*t*np.sqrt(E/Fy)*(1 - 0.38/b_t * np.sqrt(E/Fy)))
     #            Be=be+2*t
     #            Se = ((Be**4)/12.0 - (be**4)/12.0)*2/Be
     #            mcx = mcy = 0.9 *Fy*Se
     #        if sec_info["Square"].FW_s[sec] == "NC":    # Check Web after flange
     #            mp = Z*Fy
     #            mcx = mcy = 0.9 * min(mcx,mp-(mp-Fy*S)*(0.305*h_t*np.sqrt(Fy/E)-0.738))

        # if section_type == 'Round':
        #   Fy = 315.0
     #        A=sec_info["Round"].A[sec]
     #        I=sec_info["Round"].I[sec]
     #        Z=sec_info["Round"].Z[sec]
     #        D=sec_info["Round"].D[sec]
     #        S=sec_info["Round"].S[sec]
     #        d_t=sec_info["Round"].d_t[sec]
     #        r = sec_info["Round"].r[sec]
        
     #        # Pc
     #        Q = 1.0
     #        if sec_info["Round"].C_s[sec] == "S":
     #            Q = 0.038*E/(Fy*d_t) + 2/3
                
     #        fe=np.pi**2*E/(K*l/r)**2
     #        if (Q*Fy/fe)<=2.25:
     #            fcr=Q*Fy*[0.658]**(Q*Fy/fe)
     #        else:
     #            fcr=0.877*fe
     #        pc=(fcr*A)*0.9
        
     #        # Mcx & Mcy
     #        if sec_info["Round"].F_s[sec] == "C":
     #            mcx = mcy = 0.9*Z*Fy
     #        elif sec_info["Round"].F_s[sec] == "NC":
     #            mcx = mcy = 0.9* min(Z*Fy, (0.021*E/d_t + Fy)*S)      
     #        elif sec_info["Round"].F_s[sec] == "S":
     #            mcx = mcy = 0.9* min (Z*Fy, (0.021*E/d_t + Fy)*S, 0.33*E*S/d_t)
                
     #        # V        
     #        Vd = 0.9 * 0.78*E/(d_t**(3/2)) * A/2


        if section_type == 'W':
            Fy = 345.0
            A=sec_info["W"].A[sec]
            Ix=sec_info["W"].Ix[sec]
            Iy=sec_info["W"].Iy[sec]
            Zx=sec_info["W"].Zx[sec]
            Zy=sec_info["W"].Zy[sec]
            D=sec_info["W"].d[sec]
            Sx=sec_info["W"].Sx[sec]
            Sy=sec_info["W"].Sy[sec]
            d=sec_info["W"].d[sec]
            b_t=sec_info["W"].b_t[sec] #bf/2tf (half the flange width over the flange thickness)
            tf=sec_info["W"].tf[sec]
            tw=sec_info["W"].tw[sec]
            bf=sec_info["W"].bf[sec]
            h_tw=sec_info["W"].h_tw[sec]
            J=sec_info["W"].J[sec]
            rts=sec_info["W"].rts[sec]
            rx=sec_info["W"].rx[sec]
            ry=sec_info["W"].ry[sec]        
        
            #Pc
            Q = 1.0
            if sec_info["W"].CF_s[sec] == "NS":
                if b_t <= 0.56*np.sqrt(E/Fy):
                    Q = 1.0
                elif b_t > 0.56*np.sqrt(E/Fy) and b_t < 1.03*np.sqrt(E/Fy):
                    Q = 1.415 - 0.74*b_t*np.sqrt(Fy/E)
                else:
                    Q = 0.69*E/(Fy*b_t**2)
                    
            fe=np.pi**2*E/(K*l/ry)**2
            if (Q*Fy/fe)<=2.25:
                fcr=Q*Fy*[0.658]**(Q*Fy/fe)
            else:
                fcr=0.877*fe

            if pr < 0:
                pc= -(fcr*A)*0.9 # Compression strength (minus sign eliminates sign of pr)
            else:
                pc=(A*Fy)*0.9 # Tension strength
        
            # Mcx - strong axis
            Cb = 1.0
            Lp = 1.76*ry*np.sqrt(E/Fy)
            h0 = d-tf
            Lr = 1.95*rts*E/(0.7*Fy)*np.sqrt( J/(Sx*h0) + np.sqrt( (J/(Sx*h0))**2 + 6.76*(0.7*Fy/E)**2) )
            Sxc = Sxt = 4*bf*(tf**3)/(3*d) + 2*bf*tf*d - 2*bf*(tf**2)
            lamda = bf/(2*tf)
            lamda_pf = 0.38*np.sqrt(E/Fy)
            lamda_rf = np.sqrt(E/Fy)
            
            # AISC Section F2 - compact flanges and web
            if sec_info["W"].FF_s[sec] == "C" and sec_info["W"].FW_s[sec] == "C":       
                if (l <= Lp):
                    mcx = 0.9*Fy*Zx
                elif (l > Lp and l <= Lr):
                    mp = Fy*Zx
                    mcx = 0.9* min(mp, Cb*(mp - (mp-0.7*Fy*Sx)*(l-Lp)/(Lr-Lp)))
                else:
                    mp = Fy*Zx
                    mcx = 0.9* min(mp, Sx* Cb*np.pi**2*E/((l/rts)**2) * np.sqrt(1+ 0.078*J/(Sx*h0)*(l/rts)**2))
            
            # AISC Section F3 - non-compact or slender flanges, and compact web
            elif (sec_info["W"].FF_s[sec] == "NC" or sec_info["W"].FF_s[sec] == "S") and sec_info["W"].FW_s[sec] == "C":
                if (sec_info["W"].FF_s[sec] == "NC"):
                    mp = Fy*Zx
                    mcx = 0.9* (mp-(mp-0.7*Fy*Sx)*(lamda-lamda_pf)/(lamda_rf-lamda_pf))       
                elif (sec_info["W"].FF_s[sec] == "S"):
                    kc = min(max(4/np.sqrt(h_tw), 0.35),0.76)
                    mcx = 0.9* 0.9*E*kc*Sx/(lamda**2)
            
            # AISC Section F4 - non-compact web
            elif sec_info["W"].FW_s[sec] == "NC":
                return
                # No such case in our design - will handle if needed
            
            # AISC Section F5 - slender web
            elif sec_info["W"].FW_s[sec] == "S":
                return
                # No such case in our design - will handle if needed
                
            # Mcy - weak axis
            if sec_info["W"].FW_s[sec] == "C":
                mcy = 0.9* min(Fy*Zy, 1.6*Fy*Sy)
            elif sec_info["W"].FW_s[sec] == "NC":
                mp = min(Fy*Zy, 1.6*Fy*Sy)
                mcy = 0.9 * (mp - (mp - 0.7*Fy*Sy)*(lamda-lamda_pf)/(lamda_rf-lamda_pf))
            elif sec_info["W"].FW_s[sec] == "S":
                mcy = 0.9* 0.69*E/(b_t**2)*Sy
                        
            # V
            Cv = 1
            kv = 5
            if h_tw > 1.1*np.sqrt(kv*E/Fy) and h_tw <= 1.37*np.sqrt(kv*E/Fy):
                Cv = 1.1*np.sqrt(kv*E/Fy)/h_tw
            elif h_tw > 1.37*np.sqrt(kv*E/Fy):
                Cv = 1.51*kv*E/((h_tw**2)*Fy)
                
            Vd = 0.9* 0.6*Fy*(d-2*tf)*tw*Cv

        # strength check
        if (pr/pc)>=0.2:
            util=pr/pc+8.0/9.0*(abs(mrx/mcx)+abs(mry/mcy)) 
        else:
            util=pr/(2.0*pc)+(abs(mrx/mcx)+abs(mry/mcy))

        # if mem_id == 'BE07':
        #   print sec, util
        #   print pr,mrx,mry
        #   print pc,mcx,mcy
        #   assert False

        if util<=0.95:
            return 1,pc,mcx,mcy,util
        else:
            return 0,pc,mcx,mcy,util
    def read_section_data(sec_info_W,sec_info_Round,sec_info_Square):
        # sec_info_W = pd.read_csv('W.txt',sep=',',encoding='utf-8')
        # sec_info_Round = pd.read_csv('Round.txt',sep=',',encoding='utf-8')
        # sec_info_Square = pd.read_csv('Square.txt',sep=',',encoding='utf-8')
        
        Type=(sec_info_W["Type"])
        A=np.array(sec_info_W.A)
        tnom=np.array(sec_info_W['tnom'])
        W=np.array(sec_info_W['W'])
        Ix=np.array(sec_info_W['Ix / 106'])*10**6
        Iy=np.array(sec_info_W['Iy / 106'])*10**6
        Zx=np.array(sec_info_W['Zx / 103'])*10**3
        Zy=np.array(sec_info_W['Zy / 103'])*10**3
        d = np.array(sec_info_W['d'])
        b_t = np.array(sec_info_W['b/t'])
        bf = np.array(sec_info_W['bf'])
        tf = np.array(sec_info_W['tf'])
        tw = np.array(sec_info_W['tw'])
        h_tw = np.array(sec_info_W['h/tw'])
        Sx = np.array(sec_info_W['Sx / 103'])*10**3
        Sy = np.array(sec_info_W['Sy / 103'])*10**3
        J = np.array(sec_info_W['J / 103'])*10**3
        rts = np.array(sec_info_W['rts'])
        rx = np.array(sec_info_W['rx'])
        ry = np.array(sec_info_W['ry'])
        C_flange = np.array(sec_info_W['Flange Compression'])
        C_web = np.array(sec_info_W['Web Compression'])
        F_flange = np.array(sec_info_W['Flange Flexure'])
        F_web = np.array(sec_info_W['Web Flexure'])
        sec_info_W=pd.DataFrame({'AISC_W_Data [metric]':sec_info_W["AISC_Manual_Label"],"Type":Type,'tnom':tnom,'W':W,'FW_s':F_web,'FF_s':F_flange,'CW_s':C_web,'CF_s':C_flange,'tw':tw,'tf':tf,'rx':rx, "ry":ry, 'd':d,'b_t':b_t,'bf':bf,'h_tw':h_tw,'Sx':Sx,'Sy':Sy,'J':J,'rts':rts,'Ix':Ix,'Iy':Iy,'Zx':Zx,'Zy':Zy,'A':A,'W':sec_info_W.W,'unit_price':sec_info_W.unit_price})
        
        Type=(sec_info_Round["Type"])
        A=np.array(sec_info_Round.A)
        tnom=np.array(sec_info_Round['tnom'])
        W=np.array(sec_info_Round['W'])
        I=np.array(sec_info_Round['Ix / 106'])*10**6
        Z=np.array(sec_info_Round['Zx / 103'])*10**3
        S = np.array(sec_info_Round['Sx / 103'])*10**3
        D = np.array(sec_info_Round['OD'])
        d_t = np.array(sec_info_Round['D/t'])
        r = np.array(sec_info_Round['rx'])
        C_slenderness = np.array(sec_info_Round['Compression'])
        F_slenderness = np.array(sec_info_Round['Flexure'])      
        sec_info_Round=pd.DataFrame({'AISC_Round_Data [metric]':sec_info_Round["AISC_Manual_Label"],"Type":Type,'tnom':tnom,'W':W,'F_s':F_slenderness,'C_s':C_slenderness,'r':r, 'd_t':d_t,'D':D,'I':I,'S':S,'Z':Z,'A':A,'W':sec_info_Round.W,'unit_price':sec_info_Round.unit_price})
        
        Type=(sec_info_Square["Type"])
        A=np.array(sec_info_Square.A)
        tnom=np.array(sec_info_Square['tnom'])
        B=np.array(sec_info_Square.B) #outer base width
        b=np.array(sec_info_Square.b) #inner base width
        t=np.array(sec_info_Square.tdes)
        r=np.array(sec_info_Square.rx)
        b_t=np.array(sec_info_Square['b/tdes'])
        h_t=np.array(sec_info_Square['h/tdes'])
        W=np.array(sec_info_Square['W'])
        I=np.array(sec_info_Square['Ix / 106'])*10**6
        Z=np.array(sec_info_Square['Zx / 103'])*10**3
        S = np.array(sec_info_Square['Sx / 103'])*10**3
        C_s = np.array(sec_info_Square['Compression'])
        F_flange = np.array(sec_info_Square['Flange Flexure'])
        F_web = np.array(sec_info_Square['Web Flexure'])
        sec_info_Square=pd.DataFrame({'AISC_Square_Data [metric]':sec_info_Square["AISC_Manual_Label [metric]"],"Type":Type,'tnom':tnom,'A':A,'b':b,'B':B,'W':W,'I':I,'S':S,'Z':Z,'r':r,'t':t,'b_t':b_t,'h_t':h_t,'FW_s':F_web,'FF_s':F_flange,'C_s':C_s,'unit_price':sec_info_Square.unit_price})

        sec_info= dict({"W":sec_info_W, "Round":sec_info_Round, "Square":sec_info_Square})
        sec_info_indexed= dict({"W":sec_info_W, "Round":sec_info_Round, "Square":sec_info_Square})
        # sec_info["W"]=sec_info_W.sort_values(['d','W'],ascending=True).set_index('AISC_W_Data [metric]')
        # sec_info["Round"]=sec_info_Round.sort_values(['D','W'],ascending=True).set_index('AISC_Round_Data [metric]')
        # sec_info["Square"]=sec_info_Square.sort_values(['B','W'],ascending=True).set_index('AISC_Square_Data [metric]')
        # sec_info["W"]=sec_info_W.sort_values(['A','d'],ascending=True).reset_index(drop = True)
        # sec_info["Round"]=sec_info_Round.sort_values(['A','D'],ascending=True).reset_index(drop = True)
        # sec_info["Square"]=sec_info_Square.sort_values(['A','B'],ascending=True).reset_index(drop = True)
        sec_info["W"]=sec_info_W.sort_values(['A','d'],ascending=True).set_index('AISC_W_Data [metric]')
        sec_info["Round"]=sec_info_Round.sort_values(['A','D'],ascending=True).set_index('AISC_Round_Data [metric]')
        sec_info["Square"]=sec_info_Square.sort_values(['A','B'],ascending=True).set_index('AISC_Square_Data [metric]')
        sec_info_indexed["W"]=sec_info_W.sort_values(['A','d'],ascending=True).reset_index(drop=True)
        sec_info_indexed["Round"]=sec_info_Round.sort_values(['A','D'],ascending=True).reset_index(drop=True)
        sec_info_indexed["Square"]=sec_info_Square.sort_values(['A','B'],ascending=True).reset_index(drop=True)

        #num_all_sec = dict({"W":sec_info_W.count().A, "Round":sec_info_Round.count().A, "Square":sec_info_Square.count().A})

        return sec_info,sec_info_indexed
    def size_for_strength(mem_info,sec_info):
        strength_design = pd.DataFrame({'member_ID': mem_info['member_ID'],'cross_section': mem_info['size']}).set_index('member_ID')
        member_list = mem_info['member_ID']
        strength_unsizable = False

        for member in member_list:
            #strength_design['cross_section'][member] = 'test'

            section_type = mem_info[mem_info['member_ID'] == member].iloc[0]['member_type']
            # if section_type == 'W':
            #     sections = sec_info['W'].sort_values(by=['A'],ascending=[True])
            # elif section_type == 'Round':
            #     sections = sec_info['Round'].sort_values(by=['A'],ascending=[True])
            # else:
            #     sections = sec_info['Square'].sort_values(by=['A'],ascending=[True])
            if section_type == 'W':
                sections = sec_info['W']
            elif section_type == 'Round':
                sections = sec_info['Round']
            else:
                sections = sec_info['Square']

            i = 0
            proposed_sec = sections.iloc[i].name

            num_sections = len(sections['A'])
            # while check_feasibility(member,sec_info,proposed_sec,section_type,mem_info)[0] != 1 or i == num_sections-1:
            while check_feasibility(member,sec_info,proposed_sec,section_type,mem_info)[0] != 1:
                i += 1
                if i == num_sections:
                    strength_unsizable = True
                    proposed_sec = sections.iloc[i-1].name
                    break
                proposed_sec = sections.iloc[i].name



            strength_design['cross_section'][member] = proposed_sec
        return strength_design,strength_unsizable
    def drift_check(max_drift,node_displacements,critical_nodes):
        u1 = {}
        u2 = {}
        u3 = {}
        for node in node_displacements['Node_ID']:
            if abs(node_displacements[node_displacements['Node_ID'] == node].iloc[0]['u1']) > max_drift[0] and node in critical_nodes:
                u1[node] = (node_displacements[node_displacements['Node_ID'] == node].iloc[0]['u1']/max_drift[0],node_displacements[node_displacements['Node_ID'] == node].iloc[0]['Load_Combination'])
            # if abs(node_displacements[node_displacements['Node_ID'] == node].iloc[0]['u2']) > max_drift[1] and node in critical_nodes:
            #     u2[node] = (abs(node_displacements[node_displacements['Node_ID'] == node].iloc[0]['u2'])/max_drift[1],node_displacements[node_displacements['Node_ID'] == node].iloc[0]['Load_Combination'])
            if abs(node_displacements[node_displacements['Node_ID'] == node].iloc[0]['u3']) > max_drift[2] and node in critical_nodes:
                u3[node] = (node_displacements[node_displacements['Node_ID'] == node].iloc[0]['u3']/max_drift[2],node_displacements[node_displacements['Node_ID'] == node].iloc[0]['Load_Combination'])
        
        governing_combos = []
        if len(u3.keys()) != 0:
            for node in u3.keys():
                if u3[node][0] < 0:
                    governing_combos.append([u3[node][1],node + '_z',abs(u3[node][0])])
        # if len(u2.keys()) != 0:
        #   for node in u2.keys():
        #       govening_combos.append([u2[node][1],node + '_y',u2[node][0]])
        if len(u1.keys()) != 0:
          for node in u1.keys():
            if u1[node][0] > 0:
                governing_combos.append([u1[node][1],node + '_x_pos',abs(u1[node][0])])
            else:
                governing_combos.append([u1[node][1],node + '_x_neg',abs(u1[node][0])])
        return governing_combos
    def envelope_upwards(design1,design2,sec_info,mem_info):
        info = sec_info['W']
        for i in range(0,len(design1['cross_section'])):
            ix1 = design1.index[i]
            section1 = design1['cross_section'].iloc[i]
            section2 = design2.get_value(ix1,'cross_section')
            order1 = info[info['AISC_W_Data [metric]'] == section1].iloc[0].name
            order2 = info[info['AISC_W_Data [metric]'] == section2].iloc[0].name
            if order2>order1:
                design1['cross_section'].iloc[i] = section2
        return design1
    def envelope_downwards(design1,design2,sec_info,mem_info):
        info = sec_info['W']
        for i in range(0,len(design1['cross_section'])):
            ix1 = design1.index[i]
            section1 = design1['cross_section'].iloc[i]
            section2 = design2.get_value(ix1,'cross_section')
            order1 = info[info['AISC_W_Data [metric]'] == section1].iloc[0].name
            order2 = info[info['AISC_W_Data [metric]'] == section2].iloc[0].name
            if order2<order1:
                design1['cross_section'].iloc[i] = section2
        return design1
    def size_for_stiffness(current_design,mem_info,sec_info,VW,drift_factor):
        stiffness_design = pd.DataFrame({'member_ID':mem_info['member_ID'], 'cross_section':mem_info['size']}).set_index('member_ID')
        info = sec_info['W']

        # print stiffness_design
        #print VW

        for i in range(0,len(VW)):

            member = VW.get_value(i,'member_ID')
            top_member_ratio = VW.get_value(i,'VW')
            current_section = stiffness_design['cross_section'][member]
            current_section_ix = info[info['AISC_W_Data [metric]'] == current_section].iloc[0].name

            #new_section_ix = np.ceil(current_section_ix * abs(drift_factor) *top_member_ratio)  #option 1

            A = 0.10
            B = 0.10
            new_section_ix = np.ceil(max(current_section_ix,1.0) * abs(drift_factor + B) * (top_member_ratio + A)) # option 2

            new_section_ix = max(new_section_ix,current_section_ix)
            new_section_ix = min(new_section_ix,len(info['A']) - 1)
            new_section_ix = int(new_section_ix)
            stiffness_design['cross_section'][member] = info.iloc[new_section_ix]['AISC_W_Data [metric]']
            #print 'member:',member, 'current section:', current_section, 'current_section_ix:',current_section_ix, 'new_section_ix:',new_section_ix, 'new section:', stiffness_design['cross_section'][member]

        #assert False
        return stiffness_design
    def envelope_strength_stiffness_designs(current_design,design1,design2,sec_info,mem_info):

        info = sec_info['W']
        for i in range(0,len(design1['cross_section'])):
            ix1 = design1.index[i]
            section1 = design1['cross_section'].iloc[i]
            section2 = design2.get_value(ix1,'cross_section')
            order1 = info[info['AISC_W_Data [metric]'] == section1].iloc[0].name
            order2 = info[info['AISC_W_Data [metric]'] == section2].iloc[0].name
            if order2>order1 and (design2.get_value(ix1,'cross_section') != current_design.get_value(ix1,'cross_section')):
                design1['cross_section'].iloc[i] = section2
        return design1
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

    def apply_continuity(continuity,sec_info,design):
        info = sec_info['W']
        for cont_list in continuity:
            max_order = 0
            max_section = 'W150X13'
            for member in cont_list:
                section = design['cross_section'][member]
                order = info[info['AISC_W_Data [metric]'] == section].iloc[0].name
                if order >= max_order:
                    max_section = section
                    max_order = order
            for member in cont_list:
                design['cross_section'][member] = max_section
        return design
    def computeMemberCost(design,sec_info,mem_info):
        info = sec_info['W']

        member_cost = pd.DataFrame({'member_ID':list(design.index), 'cost':design['cross_section']}).set_index('member_ID')
        for member in list(design.index):
            section = design['cross_section'][member]
            member_cost['cost'][member] = (info['unit_price'][section]/1000) * (mem_info['member_length'][member]/1000) * info['W'][section]
        return member_cost
    def memberStiffnessContribution(node_displacements,critical_nodes,max_drift,mem_info,sec_info,virtual_forces):
        relevant_combos = []
        for i in range(0,len(node_displacements['Node_ID'])):
            node = node_displacements.get_value(i,'Node_ID')
            if node in critical_nodes:
                load_combo = node_displacements.get_value(i,'Load_Combination')
                u1 = node_displacements.get_value(i,'u1') #x
                # u2 = abs(node_displacements.get_value(i,'u2')) #y
                u3 = node_displacements.get_value(i,'u3') #z

                if u1 >0:
                    relevant_combos.append([load_combo,node + '_x_pos',abs(u1)/max_drift[0]])
                else:
                    relevant_combos.append([load_combo,node + '_x_neg',abs(u1)/max_drift[0]])
                relevant_combos.append([load_combo,node + '_z',abs(u3)/max_drift[2]])

        # print (mem_info)
        all_contributions = np.zeros(len(mem_info['Governing_Load_Case']))
        for case in relevant_combos:
            real_case = case[0]
            virtual_case = case[1]
            VW = calculate_virtual_work(mem_info,sec_info,virtual_forces,real_case,virtual_case,verbose=False).sort_values(by='member_ID',ascending = True).reset_index(drop=True)
            # all_contributions.append(VW)
            importance_factor = case[2] # ratio of actual to allowable drift
            all_contributions += (np.array(VW['VW']) * importance_factor)

        regularization = 0.00001
        all_contributions = all_contributions + regularization
        all_contributions = all_contributions/np.max(all_contributions)
        average_contribution = pd.DataFrame({'member_ID':mem_info['member_ID'],'VW':list(all_contributions)})
        # average_contribution = average_contribution.set_index('member_ID')
        average_contribution = average_contribution

        return average_contribution

    def continuityConstributions(stiffness_contributions,continuity):
        for cont_list in continuity:
            average_contribution = 0
            length = len(cont_list)
            for member in cont_list:
                # average_contribution += stiffness_contributions['VW'][member]/length
                average_contribution += stiffness_contributions[stiffness_contributions['member_ID']==member].iloc[0]['VW'] / length
            for member in cont_list:
                # stiffness_contributions['VW'][member] = max(stiffness_contributions['VW'][member],average_contribution)
                stiffness_contributions['VW'][member] = max(stiffness_contributions[stiffness_contributions['member_ID']==member].iloc[0]['VW'],average_contribution)
        return stiffness_contributions



    '########## Main Executable Program ##########'
    # Assemble information
    critical_nodes = ['P02','P03','P05']
    max_drift = [100.0,100.0,100.0]
    sec_info,sec_info_indexed = read_section_data(sec_info_W,sec_info_Round,sec_info_Square)
    mem_info = read_member_data(geomData,sapImember,mem_force)
    current_design = pd.DataFrame({'member_ID':mem_info['member_ID'], 'cross_section':mem_info['size']}).set_index('member_ID')
    sizing_iteration = int(parameters['Values'][12])
    drift_factors = [0]

    if np.abs(frequency['Frequency'][0]) <= 0.01:
        print ('---- INFEASIBLE SIZING LOOP ----')
        parameters['Values'][13] = False # sizing feasible
        parameters['Values'][12] = 0 # iteration number
        return sapImember, -1, -1, current_design, -1, True, parameters

    # Size for strength
    strength_design,strength_unsizable = size_for_strength(mem_info,sec_info)
    if sizing_iteration == 0:
        proposed_design = strength_design
        sizing_convergence = False
    else:
        # Check drifts and drift-governing load combos
        # max_drift = [50.8,50.8,50.8]
        stiffness_cases = drift_check(max_drift,node_displacements,critical_nodes)
        # print (stiffness_cases)
        sizing_convergence = True
        if len(stiffness_cases) != 0:
            sizing_convergence = False
            stiffness_design = pd.DataFrame({'member_ID':mem_info['member_ID'], 'cross_section':mem_info['size']}).set_index('member_ID')
            for case in stiffness_cases:
                real_case = case[0]
                virtual_case = case[1]
                drift_factor = case[2]
                drift_factors.append(drift_factor)
                VW = calculate_virtual_work(mem_info,sec_info,virtual_forces,real_case,virtual_case)
                new_design = size_for_stiffness(current_design,mem_info,sec_info_indexed,VW,drift_factor)
                stiffness_design = envelope_upwards(stiffness_design,new_design,sec_info_indexed,mem_info)

            proposed_design = envelope_strength_stiffness_designs(current_design,strength_design,stiffness_design,sec_info_indexed,mem_info)
        else:
            proposed_design = envelope_upwards(strength_design,current_design,sec_info_indexed,mem_info)


    # Check for an unsizable iteration
    # if (max(drift_factors) > 1.5 and sizing_iteration > 4) or strength_unsizable == True or sizing_iteration > 10:
    if sizing_iteration > 10 or strength_unsizable == True or sizing_iteration > 25:
        print ('---- INFEASIBLE SIZING LOOP ----')
        parameters['Values'][13] = False # sizing feasible
        parameters['Values'][12] = 0 # iteration number
        # sizing_converged = open("sizing_converged.txt", "w")
        # sizing_converged.write(str(True))
        # sizing_converged.close()
        # sys.exit(0)
        return sapImember, -1, -1, current_design, -1, True, parameters

    # Apply continuity constraints
    sizing_continuity = eval(parameters['Values'][9])
    if sizing_continuity:
        # continuity = prioritize(nodeToElementMap)
        full_continuity = assembleFullContinuityList(continuity)
        # print ('APPLYING SIZING CONTINUITY')
        proposed_design = apply_continuity(full_continuity,sec_info_indexed,proposed_design)
    else:
        continuity = pd.DataFrame({'Member_A':[],'Member_B':[]})

    # Update SAP_I_Member
    sizes = proposed_design["cross_section"]
    memberIDs = proposed_design.index.tolist()
    for i in range(0, len(memberIDs)):
        sapImember.loc[sapImember["member_ID"] == memberIDs[i], "size"] = sizes[i]

    member_cost = computeMemberCost(proposed_design,sec_info,mem_info.set_index('member_ID'))
    stiffness_contributions = memberStiffnessContribution(node_displacements,critical_nodes,max_drift,mem_info,sec_info,virtual_forces)

    if sizing_convergence == True:
        parameters['Values'][12] = 0 # iteration number
        mem_info['size'] = proposed_design['cross_section']
        mem_info = mem_info.reset_index(drop = False)
        if sizing_continuity == True:
            stiffness_contributions = continuityConstributions(stiffness_contributions,full_continuity)
    else:
        parameters['Values'][12] = int(parameters['Values'][12]) + 1 # iteration number

    parameters['Values'][13] = True
    return sapImember, continuity, member_cost, proposed_design, stiffness_contributions, sizing_convergence, parameters

# size()
