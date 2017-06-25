# Subject to BSD-3-Clause
# Copyright (c) 2016-2017 Filippo Ranalli, Forest Flager, Pratyush Havelia
# See LICENSE.txt and NOTICE.txt

import numpy as np 
import pandas as pd

# must be extended to 3D continuity
def prioritize(nodeToElementMap,parameters):
    Member_A = []
    Member_B = []
    sizing_continuity = eval(parameters['Values'][9])
    if sizing_continuity:
        for node in nodeToElementMap.keys():
            elements = nodeToElementMap[node]
            columns = []
            beams = []
            braces = []
            for element in elements:
                if element[:2] == 'CO':
                    columns.append(element)
                elif element[:2] == 'BE':
                    beams.append(element)
                elif element[:2] == 'BR':
                    braces.append(element)

            if len(columns) == 2:
                Member_A.append(columns[0])
                Member_B.append(columns[1])

            elif len(columns) <= 1 and len(beams) == 2:
                Member_A.append(beams[0])
                Member_B.append(beams[1])

            elif len(columns) <= 1 and len(beams) <= 1 and len(braces) == 2:
                Member_A.append(braces[0])
                Member_B.append(braces[1])


    continuity = pd.DataFrame({'Member_A':Member_A,'Member_B':Member_B})
    return continuity





