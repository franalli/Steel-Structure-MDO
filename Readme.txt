Before using or redistributing, refer to LICENSE.txt and NOTICE.txt for binding copyright information.

This code was developed at Stanford and used to publish the following technical report:
A Ground Structure Method to Minimize the Total Installed Cost of Steel Frame Structures (https://cife.stanford.edu/node/1493).

The main components are:
Analysis - Structural analysis in SAP2000
Prioritizing - Member continuity rules
Sizing - Member sizing per AISC code
Detailing - Connection mapping
Optimizing - Cost or weight-based subtractive ground structure optimization
MDO - Main executable wrapper file. 

Please note comtypes.client has a bug in __init__ that can be fixed by following www.stackoverflow.com directions.

Contact: franalli@stanford.edu