import pandas as pd

# S = pd.read_csv("Square.csv",encoding='utf-8')
# W = pd.read_csv("W.csv",encoding='utf-8')
# R = pd.read_csv("Round.csv",encoding='utf-8')
SAPIMEMBER = pd.read_csv("SAP_I_Member.csv",encoding='utf-8')

# S.to_csv("Square.txt", index = False,encoding='utf-8')
# W.to_csv("W.txt", index = False,encoding='utf-8')
# R.to_csv("Round.txt", index = False,encoding='utf-8')
SAPIMEMBER.to_csv("SAP_I_Member.txt",index=False,encoding='utf-8')

