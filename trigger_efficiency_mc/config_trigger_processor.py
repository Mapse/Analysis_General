''' 
Config filed used to help TriggerProcessor script

'''

# Path to the .coffea file
path_file = 'data/'
# Files name:


files = {'25_45' : 'Monte_Carlo_2018_Jpsi_25to45_Dstar_DPS_2018_13TeV.coffea',
         '45_65' : 'Monte_Carlo_2018_Jpsi_45to65_Dstar_DPS_2018_13TeV.coffea',
         '65_100' : 'Monte_Carlo_2018_Jpsi_65to100_Dstar_DPS_2018_13TeV.coffea',}

# pT list for efficiency
pt_val_list = [24, 24.4, 24.8, 25.2, 25.6, 26.,
                        27, 28, 29, 30, 32, 34, 36, 38, 40, 45, 50, 55, 60, 70, 80, 100]

# List of triggers to be applied
hlt_filter = ['HLT_Dimuon25_Jpsi']

# Csv file to save the variables
csv_file = 'trigger_mc_efficiency_2018.csv'

