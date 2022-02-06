import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

# Harp 0
# check all the masses in combination with the above stuff
gromos_atp = pd.DataFrame(
    {'name': ['O', 'OA', 'N', 'C', 'CH1', 
            'CH2', 'CH3', 'CH2r', 'NT', 'S',
            'NR', 'OM', 'NE', 'NL', 'NZ'],
     'at.num': [8, 8, 7, 6, 6, 6, 6, 6, 7, 16, 7, 8, 7, 7, 7],
     'c12': [1e-06, 1.505529e-06, 2.319529e-06, 4.937284e-06, 9.70225e-05, # CH1
            3.3965584e-05, 2.6646244e-05, 2.8058209e-05, 5.0625e-06, 1.3075456e-05,
            3.389281e-06, 7.4149321e-07, 2.319529e-06, 2.319529e-06, 2.319529e-06]
     }
)
gromos_atp.to_dict()
gromos_atp.set_index('name', inplace=True)

gromos_atp_c6 = pd.DataFrame(
    {'name': ['O', 'CH1', 'CH2', 'CH3'],
     'c6': [0.0022619536, 0.00606841, 0.0074684164, 0.0096138025]
    }
)
gromos_atp_c6.to_dict()
gromos_atp_c6.set_index('name', inplace=True)

# TODO make it more general (without the residue number)
# native reweight for TTR and ABeta. This dictionary will rename the amber topology to gromos topology
#gro_to_amb_dict = {'OC1_11' : 'O1_11', 'OC2_11':'O2_11'}
gro_to_amb_dict = {'OT1_42' : 'O1_42', 'OT2_42':'O2_42'}
