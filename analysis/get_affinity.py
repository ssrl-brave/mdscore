import json
import os
system=os.getcwd().split('/')[-1]
f=open('boltz_results_boltz_input/predictions/boltz_input/affinity_boltz_input.json','r')
a=json.load(f)
print(system,a['affinity_pred_value'])
