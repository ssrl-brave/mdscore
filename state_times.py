from itertools import groupby

import numpy as np

def get_states(clust_ids, dt=5):
  state = None
  t = 0
  states = []
  for i,l in enumerate(clust_ids[:-1]):
    if any ([x==l for x in clust_ids[i+1:i+1+dt]]):
      t += 1
      state = l
    else:
      if state is not None and t > dt:
        states.append((int(state), t))
      state = None
      t = 0
  return states



lab = np.load("cluster_labels_457x15000.npy")
all_state_times = {}
# i is the system index
for i in range(lab.shape[0]):
    states = []
    for L in (lab[i,:5000], lab[i,5000:10000], lab[i,-5000:]):
        states_L = get_states(L, dt=5) # check if it comes back within 5 frames
        states += states_L
    
    key = lambda x: x[0]
    gb = groupby(sorted(states, key=key), key=key)
    
    state_times = [(k,[t for _,t in list(v)]) for k,v in gb]
    state_times = sorted(state_times, key=lambda x: sum(x[1]))[::-1]
    all_state_times[i] = state_times
    
    if state_times:
      top_s = state_times[0][0]
      top_s_t = sum(state_times[0][1])
      print(f"done with system {i}. Top state time=state{top_s} t={top_s_t} frames")
    else:
      print(f"not state times for system {i+1}")


sum_state_t = []
for i in all_state_times:
    state_times = all_state_times[i]
    if state_times:
        t = sum(state_times[0][1])
        state = state_times[0][0]
    else:
        t = 0
        state = -1
    print(f"System {i}: top residence time cluster `{state}` was occupied for {t} frames")

