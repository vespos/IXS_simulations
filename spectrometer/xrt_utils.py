import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append(r'C:\Users\espov\Documents\Python\xrt')


def plot_layout(bl, coord1='y', coord2='z'):
    for coord, idx in zip(['x','y','z'], [0,1,2]):
        if coord1==coord:
            coord1=idx
        if coord2==coord:
            coord2=idx
            
    oes = bl.oesDict
    coords = []
    labels = []
    for key in oes:
        if 'Screen' in key:
            continue
        else:
            coords.append(oes[key][0].center)
            labels.append(oes[key][0].name)
    last_screen = bl.screens[-1]
    coords.append(last_screen.center)
    labels.append('detector')
    
    coords = np.asarray(coords)
    print(labels)
    
    fig, ax = plt.subplots()
    ax.plot(coords[:,coord1], coords[:,coord2], '-o', markersize=20)
    for lab, c1, c2 in zip(labels, coords[:,coord1], coords[:,coord2]):
        ax.annotate(lab, (c1, c2))
    plt.show()