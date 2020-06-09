#%% #import packages

import xarray as xr
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#%% 
def get_indices():
    '''
    Funciton defines longitudes and latitudes for each part of the path.
    There are four parts: Southern Ocean, Equatorial Ocean, South America
    and Western Australia (so, eo, sa, wa). A random snapshot file is opened.
    Indices are extracted in a double loop from the ds.xu/ds.xt and 
    ds.yu/ds.yt files. Function returns these indices. 
    '''
    
    ds = xr.open_dataset('4deg_000.snapshot.nc')

    so_y = np.full(35,-46)
    so_x = np.array([148, 152, 156, 160, 164, 168, 172, 176, 180, 184, 188,
                         192, 196, 200, 204, 208, 212, 216, 220, 224, 228,
                         232, 236, 240, 244, 248, 252, 256, 260, 264, 268,
                         272, 276, 280, 284])
    i_so_y = []
    i_so_x = []
    for i in range (len(so_y)):
        for j in range (len(ds.yt)):
            if so_y[i] == ds.yt[j]:
                i_so_y.append(j)
    for i in range (len(so_y)):
        for j in range (len(ds.xu)):
            if so_x[i] == ds.xu[j]:
                i_so_x.append(j)
        
                
                
              
    eo_y = np.full(36,-2)
    eo_x = np.array([276, 272, 268, 264, 260,
                        256, 252, 248, 244, 240, 236, 232, 228, 224, 220, 216,
                         212, 208, 204, 200, 196, 192, 188, 184, 180, 176, 172,
                         168, 164, 160, 156, 152, 148, 144, 140, 136])
    i_eo_y = []
    i_eo_x = []
    for i in range (len(eo_y)):
        for j in range (len(ds.yt)):
            if eo_y[i] == ds.yt[j]:
                i_eo_y.append(j)
    for i in range (len(eo_y)):
        for j in range (len(ds.xu)):
            if eo_x[i] == ds.xu[j]:
                i_eo_x.append(j)
           
        
        
    sa_x = np.array([284, 284, 284, 284, 288, 288, 288, 284, 280, 280, 276])
    sa_y = np.array([-44, -40, -36, -32, -28, -24, -20, -16, -12, -8, -4])
    i_sa_y = []
    i_sa_x = []
    for i in range (len(sa_y)):
        for j in range (len(ds.yu)):
            if sa_y[i] == ds.yu[j]:
                i_sa_y.append(j)
    for i in range (len(sa_y)):
        for j in range (len(ds.xu)):
            if sa_x[i] == ds.xu[j]:
                i_sa_x.append(j)
             

            
    wa_x = np.array([132, 128, 128, 124, 120, 116, 112, 112, 112, 116, 120,
                         124, 128, 132, 136, 140])
    wa_y = np.array([0, -4, -8, -8, -12, -16, -20, -24, -28, -32, -36, -36,
                         -36, -36, -36, -40])
    i_wa_y = []
    i_wa_x = []
    for i in range (len(wa_y)):
        for j in range (len(ds.yu)):
            if wa_y[i] == ds.yu[j]:
                i_wa_y.append(j)
    for i in range (len(wa_y)):
        for j in range (len(ds.xu)):
            if wa_x[i] == ds.xu[j]:
                i_wa_x.append(j)

                
                
                
    return i_so_x, i_so_y, i_eo_x, i_eo_y, i_sa_x, i_sa_y, i_wa_x, i_wa_y


#%%
def Island_Rule(wind, i_so_x, i_so_y, i_eo_x, i_eo_y, i_sa_x, i_sa_y, i_wa_x, i_wa_y):
    ''' Function opens a nc snapshot file given by the wind input. NC snapshot 
    files should be named as 4deg_{'inser wind factor'}.snapshot.nc. The value of
    the throughflow from the simulation is extracted in the point 128°e, -8°s.
    Integration for each part is done in an integration loop, integrating the
    surface wind stress in the correct direction. Where the points are located
    on land a windstress of zero is added. All ds.dxu/ds.dxt values are the same,
    this also applies to all ds.dyu/ds.dyt values, so the first is just taken
    (ds.dxu[0], ds.dyu[0]). The sum of the part integrations is taken and divided
    by the constants as defined in the Island Rule. Function returns the Island Rule 
    integral given in Sverdrups, the simulated throughflow in Sverdrups, the
    contribution to the integral of each part of the path in Sverdrups and the 
    relative Error of the Island Rule compared to the simulated throughflow.
    '''

    #defining constants:
    f_2s = 2 * (2 * np.pi / (24 * 60 * 60)) * np.sin(-(1/90) * np.pi)
    f_44s = 2 * (2 * np.pi / (24 * 60 * 60)) * np.sin(-(23/90) * np.pi)
    rho = 1020
    
    #opens snapshot file for specific windstress
    ds = xr.open_dataset('4deg_{}.snapshot.nc'.format(wind))
    
    #Throughflow result of simulation is found (128°e, 8°s):
    psi = ds.psi[-1, 17, 31]

    
    
    
    #Integral along Southern Ocean: 
    SO = 0
    for i in range(len(i_so_y)):
        if i == 4 or i == 5 or i ==6: 
            SO += 0 
        else:
            SO += ds.surface_taux[-1, i_so_y[i], i_so_x[i]] * ds.dxu[0]

    #Integral along Equatorial Ocean: 
    EO = 0
    for i in range(len(i_eo_y)):
        EO += - ds.surface_taux[-1, i_eo_y[i], i_eo_x[i]] * ds.dxu[0]


    #Integral along South America: 
    SA = 0
    for i in range(len(i_sa_y)):
        if i < 3: 
            SA += ds.surface_tauy[-1, i_sa_y[i], i_sa_x[i]] * ds.dyu[0]
        if i == 3: 
            SA += (ds.surface_tauy[-1, i_sa_y[i], i_sa_x[i]] * ds.dyu[0]
                    + ds.surface_taux[-1, i_sa_y[i], i_sa_x[i]] * ds.dxu[0])
        if i == 4 or i == 5:
            SA += ds.surface_tauy[-1, i_sa_y[i], i_sa_x[i]] * ds.dyu[0]
        if i == 6 or i == 7:
            SA += (ds.surface_tauy[-1, i_sa_y[i], i_sa_x[i]] * ds.dyu[0] 
                    - ds.surface_taux[-1, i_sa_y[i], i_sa_x[i]] * ds.dxu[0])
        if i == 8: 
            SA += ds.surface_tauy[-1, i_sa_y[i], i_sa_x[i]] * ds.dyu[0]
        if i == 9: 
            SA += (ds.surface_tauy[-1, i_sa_y[i], i_sa_x[i]] * ds.dyu[0] 
                    - ds.surface_taux[-1, i_sa_y[i], i_sa_x[i]] * ds.dxu[0])
        if i == 10: 
            SA += ds.surface_tauy[-1, i_sa_y[i], i_sa_x[i]] * ds.dyu[0]            


    #Integral along western Australia:
    WA = 0
    for i in range (len(i_wa_y)):
        if i == 0: 
            WA += (- ds.surface_tauy[-1, i_wa_y[i], i_wa_x[i]] * ds.dyu[0] 
                    - ds.surface_taux[-1, i_wa_y[i], i_wa_x[i]] * ds.dxu[0])
        if i == 1:
            WA += - ds.surface_tauy[-1, i_wa_y[i], i_wa_x[i]] * ds.dyu[0]
        if i == 2:
            WA += - ds.surface_taux[-1, i_wa_y[i], i_wa_x[i]] * ds.dxu[0]
        if i == 3 or i == 4 or i == 5:
            WA += (- ds.surface_tauy[-1, i_wa_y[i], i_wa_x[i]] * ds.dyu[0] 
                    - ds.surface_taux[-1, i_sa_y[i], i_sa_x[i]] * ds.dxu[0])
        if i == 6 or i == 7:
            WA += - ds.surface_tauy[-1, i_wa_y[i], i_wa_x[i]] * ds.dyu[0]
        if i == 8 or i == 9: 
            WA += (- ds.surface_tauy[-1, i_wa_y[i], i_wa_x[i]] * ds.dyu[0]  
                    + ds.surface_taux[-1, i_wa_y[i], i_wa_x[i]] * ds.dxu[0])
        if i == 10 or i == 11 or i == 12 or i == 13: 
            WA += ds.surface_taux[-1, i_wa_y[i], i_wa_x[i]] * ds.dxu[0]
        if i == 14 or i == 15:
            WA += (- ds.surface_tauy[-1, i_wa_y[i], i_wa_x[i]] * ds.dyu[0] 
                    + ds.surface_taux[-1, i_wa_y[i], i_wa_x[i]] * ds.dxu[0])      
    
    tau = SO + EO + SA + WA
    IR = - tau / (rho * (f_2s - f_44s))

    
    SO_c = (- SO / (rho * (f_2s - f_44s))) 
    EO_c = (- EO / (rho * (f_2s - f_44s))) 
    SA_c = (- SA / (rho * (f_2s - f_44s))) 
    WA_c = (- WA / (rho * (f_2s - f_44s)))  
    
    RELERR = np.abs((1 - (IR / psi)) * 100)
    return IR / 10**6, psi / 10**6, SO_c / 10**6, EO_c / 10**6, SA_c / 10**6, \
                WA_c / 10**6, RELERR
#%%
def call_Island_Rule(winds):
    '''
    Function calls up the Island_Rule function for all winds in order to
    assemble arrays of Island Rules, simulated throughflow values, 
    contributions of paths and relative errors. Results are saved into
    a csv file.
    '''
    
    i_so_x, i_so_y, i_eo_x, i_eo_y, i_sa_x, i_sa_y, i_wa_x, i_wa_y = get_indices()
    
    Island_rules = np.zeros(len(winds))
    psis = np.zeros(len(winds))
    SO_cs = np.zeros(len(winds))
    EO_cs = np.zeros(len(winds))
    SA_cs = np.zeros(len(winds))
    WA_cs = np.zeros(len(winds))
    RELERRS = np.zeros(len(winds))
    for i in range(len(winds)):
        Island_rules[i], psis[i], SO_cs[i], EO_cs[i], SA_cs[i], WA_cs[i], RELERRS[i] = \
            Island_Rule(winds[i], i_so_x, i_so_y, i_eo_x, i_eo_y, i_sa_x, i_sa_y,
                        i_wa_x, i_wa_y)  


    Data = {'winds': winds, 'islandrules': Island_rules, 'psis': psis,
            'southernOcean': SO_cs, 'equatorialOcean': EO_cs,
           'SA': SA_cs, 'WA': WA_cs, 'RELERR': RELERRS}
    df = pd.DataFrame(Data)
    df.to_csv('solutions.csv')
    
    #return Island_rules, psis, SO_cs, EO_cs, SA_cs, WA_cs, RELERRS

#%%
if __name__ == '__main__':
    call_Island_Rule(['000', '0', '01', '015', '025', '05', '075',
                          '1', '125', '135', '145' ,'15', '2'])
