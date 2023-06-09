import ROOT
import uproot
import pandas as pd
import awkward as ak

def eventHistoryDf(filepath,treeName='eventHistory', keysDict=[], eventWeight = False):

    """
    Converts a eventHistory .root file to a pandas dataframe 

    Args:
        filepath : location of the root file 
        keys_dict (list, optional): keys (branches) of root files. Defaults to [] for all branches.
        eventWeight (bool, optional): True only allows for recorded events. Defaults to False.

    Returns:
        dataframe: All recorded events. 
    """
    file = uproot.open(filepath)
    eventHistory = file[treeName]

    if not keysDict:
        keys = eventHistory.keys()
    else:
        keys = keysDict
    df_list = []
    for i in keys:
        awks = eventHistory[i].array()
        df_loc = ak.to_dataframe(awks)
        if eventWeight:
            df_loc = df_loc[df_loc["eventWeight"]== 50]
        df_loc['loc'] = i
        df_list.append(df_loc)
    df = pd.concat(df_list)
    return df 