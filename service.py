#!/usr/bin/env python
# coding: utf-8

# In[115]:


from tkinter import *
from tkinter import filedialog
import myvariant
import pandas as pd
import numpy as np
import lightgbm

filename = None
#FUNCTIONS
def predict():
    #read input
    if filename != None:
        with open(filename, 'r') as file:
            lines = list(file)
    else:
        lines = text_query.get('1.0', 'end-1c').split('\n')
    
    #create queries    
    queries = []
    for line in lines:
        try:
            ch, pos, ref, mis = line.split()
            queries.append('chr{}:g.{}{}>{}'.format(ch, pos, ref, mis))
        except Exception as e:
            print(e)

    #retrieve dbNSFP
    mv = myvariant.MyVariantInfo()
    dbnsfp_results = {q: None for q in queries}
    results = mv.getvariants(queries)
    for result in results:
        dbnsfp_results[result['query']] = result.get('dbnsfp')
        
    #get predictors' result
    access_dict = {'bayesdel': ['add_af', 'rankscore'], #
               'bstatistic': ['rankscore'], 
               'clinpred': ['rankscore'], 
               'dann': ['rankscore'], 
               'deogen2': ['rankscore'], 
               'eigen': ['raw_coding_rankscore'], 
               'eigen-pc': ['raw_rankscore'], 
               'fathmm': ['rankscore'], 
               'fathmm-mkl': ['coding_rankscore'], 
               'fathmm-xf': ['coding_rankscore'], 
               'genocanyon': ['rankscore'], 
               'gerp++': ['rs_rankscore'], 
               'gm12878': ['fitcons_rankscore'], 
               'h1-hesc': ['fitcons_rankscore'], 
               'huvec': ['fitcons_rankscore'], 
               'integrated': ['fitcons_rankscore'], 
               'list-s2': ['rankscore'], 
               'lrt': ['converted_rankscore'], 
               'm_cap_score': ['rankscore'], 
               'metalr': ['rankscore'], 
               'metasvm': ['rankscore'], 
               'mpc': ['rankscore'], 
               'mutationassessor': ['rankscore'], 
               'mutationtaster': ['converted_rankscore'], 
               'mutpred': ['rankscore'], 
               'mvp': ['rankscore'], 
               'polyphen2_hdiv': ['hdiv', 'rankscore'], #
               'polyphen2_hvar': ['hvar', 'rankscore'], #
               'primateai': ['rankscore'], 
               'provean': ['rankscore'], 
               'revel': ['rankscore'], 
               'sift': ['converted_rankscore'], 
               'sift4g': ['converted_rankscore'], 
               'siphy_29way': ['logodds_rankscore'], 
               'vest4': ['rankscore']}
    
    def safe_recursive_get(d, query, keys):
        try:
            v = d[query]
            for key in keys:
                if key[:8] == 'polyphen':
                    key = 'polyphen2'
                v = v[key]
            return v
        except Exception:
            return np.nan
    
    predictor_values = {}
    for predictor in access_dict.keys():
        rank_name = access_dict[predictor]
        predictor_values[predictor] = [safe_recursive_get(dbnsfp_results, query, [predictor] + rank_name) for query in queries]
    
    df_all_predictions = pd.DataFrame(predictor_values, index=queries)
        
    #classification
    clf = lightgbm.Booster(model_file='clf.txt')
    answer = clf.predict(df_all_predictions)
    
    #output results
    pd.DataFrame({'score': answer}, index=queries).to_csv('predictions.csv', sep='\t')
    print('predictions are saved in predictions.csv')
    
def browseFiles():
    global filename 
    filename = filedialog.askopenfilename(initialdir = "/",
                                          title = "Select a File",
                                          filetypes = (("Text files",
                                                        "*.txt*"),
                                                       ("all files",
                                                        "*.*")))
    print('opened', filename)
    label_file_upload['text'] = filename.split('/')[-1]

            
#WINDOW
window = Tk()
window.title('Predictor')
window.geometry("490x350")
window.config(background = "white")
window.resizable(0,0)

#WIDGETS


#label_file_explorer = Label(window,text = "File Explorer using Tkinter",width = 100, height = 4,fg = "blue")
  
label_query = Label(window, text='Write in mutations in the following format:', bg='white', width=60)
label_query2 = Label(window, text = 'Chromosome\tGenomicCoordinates\tReferenceAllele\tMissenseAllele\n', bg='white', width=60)
label_query3 = Label(window, text = ' Example:', bg='white', anchor='w', width=60)
label_query4 = Label(window, text = ' 10\t100177368\tG\tA\n', bg='white', anchor='w', width=60)

text_query = Text(window, width=60, height=10, bd=3)

label_file_upload = Label(window, text='Or upload variations as a file (same format):', bg='white', width=45)

button_explore = Button(window, text = "Browse Files", command = browseFiles)
  
button_compute = Button(window, text = 'Submit', command = predict)

# Grid method is chosen for placing
# the widgets at respective positions
# in a table like structure by
# specifying rows and columns
label_query.grid(column=1, row=1, columnspan=4, rowspan=1)
label_query2.grid(column=1, row=2, columnspan=4, rowspan=1)
label_query3.grid(column=1, row=3, columnspan=4, rowspan=1)
label_query4.grid(column=1, row=4, columnspan=4, rowspan=1)

text_query.grid(column=1, row=5, columnspan=4, rowspan=4)

label_file_upload.grid(column=1, row=9, columnspan=3, rowspan=1)
  
button_compute.grid(column = 4, row = 10)
  
button_explore.grid(column = 2, row = 10)

# Let the window wait for any events
window.mainloop()

