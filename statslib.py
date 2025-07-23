#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Monica Keith"
__status__ = "Production"

import scipy.stats as stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pingouin as pg
import mylib
import fileslib
import math
import os
from matplotlib.colors import LinearSegmentedColormap

# Have to do this otherwise gives horrible errors
import  scipy.signal.signaltools
def _centered(arr, newsize):
    # Return the center newsize portion of the array.
    newsize = np.asarray(newsize)
    currsize = np.array(arr.shape)
    startind = (currsize - newsize) // 2
    endind = startind + newsize
    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]
    return arr[tuple(myslice)]

scipy.signal.signaltools._centered = _centered

# Performs a ttest between two groups for one or more continuous variables
# INPUTS:
# grp1: list of sbjs or sess that belong to grp1
# grp2: list of sbjs or sess that belong to grp2
# dic_vars: key: variable name, value: variable dictionary
# variable dictionary in dic_vars: key: sbj or sess, value: value of that variable for the corresponding sbj or sess
# OUTPUTS:
# dictionary: key: variable name, value: array[0]: t-statistic, array[1]: p-value
def ttest(grp1,grp2,dic_vars):
    # For each variable, create the two arrays to perform the ttest (array1,array2)
    # Perform the t-test and save the output [t,p] in the output dictionary
    dic_out = {}
    for var,dic in dic_vars.items():
        array1 = []
        for sbj in grp1:
            if sbj not in dic.keys():
                mylib.printError(f"{sbj} not found in {var} dictionary")
                return
            array1+=[float(dic[sbj])]
        
        array2 = []
        for sbj in grp2:
            if sbj not in dic.keys():
                mylib.printError(f"{sbj} not found in {var} dictionary")
                return
            array2+=[float(dic[sbj])]
        
        res = stats.ttest_ind(array1,array2)
        dic_out[var] = [res[0],res[1]]
        
    return dic_out

# Performs and ancova between two or more groups for one continous variable controling for zero or more covariates (confounders)
# If there are no confounders it is an anova
# INPUTS:
# dic_groups: key: group name, value: list of sbjs or sess that belong to that group
# var_name: continous (dependent) variable name
# dic_var: key: sess or sbj, value: value of the dependent variable for that sess or sbj
# dic_covars: key: covariate name, value: covariate dictionary
# covariate dictionary in dic_covats: key: sess or sbj, value: value of the covariate for that sess or sbj
# OUTPUTS: p-value
def ancova(dic_groups,var_name,dic_var,dic_covars={}):
    # Create the data frame
    # Create the columns for sbj, group and dependent variable
    sbj_col = []
    grp_col = []
    dv_col = []
    for grp,array in dic_groups.items():
        for sbj in array:
            sbj_col+=[sbj]
            grp_col+=[grp]
            dv_col+=[float(dic_var[sbj])]
    df = pd.DataFrame({"sbj":sbj_col,"grp":grp_col,var_name:dv_col})
    
    # Create one additional column for each covariate
    for covar,dic in dic_covars.items():
        array = []
        for sbj in sbj_col:
            array+=[float(dic[sbj])]
        df[covar] = array
        
    # Run the ancova or anova (if no covars)
    res = pg.ancova(data=df,dv=var_name,covar=list(dic_covars.keys()),between="grp")
    return res["p-unc"][0]

# Generates plots for comparing a continuous variable between groups
# INPUTS:
# dic_grps: key: group name, value: list of sbjs or sessions that belong to that group
# dic_var: key: sbj or sess, value: value of the variable for that sbj or sess
# var_name: variable name
# prefix: empty if not saving the plot, output prefix otherwise
# OUTPUTS:
# prefix_bars.png
def groupsPlots(dic_grps,dic_var,var_name,prefix="",title="",bars=True,box=False,violin=False,pairplot=False):
    # Create the two columns of the dataframe
    grouping = []
    var = []
    for grp,array in dic_grps.items():
        for sbj in array:
            if sbj not in dic_var.keys():
                mylib.printError(f"{sbj} not found in variable dictionary")
                return
            grouping+=[grp]
            var+=[float(dic_var[sbj])]
    
    # Create the dataframe        
    df = pd.DataFrame(data={"grouping":grouping,var_name:var})
    
    # Create the plot
    if bars:
        sns.set(style="whitegrid")
        sns.barplot(x="grouping",y=var_name,data=df,capsize=.1,hue="grouping")
        sns.swarmplot(x="grouping",y=var_name,data=df,color='0',alpha=.35,size=3)
        output = prefix+"_bars.png"
    elif box:
        sns.set(style="whitegrid")
        sns.boxplot(x="grouping",y=var_name,data=df)
        sns.swarmplot(x="grouping",y=var_name,data=df,color='0',alpha=.35,size=3)
        output = prefix+"_box.png"
    elif violin:
        sns.set(style="whitegrid")
        sns.violinplot(x="grouping",y=var_name,data=df)
        sns.swarmplot(x="grouping",y=var_name,data=df,color='0',alpha=.35,size=3)
        output = prefix+"_violin.png"
    else:
        sns.pairplot(df,vars=[var_name],kind='reg',hue="grouping")
        output = prefix+"_pairplot.png"
    plt.title(title)
    
    if prefix!="":
        plt.savefig(output, dpi=150)
        plt.close()
    else:
        plt.show()
 
# Perform correlations between two or more continuous variables controlling for zero or more confounders
# If there are zero confounders, it is a pearson correlation. Otherwise, it is a partial correlation.
# It is better than using dictionaries with just arrays of values because this way I can make sure that the value of one subject is present for all variables     
# INPUTS:
# sbj_list: list of subjects or sessions
# dic_cont_vars: key: variable name, value: variable dictionary
# variable dictionary in dic_cont_vars: key: sbj or sess, value: value of that variable for the corresponding sbj or sess
# dic_covars: key: confounder name, value: confounder dictionary
# confounder dictionary in dic_covars: key: sbj or sess, value: value of that confoudner for the corresponding sbj or sess
# OUTPUTS:
# dic["pvals"] = matrix of p-values
# dic["rvals"] = matrix of r-values
def pearson_partial(sbj_list,dic_cont_vars,dic_conf={}):
    # Initialize the dataframe
    df = pd.DataFrame({"sbj":[]})
    df.set_index("sbj",inplace=True)
    variables = list(dic_cont_vars.keys())
    confounders = list(dic_conf.keys())
    for col in variables+confounders:
        df[col] = []
        
    # Fill the DF
    for sbj in sbj_list:
        dic = {}
        complete = True
        
        # Add the continous variables
        for var,dic1 in dic_cont_vars.items():
            if sbj not in dic1.keys():
                complete = False
                break
            dic[var] = float(dic1[sbj])
        if not complete:
            continue
        
        # Add the confounders
        for var,dic2 in dic_conf.items():
            if sbj not in dic2.keys():
                complete = False
                break
            dic[var] = float(dic2[sbj])
        if not complete:
            continue
        
        # Add sbj to the DF
        for var,val in dic.items():
            df.loc[sbj,var] = val
    
    pvals = np.zeros((len(dic_cont_vars),len(dic_cont_vars)))
    rvals = np.zeros((len(dic_cont_vars),len(dic_cont_vars)))
    # For each variable, create one line of the matrix
    for i in range(len(variables)):
        # For each variable, create one column of the matrix
        for j in range(len(variables)):
            if i==j:
                pvals[i,j] = float("nan")
                rvals[i,j] = float("nan")
                continue
            elif pvals[i,j]!=0:
                continue
            
            # Perform the partial correlation
            res = pg.partial_corr(data=df, x=variables[i], y=variables[j], covar=confounders)
            r = res.loc["pearson","r"]
            p = res.loc["pearson","p-val"]
            pvals[i,j] = p
            pvals[j,i] = p
            rvals[i,j] = r
            rvals[j,i] = r
            
    return {"pvals":pvals,"rvals":rvals}

# Create a scatter plot for continuous variables
# OUTPUT: scatter plot if output!=""
def plotCont1(var1_name,var1_values,var2_name,var2_values,legend="",output="",noTics=False,title="",color_dots="blue"):
    plt.plot(var1_values,var2_values,'o',markersize=4,color=color_dots)
    plt.xlabel(var1_name)
    plt.ylabel(var2_name)
    plt.title(var1_name+" vs "+var2_name) if title=="" else plt.title(title)
    plt.grid(True)
    
    if noTics:
        plt.gca().axes.get_yaxis().set_ticklabels([])
        plt.gca().axes.get_xaxis().set_ticklabels([])
    
    # The legend doesn't fit the image if saved as png
    if legend!="" and output!="":
        # transform=plt.gca().transAxes makes it be independent of the size of the plot
        # If x and y coords are less or equal to 1 it will be inside the plot. Greater than 1 will make it outside the plot.
        plt.text(1.05 ,0.85, legend, transform=plt.gca().transAxes, fontsize=10)
    
    if output!="":
        plt.savefig(output, dpi=150)
    plt.show()
    plt.close()

# Create a scatter plot of contious variables: values of one dictionary vs the other
def plotCont1_dics(key_list,var1_name,dic1,var2_name,dic2,color="blue"):
    array1 = []
    array2 = []
    for key in key_list:
        if (key in dic1.keys()) and (key in dic2.keys()):
            array1+=[dic1[key]]
            array2+=[dic2[key]]
    plotCont1(var1_name,array1,var2_name,array2,color_dots=color)

# Create scatter plot for continous variables with the option of differenciating between groups
# INPUTS:
# sbj_list: list of subjects or sessions
# dic_vars: key: variable name, value: variable dictionary
# variable dictionary in dic_vars: key: sbj or sess, value: variable value for that sbj or sess
# output="" if not saving the plot
# grp_var: grouping variable to diffetenciate the dots of one group and the other(s)    
def plotCont2(sbj_list,dic_vars,output="",grp_var=""):
    # Create the columns of the dataFrame (one col per continuous variable)
    df = pd.DataFrame({"sbj":sbj_list})
    for var,dic in dic_vars.items():
        array = []
        for sbj in sbj_list:
            if sbj not in dic.keys():
                mylib.printError(f"{sbj} not in {var} dictionary")
                return
            array+=[float(dic[sbj])]
        df[var] = array
    
    # Create the plot
    sns.pairplot(df,vars=list(dic_vars.keys()),kind='reg') if grp_var=="" else sns.pairplot(df,vars=list(dic_vars.keys()),kind='reg', hue=grp_var)
    
    if output!="":
        plt.savefig(output, dpi=150)
    plt.show()
    plt.close()
    
# Create a scatter plot for continuous variables
# OUTPUT: scatter plot if output!=""
# var1_values[grp_names[i]] = list of values for variable 1 in sbjs of the corresponding grp    
# var2_values[grp_names[i]] = list of values for variable 2 in sbjs of the corresponding grp
def plotCont3(var1_name,var1_values,var2_name,var2_values,grp_names,output="",noTics=False):
    for i in range(len(grp_names)):
        plt.plot(var1_values[grp_names[i]],var2_values[grp_names[i]],'o',markersize=4,label=grp_names[i])
    
    plt.xlabel(var1_name)
    plt.ylabel(var2_name)
    plt.title(var1_name+" vs "+var2_name)
    plt.grid(True)
    plt.legend(bbox_to_anchor=(1.05,1.0),loc='upper left')
    
    if noTics:
        plt.gca().axes.get_yaxis().set_ticklabels([])
        plt.gca().axes.get_xaxis().set_ticklabels([])
    
    if output!="":
        plt.savefig(output, dpi=150)
    plt.show()
    plt.close()

# Create multiple boxplots for continous variables
# dics: list of dictionaries with the values to be plotted
# titles: title for each plot in the same order as the dicitonaries
def multipleBoxPlots(dics,titles,n_rows,n_cols,x_size,y_size):
    fig = plt.figure(figsize=(x_size,y_size))
    axes = []
    for i in range(1,len(dics)+1):
        ax = fig.add_subplot(n_rows,n_cols,i)
        ax.set_title(titles[i-1],fontsize='large')
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        axes+=[ax]
        
    for i in range(0,len(dics)):
        axes[i].boxplot(list(dics[i].values()));
    plt.show()
    plt.close()
        
def createCorrMatrix(tsfiles,tsfolder,motionfile='',censorfile=''):
    print("Creating correlation matrix...")
    
    premot = np.loadtxt(motionfile,delimiter=' ',usecols=0) if fileslib.linesCols(motionfile,' ')!=[0,0] else np.array([])
    censor = np.loadtxt(censorfile) if fileslib.linesCols(censorfile,' ')!=[0,0] else np.array([])
    if len(premot)>0 and len(censor)>0:
        mot = np.array([premot[i] for i in range(len(premot)) if censor[i]!=0])
    elif len(premot)>0 and len(censor)==0:
        mot = premot
    else:
        mot = np.array([])

    roi = []
    rmat = []
    pmat = []
    for tsfile1 in tsfiles:
        [nlines,ncols] = fileslib.linesCols(tsfile1,' ')
        if nlines>0 and ncols>0:
            x = np.loadtxt(tsfile1)
            fname = os.path.basename(tsfile1).split('.')[0]
            roi+=[fname]
            r = []
            p = []
            for tsfile2 in tsfiles:
                [nlines,ncols] = fileslib.linesCols(tsfile2,' ')
                if nlines>0 and ncols>0:
                    if tsfile1!=tsfile2:
                        y = np.loadtxt(tsfile2)
                        ans = stats.pearsonr(x, y) if len(x)==len(y) else (0,1)
                    else:
                        ans = (0,1)
                    r+=[ans[0]]
                    p+=[ans[1]]
            rmat+=[r]
            pmat+=[p]

        if len(mot)>0:
            if len(x)==len(mot):
                slope,intercept,motr,motp,stderr = stats.linregress(x, mot)
                print(f"Creating correlation graph with motion for {fname}...")
                line = "y="+str(round(intercept,4))+"+"+str(round(slope,4))+"x, r="+str(round(motr,4))
                plt.style.use('ggplot')
                fig, ax = plt.subplots()
                ax.plot(x, mot, linewidth=0, marker='s', label='Data points')
                ax.plot(x, intercept + slope * x, label=line)
                ax.legend(facecolor='white')
                outfig = tsfolder+fname+"_movCorr.png"
                plt.savefig(outfig, dpi=150)
                plt.close()
                print(f"Figure saved in {outfig}")
            else:
                mylib.printError("ts dont have the same size, could not perform correlation")
    
    rout = open(tsfolder+"rmat.csv",'w')
    pout = open(tsfolder+"pmat.csv",'w')

    # Write the first line of the file
    # Except in the group matrix for image files
    print("Saving the csv files...")
    rout.write(',')
    pout.write(',')
    for i in range(len(roi)):
        if i<len(roi)-1:
            rout.write(roi[i]+',')
            pout.write(roi[i]+',')
        else:
            rout.write(roi[i]+'\n')
            pout.write(roi[i]+'\n')

    # Write one line per roi
    for i in range(len(roi)):
        rout.write(roi[i])
        pout.write(roi[i])
        for j in range(len(roi)):
            rout.write(','+str(rmat[i][j]))
            pval = pmat[i][j]
            pout.write(','+str(pval))
            
        rout.write('\n')
        pout.write('\n')

    rout.close()
    pout.close()

def plotTS(ts):
    png = ts.replace('.txt','.png')
    while os.path.isfile(png):
        png = png.replace(".png","+.png")
    
    print(f"Generating graph for {ts}...")
    y = fileslib.file2array(ts)
    npoints = len(y)

    title = os.path.basename(ts).replace('.txt','')
    fig, ax = plt.subplots(1, 1, figsize=(15, 10))
    fig.suptitle(title, fontsize=32)
    
    myxticks = np.arange(1, npoints+1, step=1)
    ax.plot(myxticks, y,'g')
    ax.ticklabel_format(useMathText=True)
    plt.xlabel("Time points")
    plt.ylabel("Bold signal")
    myxticks = np.arange(1, npoints+1, step=10)
    myxticks = np.append(myxticks,npoints)
    plt.xticks(ticks=myxticks, rotation=90)
    plt.grid(axis='x', linestyle='--')
    plt.savefig(png, dpi=150)
    plt.close()
    
    print("done: "+png)

# movementFile can have 6 or 1 column(s): either the eucledian norm or [roll pitch yaw dS dL dP]
# lim is an arbitrary number where 0.2 is strict and 0.5 is lenient
# outliers[i]='1' if volume i should be included, outliers[i]='0' if volume should be censored
def plotmotion(outtitle,movementFile,outfile,screen=True,lim=0.2,outliers=[]):
    npoints = int(mylib.systemOut("cat "+movementFile+" | wc -l"))
    mov = np.zeros((npoints)) # all motion
    mov_cens = [] # only censored motion
    ln = np.zeros((npoints)) # limit
    shade1 = [False] * npoints # shaded areas
    
    i=0
    fin = open(movementFile,'r')
    for line in fin:
        array = list(filter(None,line.replace("\n","").split(' ')))
        
        # File has a single column with the eucledian norm
        if len(array)==1:
            enorm = float(array[0])
        
        # Calculate the eucledian norm if what we have is the raw movement values
        elif len(array)==6:
            enorm = math.sqrt((float(array[0])**2)+(float(array[1])**2)+(float(array[2])**2)+(float(array[3])**2)+(float(array[4])**2)+(float(array[5])**2))
        
        # The file has wrong format
        else:
            return
        
        mov[i] = enorm
        ln[i] = float(lim)
        
        # Define the outliers from the limit
        if len(outliers)==0: 
            if enorm>=float(lim):
                shade1[i] = True
                mov_cens+=[enorm]
        # Define the outliers from a vector of outiers
        elif outliers[i]=='0':
            shade1[i] = True
            mov_cens+=[enorm]
        
        i+=1
    fin.close()

    myxticks = np.arange(1, npoints+1, step=1)
    plt.plot(myxticks, mov,'b', label='enorm')
    plt.plot(myxticks, ln, 'r', label=">"+str(lim))
    plt.fill_between(myxticks, 0, max([max(mov),lim]), where=shade1, color='red', alpha=0.5)
    newxticks = np.arange(1, npoints+1, step=10)
    newxticks = np.append(newxticks,npoints)
    plt.xticks(ticks=newxticks,labels=newxticks,rotation=90,fontsize=8)
    plt.grid(axis='x', linestyle='--')

    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"Motion plot saved in {outfile}")
    if screen:
        os.system("gio open "+outfile)
    return [mov_cens,mov]

def simplepie(cursor,pipeline,n_exlc,grplist,steps,outfig):
    n=1
    for step in steps:
        if n==1:
            first = step
        elif n==len(steps):
            last = step
        n+=1
    
    n_tot = n_exlc
    n_fin = 0
    n_run = 0
    n_red = 0
    print("Hasn't start: ")
    for grp in grplist:
        n_tot+=len(grp)
        for sess in grp:
            if pipeline=="TBSS":
                cursor.execute("select DIFFsbj from subjects where sbjID='"+mylib.getSbjID(cursor,sess)+"'")
                entry = cursor.fetchone()
                if not entry==None or entry[0]!=1:
                    continue
            
            if mylib.getStatus(cursor,sess,first)=="checked":
                n_fin+=1
            elif mylib.getStatus(cursor,sess,last)!="None":
                n_run+=1
            else:
                print(sess)
                n_red+=1
    
    labels = []
    sizes = []
    if n_exlc!=0:
        labels+=['Excluded']
        sizes+=[n_exlc/n_tot]
    if n_fin!=0:
        labels+=['Completed']
        sizes+=[n_fin/n_tot]
    if n_run!=0:
        labels+=['Running']
        sizes+=[n_run/n_tot]
    if n_red!=0:
        labels+=['Ready']
        sizes+=[n_red/n_tot]
    
    fig1, ax1 = plt.subplots()
    plt.title(pipeline+" pipeline", weight="bold")
    ax1.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90, normalize=False)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig(outfig, dpi=150)
    plt.close()

def grpimg(cursor,grp,steps,outfig):
    mx = []   
    for sess in grp:
        sbjarray = []
        for step in steps:
            status = mylib.getStatus(cursor,sess,step)
            if status=='running' or status=='error' or status=='ready':
                sbjarray+=[1]
            elif status=='done':
                sbjarray+=[2]
            elif status=='checked':
                sbjarray+=[3]
            else:
                sbjarray+=[0]
        sbjarray+=[0,1,2,3]
        mx+=[sbjarray]
    
    mx = np.array(mx)
    fig, ax1 = plt.subplots(figsize=(15, 10))
    ticksx = steps
    ticksx+=["Not run", "Running", "Done", "Checked"]
    ax1.xaxis.set_ticks(range(len(ticksx)))
    ax1.xaxis.set_ticklabels(ticksx)
    ax1.yaxis.set_ticks(range(len(grp)))
    ax1.yaxis.set_ticklabels(grp)
    ax1.grid(True, linestyle=':')
    plt.setp(ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=10)
    plt.setp(ax1.get_yticklabels(), fontsize=10) 
    cm1_colors = [(0.52,0.32,0.5), (0.05,0.3,0.57), (0,0.56,0.8), (0,0.5,0.5)]
    cm1 = LinearSegmentedColormap.from_list("grp_map",cm1_colors,N=mx.max()+1)
    ax1.imshow(mx,cmap=cm1)
    fig.tight_layout()
    plt.axvline(x=len(steps)-0.5, color='black', linewidth=3)
    plt.savefig(outfig, dpi=150)
    #plt.show()
    plt.close()

def convert2zscores(var,cnx,table="demographicData",col="measure"):
    cursor = cnx.cursor()

    # Get the list of values for each of the two measures. Save in a dic
    # Get the mean & std of the values
    # Calculate the zscore of each value
    # Save in the DB
    dic = {}
    
    cursor.execute(f"select sess,value from {table} where {col}='{var}'")
    for (sess,value) in cursor.fetchall():
        dic[sess] = float(value)
        
    vals = list(dic.values())
    M = np.mean(vals)
    STD = np.std(vals)

    for (sess,value) in dic.items():
        new_val = (value-M)/STD
        mylib.insertEntryToTable(cursor,table,["sess",col,value],[sess,f'{var}_ZScore',new_val])
        mylib.updateDB(cnx,cursor)

    cursor.close()
    
# Remove pairs from dictionary based on a threshold for the values
def rmFromDic_lessThan(dic,thr):
    list_del = []
    for key,value in dic.items():
        if value<thr:
            list_del+=[key]
    for key in list_del:
        del dic[key]
        
# Remove pairs from dictionary based on a threshold for the values
def rmFromDic_moreThan(dic,thr):
    list_del = []
    for key,value in dic.items():
        if value>thr:
            list_del+=[key]
    for key in list_del:
        del dic[key]

# df from mylib.getSummaryTable or a subframe of it
def summaryNumbers(df,grp_cols,stats_cols=[]):
    # Merge group columns if necessary
    res = pd.DataFrame()
    if len(grp_cols)==0 or not set(grp_cols)<=set(df.columns.tolist()):
        return res
    if len(grp_cols)==1:
        grp_col = grp_cols[0]
    else:
        grp_col = "_".join(grp_cols)
        # Add the new grp column at the beginning of the df
        df = mylib.mergeColsDF(df,grp_cols)
    
    # Get list of unique groups
    groups = [grp for grp in list(set(df[grp_col].astype(str).values.tolist()))]
    res["grp"] = sorted(groups)
    res.set_index("grp",inplace=True)
    
    for grp in res.index.values.tolist():
        if grp=="None":
            grp_df = df[df[grp_col].isnull()]
        else:
            grp_df = df[df[grp_col]==grp]
        res.loc[grp,"N"] = len(grp_df)
            
        for col in stats_cols:
            col_vals = grp_df[col].astype(float).values.tolist()
            col_vals = np.delete(col_vals,np.isnan(col_vals)).tolist()
            m = round(np.mean(col_vals),3)
            std = round(np.std(col_vals),3)
            res.loc[grp,col] = f"{m},{std} [{min(col_vals)}-{max(col_vals)}]"
            
    for col in stats_cols:
        col_vals = df[col].astype(float).values.tolist()
        col_vals = np.delete(col_vals,np.isnan(col_vals)).tolist()
        m = round(np.mean(col_vals),3)
        std = round(np.std(col_vals),3)
        print(f"{col}: {m},{std} [{min(col_vals)}-{max(col_vals)}]")
    
    print(res)
    return res

# df from mylib.getSummaryTable or a subframe of it
def summaryStats(df,grp_cols,stats_cols):
    # Merge group columns if necessary
    res = pd.DataFrame()
    if len(grp_cols)==0 or not set(grp_cols)<=set(df.columns.tolist()) or len(stats_cols)==0:
        return res
    if len(grp_cols)==1:
        grp_col = grp_cols[0]
    else:
        grp_col = "_".join(grp_cols)
        # Add the new grp column at the beginning of the df
        df = mylib.mergeColsDF(df,grp_cols)
    
    # Get list of unique groups
    groups = [grp for grp in list(set(df[grp_col].astype(str).values.tolist()))]
    for col in stats_cols:
        # Remove all entries with empty values in col
        df_col = df.loc[:,[grp_col,col]]
        df_col = df_col[~df_col[col].isnull()]
        
        # Make sure the values in col are float
        df_col[col] = df_col[col].astype(float)
        
        # Performs and anova between the groups in grp_col for the continous variable col
        if len(np.unique(df_col[col].values.tolist()))==1:
            mylib.printWarning("Cant do anova for "+col+" because all lines have the same value")
            continue
        av = df_col.anova(dv=col,between=grp_col)
        pval = round(av["p-unc"].values[0],3)
        res.loc[col,"pvalue"] = pval
        
        for grp in groups:
            # Select only the entries that belong to grp
            df_grp = df_col[df_col[grp_col]==grp]
            grp_vals = df_grp[col].values.tolist()
            res.loc[col,grp] = round(np.mean(grp_vals),3)
    
    print(res)    
    return res
