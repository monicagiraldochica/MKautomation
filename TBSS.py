#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import mylib
import os
import datetime
import fileslib
import mrilib
import numpy as np
import pandas as pd
import importlib
from mysql.connector import Error
from mysql.connector import errorcode
import glob
import time
import ECP

# Examples of covariates: age, snr
# Examples of independent variables: grp, AEDgrp
# random_factors must be a subset of independent_vars (only needed in 3way & 2way)
# img_path inside the sess folder, i.e. "data/camino/wlf/ants/fa.nii.gz"
def createTBSS(cnx,cursor,name,test,img_path,project,covariates=[],independent_vars=[],effects="fixed",random_factors=[],dependent_variable="fa",rec_method="camino_wlf",description="",location="petrov",parent_pipe=""):
    ## Check that TBSS doesnt already exist
    if mylib.getNumberEntries(cursor,"tbss",["name"],[name])>0:
        mylib.printError(name+" already exist")
        return False
    
    # Check that it is a valid test
    valid_test = ["ttest","1grp1cov","1grp_ttest","2grpCovInteraction","anova_1way","anova_2way","anova_2FxNL","anova_3way"]
    if test not in valid_test:
        mylib.printError(f"{test} is not a valid test from {','.join(valid_test)}")
        return False
    
    # For 1grp1cov the main covariate must go first
    if test=="1grp1cov":
        if len(covariates)==0:
            mylib.printError("a 1grp1cov must have at least one covariate")
        
        if len(covariates)>1:
            print(f"list of covariates: {','.join(covariates)}")
            main = input(f"Choose the MAIN covariate ([Enter] for {covariates[0]}): ")
            if main!="":
                if main not in covariates:
                    mylib.printError(f"{main} is not in the covariates list: {','.join(covariates)}")
                    return False
                covariates.remove(main)
                covariates = [main]+covariates
            
    # Check the independent variables
    if ("3way" in test) and len(independent_vars)!=3:
        mylib.printError(f"{test} must have 3 independent variables but {len(independent_vars)} were given")
        return False
        
    if (("2way" in test) or ("2FxNL" in test)) and len(independent_vars)!=2:
        mylib.printError(f"{test} must have 2 independent variables but {len(independent_vars)} were given")
        return False
        
    if test!="1grp1cov" and test!="1grp_ttest" and len(independent_vars)<1:
        mylib.printError(f"{test} must have at least 1 independent variable but {len(independent_vars)} were given")
        return False
    
    # check effects and random_factors (only matters for 2way and 3way)
    if effects not in ["fixed","random","mixed"]:
        mylib.printError(f"{effects} effects not valid. must be one of these: fixed,random,mixed")
        return False
    if effects!="fixed" and not set(random_factors)<=set(independent_vars):
        mylib.printError(f"random_factors ({random_factors}) is not a subset of the independent variables ({independent_vars})")
        return False
    if (not "2way" in test) and (not "3way" in test):
        effects = ""
        random_factors = []
    
    # Insert new TBSS in tbss table
    cols = ["name","description","test","covariates","independent_variable","dependent_variable","rec_method","img_path","effects","random_factors"]
    values = [name,description,test,",".join(covariates),",".join(independent_vars),dependent_variable,rec_method,img_path,effects,",".join(random_factors)]
    mylib.insertEntryToTable(cursor,"tbss",cols,values)
    mylib.insertEntryToTable(cursor,"tbss_results",["name"],[name])
    
    # Insert new TBSS in sessions table
    cursor.execute("alter table sessions add column "+name+" tinyint(1) not null default 0")
    
    # Create folder in scratch
    pipes = mylib.getColumnFromTable(cursor,"pipename","pipelines")
    if parent_pipe=="":
        parent_pipe = input("Select one pipe from ["+",".join(pipes)+"]: ")
    if parent_pipe not in pipes:
        mylib.printError(parent_pipe+" is not a valid pipeline")
        return False
    
    scratchdir = mylib.getScratchDirPipe(cursor,location,parent_pipe)
    outdir = f"{scratchdir}/{name}"
    os.system("rm -rf "+outdir)
    os.system("mkdir "+outdir)
    
    mylib.updateDB(cnx,cursor)
    return True

def renameTBSS(cnx,cursor,name,new_name,location,parent_pipe):
    # Check that the TBSS exists
    if mylib.getNumberEntries(cursor,"tbss",["name"],[name])==0:
        mylib.printError(name+" doesnt exist")
        return False
    
    # Check that there's not another tbss with new_name
    if mylib.getNumberEntries(cursor,"tbss",["name"],[new_name])>0:
        mylib.printError(f"{new_name} already exists")
        return False
    
    # Remove the constraint
    cursor.execute("alter table tbss_results drop foreign key tbss_results_ibfk_1;")
    
    # Rename on the tbss_results table
    if mylib.updateEntryFromTable(cursor,"tbss_results",["name"],[new_name],["name"],[name])==-1:
        return False
    
    # Rename on the tbss table
    if mylib.updateEntryFromTable(cursor,"tbss",["name"],[new_name],["name"],[name])==-1:
        return False
    
    # Re-add the constraint
    cursor.execute("alter table tbss_results add constraint tbss_results_ibfk_1 foreign key(name) references tbss(name);")
    
    # Rename folder on scratch
    scratch = mylib.getScratchDirPipe(cursor,location,parent_pipe)
    if os.path.isdir(scratch+'/'+name):
        os.system(f"rm -rf {scratch}/{new_name}")
        os.system(f"mv {scratch}/{name} {scratch}/{new_name}")
    else:
        other_location = input(f"{name} folder to rename (full path): ")
        if os.path.isdir(other_location):
            os.system(f"mv {other_location} {other_location.replace(name,new_name)}")
        else:
            mylib.printError(other_location+" not found")
            return False
        
    mylib.updateDB(cnx,cursor)
    return True

def delTBSS(cnx,cursor,name,location,parent_pipe):
    # Check that the TBSS exists
    # The checking name is just being extra cautious that it doesnt remove the whole scratchdir in the last step of this function
    if name=="" or mylib.getNumberEntries(cursor,"tbss",["name"],[name])==0:
        mylib.printError(name+" doesnt exist")
        return False
    
    # Delete tbss_results
    if mylib.deleteEntryFromTable(cursor,"tbss_results",["name"],[name])==-1:
        return False
    
    # Delete from tbss table
    if mylib.deleteEntryFromTable(cursor,"tbss",["name"],[name])==-1:
        return False
    
    # Delete from scratch
    scratch = mylib.getScratchDirPipe(cursor,location,parent_pipe)
    if os.path.isdir(scratch+'/'+name):
        os.system(f"rm -r {scratch}/{name}")
    else:
        other_location = input(f"{name} folder to delete (full path): ")
        os.system(f"rm -rf {other_location}")
    
    # Delete any associated results     
    if input("Delete anything from results? [Y]: ")=="Y":
        print("Select results to delete divided by comma:")
        for result in mylib.getColumnFromTable(cursor,"result","results",distinct=True):
            print(result)
        delete = input(">> ").split(",")
        if delete!="":
            delete = "("+",".join(["'"+result+"'" for result in delete])+")"
            mylib.deleteEntriesFromTable(cursor,"results",["result"],[delete])
        
    mylib.updateDB(cnx,cursor)
    
def clearTBSS(cnx,cursor,name,location,parent_pipe,del_results=True):
    # Check that the TBSS exists
    # The checking name is just being extra cautious that it doesnt remove the whole scratchdir in the last step of this function
    if name=="" or mylib.getNumberEntries(cursor,"tbss",["name"],[name])==0:
        mylib.printError(name+" doesnt exist")
        return False
    
    # Delete tbss_results
    if del_results and mylib.deleteEntryFromTable(cursor,"tbss_results",["name"],[name])==-1:
        return False
    
    # Delete randomise, preStats_part2, preStats_part1
    if mylib.updateEntryFromTable(cursor,"tbss",["randomise","preStats_part2","preStats_part1"],["","",""],["name"],[name])==-1:
        return False
    
    # Mark tbss as 0 in the sessions table for everybody
    if mylib.updateEntriesFromTable(cursor,"sessions",[name],[0])==-1:
        return False
    
    mylib.updateDB(cnx,cursor)
    return True
    
def startTBSS(cnx,cursor,name,df,location,parent_pipe,test_mode=False):
    # Clear TBSS
    print("Clearing "+name+"...")
    if not clearTBSS(cnx,cursor,name,location,parent_pipe,False):
        return False
    print("done")
    
    # Mark as 1 in sessions table
    print("\nAdding sessions to "+name+"...")
    sess_list = "("+",".join(["'"+sess+"'" for sess in df.index.values.tolist()])+")"
    n = mylib.updateEntriesFromTable(cursor,"sessions",[name],[1],["sess"],[sess_list],test_mode)
    if not test_mode and n!=len(df):
        if n==-1:
            mylib.printError(f"no entries where updated")
        else:
            mylib.printError(f"{n} entries where updated instead of {len(df)}")
        cnx.rollback()
        return False
    print(f"done: {n} entries updated in sessions")
    
    # Create design files and mark preStats_part1 as ready
    createDesignFiles(cursor,name,df,location,parent_pipe,test_mode)
    if not createDesignFiles(cursor,name,df,location,parent_pipe,test_mode):
        cnx.rollback()
        return False
    
    if mylib.updateEntryFromTable(cursor,"tbss",["preStats_part1"],["ready"],["name"],[name],test_mode)==-1:
        cnx.rollback()
        return False
    
    mylib.updateDB(cnx,cursor)
    return True

# Demean values for each covariate, if any
# If covariate is already demeaned, it won't do anything because vals.mean will be zero
# This must be done across ALL sbjs and NOT within grp
# It should be done also for categorical covariates represented as 0/1
def demeanCovars(covars,df):
    try:
        for covar in covars:
            vals = df[covar].astype(float).values.tolist()
            demeaned = vals-np.mean(vals)
            
            if np.isnan(np.array(demeaned)).any():
                mylib.printError(f"{covar} demeaned vals contain nan")
                return pd.DataFrame()
            
            df.insert(len(df.columns), covar+"Demeaned", demeaned, True)
        
        return df
    
    except Error as err:
        mylib.printError(f"something went wrong deamining {','.join(covars)}\n{err}")
        return pd.DataFrame()

# Since we are including at least one covariate, we need to model the global mean.
# Since we are not including a column of 1s in the model, the randomise command MUST have the -D flag
def createDesignFile_1grp1cov(outdir,main_covar,covars,df,test_mode=False):
    # Demean the data
    df = demeanCovars([main_covar]+covars,df)
    if len(df)==0:
        return ''
    
    if test_mode:
        print(df)
    
    # Create the design matrix file design.mat 
    design = outdir+"/design.mat" if not outdir.endswith("/") else outdir+"design.mat"
    try:
        fout = open(design,'w')
        # NumWaves: # cols design.con, # cols design.mat: # EVs
        nEVs = 1+len(covars)
        fout.write(f"/NumWaves {nEVs}\n")
        fout.write(f"/NumPoints {len(df)}\n")
        fout.write(f"/EV1 {main_covar}Demeaned\n")
        for i in range(len(covars)):
            fout.write(f"/EV{i+2} {covars[i]}Demeaned\n")
        fout.write("/Matrix")
        
        for sess in df.index.values.tolist():
            fout.write('\n'+str(df.loc[sess,main_covar+"Demeaned"]))
            for covar in covars:
                fout.write(" "+str(df.loc[sess,covar+"Demeaned"]))
        
    except Error as err:
        mylib.printError(f"an exception happened while creating {design}\n{err}")
        design = ""
    
    finally:
        fout.close()
        if design!="":
            input(f"Design matrix created: {design}")
        return design

def createContrastFile_1grp1cov(outdir,main_covar,covars):
    contrast = outdir+"/design.con"
    try:
        fout = open(contrast,'w')
        nEVs = 1+len(covars)
        
        # NumWaves: # cols design.con, # cols design.mat: # EVs
        fout.write(f"/NumWaves {nEVs}\n")
        # NumContrasts: # lines design.con
        fout.write("/NumContrasts 2\n")   
        fout.write(f"/ContrastName1 pos_{main_covar}\n")
        fout.write(f"/ContrastName2 neg_{main_covar}\n")
        fout.write("/Matrix\n")
        
        C1 = "1"
        for i in range(len(covars)):
            C1+=" 0"
        fout.write(C1+"\n")
        C2 = "-1"
        for i in range(len(covars)):
            C2+=" 0"
        fout.write(C2)
        
    except Error as err:
        mylib.printError(f"an exception happened while creating {contrast}\n{err}")
        contrast = ""
    
    finally:
        fout.close()
        
        if contrast!="":
            input(f"Contrast file created: {contrast}")
        
        return contrast

def singleGrpSingleCov(outdir,df,main_covar,other_covars=[],test_mode=False):
    # Create the design matrix file design.mat
    design = createDesignFile_1grp1cov(outdir,main_covar,other_covars,df,test_mode)
    if not os.path.isfile(design):
        return False
    
    # Create the t-contrasts file design.con
    contrast = createContrastFile_1grp1cov(outdir,main_covar,other_covars)
    if not os.path.isfile(contrast):
        return False
        
    if test_mode:
        input("Check {contrast} before it's deleted [Enter]")
        input("Check {design} before it's deleted [Enter]")
        
        os.system("rm "+design+" "+contrast)
    
    return True

# 1-way, 2-way, 2xN between-subjects ANOVA (cell means models)
# For more than 2 factors one should use the factor effects approach
# factor_levels: dictionary factor_levels[factor] = list of levels (groups) for that factor (group_name)
# The global mean is not modeled because it is already included in the GLM. It is represented jointly by the regressors of each grp.
def createDesignFile_ftest(cursor,outdir,df,factor_levels,covariates=[],test_mode=False):
    # This function only works with 1 or 2 factors
    factors = list(factor_levels.keys())
    n_factors = len(factors)
    if n_factors not in [1,2]:
        return ['',[]]
    
    if test_mode:
        print("factors and levels:")
        print(factor_levels)
    
    # Demean values for each covariate, if any
    df = demeanCovars(covariates,df)
    if len(df)==0:
        return ['',[]]
    
    if test_mode and len(covariates)>0:
        print(df.loc[:,covariates+[covar+"Demeaned" for covar in covariates]])
    
    design = outdir+"/design.mat"
    EVs = []
    try:
        fout = open(design,'w')
        
        # Levels of the first factor
        levels_factor1 = factor_levels[factors[0]]
        # Levels of the second factor (or an empty array if there's only one factor)
        levels_factor2 = [] if n_factors==1 else factor_levels[factors[1]]
        
        # If there's 1 factor, nEVs = number of levels (groups) for that factor
        # If there are 2 factors, nEVs = number_levels_factor1 x number_levels_factor2
        nEVs = len(levels_factor1) if n_factors==1 else len(levels_factor1)*len(levels_factor2)
        nEVs+=len(covariates)
        
        # NumWaves: # cols design.con, # cols design.mat: nEVs
        fout.write(f"/NumWaves {nEVs}\n")
        fout.write(f"/NumPoints {len(df)}\n")
        matrix = "/Matrix"
        
        for sess in df.index.values.tolist():
            matrix+="\n"
            
            for level1 in levels_factor1:
                l1_name = mylib.getStringVal(cursor,factors[0],level1).replace("_","") if not isinstance(level1,str) else level1
                
                if n_factors==1:
                    matrix+=str(int(df.loc[sess,factors[0]]==level1))+' '
                    if l1_name not in EVs:
                        EVs+=[l1_name]
                        fout.write("/EV"+str(len(EVs))+" "+l1_name+"\n")
                for level2 in levels_factor2:
                    l2_name = mylib.getStringVal(cursor,factors[1],level2).replace("_","") if not isinstance(level2,str) else level2
                    matrix+=str(int(df.loc[sess,factors[0]]==level1 and df.loc[sess,factors[1]]==level2))+' '
                    if l1_name+"_"+l2_name not in EVs:
                        EVs+=[l1_name+"_"+l2_name]
                        fout.write("/EV"+str(len(EVs))+" "+l1_name+"_"+l2_name+"\n")
            
            for covar in covariates:
                matrix+=str(df.loc[sess,covar+"Demeaned"])+' '
                if covar.replace("_","") not in EVs:
                    EVs+=[covar.replace("_","")]
                    fout.write(f"/EV{len(EVs)} {covar}\n")
                
        fout.write(matrix)
    
    except Error as err:
        mylib.printError(f"an exception happened while creating {design}\n{err}")
        design = ""
    
    finally:
        fout.close()
        
        if design!="":
            input(f"Design file created: {design}")
        
        if test_mode and len(EVs)>0:
            print("EVs: "+",".join(EVs))
        
        return [design,EVs]

# pairs_array must be initialized before calling, it'll be modified (new pairs added)
def getAllPairs(groups,pairs_array,join_char):
    for grp1 in groups:
        for grp2 in groups:
            if grp1!=grp2 and (grp1+join_char+grp2 not in pairs_array):
                pairs_array+=[grp1+join_char+grp2]
                
# This function only works with a 1-way between subjects ANOVA (1 factor, k levels)
# If the number of factors is different, the contrast file has a different format
def createContrastFile_1wayANOVAcell(cursor,outdir,group_name,levels,EVs,test_mode=False):
    groups = []
    for level in levels:
        if not isinstance(level,str):
            groups+=[mylib.getStringVal(cursor,group_name,level).replace("_","")]
        else:
            groups+=[level]
    
    # Get the main contrasts
    n_main_contrast = len(groups)-1
    contrasts = []
    for cont in range(1,n_main_contrast+1):
        contrasts+=[groups[cont-1]+"-"+groups[-1]]
        
    if test_mode:
        print("main contrast:")
        print(contrasts)
        
    # Get additional contrasts
    getAllPairs(groups,contrasts,"-")
    
    if test_mode:
        print("\nall contrast:")
        print(contrasts)
    
    # NumWaves: # cols design.con, # cols design.mat: # EVs: number of groups (levels)
    contrast = outdir+"/design.con"
    ftest_array = []
    try:
        fout = open(contrast,'w')
        fout.write(f"/NumWaves {len(EVs)}\n")
        
        # NumContrasts: # lines design.con
        fout.write(f"/NumContrasts {len(contrasts)}\n")
        
        n=1
        for cname in contrasts:
            fout.write(f"/ContrastName{n} {cname}\n")
            n+=1
        
        # One line per contrast
        fout.write("/Matrix")
        for cname in contrasts:
            fout.write('\n')
            info = cname.split("-")
            
            # One column per EV
            for ev in EVs:
                fout.write("-1") if ev==info[1] else fout.write(str(int(ev==info[0])))
                if ev!=EVs[-1]:
                    fout.write(" ")
                    
        # The fts file will have only one ftest formed by the main contrasts, the other contrasts will have value 0
        ftest_array = [int(i) for i in list(np.ones(n_main_contrast))+list(np.zeros(len(contrasts)-n_main_contrast))]
    
    except Error as err:
        mylib.printError(f"an exception happened while creating {contrast}\n{err}")
        contrast = ""
    
    finally:
        fout.close()
        
        if contrast!="":
            input(f"Contrast file created: {contrast}")             
        
        return [contrast,[ftest_array]]

# arrays: one array per f-test. All arrays should have the same length.
# Only the first N contrast need to be included in the F test
# The first N contrast that cover all the groups
# Once the difference between those groups is known, all other differences are fully determined. So, they don't need to go explicity in the f-test.
def createFtestFile(outdir,arrays,ftest_names=[]):
    ftest = outdir+"/design.fts"
    try:
        fout = open(ftest,'w')
        
        # NumWaves: #contrasts/#lines in design.con file, # columns in F-test (different than # waves in .con file, instead # contrast in .con file)
        fout.write(f"/NumWaves {len(arrays[0])}\n")
        # One contrast per f-test
        fout.write(f"/NumContrasts {len(arrays)}\n")
        
        for i in range(len(ftest_names)):
            fout.write(f"/ContrastName{i+1} {ftest_names[i]}\n")
        
        # One line per f-test
        fout.write("/Matrix")
        for array in arrays:
            fout.write("\n")
            
            # One column per contrast in the .con file
            for n in range(len(array)):
                fout.write(str(int(array[n])))
                if n<len(array)-1:
                    fout.write(" ")
        
    except Error as err:
        mylib.printError(f"an exception happened while creating {ftest}\n{err}")
        ftest = ""
        
    finally:
        fout.close()
        
        if ftest!="":
            input(f"F-test file created: {ftest}")
            
        return ftest

# 1-way between sbjs anova: 1 factor with k levels
# Create k dummy variables, one for each level
# k-th dummy variable has value 1 for observations in the in the k-th level of the factor
# The k dummy variables also model the grand mean, so randomaise command should be used without the -D flag
# This can lead to problems when using more than one factor & lead to rank-deficiency matrices. So, do not use with more than one factor
def ANOVAcell_1way(cursor,outdir,group_name,tbss_name,df,covariates=[],test_mode=False):
    groups = np.unique(df[group_name].astype(str).values).tolist()
    if len(groups)<2:
        return False
    
    # Create the design matrix file design.mat
    [design,EVs] = createDesignFile_ftest(cursor,outdir,df,{group_name:groups},covariates,test_mode)
    if not os.path.isfile(design):
        return False
        
    # Create the t-contrasts file design.con
    [contrast,ftest_array] = createContrastFile_1wayANOVAcell(cursor,outdir,group_name,groups,EVs,test_mode)
    if not os.path.isfile(contrast):
        return False
    
    # If there's only two groups, don't do an f-test 
    # And change the test to a ttest because in that case it's the same as a ttest & the script is the one without fts file
    if len(groups)>2:
        ftest = createFtestFile(outdir,ftest_array,[group_name])
        if not os.path.isfile(ftest):
            return False
    else:
        ftest = ""
        mylib.updateEntriesFromTable(cursor,"tbss",["test"],["ttest"],["name"],[tbss_name])
        
    if test_mode:
        input("Check {contrast} before it's deleted [Enter]")
        input("Check {design} before it's deleted [Enter]")
        if os.path.isfile(ftest):
            input("Check {ftest} before it's deleted [Enter]")
        
        os.system("rm -f "+design+" "+contrast+" "+ftest)
    
    return True

def mergeFiles(cursor,df,tbss_dir,tbss_name,parent_pipe,location,test_mode=False):
    [img,img_path] = mylib.getEntryFromTable(cursor,["dependent_variable","img_path"],"tbss",["name"],[tbss_name])
    
    # Dir of symlinks for all the maps for the current tbss
    links_dir = tbss_dir+"/"+img.upper()
    if not test_mode:
        os.system("rm -rf "+links_dir)
        os.system("mkdir "+links_dir)
    
    # Create symlinks to the img files in the tbss folder
    merged = tbss_dir+"/all_"+img.upper()+".nii.gz"
    if not test_mode:
        os.system("rm -f "+merged)
        
    print(f"\nCopying sessions {img} files...")
    cmd = "fslmerge -t "+merged
    for sess in df.index.values.tolist():
        # Link to the original in the tbss folder
        copy = links_dir+"/"+sess+".nii.gz"
        
        # Original sbj map in the project directory should not be missing
        orig = mylib.getSessDir(cursor,parent_pipe,location,sess)+"/"+img_path.replace("$sess",sess).replace("$sbj",df.loc[sess,"sbjID"])
        if not os.path.isfile(orig):
            mylib.printError(f"{orig} missing")
            return False
        
        if not test_mode:
            os.system("ln -s "+orig+" "+copy)
        else:
            print("ln -s "+orig+" "+copy)
        cmd+=" "+copy
    print("done")
        
    # Create the merged file
    print("Creating the merged file...")
    print(cmd)
    if not test_mode:
        os.system(cmd)
    print("done")
    
    # Is better to save the exact order at which subjects were concatenated
    fileslib.array2file(df.index.values.tolist(),tbss_dir+"/merged_order.txt","\n")
    
    # Check the merged file
    if not test_mode:
        n_vols = mrilib.getDims(merged)[3]
        if n_vols!=len(df):
            mylib.printError("could not merge the files")
            return False
        
        print(f"{n_vols} images successfully merged")
        
    return True

def arrayByElement(array,element,join_char):
    new_array = []
    for item in array:
        new_array+=[item+join_char+element]
    return new_array

def negArray(array):
    new_array = []
    for item in array:
        if item.startswith("-"):
            new_array+=[item[1:]]
        else:
            new_array+=["-"+item]
    return new_array

def invertEV(ev):
    array = ev.split("_")
    new_ev = array[-1]
    i = len(array)-2
    while i>=0:
        new_ev+="_"+array[i]
        i-=1
    return new_ev

# Does NOT work for the 1-way nor 3-way
# ONLY WORKS for ANOVA with 2 factors
def createContrastFile_nwayANOVAcell_2factors(cursor,outdir,factor_levels,EVs,test_mode):
    ## Get the contrasts for the main effects of each factor
    all_contrasts = {} # All the contrasts of each factor
    main_contrasts = {} # The main contrasts of each factor
    num_interactions = 0
    n_contrast = 0
    factor_groups = {}
    for factor,levels in factor_levels.items():
        # Get the grp names for the factor
        groups = []
        for level in levels:
            if not isinstance(level,str):
                groups+=[mylib.getStringVal(cursor,factor,level).replace("_","")]
            else:
                groups+=[level]
        factor_groups[factor] = groups
            
        # Get the main interactions for the levels of this factor
        main_interact = []
        for i in range(1,len(groups)):
            main_interact+=["-".join([groups[i],groups[0]])]
        main_contrasts[factor] = main_interact
        if len(main_interact)>num_interactions:
            num_interactions = len(main_interact)
        
        # Get all interactions for the levels of this factor
        # getAllPairs modifies all_interact
        all_interact = []
        getAllPairs(groups,all_interact,"-")
        all_contrasts[factor] = all_interact
        n_contrast+=len(all_interact)
        
    ## Get the contrast names for the interactions
    factors = list(factor_levels.keys())
    factor0 = factors[0]
    groups_factor0 = factor_groups[factor0]
    factor1 = factors[1]
    groups_factor1 = factor_groups[factor1]
    all_contrasts["interaction"] = []
    interactions = []
    
    if test_mode:
        print("groups for "+factor0+": "+",".join(groups_factor0))
        print("groups for "+factor1+": "+",".join(groups_factor1))
    
    # 2-way ANOVA has only one interaction (between the two factors)
    if len(groups_factor1)<=2 and len(groups_factor0)<=2:
        all_contrasts["interaction"]+=[groups_factor0[0]+"-"+groups_factor0[1]+"_across_"+factor1]
        n_contrast+=1
        
        # Interactions for the formula
        array = []
        for group1 in groups_factor1:
            array+=[arrayByElement([groups_factor0[0],"-"+groups_factor0[1]],group1,"_")]
        interactions+=[array]
        
    # One or the two factors have more than 2 levels
    else:
        # Get the interactions of factor0 across factor1
        # groups_factor0[0] minus each other factor0_group across factor1
        if len(groups_factor1)>2:
            for group in groups_factor0[1:]:
                for i in range(1,len(groups_factor1)):
                    all_contrasts["interaction"]+=[groups_factor0[0]+"-"+group+"_across_"+factor1+"-cont"+str(i)]
                    n_contrast+=1
                    
                # Interactions for the formula
                array = []
                for group1 in groups_factor1:
                    array+=[arrayByElement([groups_factor0[0],"-"+group],group1,"_")]
                interactions+=[array]
        
        # Get the interactions of factor1 across factor0
        # groups_factor1[0] minus each other factor1_group across factor0
        if len(groups_factor0)>2:
            for group in groups_factor1[1:]:
                for i in range(1,len(groups_factor0)):
                    all_contrasts["interaction"]+=[groups_factor1[0]+"-"+group+"_across_"+factor0+"-cont"+str(i)]
                    n_contrast+=1
                    
                # Interactions for the formula
                array = []
                for group0 in groups_factor0:
                    array+=[arrayByElement([groups_factor1[0],"-"+group],group0,"_")]
                interactions+=[array]
    
    # Simplify the interactions formula
    simplified_interactions = []
    for array in interactions:
        neg = negArray(array[-1])
        for interact in array[:-1]:
            simplified_interactions+=[interact+neg]
    
    ## Get the interaction lines for the matrix
    interaction_lines = []
    for interact in simplified_interactions:
        array_line = []
        for ev in EVs:
            if (ev in interact) or (invertEV(ev) in interact):
                array_line+=[1]
            elif ("-"+ev in interact) or ("-"+invertEV(ev) in interact):
                array_line+=[-1]
            else:
                array_line+=[0]
        interaction_lines+=[array_line]
            
    # All the interactions will be included in the f-tests
    main_contrasts["interaction"] = all_contrasts["interaction"]
        
    if test_mode:
        print("all contrasts for each facotr:")
        print(all_contrasts)
        
        print("\nmain contrasts for the f-test:")
        print(main_contrasts)
    
    ## Create the contrast file   
    design = outdir+"/design.con"
    try:
        fout = open(design,'w')
        
        # NumWaves: # cols design.con, # cols design.mat: # EVs
        fout.write(f"/NumWaves {len(EVs)}\n")
        
        # NumContrasts: # lines design.con
        fout.write(f"/NumContrasts {n_contrast}\n")
        
        ## Write the contrast names and initialize the matrix
        i = 1
        j = 0
        matrix = []
        ftest_array = []
        ftest_names = []
        for factor,cont_array in all_contrasts.items():
            for cont_name in cont_array:
                fout.write(f"/ContrastName{i} {cont_name}\n")
                i+=1
                
                if factor!="interaction":
                    pos = cont_name.split("-")[0]
                    neg = cont_name.split("-")[1]
                    
                    line = []
                    for ev in EVs:
                        array = ev.split("_")
                        if pos in array:
                            line+=[1]
                        elif neg in array:
                            line+=[-1]
                        else:
                            line+=[0]
                else:
                    line = interaction_lines[j]
                    j+=1
                    for k in range(len(line),len(EVs)):
                        line+=[0]
                
                matrix+=[line]
            
            farray = []
            for factor2,cont_array2 in all_contrasts.items():
                for cont_name2 in cont_array2:
                    farray+=[int(factor2==factor and (cont_name2 in main_contrasts[factor2]))]
            ftest_array+=[farray]
            ftest_names+=[factor]
        
        ## Write the matrix        
        fout.write("/Matrix")
        for line in matrix:
            fout.write("\n")
            for col in line:
                fout.write(str(col)+" ")
    
    except Error as err:
        mylib.printError(f"an exception happened wile creating {design}\n{err}")
        design = ""
        ftest_array = []
        ftest_names = []
    
    finally:
        fout.close()
        
        if design!="":
            input(f"Contrast file created: {outdir}/design.con")
        
        return [design,ftest_array,ftest_names]

# 2 factors ANOVA
# The k dummy variables also model the grand mean, so randomaise command should be used without the -D flag
# This can lead to problems when using more than one factor & lead to rank-deficiency matrices. So, do not use with more than one factor
def ANOVAcell_2way_2FxNL(cursor,outdir,tbss_name,df,factors,covariates,test_mode=False):
    if len(factors)!=2:
        return False
    
    # Get the existing levels for each factor
    factor_levels = {}
    for factor in factors:
        factor_levels[factor] = np.unique(df[factor].values.tolist())
    
    # Create the design matrix file design.mat
    [design,EVs] = createDesignFile_ftest(cursor,outdir,df,factor_levels,covariates,test_mode)
    if not os.path.isfile(design):
        return False
    
    # Create the t-contrasts file design.con
    [contrast,ftest_array,ftest_names] = createContrastFile_nwayANOVAcell_2factors(cursor,outdir,factor_levels,EVs,test_mode)
    if not os.path.isfile(contrast):
        return False
    
    # Create the F-test file
    ftest = createFtestFile(outdir,ftest_array,ftest_names)
    if not os.path.isfile(ftest):
        return False
    
    if test_mode:
        input("Check {contrast} before it's deleted [Enter]")
        input("Check {design} before it's deleted [Enter]")
        input("Check {ftest} before it's deleted [Enter]")
        
        os.system("rm "+contrast+" "+design+" "+ftest)
    
    return True

def getUnpermuttedPairs(groups,pairs_array,join_char):
    for grp1 in groups:
        for grp2 in groups:
            if grp1!=grp2 and (grp1+join_char+grp2 not in pairs_array) and (grp2+join_char+grp1 not in pairs_array):
                pairs_array+=[grp1+join_char+grp2]
                
# pairs_array must be initialized before calling, it'll be modified (new pairs added)
def getUnpermuttedTrios(groups,trios_array,join_char):
    for grp1 in groups:
        for grp2 in groups:
            for grp3 in groups:
                if grp1!=grp2 and grp1!=grp3 and grp2!=grp3 and (grp1+join_char+grp2+join_char+grp3 not in trios_array) and (grp1+join_char+grp3+join_char+grp2 not in trios_array) and (grp2+join_char+grp1+join_char+grp3 not in trios_array) and (grp2+join_char+grp3+join_char+grp1 not in trios_array) and (grp3+join_char+grp1+join_char+grp2 not in trios_array) and (grp3+join_char+grp2+join_char+grp1 not in trios_array):
                    trios_array+=[grp1+join_char+grp2+join_char+grp3]
                
# Works for at most 3 levels for each factor. If a factor has more than 3 levels, it won't work
def createDesignFile_nwayANOVAfactor(outdir,df,factor_levels,covariates,test_mode=False):
    ######################
    ## Get the list EVs ##
    ######################    
    ## Find the main factors effects
    factors = list(factor_levels.keys())
    EVs = []
    for factor in factors:
        # Each factor will have nLevels-1 EVs
        nEVs_factor = len(factor_levels[factor])-1
        factor = factor.replace("_","")
        if nEVs_factor==1:
            EVs+=[factor]
        else:
            for i in range(1,nEVs_factor+1):
                EVs+=[factor+"-ev"+str(i)]
                
    ## Find all possible interactions
    two_way_interact = []
    getUnpermuttedPairs(EVs,two_way_interact,"_")
    three_way_interact = []
    getUnpermuttedTrios(EVs,three_way_interact,"_")
                    
    ## Add all the interactions to the EVs
    # Print first the 2-way interactions, then the 3-way interactions, etc
    # Last EV is the grand mean or intercept
    EVs = EVs+two_way_interact+three_way_interact+["mean"]+covariates
    
    if test_mode:
        print("EVs: "+",".join(EVs))
    
    ############################
    ## Create the design file ##
    ###################### #####
    design = outdir+"/design.mat"
    try:
        fout = open(design,'w')
        
        # NumWaves: # cols design.con, # cols design.mat: # EVs
        fout.write(f"/NumWaves {len(EVs)}\n")
        fout.write(f"/NumPoints {len(df)}\n")
        
        for i in range(len(EVs)):
            fout.write(f"/EV{i+1} {EVs[i]}\n")
        
        fout.write("/Matrix")
        ok = True
        for sess in df.index.values.tolist():
            fout.write("\n")
            dic = {}
            
            # for each EV except the last one which is all 1 and except the covariates
            for i in range(len(EVs)-1-len(covariates)):
                ev = EVs[i]
                factors = ev.split("_")
                
                # The EV is composed by a single factor
                if len(factors)==1:
                    factor = factors[0]
                    group_name = factor if ("-ev" not in factor) else factor.split("-ev")[0]
                    grp = df.loc[sess,group_name]
                    levels = factor_levels[group_name]
                    if len(levels)>3:
                        design = ""
                        EVs = []
                        ok = False
                        break
                    
                    # Levels[0] will be the reference level. EVs=-1 for sess in that group
                    if grp==levels[0]:
                        dic[ev] = -1
                    # There's only on EV for this factor, so if it's not the reference level, value is 1
                    elif "-ev" not in factor:
                        dic[ev] = 1
                    # There's more than one EV for this factor, it will have value 1 for the level of interest, 0 otherwise
                    else:
                        level_interest = int(factor.split("-ev")[1])-1
                        dic[ev] = int(grp==levels[level_interest])
                    fout.write(str(dic[ev])+" ")
                
                # The EV is a 2-way interaction effect
                elif len(factors)==2:
                    fout.write(str(dic[factors[0]]*dic[factors[1]])+" ")
                
                # The EV is a 3-way interaction effect
                elif len(factors)==3:
                    fout.write(str(dic[factors[0]]*dic[factors[1]]*dic[factors[2]])+" ")
                
                # Since this function only works with a max of 3 levels, there can't be any larger interactions
                else:
                    design = ""
                    EVs = []
                    ok = False
                    break
                
            if not ok:
                break
                           
            # Last column is the grand mean
            fout.write("1")
            
            # Insert the columns for covariates
            for covar in covariates:
                fout.write(" "+str(df.loc[sess,covar]))
    except Error as err:
        print(f"an exception happened creating {design}\n{err}")
        design = ""
        EVs = []
    
    finally:                                                                                    
        fout.close()
        
        if ok:
            input(f"Design file created: {design}")
        else:
            os.system("rm "+design)
            design = ""
            EVs = []
            
        return [design,EVs]

def createFarrays(ncontrast_array):
    num_ftest = len(ncontrast_array)
    ftest_array = []
    
    # One line per f-test
    for i in range(num_ftest):
        array = []
        # For each ftest write 1s in the corresponding column if i is the corresponding ftest, 0 otherwise
        for j in range(num_ftest):
            for k in range(ncontrast_array[j]):
                array+=[int(i==j)]
        ftest_array+=[array]
        
    return ftest_array

# Does NOT work for the 1-way
# ONLY WORKS for 2-way, 3-way & 2x3 between sbjs ANOVA
def createContrastFile_nwayANOVAfactor(cursor,outdir,EVs,covariates):
    contrast = outdir+"/design.con"
    try:
        fout = open(contrast,'w')
        
        # NumWaves: # cols design.con, # cols design.mat: # EVs
        fout.write(f"/NumWaves {len(EVs)}\n")
        
        # NumContrasts: # lines design.con: nEVs-1 (subtracting the grand mean), then subtracting the covariates
        fout.write(f"/NumContrasts {len(EVs)-1-len(covariates)}\n")
        
        # ncontrast_dic contains the number of contrasts for each ftest
        ncontrast_dic = {}
        contrasts = []
        for i in range(len(EVs)):
            cname = EVs[i-1]
            if (cname in covariates) or cname=='mean':
                continue
            contrasts+=[cname]
            fout.write(f"/ContrastName{i} {cname}\n")
            
            ftest = cname.split("-ev")[0]
            if ftest not in ncontrast_dic.keys():
                ncontrast_dic[ftest] = 1
            else:
                ncontrast_dic[ftest]+=1
        
        fout.write("/Matrix")
        for line in range(len(contrasts)):
            fout.write("\n")
            
            # One col per EV (except mean & covariates)
            for col in range(len(contrasts)):
                fout.write(str(int(line==col))+" ")
                
            # The col corresponsing to mean EV is always 0
            fout.write("0")
            
            # One zero col per covar
            for covar in covariates:
                fout.write(" 0")
                
        # The fts file will have len(ncontrast_dic) f-test formed by ncontrast_dic[f-test] contrasts
        ftest_array = createFarrays(list(ncontrast_dic.values()))
    
    except Error as err:
        mylib.printError(f"an exception happened while creating {contrast}\n{err}")
        contrast = ""
        ftest_array = []
    
    finally:
        fout.close()
        
        if contrast!="":
            input(f"Contrast file created: {contrast}")
        
        return [contrast,ftest_array]

# 3-way between sbjs anova: 3 factors, 2 levels
# The design includes the grand mean, so randomaise command should be used without the -D flag
# The statistics generated by randomise using this design are the final results if it is a fixed-effects model
# However, if it's a mixed or random effects model it will need some post processing after randomise finishes running. See checkDoneRandomise function which deals with this.
def ANOVAfactor_3way(cursor,outdir,tbss_name,df,factors,covariates,test_mode=False):
    if len(factors)!=3:
        return False
    
    factor_levels = {}
    for factor in factors:
        factor_levels[factor] = np.unique(df[factor].values).tolist()
    
    # Create the design matrix file design.mat 
    [design,EVs] = createDesignFile_nwayANOVAfactor(outdir,df,factor_levels,covariates,test_mode)
    if not os.path.isfile(design):
        return False
        
    # Create the t-contrasts file design.con
    [contrast,ftest_array] = createContrastFile_nwayANOVAfactor(cursor,outdir,EVs,covariates)
    if not os.path.isfile(contrast):
        return False
    
    # Create the F-test file
    ftest = createFtestFile(outdir,ftest_array,EVs[:-(1+len(covariates))],test_mode)
    if not os.path.isfile(ftest):
        return False
    
    if test_mode:
        input("Check {contrast} before it's deleted [Enter]")
        input("Check {design} before it's deleted [Enter]")
        input("Check {ftest} before it's deleted [Enter]")
        
        os.system("rm "+contrast+" "+design+" "+ftest)
    
    return True

# pairs_array must be initialized before calling, it'll be modified (new pairs added)
# The global mean is not modeled because it is already included in the GLM. It is represented jointly by the regressors of each grp.
def createDesignFile_2grpCovarInteraction(outdir,df,groups,group_name,covar,test_mode=False):
    # Demean values of the covariate
    df = demeanCovars([covar],df)
    
    # outdir would not be created so the design file cant be created
    if test_mode:
        print(df.loc[:,[group_name,covar,covar+"Demeaned"]])
    
    # Generate the design.mat file
    design = outdir+"/design.mat"
    try:
        fout = open(design,'w')
        # NumWaves: # cols design.con, # cols design.mat: # EVs
        nEVs = 4
        fout.write(f"/NumWaves {nEVs}\n")
        fout.write(f"/NumPoints {len(df)}\n")
        fout.write("/Matrix")
        
        for sess in df.index.values.tolist():
            # Get the grp membership of the sbj
            grp = df.loc[sess,group_name]
            # Split the covariate in two columns
            if grp==groups[0]:
                fout.write("\n1 0 "+str(df.loc[sess,covar+"Demeaned"])+" 0.0")
            else:
                fout.write("\n0 1 0.0 "+str(df.loc[sess,covar+"Demeaned"]))
                
    except Error as err:
        mylib.printError(f"an exception happened while creating {design}\n{err}")
        design = ""
    
    finally:        
        fout.close()
        
        if design!="":
            input(f"Design file created: {design}")
        
        return design

def createContrastFile_2grpCovarInteraction(cursor,outdir,groups,group_name):
    contrast = outdir+"/design.con"
    try:
        fout = open(contrast,'w')
        
        # NumWaves: # cols design.con, # cols design.mat: # EVs
        fout.write("/NumWaves 4\n")
        # NumContrasts: # lines design.con
        fout.write("/NumContrasts 2\n")    
        grp0 = mylib.getStringVal(cursor,group_name,groups[0]) if not isinstance(groups[0],str) else groups[0]
        grp1 = mylib.getStringVal(cursor,group_name,groups[1]) if not isinstance(groups[1],str) else groups[1]
        fout.write(f"/ContrastName1 slope_{grp0} > slope_{grp1}")
        fout.write(f"/ContrastName2 slope_{grp1} > slope_{grp0}")
        fout.write("/Matrix\n")
        
        fout.write("0 0 1 -1\n")
        fout.write("0 0 -1 1")
        
    except Error as err:
        mylib.printError(f"an exception happened while creating {contrast}\n{err}")
        contrast = ""
    
    finally:
        fout.close()
        
        if contrast!="":
            input(f"Contrast file created: {contrast}")
        
        return contrast

# Test if the linear relationship between the dependent var & covar differs between the 2 grps
# A significant interaction indicates that the grp difference varies as a function of the covar (i.e. age). In this case, only the inferences for the interaction are of interest.
# If the interaction is not significant, the grp mean differences can be obtained with unpairedTtest
def twoGrpCovarInteraction(cursor,outdir,group_name,tbss_name,df,covar,test_mode=False):
    groups = np.unique(df[group_name].astype(str).values)
    if len(groups)!=2:
        mylib.printError(f'A ttest must have two levels. There are {len(groups)} levels.')
        return False
    
    # Create the design matrix file design.mat
    design = createDesignFile_2grpCovarInteraction(outdir,df,groups,group_name,covar,test_mode)
    if not os.path.isfile(design):
        return False
        
    # Create the t-contrasts file design.con
    contrast = createContrastFile_2grpCovarInteraction(cursor,outdir,groups,group_name)
    if not os.path.isfile(contrast):
        return False
    
    if test_mode:
        input("Check {contrast} before it's deleted [Enter]")
        input("Check {design} before it's deleted [Enter]")
        
        os.system("rm "+design+" "+contrast)
    
    return True

# To create the df use the function mylib.getSummaryTable
def createDesignFiles(cursor,name,df,location,parent_pipe,test_mode=False):
    [test,covariates,independent_var] = mylib.getEntryFromTable(cursor,["test","covariates","independent_variable"],"tbss",["name"],[name])
    scratchdir = mylib.getScratchDirPipe(cursor,location,parent_pipe)
    outdir = f"{scratchdir}/{name}"
    covariates = [] if covariates==None else covariates.split(",")
    independent_var = [] if independent_var==None else independent_var.split(",")
    
    print("\n*** Creating "+test+" design files for "+name+" ***")
    
    if test=="1grp1cov":
        ok = singleGrpSingleCov(outdir,df,covariates[0],covariates[1:],test_mode)
    elif test=="anova_1way" or test=="ttest":
        ok = ANOVAcell_1way(cursor,outdir,independent_var[0],name,df,covariates,test_mode)
    elif test=="anova_2way" or test=="anova_2FxNL":
        ok = ANOVAcell_2way_2FxNL(cursor,outdir,name,df,independent_var,covariates,test_mode)
    elif test=="anova_3way":
        ok = ANOVAfactor_3way(cursor,outdir,name,df,independent_var,covariates,test_mode)
    elif test=="2grpCovInteraction":
        ok = twoGrpCovarInteraction(cursor,outdir,independent_var[0],name,df,covariates[0],test_mode)    
    else:
        ok = False
    
    # mergeFiles
    if ok:
        return mergeFiles(cursor,df,outdir,name,parent_pipe,location,test_mode)
        
    return False
        
def copyTBSSscript(infile,scratch,step,tbssname,img="FA",ext="sh"):
    if not os.path.isfile(infile):
        mylib.printError(infile+" doesnt exist")
    
    fname = step+'_'+tbssname
    try:
        fin = open(infile,'r')
        fout = open(scratch+'/'+fname+"."+ext,'w')
        
        for line in fin:
            if "job-name" in line:
                line = "#SBATCH --job-name=\""+fname+"\"\n"
            elif "img=FA" in line:
                line = "img="+img+'\n'
            elif "outDir=TBSS" in line:
                line = "outDir="+tbssname+'\n'
            elif "proj=proj" in line:
                line = "proj="+scratch.split("/")[-1]+'\n'
            elif "DONE "+step in line:
                line = "echo \"DONE "+fname+"\"\n"
            fout.write(line)
        
        fout.close()
        fin.close()
        
        return True
    
    except:
        mylib.printError("an exception happened creating "+fname+" in scratch")
        return False

def checkReadyPreStats(cnx,cursor,location,step,parent_pipe):
    print(f'\n*** {step} READY ***')
    
    scriptsdir = mylib.getScriptsDir(cursor,parent_pipe,location)
    scratch = mylib.getScratchDirPipe(cursor,location,parent_pipe)
    
    # Get the list of tbss that are ready to run preStats_part1/part2
    for name,line in mylib.getEntriesFromTable(cursor,["name","description","dependent_variable"],"tbss",[step],["ready"],index="name").iterrows():  
        print(f"\n{name}: {str(line['description'])}")
        running_dir = scratch+"/"+name
        print(f'running dir: {running_dir}\n')
        
        #### REMOVE FOLLOWING LINES AFTER preStats_part2.sh IS READY TO EXECUTE FOR MD ####
        if line["dependent_variable"].upper()!="FA" and step=="preStats_part2":
            mylib.updateEntriesFromTable(cursor,"tbss",["preStats_part2"],["error.preStatsPart2 onlyt runs on FA"],["name"],[name])
            mylib.updateDB(cnx,cursor)
            continue
        
        # Copy the script
        if not copyTBSSscript(scriptsdir+"/"+step+".py",running_dir,step,name,line["dependent_variable"],"py"):
            mylib.updateEntriesFromTable(cursor,"tbss",[step],["error.exception creating script"],["name"],[name])
            mylib.updateDB(cnx,cursor)
            continue
        
        # Run the job
        os.chdir(running_dir)
        
        if line["dependent_variable"]!="fa":
            input("Copy mean_FA from FA tbss into "+running_dir+" [Enter]")
            input("Copy mean_FA_skeleton_mask_dst from FA tbss into "+running_dir+" [Enter]")
            input("Copy all_FA from FA tbss into "+running_dir+" [Enter]")
        
        fname = step+"_"+name+".py"
        project = mylib.getValueFromTable(cursor,"project","pipelines",["pipename"],[parent_pipe])
        if step=="preStats_part1":
            cmd = "nohup python3 ./"+fname+" "+name+" "+project+" "+line["dependent_variable"].upper()+" "+step
        else:
            cmd = "nohup python3 ./"+fname+" "+name+" "+project+" "+step
        print(cmd)
        os.system(cmd)
        
        # Update the DB
        if os.path.isfile("nohup.out"):
            os.system(f"rm -f {fname}.log")
            print(f"mv nohup.out {fname}.log")
            os.system(f"mv nohup.out {fname}.log")
            print("done")
            mylib.updateEntriesFromTable(cursor,"tbss",[step],["done"],["name"],[name])
        else:
            mylib.updateEntriesFromTable(cursor,"tbss",[step],["error.nohup file not generated during local run"],["name"],[name])
        mylib.updateDB(cnx,cursor)
 
# It cannot run in automatic mode cause it has an input() to copy files
def checkReadyRandomise(cnx,cursor,location,parent_pipe):
    print("*** randomise READY ***")
    
    scratch = mylib.getScratchDirPipe(cursor,location,parent_pipe)
    scriptsdir = mylib.getScriptsDir(cursor,parent_pipe,location)
    
    # Get the list of tbss that are ready to run randomise
    for name,line in mylib.getEntriesFromTable(cursor,["name","description","test","dependent_variable"],"tbss",["randomise"],["ready"],["test"],"name").iterrows():
        print(f"\n{name}: {str(line['description'])}")
        
        if line["test"] in ["ttest","2grpCovInteraction"]:
            script = scriptsdir+"/27_randomise_uTtest.sh"
        elif line["test"]=="ftest" or ("anova" in line["test"]):
            script = scriptsdir+"/27_randomise_fTest.sh"
        elif line["test"]=="1grp1cov":
            script = scriptsdir+"/27_randomise_1GRP1COV.sh"
        else:
            mylib.printError(f"unrecognized test: {line['test']}")
            continue
            
        # Copy the script
        print(f"script: {script}")
        copyTBSSscript(script,scratch,"randomise",name,line["dependent_variable"].upper())
        
        if line["dependent_variable"]!="fa":
            input("Copy mean_FA_skeleton_mask from FA tbss folder [Enter]")
        
        # Submit the job
        fname = "randomise_"+name+".sh"
        jobID = ""
        if location=="hpc":
            os.chdir(scratch)
            print(f"qsub {fname}")
            jobID = mylib.systemOut("qsub "+fname)
            print(jobID)
        else:
            proj = scratch.split("/")[-1]
            print(f"cd /scratch/g/jbinder/mkeith/{proj}")
            print(f"qsub {fname}")
            jobID = input("jobID [Enter to skip submission]: ")
        if jobID=="":
            continue
        
        # Get the current time
        # Can be useful to know when I submitted in case is stuck in queue
        now_date = str(datetime.date.today().strftime("%Y_%m_%d"))
        now_time = str(datetime.datetime.today().strftime("%H:%M:%S"))
        now = now_date+'_'+now_time
            
        # Update the DB
        mylib.updateEntriesFromTable(cursor,"tbss",["randomise"],[f"job.{jobID}.{now}"],["name"],[name])
        mylib.updateDB(cnx,cursor)
    
def checkRunningTBSS(cnx,cursor,step,location,parent_pipe):
    print(f"*** {step} RUNNING ***")
    
    scratch = mylib.getScratchDirPipe(cursor,location,parent_pipe)
    scriptsdir = mylib.getScriptsDir(cursor,parent_pipe,location)
    
    # Get all the tbss where the step is running
    tbss = {}
    cursor.execute(f"select name,{step} from tbss where {step} like 'job.%'")
    for (name,info) in cursor:
        tbss[name] = info
    
    for name,info in tbss.items():
        jobID = info.split('.')[1]
        script = scratch+'/'+step+'_'+name+".sh"
        elog = script+".e"+jobID
        olog = script+".o"+jobID
        print(name)
        
        # Check if the script ever start running
        if not os.path.isfile(olog) or not os.path.isfile(elog):
            mylib.printError("Log files not found:")
            print(f"{elog}\n{olog}")
            continue
        
        # Check if the script finished running
        if fileslib.lineInFile(olog,"DONE "+step+'_'+name):
            # Take note of any warnings
            warning1 = fileslib.warningInLog(elog)
            warning2 = fileslib.warningInLog(olog)
            if warning1!="":
                mylib.updateNotesTBSS(cursor,name,warning1)
            elif warning2!="":
                mylib.updateNotesTBSS(cursor,name,warning2)
                
            # Move the logs to the project folder and the script to the TBSS folder
            os.system(f"mv {elog} {olog} {scriptsdir}/logs")
            os.system(f"mv {script} {scratch}/{name}")
                
            # Mark as done
            mylib.updateEntryFromTable(cursor,"tbss",[step],[f"done.{jobID}"],["name"],[name])
            print("DONE")          
                
        # Check if the log file has an error
        else:
            err1 = fileslib.errorInLog(elog).replace("'","")
            err2 = fileslib.errorInLog(olog).replace("'","")
            if err1!='':
                mylib.updateEntryFromTable(cursor,"tbss",[step],[f"error.{jobID}.{err1}"],["name"],[name])
            elif err2!='':
                mylib.updateEntryFromTable(cursor,"tbss",[step],[f"error.{jobID}.{err2}"],["name"],[name])
                    
        mylib.updateDB(cnx,cursor)
    
def checkDonepreStatsp1(cnx,cursor,location,parent_pipe):
    print("*** preStats_part1 DONE ***")
    
    # Get all the tbss where the preStats_part1 is done
    scratch = mylib.getScratchDirPipe(cursor,location,parent_pipe)
    for name,line in mylib.getEntriesFromTable(cursor,["name","dependent_variable"],"tbss",["preStats_part1"],["done"],index="name").iterrows():
        outdir = scratch+'/'+name
        line['dependent_variable'] = line['dependent_variable'].upper()
        print(f"{name} {line['dependent_variable']}")
        
        # Check that all the outputs are present
        if line['dependent_variable']=="FA":
            if not os.path.isfile(outdir+"/mean_FA.nii.gz") or not os.path.isfile(outdir+"/mean_FA_mask.nii.gz") or not os.path.isfile(outdir+"/mean_FA_skeleton.nii.gz"):
                mylib.updateEntriesFromTable(cursor,"tbss",["preStats_part1"],["error.output file missing"],["name"],[name])
                mylib.updateDB(cnx,cursor)
                continue
        elif not os.path.isfile(outdir+"/all_"+line['dependent_variable']+"_skeletonised.nii.gz"):
            mylib.updateEntriesFromTable(cursor,"tbss",["preStats_part1"],["error.output file missing"],["name"],[name])
            mylib.updateDB(cnx,cursor)
            continue
        
        # Check that the number of subjects is the same in the list and merged file
        n_list = mylib.getNumberEntries(cursor,"sessions",[name],[1])
        n_merged = mrilib.getDims(outdir+"/all_"+line['dependent_variable']+".nii.gz")[3]
        if n_merged!=n_list:
            mylib.updateEntriesFromTable(cursor,"tbss",["preStats_part1"],["error.wrong number of subjects"],["name"],[name])
            mylib.updateDB(cnx,cursor)
            continue
        
        if line['dependent_variable']!="FA":
            n_merged2 = mrilib.getDims(outdir+"/all_"+line['dependent_variable']+"_skeletonised.nii.gz")[3]
            if n_list!=n_merged2:
                mylib.updateEntriesFromTable(cursor,"tbss",["preStats_part1"],["error.wrong number of subjects"],["name"],[name])
                mylib.updateDB(cnx,cursor)
                continue
        
        # Get the execution time & update DB
        log = scratch+"/"+name+"/preStats_part1_"+name+".py.log"
        if os.path.isfile(log):
            # I could read days too although this should never take even a day
            exec_time = mylib.systemOut("cat "+log+" | tail -n1 | awk '{print $5,$7,$9}'",False).split(" ")
            sec = str(int(exec_time[2])+int(exec_time[1])*60+int(exec_time[0])*120)
            mylib.updateEntriesFromTable(cursor,"tbss",["preStats_part1"],[f"checked.{sec}sec"],["name"],[name])
        else:
            mylib.updateEntriesFromTable(cursor,"tbss",["preStats_part1"],["checked"],["name"],[name])
            
        # Mark preStats2 as ready (for FA) and randomise as ready (for all other imgs)
        if line['dependent_variable']=="FA":
            mylib.updateEntriesFromTable(cursor,"tbss",["preStats_part2"],["ready"],["name"],[name])
        else:
            mylib.updateEntriesFromTable(cursor,"tbss",["randomise"],["ready"],["name"],[name])
        mylib.updateDB(cnx,cursor)
        
def checkDonepreStatsp2(cnx,cursor,location,parent_pipe):
    print("*** preStats_part2 DONE ***")
    
    # Get all the tbss where the preStats_part2 is done
    scratch = mylib.getScratchDirPipe(cursor,location,parent_pipe)
    for name,line in mylib.getEntriesFromTable(cursor,["name","dependent_variable"],"tbss",["preStats_part2"],["done"],index="name").iterrows():
        line["dependent_variable"] = line["dependent_variable"].upper()
        print(f"{name} {line['dependent_variable']}")
        
        if line["dependent_variable"]!="FA":
            mylib.printError(f"{name} shouldnt run preStats_part2. This is only for FA!")
            mylib.updateEntriesFromTable(cursor,"tbss",["preStats_part2"],["error.shouldnt run preStats_part2. This is only for FA!"],["name"],[name])
            mylib.updateDB(cnx,cursor)
            continue
    
        # Check that all the outputs are present
        outdir = scratch+'/'+name
        outputs = [outdir+"/mean_FA_skeleton_mask.nii.gz",outdir+"/mean_FA_skeleton_mask_dst.nii.gz",outdir+"/all_FA_skeletonised.nii.gz"]
        OK = True
        for fout in outputs:
            if not os.path.isfile(fout):
                mylib.updateEntriesFromTable(cursor,"tbss",["preStats_part2"],["error.output file missing"],["name"],[name])
                OK = False
                break
        if not OK:
            continue
        
        # Get the execution time & update the DB
        log = scratch+"/"+name+"/preStats_part2_"+name+".py.log"
        if os.path.isfile(log):
            # I could read days too although this should never take even a day
            exec_time = mylib.systemOut("cat "+log+" | tail -n1 | awk '{print $7,$9,$11}'",False).split(" ")
            sec = str(int(exec_time[2])+int(exec_time[1])*60+int(exec_time[0])*120)
            mylib.updateEntriesFromTable(cursor,"tbss",["preStats_part2"],[f"checked.{sec}sec"],["name"],[name])
        else:
            mylib.updateEntriesFromTable(cursor,"tbss",["preStats_part2"],["checked"],["name"],[name])
        
        # Mark as ready randomise
        mylib.updateEntriesFromTable(cursor,"tbss",["randomise"],["ready"],["name"],[name])
        mylib.updateDB(cnx,cursor)
        
# num_contrast: total number of contrasts in the fts file (# of f-test)
# n_contrast: current contrast number
def getRawRandomEffects(outdir,prefix,stat,n_contrast,num_contrast):
    raw1 = outdir+'/'+prefix+'_'+stat+str(n_contrast)+".nii.gz"
    if not os.path.isfile(raw1) or mrilib.getnvoxels(raw1)==0:
        mylib.printError(f'raw statistic {raw1} missing or empty')
        return ""
    
    # The last contrast image doesn't change the denominator
    # If there's more than 3 contrast, the main contrast (non-interactions, n_contrast<=3) don't change the denominator
    if n_contrast==num_contrast or (num_contrast>3 and n_contrast<=3):
        return raw1
    
    # The denominator would be the last contrast which is the highest interaction
    raw2 = outdir+'/'+prefix+'_'+stat+str(num_contrast)+".nii.gz"
    if not os.path.isfile(raw2) or mrilib.getnvoxels(raw2)==0:
        mylib.printError(f'raw statistic {raw2} missing or empty')
        return ""
    
    raw = outdir+'/'+prefix+'_'+stat+str(num_contrast+1)+".nii.gz"
    os.system(f"fslmaths {raw1} -div {raw2} {raw}")
    return raw

# num_contrast: total number of contrasts in the fts file (# of f-test)
# n_contrast: current contrast number
def getRawMixedEffects(outdir,prefix,stat,n_contrast,num_contrast,contrast_name,random_factors,ind_var):
    # In the mixed effects model, the random effects stay the same & only the fixed factors change the denominator
    # So, for those contrast that contain all the random effects, just return the corresponding raw image
    random_factors = random_factors.split(",")
    if (set(random_factors).issubset(set(contrast_name.split("_")))):
        return outdir+'/'+prefix+'_'+stat+str(n_contrast)+".nii.gz"
    else:
        # The exception to the above is for the fixed factor when there's only one fix factor & all the others are random
        nFixed = len(ind_var.split(","))-len(random_factors)
        if nFixed==1 and ("_" not in contrast_name) and (contrast_name not in random_factors):
            return outdir+'/'+prefix+'_'+stat+str(n_contrast)+".nii.gz"
        
        if num_contrast<4 or n_contrast>3:
            return getRawRandomEffects(outdir,prefix,stat,n_contrast,num_contrast)
        
        # In these cases the denominator is fstat6 instead of fstat7
        if n_contrast>1:
            return getRawRandomEffects(outdir,prefix,stat,n_contrast,num_contrast-1)
        
        # In this case the denominator is fstat5
        return getRawRandomEffects(outdir,prefix,stat,n_contrast,num_contrast-2)
    
def getResultsCols(cursor):
    n_results = len([x for x in mylib.getColsFromTable(mylib.getTableInfo(cursor,"tbss_results")) if x.startswith("results")])
    res_cols = []
    for i in range(1,n_results+1):
        res_cols+=[f"results{i}"]
    return res_cols

def checkDoneRandomise(cnx,cursor,automatic,location,parent_pipe):
    print("*** randomise DONE ***")
    
    scratch = mylib.getScratchDirPipe(cursor,location,parent_pipe)
    scriptsdir = mylib.getScriptsDir(cursor,parent_pipe,location)
    
    ## Get all the tbss where the randomise is done
    for name,line in mylib.getEntriesFromTable(cursor,["name","randomise","effects","random_factors","independent_variable"],"tbss",["randomise"],["done.%"],index="name").iterrows():
        jobID = str(line["randomise"]).split('.')[1]
        outdir = f"{scratch}/{name}"
        script = f"{outdir}/randomise_{name}.sh"
        confile = f"{outdir}/design.con"
        ftsfile = f"{outdir}/design.fts"
        stats = ["fstat","tstat"] if os.path.isfile(ftsfile) else ["tstat"]
        
        if os.path.isfile(script):
            prefix = mylib.systemOut("grep 'randomise -i' "+script+" | awk '{print $5}'")
        else:
            mylib.printError("Could not find script: "+script)
            continue
        
        print(f'\n{name}')
        OK = True
        for stat in stats:
            # Number of contrasts
            n_con = int(mylib.systemOut("grep /NumContrasts "+ftsfile+" | awk '{print $2}'")) if stat=="fstat" else int(mylib.systemOut("grep /NumContrasts "+confile+" | awk '{print $2}'"))
            print(f'\n{n_con} contrasts for the {stat}')
            
            ## Check each contrast
            for i in range(1,n_con+1):
                contrast_name = mylib.systemOut("grep \"/ContrastName"+str(i)+" \" "+confile+" | awk '{print $2}'") if stat=="tstat" else mylib.systemOut("grep \"/ContrastName"+str(i)+" \" "+ftsfile+" | awk '{print $2}'")
                print(f'\nContrast: {contrast_name}')
                
                ## Raw statistic
                
                # It is not an ANOVA or it is an ANOVA with fixed effects
                # Fixed effects assumes that the explanatory (independent) variable has a fixed or constant relationship with the response (dependent) variable across all observations.
                # All the EVs have pre-determined categories & the response is made for the categories of the EVs used in the model (even if the experiments are repeated multiple times, each EV will have the same categories in each experiment)
                # Examples of fixed effects would be two treatments or diseases
                if stat=="tstat" or line["effects"]==None or line["effects"]=="fixed":
                    raw = outdir+'/'+prefix+'_'+stat+str(i)+".nii.gz"
                
                else:
                    # ANOVA random effects: variables that vary over time or from one sbj to another, or there could be more levels other than the ones included in the study & the outcome of the model is supposed to be applied to all levels of the variable (even the ones not included in the study)
                    # Example of random effects would be disease scores
                    if line["effects"]=="random":
                        raw = getRawRandomEffects(outdir,prefix,stat,i,n_con)
                            
                    # ANOVA mixed effects
                    else:
                        # I use the con file instead of the fts because the fts will not have the contrast names
                        # contrast_name will contain the name of the factor
                        raw = getRawMixedEffects(outdir,prefix,stat,i,n_con,contrast_name,str(line["random_factors"]),line["independent_variable"])
                            
                if not os.path.isfile(raw) or mrilib.getnvoxels(raw)==0:
                    mylib.printError(f'raw statistic {raw} missing or empty')
                    OK = False
                    break
                print(f'raw statistic: {raw}')
                
                ## tfce output: map with corrected p-values
                tfce = outdir+'/'+prefix+"_tfce_corrp_"+stat+str(i)+".nii.gz"
                if not os.path.isfile(tfce):
                    mylib.printError(f'tfce statistic {tfce} missing')
                    OK = False
                    break
                print(f'corrected p-values: {tfce}')
                    
                ## Get the column where the result should be saved
                for result in getResultsCols(cursor):
                    if mylib.getNumberEntries(cursor,"tbss_results",["name",result],[name,""])>0:
                        break
                    
                ## If tfce corrected has no significant voxels, there is nothing to check
                max_val = float(mylib.systemOut("fslstats "+tfce+" -R | awk '{print $2}'"))
                if max_val<0.95:
                    if contrast_name!="":
                        OK = mylib.updateEntryFromTable(cursor,"tbss_results",[result],[f"no significance, min p={round(1-max_val,3)} ({stat}{i}:{contrast_name})"],["name"],[name])==1
                    else: 
                        OK = mylib.updateEntryFromTable(cursor,"tbss_results",[result],[f"no significance, min p={round(1-max_val,3)} ({stat}{i})"],["name"],[name])==1
                    if OK:
                        continue
                    else:
                        break
                    
                ## Threshold and add to array of images to visualize if OK
                # The raw statistic (tstat or fstat) is thresholded with the significant voxels from the tfce corrected (corrected p-values map)
                print("Thresholding raw statistic with tfce significant voxels...")
                os.system(f"fslmaths {tfce} -thr 0.95 -bin -mul {raw} "+raw.replace(".nii.gz","_sign5.nii.gz"))
                    
                if max_val>=0.99:
                    os.system(f"fslmaths {tfce} -thr 0.99 -bin -mul {raw} "+raw.replace(".nii.gz","_sign1.nii.gz"))
                    n = mrilib.getnvoxels(raw.replace(".nii.gz","_sign1.nii.gz"))
                    if contrast_name!="":
                        mylib.updateEntriesFromTable(cursor,"tbss_results",[result],[f"{n} significant voxels at .001 thr ({stat}{i}:{contrast_name})"],["name"],[name])
                    else:
                        mylib.updateEntriesFromTable(cursor,"tbss_results",[result],[f"{n} significant voxels at .001 thr ({stat}{i})"],["name"],[name])
                else:
                    n = mrilib.getnvoxels(raw.replace(".nii.gz","_sign5.nii.gz"))
                    if contrast_name!="":
                        mylib.updateEntriesFromTable(cursor,"tbss_results",[result],[f"{n} significant voxels at .05 thr ({stat}{i}:{contrast_name})"],["name"],[name])
                    else:
                        mylib.updateEntriesFromTable(cursor,"tbss_results",[result],[f"{n} significant voxels at .05 thr ({stat}{i})"],["name"],[name])
            
            if not OK:
                break
            
        if not OK:
            mylib.updateEntriesFromTable(cursor,"tbss",["randomise"],[f"error.{jobID}.output file missing"],["name"],[name])
            mylib.updateDB(cnx,cursor)
            continue
        
        ## Get the execution time & update the DB
        log = scriptsdir+"/logs/randomise_"+name+".sh.o"+jobID
        if os.path.isfile(log):
            exec_time = mylib.systemOut("cat "+log+" | tail -n1 | awk '{print $5,$7,$9}'",False).split(" ")
            sec = str(int(exec_time[2])+int(exec_time[1])*60+int(exec_time[0])*120)
            mylib.updateEntriesFromTable(cursor,"tbss",["randomise"],[f"checked.{jobID}.{sec}sec"],["name"],[name])
        else:
            mylib.updateEntriesFromTable(cursor,"tbss",["randomise"],[f"checked.{jobID}"],["name"],[name])
        
        mylib.updateDB(cnx,cursor)
        
def getMeanFromCluster(cursor,location,tbss_name,parent_pipe,sess):
    # Get the TBSS folder
    tbss_scratch = mylib.getScratchDirPipe(cursor,location,parent_pipe)+'/'+tbss_name
    if not os.path.isdir(tbss_scratch):
        tbss_scratch = input("Location of tbss folder for "+tbss_name+": ")
    if not os.path.isdir(tbss_scratch):
        return -1
    
    # Get the cluster (choose from sign files)
    sign_clusters = glob.glob(tbss_scratch+"/*sign*")
    print("\nSelect one cluster:")
    for cluster in sign_clusters:
        print(cluster)
    cluster = input(">> ")
    if not os.path.isfile(cluster):
        return -1
    
    # Get the fa img
    img = mylib.getValueFromTable(cursor,"dependent_variable","tbss",["name"],[tbss_name])
    fa = f"{tbss_scratch}/{img.upper()}/{sess}.nii.gz"
    if not os.path.islink(fa) and not os.path.isfile(fa):
        return -1
    
    # Get the mean value
    os.system("fslmaths "+cluster+" -bin -mul "+fa+" tmp.nii.gz")
    M = mrilib.getMean("tmp.nii.gz")
    os.system("rm tmp.nii.gz")
    
    return M

def printTBSSresults(cursor,name):
    df = mylib.getColumnsFromTable(cursor,getResultsCols(cursor),"tbss_results",[name])
    for idx,line in df.iterrows():
        result = line[name]
        if result==None:
            break
        print(idx+": "+result)

def printTBSSinfo(cursor,tbss_name):
    [preStats_part1,preStats_part2,randomise] = mylib.getEntryFromTable(cursor,["preStats_part1","preStats_part2","randomise"],"tbss",["name"],[tbss_name])
    if preStats_part1==None:
        status = "hasn't start"
    elif not preStats_part1.startswith("checked"):
        status = "preStats_part1:"+preStats_part1
    elif not preStats_part2.startswith("checked"):
        status = "preStats_part2:"+preStats_part2
    elif not randomise.startswith("checked"):
        status = "randomise:"+randomise
    else:
        status = "finished running!\n"
    print(f'{tbss_name} {status}\n')
    
    cols = [col for col in mylib.getColsFromTable(mylib.getTableInfo(cursor,"tbss")) if col not in ["name","preStats_part1","preStats_part2","randomise"]]
    df = mylib.getColumnsFromTable(cursor,cols,"tbss",[tbss_name])
    if len(df)>0:
        print(df)
        
        if str(randomise).startswith("checked"):
            print("\nResults:")
            printTBSSresults(cursor,tbss_name)
    
def menu():
    print("\nChoose an option:\n"
          "0. Run all (automatic mode)\n"
          "1. Exit\n"
          "2. Check ready\n"
          "3. Check running\n"
          "4. Check done\n"
          "5. Rename TBSS\n"
          "6. Delete TBSS\n"
          "7. Clear TBSS\n"
          "8. Get TBSS info")
    
    return input(">> ")

def runPipeline(location):
    print ("### TBSS processing ###")
           
    database = input("Database: ")
    [username,password,hostdir] = ECP.getUserPassHost(location)
    cnx = mylib.connect(username,password,hostdir,database)
    cursor = cnx.cursor()
    pipes = mylib.getColumnFromTable(cursor,"pipename","pipelines")
    cursor.close()
    cnx.close()
    [username,password,hostdir] = ECP.getUserPassHost(location)
    
    if len(pipes)==0:
        return
    if len(pipes)==1:
        parent_pipe = pipes[0]
    else:
        parent_pipe = input("Select one pipe from ["+",".join(pipes)+"]: ")
    if parent_pipe not in pipes:
        mylib.printError(parent_pipe+" is not a valid pipeline")
        return
           
    opt = menu()
    while opt!="1":
        try:
            cnx = mylib.connect(username,password,hostdir,database)
            cursor = cnx.cursor()
            
            if opt=='0':
                print("*** Running in automatic mode, use Ctrl+C to end ***")
                importlib.reload(mylib)
                importlib.reload(mrilib)
                
            if opt in ['0','2']:
                checkReadyPreStats(cnx,cursor,location,"preStats_part1",parent_pipe)
                checkReadyPreStats(cnx,cursor,location,"preStats_part2",parent_pipe)
                checkReadyRandomise(cnx,cursor,location,parent_pipe,opt=='0')
                
            if opt in ['0','3']:
                for step in ["preStats_part1","preStats_part2","randomise"]:
                    checkRunningTBSS(cnx,cursor,step,location,parent_pipe)
                    
            if opt in ['0','4']:
                checkDonepreStatsp1(cnx,cursor,location,parent_pipe)
                checkDonepreStatsp2(cnx,cursor,location,parent_pipe)
                checkDoneRandomise(cnx,cursor,opt=='0',location,parent_pipe)
                
            elif opt=='5':
                name = input("Current name: ")
                new_name = input("New name: ")
                renameTBSS(cnx,cursor,name,new_name,location,parent_pipe)
                
            elif opt=='6':
                name = input("Name of the tbss you whish to delete: ")
                delTBSS(cnx,cursor,name,location,parent_pipe)
                
            elif opt=='7':
                name = input("Name of the tbss you whish to clear: ")
                clearTBSS(cnx,cursor,name,location,parent_pipe)
                
            elif opt=='8':
                tbss_name = input("Name of the tbss to print info from: ")
                printTBSSinfo(cursor,tbss_name)
        
        except Error as err:
            if err.errno == errorcode.ER_LOCK_WAIT_TIMEOUT or err.errno == errorcode.ER_LOCK_DEADLOCK:
                mylib.printError("Database deadlock exception. Some changes might have not been saved in the DB.")
                continue
            else:
                mylib.printError(err)
                break
        
        finally:
            cursor.close()
            cnx.close()
            
        if opt=='0':
            # Wait one hour between automatic loops
            print(f'Went to sleep at {datetime.datetime.today().strftime("%H:%M:%S")}')
            print("ZZZZZZ...")
            time.sleep(3600)
            print("Wake up!")
            
    print("Good bye")
    
def main():
        try:
            location = input("Location? [hpc/petrov/petrov_damaged*/oneDrive]: ")
            if location=='':
                location = "petrov_damaged"
                
            runPipeline(location)
            
        except Error as err:
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                mylib.printError("Wrong username or password")
            else:
                mylib.printError(err)

if __name__ == "__main__":
        main()
