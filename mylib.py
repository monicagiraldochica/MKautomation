#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Monica Keith"
__status__ = "Production"

import mysql.connector
from mysql.connector import Error
from mysql.connector import errorcode
import os
from PyPDF2 import PdfFileMerger, PdfFileReader
import re
import glob
import time
import datetime
import pandas as pd
import fileslib
import statslib
import zipfile
import numpy as np
import subprocess
import ECP

empty_values = ["nan","NaT","None","\\N","NULL","NaN","none","nat",""]

def connect(username,password,hostdir="localhost",DB=""):
    if DB!="":
        cnx = mysql.connector.connect(user=username, password=password, host=hostdir, database=DB, allow_local_infile=True)
    else:
        cnx = mysql.connector.connect(user=username, password=password, host=hostdir, allow_local_infile=True)
    
    cursor = cnx.cursor()
    # Set the number of seconds the server waits for activity on a connection before closing it to 86400 (24 hours)
    cursor.execute("SET GLOBAL wait_timeout=86400")
    # Set the number of seconds to wait for a lock timeout to 1 minutes
    cursor.execute("SET GLOBAL innodb_lock_wait_timeout=60")
    cursor.close()
    
    return cnx

def getTablesDB(cursor,database=""):
    cursor.execute("show tables") if database=="" else cursor.execute(f"show tables from {database}")
    return [item[0] for item in cursor.fetchall()]

# Table can have only the table name if it is from the db the cursor is connected to
# Otherwise include db_name (db_name.table)
def getTableInfo(cursor,table):
    cursor.execute("show columns from "+table)
    df = pd.DataFrame(cursor.fetchall())
    if len(df)==0:
        return df
    
    df.columns = ["Field","Type","Null","Key","Default","Extra"]
    df.set_index("Field",inplace=True)
    return df

def getPrimaryKeys(df_table_info):
    if len(df_table_info)==0:
        return []
    
    return df_table_info[df_table_info["Key"]=="PRI"].index.values.tolist()
    
def getColsFromTable(df_table_info):
    if len(df_table_info)==0:
        return []
    
    return df_table_info["Key"].index.values.tolist()
    
def getStringFromTable(cursor,column_out,table,columns_in=[],values_in=[]):
    return str(getValueFromTable(cursor,column_out,table,columns_in,values_in))

def getValueFromTable(cursor,column_out,table,columns_in=[],values_in=[]):
    entries = getEntryFromTable(cursor,[column_out],table,columns_in,values_in)
    return entries[0] if len(entries)>0 else None

def getEntryFromTable(cursor,columns_out,table,columns_in=[],values_in=[]):
    entries = getEntriesFromTable(cursor,columns_out,table,columns_in,values_in)
    return entries.iloc[0,:].values.tolist() if len(entries)>0 else [None]*len(columns_out)

def getColumnFromTable(cursor,column_out,table,columns_in=[],values_in=[],order_by=[],distinct=False):
    entries = getEntriesFromTable(cursor,[column_out],table,columns_in,values_in,order_by,distinct=distinct)
    return entries.iloc[:,0].values.tolist() if len(entries)>0 else []

def whereStatement(columns_in,values_in):
    if len(columns_in)==0 or (len(columns_in)!=len(values_in)):
        return ""
    
    values_in = [str(val) if str(val) not in empty_values else "" for val in values_in]
    
    if values_in[0]=="":
        cmd = " where "+columns_in[0]+" is null"
    elif values_in[0]=="*":
        cmd = " where "+columns_in[0]+" is not null"
    elif "%" in values_in[0]:
        if not values_in[0].startswith("!"):
            cmd = " where "+columns_in[0]+" like '"+values_in[0]+"'"
        else:
            cmd = " where "+columns_in[0]+" not like '"+values_in[0][1:]+"'"
    elif values_in[0].startswith("(") and values_in[0].endswith(")"):
        cmd = " where "+columns_in[0]+" in "+values_in[0]
    elif values_in[0].startswith("!"):
        cmd = " where "+columns_in[0]+"!='"+values_in[0][1:]+"'"
    else:
        cmd = " where "+columns_in[0]+"='"+values_in[0]+"'"
        
    for i in range(1,len(columns_in)):
        if values_in[i]=="":
            cmd+=" and "+columns_in[i]+" is null"
        elif values_in[i]=="*":
            cmd+=" and "+columns_in[i]+" is not null"
        elif "%" in values_in[i]:
            if not values_in[i].startswith("!"):
                cmd+=" and "+columns_in[i]+" like '"+values_in[i]+"'"
            else:
                cmd+=" and "+columns_in[i]+" not like '"+values_in[i][1:]+"'"
        elif values_in[i].startswith("(") and values_in[i].endswith(")"):
            cmd+=" and "+columns_in[i]+" in "+values_in[i]
        elif values_in[i].startswith("!"):
            cmd+=" and "+columns_in[i]+"!='"+values_in[i][1:]+"'"
        else:
            cmd+=" and "+columns_in[i]+"='"+values_in[i]+"'"
            
    return cmd

# key_values: the primary key value for the entry for which I want to obtain the values of the columns
# since it is possible that a table has more than one key, this is an array that most times will be of length 1
def getColumnsFromTable(cursor,columns,table,key_values):
    prim_keys = getPrimaryKeys(getTableInfo(cursor,table))
    if len(prim_keys)!=len(key_values):
        printError(f"{table} has {len(prim_keys)} but {len(key_values)} were entered")
        return pd.DataFrame()
    
    df = getEntriesFromTable(cursor,columns,table,prim_keys,key_values)
    if len(df)==0:
        printError("_".join(key_values)+" is not a valid key for "+table)
        return df
    
    df = df.T
    df.columns = ["_".join(key_values)]
    return df

def getEntriesFromTable(cursor,columns_out,table,columns_in=[],values_in=[],order_by=[],index="",distinct=False):
    df = pd.DataFrame()
    if len(columns_in)!=len(values_in):
        printError('size of columns_in and values_in dont match')
        return df
    
    if len(columns_out)>0:
        cmd = "select "+",".join(columns_out)+" from "+table if not distinct else "select distinct "+",".join(columns_out)+" from "+table
    else:
        cmd = "select count(*) from "+table
    cmd+=whereStatement(columns_in,values_in)
    
    try:
        cursor.execute(cmd)
        df = pd.DataFrame(cursor.fetchall())
        if len(df)>0:
            if len(columns_out)>0:
                df.columns = columns_out
            if len(order_by)>0:
                df.sort_values(order_by,inplace=True,ignore_index=True)
            if index!="":
                df.set_index(index,inplace=True)
        return df

    except Error as err:
        printError(err)
        return df

def getNumberEntries(cursor,table,columns_in=[],values_in=[]):
    return getEntriesFromTable(cursor,[],table,columns_in,values_in).iloc[0,0]
    
def recoveryFile2dataFrame(cursor,table,location="petrov"):
    cursor.execute("select database()")
    database = cursor.fetchone()[0]
    
    db_dir = getValueFromTable(cursor,"db_dir","generalInfo.locations",["location"],[location])
    recov_file = f"{db_dir}/recovery/{database}.{table}.txt"
    if not os.path.isfile(recov_file):
        return pd.DataFrame()
    
    # Create the data frame
    df = fileslib.fileToDataFrame(recov_file,True,"\t")
    
    # Update the titles of the columns
    if len(df.columns)!=len(getColsFromTable(getTableInfo(cursor,table))):
        return pd.DataFrame()
    
    return df

def createRecoveryFiles(cursor,location):
    cursor.execute("select database()")
    db = cursor.fetchone()[0]
    
    table_cols = getColumnsDB(cursor)
    db_dir =  getValueFromTable(cursor,"db_dir","generalInfo.locations",["location"],[location])
    recovery = f"{db_dir}/recovery/"
    
    if os.path.isdir(recovery):    
        fout = open(recovery+"recreateTables.sql",'a')
        fout.write("\n")
        for cmd in commandsNewTables(db):
            fout.write("\n"+cmd+";")
        fout.close()
        
        fout = open(recovery+"reloadTables.sql",'a')
        for table,cols in table_cols.items():
            load = "\nload data local infile '~/Insync/OneDrive/DB/recovery/"+db+"."+table+".txt' into table "+db+"."+table+" fields optionally enclosed by '\"' ("+cols[0]
            for i in range(1,len(cols)):
                load+=","+cols[i]
            load+=");"
            fout.write(load)
        fout.close()
    else:
        printError("recreateTables.sql and reloadTables.sql cannot be updated. Recovery folder not found.")

def updateEntryFromTable(cursor,table,columns_edit,values_edit,columns_in=[],values_in=[],test_mode=False):
    n = updateEntriesFromTable(cursor,table,columns_edit,values_edit,columns_in,values_in,test_mode=test_mode)
    if n>1:
        printError(str(n)+" records where updated. rolling back.")
        return -1
    return n

def updateEntriesFromTable(cursor,table,columns_edit,values_edit,columns_in=[],values_in=[],test_mode=False):
    if len(columns_edit)!=len(values_edit) or len(columns_edit)==0:
        printError("size of columns_edit and values_edit dont match or are empty")
        return -1
    
    # Make sure the table exists
    table_info = getTableInfo(cursor,table)
    if len(table_info)==0:
        printError(table+" table doesnt exist")
        return -1
    
    # Make sure all values in column_in are a column in the table
    if not set(columns_edit)<=set(getColsFromTable(table_info)):
        printError("one or more values in ["+",".join(columns_edit)+"] are not valid fields for table "+table)
        return -1
    
    # Make sure its not trying to erase values that cant be None and that values are in string mode
    values_edit = [str(val) if str(val) not in empty_values else "" for val in values_edit]
    none_problem = [columns_edit[i] for i in range(len(values_edit)) if values_edit[i]=="" and table_info.loc[columns_edit[i],"Null"]=="NO"]
    if len(none_problem)>0:
        printError("trying to assign an empty value to a column that cannot be None: "+",".join(none_problem))
        return -1
    
    # Create the command
    if values_edit[0]!="":
        cmd = "update "+table+" set "+columns_edit[0]+"='"+values_edit[0]+"'"
    else:
        cmd = "update "+table+" set "+columns_edit[0]+"=NULL"
        
    for i in range(1,len(columns_edit)):
        if values_edit[i]!="":
            cmd+=", "+columns_edit[i]+"='"+values_edit[i]+"'"
        else:
            cmd+=", "+columns_edit[i]+"=NULL"
            
    cmd+=whereStatement(columns_in,values_in)
    
    # Print if in test mode
    if test_mode:
        print(cmd)
        return 0
    
    # Execute if not in test mode
    try:
        cursor.execute(cmd)
        return cursor.rowcount
    except Error as err:
        printError(err)
        return -1

def insertEntryToTable(cursor,table,columns_add,values_add,test_mode=False):
    if len(columns_add)!=len(values_add) or len(columns_add)==0:
        printError("size of columns_edit and values_edit dont match or are empty")
        return -1
    
    # Make sure the table exists
    table_info = getTableInfo(cursor,table)
    if len(table_info)==0:
        printError(table+" table doesnt exist")
        return -1
    
    # Make sure all values in columns_add are a column in the table
    if not set(columns_add)<=set(getColsFromTable(table_info)):
        printError("one or more values in ["+",".join(columns_add)+"] are not valid fields for table "+table)
        return -1
    
    # Make sure its not trying to add NULL where that cant happen
    values_add = ["'"+str(val)+"'" if str(val) not in empty_values else "NULL" for val in values_add ]
    none_problem = [columns_add[i] for i in range(len(values_add)) if values_add[i]=="NULL" and table_info.loc[columns_add[i],"Null"]=="NO"]
    if len(none_problem)>0:
        printError("trying to assign an empty value to a column that cannot be NULL: "+",".join(none_problem))
        return -1
    
    # Create the command
    cmd = "insert into "+table+"("+",".join(columns_add)+") values("+",".join(values_add)+")"
    
    # Print if in test mode
    if test_mode:
        print(cmd)
        return 0
    
    # Execute if not in test mode
    try:
        cursor.execute(cmd)
        return cursor.rowcount
    except Error as err:
        printError(err)
        return -1

def commandsNewTables(db):
    # not creating excluded table by default
    # excluded should only be used when there's more than one img per sbj and one can exclude one of them without excluding the session or sbj
    # for example on the ECP one can exclude 75_AP without excluding the session
    # if excluding a session, can add the explanation in the notes column of sessions table
    # if excluding a subject, can add the explanation in the notes column of subjects table
    return [
            "CREATE TABLE IF NOT EXISTS "+db+".projects (title varchar(255) NOT NULL, PRIMARY KEY (title)) ENGINE=InnoDB",
            "CREATE TABLE IF NOT EXISTS "+db+".subjects (sbjID varchar(255) NOT NULL, grp enum('CON','PAT'), site varchar(255), notes varchar(255), MRIsbj tinyint(1), PRIMARY KEY (sbjID)) ENGINE=InnoDB",
            "CREATE TABLE IF NOT EXISTS "+db+".sessions (sess varchar(255) NOT NULL, sbjID varchar(255) NOT NULL, sessID varchar(255) NOT NULL, notes varchar(255), PRIMARY KEY (sess), KEY sbjID (sbjID), CONSTRAINT sessions_ibfk_1 FOREIGN KEY (sbjID) REFERENCES subjects (sbjID) ON DELETE CASCADE) ENGINE=InnoDB",
            "CREATE TABLE IF NOT EXISTS "+db+".pipelines (pipename varchar(255) NOT NULL, project varchar(255) NOT NULL, data_dir varchar(255) NOT NULL, sess_dir varchar(255) NOT NULL, scripts_dir varchar(255) NOT NULL, first_step varchar(255), PRIMARY KEY (pipename), KEY project (project), CONSTRAINT pipelines_ibfk_1 FOREIGN KEY (project) REFERENCES projects (title) ON DELETE CASCADE) ENGINE=InnoDB",
            "CREATE TABLE IF NOT EXISTS "+db+".pipesteps (step varchar(255) NOT NULL, description varchar(255), pipeline varchar(255), local_script varchar(255), script varchar(255), input_folder varchar(255), input_files varchar(550), input_files_opt varchar(550), output_folder varchar(255), output_files varchar(550), output_files_opt varchar(550), next_step varchar(1142), move_files tinyint(1) DEFAULT 1 NOT NULL, PRIMARY KEY (step), KEY pipeline (pipeline), CONSTRAINT pipesteps_ibfk_1 FOREIGN KEY (pipeline) REFERENCES pipelines (pipename) ON DELETE CASCADE) ENGINE=InnoDB",
            "CREATE TABLE IF NOT EXISTS "+db+".demographicData (sbjID varchar(255) NOT NULL, measure varchar(255) NOT NULL, value varchar(255) NOT NULL, PRIMARY KEY (sbjID,measure), CONSTRAINT demographicData_ibfk_2 FOREIGN KEY (sbjID) REFERENCES subjects (sbjID) ON DELETE CASCADE) ENGINE=InnoDB",
            "CREATE TABLE IF NOT EXISTS "+db+".behavioralData (sess varchar(255) NOT NULL, measure varchar(255) NOT NULL, value varchar(255) NOT NULL, PRIMARY KEY (sess,measure), CONSTRAINT behavioralData_ibfk_2 FOREIGN KEY (sess) REFERENCES sessions (sess) ON DELETE CASCADE) ENGINE=InnoDB",
            "CREATE TABLE IF NOT EXISTS "+db+".procs (sess varchar(255) NOT NULL, step varchar(255) NOT NULL, status enum('ready','hold','running','done','checked','error') NOT NULL DEFAULT 'ready', jobID INT, exectimesec INT, notes varchar(255), qsub_time varchar(255), cluster enum('old','hpc'), PRIMARY KEY (sess,step), KEY sess (sess), KEY step (step), CONSTRAINT procs_ibfk_1 FOREIGN KEY (sess) REFERENCES sessions (sess) ON DELETE CASCADE, CONSTRAINT procs_ibfk_2 FOREIGN KEY (step) REFERENCES pipesteps (step), unique key jobID (jobID)) ENGINE=InnoDB",
            "CREATE TABLE IF NOT EXISTS "+db+".results (sess varchar(255) NOT NULL, result varchar(255) NOT NULL, value varchar(800) NOT NULL, PRIMARY KEY (sess,result), KEY sess(sess), CONSTRAINT results_ibfk_1 FOREIGN KEY (sess) REFERENCES sessions (sess) ON DELETE CASCADE) ENGINE=InnoDB",
            "CREATE TABLE IF NOT EXISTS "+db+".string2num (variable varchar(255) NOT NULL, string_val varchar(255) NOT NULL, numerical_val int NOT NULL, notes varchar(255), PRIMARY KEY (variable,string_val)) ENGINE=InnoDB",
            "CREATE TABLE IF NOT EXISTS "+db+".measures_description (measure varchar(255) NOT NULL, description varchar(255) NOT NULL, from_table varchar(50) NOT NULL, PRIMARY KEY (measure)) ENGINE=InnoDB"
            ]

def updateDB(cnx,cursor):
    try:
        cnx.commit()
    except Error as err:
        if err.errno == errorcode.ER_LOCK_WAIT_TIMEOUT or err.errno == errorcode.ER_LOCK_DEADLOCK:
            printWarning("Database deadlock exception during commit. Will re-try in five minutes.")
            time.sleep(300)
            print("Re-trying commit...")
            try:
                cnx.commit()
            except Error:
                printError("Could not complete transaction")
                raise
        elif err.errno == errorcode.CR_SERVER_GONE_ERROR or err.errno == errorcode.CR_SERVER_LOST:
            printError("Connection lost before commiting queries. Will try to re-connect and commit.")
            try:
                getNumberEntries(cursor,"pipesteps")
                cnx.commit()
                print("Successfully re connected and commited!")
            except Error:
                printError("could not re-connect")
                raise
        else:
            printError("An exception happend in updateDB: "+str(err))
            raise

def insertBehavioralData(cursor,sess,measure,value):
    if inBehavioralData(cursor,sess,measure):
        updateEntriesFromTable(cursor,"behavioralData",["value"],[value],["sess","measure"],[sess,measure])
    else:
        insertEntryToTable(cursor,"behavioralData",["sess","measure","value"],[sess,measure,value])
 
def insertDemographicData(cursor,sbjID,measure,value):
    if inDemographicData(cursor,sbjID,measure):
        updateEntriesFromTable(cursor,"demographicData",["value"],[value],["sbjID","measure"],[sbjID,measure])
    else:
        insertEntryToTable(cursor,"demographicData",["sbjID","measure","value"],[sbjID,measure,value])

def insertResult(cursor,sess,result,value,test_mode=False):
    if inResults(cursor,sess,result):
        updateEntriesFromTable(cursor,"results",["value"],[value],["sess","result"],[sess,result],test_mode=test_mode)
    else:
        insertEntryToTable(cursor,"results",["sess","result","value"],[sess,result,value],test_mode=test_mode)

# Recommended output files for each default step
def getOutputDic(pipe_type):
    diff_dic = {
                "postEddy":"bvals,bvecs,data.nii.gz,nodif_brain_mask.nii.gz,nodif_brain.nii.gz",
                "denoise":"DTI_masked_ds.nii.gz",
                "gibbs":"DTI_masked_ds_gs.nii.gz",
                "dtifit":"dti_V1.nii.gz,dti_MD.nii.gz",
                "camino_wlf":"bvector_$diffdir.scheme,dwi_$diffdir.Bfloat,WM_ROI_$diffdir.nii.gz,estimatesnr_$diffdir.txt,outliermap_$diffdir.nii.gz,dt_RESTORE_$diffdir.Bdouble,dteig_$diffdir.Bdouble,fa_$diffdir.nii.gz"
                }
    if pipe_type=="RS":
        return {
                "tcat":"pb00.$sbj.r01.tcat.nii.gz,tr_counts.txt",
                "despike":"pb01.$sbj.r01.despike.nii.gz",
                "maskEPI":"pb02.$sbj.r01.fslbet_mask.nii.gz",
                "tshift":"polort.txt,outcount.r01.1D,out.min_outlier.txt,pb03.$sbj.r01.tshift.nii.gz,vr_base_min_outlier.nii.gz,outcount_rall.1D,outcount_$sbj_censor.1D",
                "align":"brain_al.nii.gz,brain_al.mat,aseg_al.nii.gz",
                "volreg":"dfile.r01.1D,pb04.$sbj.r01.volreg.nii.gz,mat.r01.vr.aff12.1D,dfile_rall.1D,motion_$sbj_enorm.1D,motion_$sbj_censor.1D,motion_demean.1D,motion_deriv.1D,censor_$sbj_combined_2.1D",
                "segm":"wm.nii.gz,csf.nii.gz,brain_mask.nii.gz,wm_dil.nii.gz",
                "tbregs":"$sbj.r01.wm_timeCourse.1D,$sbj.r01.csf_timeCourse.1D,$sbj.r01.wm_timeCourse_deriv.1D,$sbj.r01.csf_timeCourse_deriv.1D",
                "scale":"mean_r01.nii.gz,pb05.$sbj.r01.scale.nii.gz,bandpass_rall.1D",
                "nuissReg_mot":"stats.REML_cmd,X.xmat.1D,X.jpg,X.nocensor.xmat.1D,X.nocensor.xmat.1D,errts.${sbj}.tproject.nii.gz,out.cormat_warn.txt,TSNR.$sbj.nii.gz"
                }
    elif pipe_type=="ECP":
        diff_dic["pre_mask"] = "[AP_1_mask.nii.gz],[AP_2_mask.nii.gz],[PA_1_mask.nii.gz],[PA_2_mask.nii.gz]",
        diff_dic["pre_dtifit"] = "[AP_1_V1.nii.gz],[AP_2_V1.nii.gz],[PA_1_V1.nii.gz],[PA_2_V1.nii.gz]",
        diff_dic["int_norm"] = "[AP_1.nii.gz],[AP_2.nii.gz],[PA_1.nii.gz],[PA_2.nii.gz]"
        diff_dic["preTopupEddy"] = "eddy/index.txt,eddy/acqparams.txt,eddy/Pos_Neg.nii.gz,eddy/Pos_Neg.bval,eddy/Pos_Neg.bvec,eddy/Pos.nii.gz,eddy/Neg.nii.gz,eddy/Pos.bval,eddy/Neg.bval,eddy/Pos.bvec,eddy/Pos.bval,topup/Pos_Neg_b0.nii.gz,topup/Pos_b0.nii.gz,topup/Neg_b0.nii.gz,topup/acqparams.txt"
        diff_dic["topup"] = "hifib0.nii.gz,nodif_brain_mask.nii.gz,topup_Pos_Neg_b0_movpar.txt,topup_Pos_Neg_b0_fieldcoef.nii.gz,Pos_Neg_b0.topup_log"
        diff_dic["eddyCuda"] = "acqparams.txt,eddy_unwarped_images.nii.gz,eddy_unwarped_images.eddy_rotated_bvecs,index.txt,eddy_unwarped_images.eddy_outlier_report,eddy_unwarped_images.eddy_movement_over_time"
        return diff_dic
    elif pipe_type=="DTI":
        diff_dic["pre_mask"] = "nodif_brain_mask.nii.gz"
        diff_dic["pre_dtifit"] = "dti_V1.nii.gz",
        diff_dic["eddy"] = "acqparams.txt,eddy_unwarped_images.nii.gz,eddy_unwarped_images.eddy_rotated_bvecs,index.txt,eddy_unwarped_images.eddy_outlier_report,eddy_unwarped_images.eddy_movement_over_time"
        return diff_dic
    return {}

# Recommended input files for each default step
def getInputDic(pipe_type):
    diff_dic = {
                "postEddy":"eddyCuda/eddy_unwarped_images.nii.gz,bvals,eddyCuda/eddy_unwarped_images.eddy_rotated_bvecs,nodif_brain_mask.nii.gz",
                "denoise":"DTI_masked.nii.gz,nodif_brain_mask.nii.gz,nodif_brain_4dmask.nii.gz",
                "gibbs":"DTI_masked_ds.nii.gz",
                "dtifit":"DTI.nii.gz,nodif_brain_mask.nii.gz,bvals,bvecs",
                "camino_wlf":"data/bvecs,data/bvals,data/data.nii.gz,camino/WM_ROI.nii.gz,data/nodif_brain_mask.nii.gz",
                }
    if pipe_type=="RS":
        return {
                "tcat":"$sbj_r01.nii.gz",
                "despike":"pb00.$sbj.r01.tcat.nii.gz",
                "maskEPI":"pb01.$sbj.r01.despike.nii.gz",
                "tshift":"pb02.$sbj.r01.masked.nii.gz,pb02.$sbj.r01.fslbet_mask.nii.gz",
                "align":"brain.nii.gz,vr_base_min_outlier.nii.gz,aseg.nii.gz",
                "volreg":"vr_base_min_outlier.nii.gz,pb03.$sbj.r01.tshift.nii.gz,pb02.$sbj.r01.fslbet_mask.nii.gz",
                "segm":"aseg_al.nii.gz",
                "tbregs":"wm.nii.gz,csf.nii.gz",
                "scale":"pb04.$sbj.r01.volreg.nii.gz,pb02.$sbj.r01.fslbet_mask.nii.gz",
                "nuissReg_mot":"polort.txt,pb05.$sbj.r01.scale.nii.gz,bandpass_rall.1D,motion_demean.1D,motion_deriv.1D,censor_$sbj_combined_2.1D,pb02.$sbj.r01.fslbet_mask.nii.gz"
                }
    elif pipe_type=="ECP":
        diff_dic["pre_mask"] = "[AP_1.nii.gz],[AP_2.nii.gz],[PA_1.nii.gz],[PA_2.nii.gz]",
        diff_dic["pre_dtifit"] = "[AP_1_masked.nii.gz],[AP_2_masked.nii.gz],[PA_1_masked.nii.gz],[PA_2_masked.nii.gz],[AP_1_mask.nii.gz],[AP_2_mask.nii.gz],[PA_1_mask.nii.gz],[PA_2_mask.nii.gz]",
        diff_dic["int_norm"] = "DTI_$diffdir_masked_ds_gs.nii.gz,bvals_$diffdir,bvecs_$diffdir"
        diff_dic["preTopupEddy"] = "[AP_1_masked.nii.gz],[AP_2_masked.nii.gz],[PA_1_masked.nii.gz],[PA_2_masked.nii.gz],[AP_1.bval],[AP_2.bval],[PA_1.bval],[PA_2.bval],[AP_1.bvec],[AP_2.bvec],[PA_1.bvec],[PA_2.bvec]"
        diff_dic["topup"] = "Pos_Neg_b0.nii.gz,acqparams.txt,Pos_b0.nii.gz,Neg_b0.nii.gz"
        diff_dic["eddyCuda"] = "eddy/Pos_Neg.nii.gz,eddy/index.txt,eddy/acqparams.txt,eddy/Pos_Neg.bvec,eddy/Pos_Neg.bval,eddy/slspecFile.txt,topup/nodif_brain_mask.nii.gz,topup/topup_Pos_Neg_b0_fieldcoef.nii.gz,topup/topup_Pos_Neg_b0_movpar.txt"
        return diff_dic
    elif pipe_type=="DTI":
        diff_dic["pre_mask"] = "DTI.nii.gz",
        diff_dic["pre_dtifit"] = "DTI_masked.nii.gz,nodif_brain_mask.nii.gz",
        diff_dic["eddy"] = "DTI_masked_ds_gs.nii.gz,nodif_brain_mask.nii.gz,bvecs,bvals"
        return diff_dic
    return {}

# Recommended description for each default step
def getDescDic(pipe_type):
    if pipe_type=="RS":
        return {
                "tcat":"Remove the first 3 TRs, save the number of TRs per run.",
                "despike":"Remove spikes from the dataset, remove negative values from output.",
                "maskEPI":"Mask despike data using BET",
                "tshift":"Align slices to the beginning of the TR (slice-timing correction), compute polort and outlier fraction per volume, censor outlier TRs per run (>5% voxels), extract base img for registration.",
                "align":"Align anatomical to RS using Flirt, apply registration to ROI and brain mask.",
                "volreg":"Register each volume to the base image, do motion correction, obtain motion parameters, combine multiple censor files.",
                "segm":"Perform tissue based segmentation, dilatate WM mask to make sure we are not including WM in the correltions, rm WM from the ROI.",
                "tbregs":"Extract avg time series from WM and CSF (must be done BEFORE blurring).",
                "scale":"Compute the mean of each voxel time series, scale each voxel time series to have a mean of 100.",
                "nuissReg_mot":"Perform nuissance regression and bandpass filtering (run 3dDeconvolve to regress out signals of no interest), generate clean data and tSNR dataset, extract time series.",
                "func_graph":"Generate connectivity matrices and graphs"
                }
    return {}

# Recommended steps for each standard pipeline
def getDefaultSteps(pipe_type):
    if pipe_type=="ECP":
        return ["mask","denoise","pre_ecc","dtifitPreTopUp","int_norm","preTopupEddy","topup","eddyCuda","postEddy","camino_wlf","ants_wlf","transform_sklt","camino_snr"]
    elif pipe_type=="DTI":
        return ["mask","denoise","dtifit","eddyCuda","postEddy","camino_wlf","ants_wlf","camino_snr"]
    elif pipe_type=="RS":
        return ["tcat","despike","maskEPI","tshift","bet_anat","align","volreg","scale","segm","tbregs","nuissReg_mot","registration","roi_segment","extract_ts","con_mat","con_map"]
    else:
        return []

def getStringVal(cursor,variable,numerical_val):
    return getStringFromTable(cursor,"string_val","string2num",["variable","numerical_val"],[variable,numerical_val])
    
def getNotesTBSS(cursor,tbss_name):
    return getStringFromTable(cursor,"notes","tbss",["name"],[tbss_name])

def getNotesProc(cursor,sess,step):
    return getStringFromTable(cursor,"notes","procs",["sess","step"],[sess,step])

def getNotesSess(cursor,sess):
    return getStringFromTable(cursor,"notes","sessions",["sess"],[sess])

def getNotesSbj(cursor,sbjID):
    return getStringFromTable(cursor,"notes","subjects",["sbjID"],[sbjID])

def getDatabases(cursor):
    exclude = ["information_schema","mysql","performance_schema","sys","generalInfo"]
    cursor.execute("show databases")
    return [item[0] for item in cursor.fetchall() if item[0] not in exclude]

def getDB(cursor):
    cursor.execute("select database()");
    return cursor.fetchone()[0]

def getColumnsDB(cursor,database=""):
    databases = getDatabases(cursor) if database=="" else [database]
    tables = {}
    for db in databases:
        for table in getTablesDB(cursor,db):
            tables[db+"."+table] = getColsFromTable(getTableInfo(cursor,db+"."+table))
    return tables

def colvalsTable(cursor,table,column,toPrint=False):
    values = np.unique(getColumnFromTable(cursor,column,table,order_by=[column]))
        
    if toPrint:
        print("\n---------------------------")
        print("# Current "+table+":")
        print("\n".join(values))
        print("---------------------------")
    
    return values

def getPipelines(cursor,project):
    return getColumnFromTable(cursor,"pipename","pipelines",["project"],[project])

def getSbjID(cursor,sess):
    return getStringFromTable(cursor,"sbjID","sessions",["sess"],[sess])

def getSessList(cursor,pipeline):
    return getColumnFromTable(cursor,"sess","sessions",[pipeline],["1"])

def getIncludedSessions(cursor,pipeline):
    return [sess for sess in getSessList(cursor,pipeline) if not isExcluded(cursor,sess)]

def getSessions(cursor,sbjID):
    return getColumnFromTable(cursor,"sess","sessions",["sbjID"],[sbjID])

def getProcsRunningLocal(cursor):
    return getNumberEntries(cursor,"procs",["status","jobID"],["running",""])

def isExcluded(cursor,sess):
    return getNumberEntries(cursor,"excluded",["sess"],[sess])!=0

def inProcs(cursor,sess,step):
    return getNumberEntries(cursor,"procs",["sess","step"],[sess,step])!=0

def inResults(cursor,sess,result):
    return getNumberEntries(cursor,"results",["sess","result"],[sess,result])!=0
    
def inDemographicData(cursor,sbjID,measure):
    return getNumberEntries(cursor,"demographicData",["sbjID","measure"],[sbjID,measure])!=0

def inBehavioralData(cursor,sess,measure):
    return getNumberEntries(cursor,"behavioralData",["sess","measure"],[sess,measure])!=0

def procFinishedRunning(step,outfile,location,jobID):
    if location=="hpc":
        n_run = int(systemOut("qstat "+jobID+" | wc -l"))
        return True if n_run==0 else False
    else:
        ans = False
        try:
            fo = open(outfile,'r')
            for lastline in fo:
                if lastline.replace('\n','')=="DONE "+step:
                    ans = True
                    break
        except ValueError as err:
            array = outfile.split('/')
            error = ("an exception occured while reading "+array[-1]+": "+str(err)).replace('\n','').replace('\'','')
            raise ValueError(error)
        finally:
            fo.close()
        return ans

def mvFiles(cursor,step):
    # If there's more than one step contained in step, they should all have the same move_files
    return getValueFromTable(cursor,"move_files","pipesteps",["step"],[step.split(",")[0]])

def getScanDateFromDICOM(cursor,pipe,location,sess):
    sess_dir = getSessDir(cursor,pipe,location,sess).replace("$sess",sess).replace("$sbj",getSbjID(cursor,sess))+"/"
    array = glob.glob(sess_dir+"*.zip")
    if len(array)==0:
        return ""
    
    # Unzip the dicom folder
    zipped = array[0]
    unzipped = zipped.replace(".zip","")
    if not os.path.isdir(unzipped):
        with zipfile.ZipFile(zipped, 'r') as zip_ref:
            zip_ref.extractall(unzipped)
            
    # Get the scan date from the first dicom
    dcm1 = glob.glob(unzipped+"/"+os.path.basename(unzipped)+"/"+sess+"/?/DICOM/*.dcm")[0]
    line = systemOut("dicom_hdr "+dcm1+" | grep 'ID Acquisition Date'")
    dcmdate = str(line.split("//")[-1])
    dcmdate = dcmdate[0:4]+"-"+dcmdate[4:6]+"-"+dcmdate[6:8]
    return dcmdate

def printDescription(cursor,step):
    print(f"\n*** Step: {step} ***")
    desc = getValueFromTable(cursor,"description","pipesteps",["step"],[step])
    if desc!=None:
        print(desc+'\n')

def systemOut(cmd,noSpace=True,noLines=True):
    out = os.popen(cmd).read()
    if noLines:
        out = out.replace("\n","")
    return out.replace(" ","").replace("\t","") if noSpace else out

def getLogs(logsDir,sess,step,jobID,cluster):
    if not logsDir.endswith("/"):
        logsDir+="/"
    
    # Job ran in the cluster
    if str(cluster) not in empty_values:
        jobID = str(jobID).replace(".0","")
        if cluster=="hpc":
            ofile = logsDir+step+'_'+sess+".sh.o"+jobID
            efile = logsDir+step+'_'+sess+".sh.e"+jobID
        else:
            ofile = logsDir+step+".o"+jobID
            efile = logsDir+step+".e"+jobID
            
    # Job ran locally
    else:
        ofile = logsDir+step+'.'+sess
        efile = ofile
        
    return [ofile,efile]

def rollBackSess(db_name,sess,pipe,location,username,password,hostdir,test_mode=False,first=""):
    cnx = connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    
    ## Select what step to roll to for sess in pipe
    steps_pipe = "("+",".join(["'"+step+"'" for step in getStepsPipeline(cursor,pipe)])+")"
    sess_steps_pipe = getEntriesFromTable(cursor,["step","jobID","cluster"],"procs",["sess","step"],[sess,steps_pipe],index="step")
    if len(sess_steps_pipe)==0:
        printError("No steps from "+pipe+" are running on "+sess)
        return False
    
    if first=="" or first not in sess_steps_pipe.index.tolist():
        first = input("To what step would you like to roll back "+sess+"? ("+", ".join(sess_steps_pipe.index.tolist())+"): ")
        if first not in sess_steps_pipe.index.tolist():
            printError("Not a valid step")
            return False
    
    ## Set that step as ready
    if updateEntryFromTable(cursor,"procs",["status","jobID","exectimesec","qsub_time","cluster"],["ready","","","",""],["sess","step"],[sess,first],test_mode=test_mode)==-1:
        cnx.rollback()
        return False
    if input("\nDelete output files for "+first+"? [Y]: ")=='Y':
        outs = getOutputFiles(cursor,sess,first,location,optional=True)
        if len(outs)>0:
            print("rm -f "+" ".join(outs)+"\n")
            if not test_mode: 
                os.system("rm -f "+" ".join(outs)+"\n")
    
    # Remove logs
    jobID = str(sess_steps_pipe.loc[first,"jobID"])
    for val in empty_values:
        jobID = jobID.replace(val,"")
    logsDir = getScriptsDir(cursor,pipe,location)+"/logs/"
    logs = list(set(getLogs(logsDir,sess,first,jobID,str(sess_steps_pipe.loc[first,"cluster"]))))
    print("\nrm -f "+" ".join(logs)+"\n")
    if not test_mode:
        os.system("rm -f "+" ".join(logs))
    
    # From all the steps running on this sess, select the ones that should be deleted
    sess_steps_all = getEntriesFromTable(cursor,["step","jobID","cluster"],"procs",["sess","step"],[sess,f"!{first}"],index="step")
    delete = sess_steps_all.loc[[step for step in getStepsPipeline(cursor,step_start=first) if step in sess_steps_all.index.tolist()[1:]]]
    
    # Delete steps
    del_steps = "("+",".join(["'"+step+"'" for step in delete.index.tolist()])+")"
    print("Deleting procs for "+sess+" "+del_steps+"...")
    deleteEntriesFromTable(cursor,"procs",["sess","step"],[sess,del_steps],test_mode=test_mode)
    print("done")
    
    # Delete logs and output files if desired
    if input("\nDelete output files "+del_steps+"? [Y]: ")=='Y':
        for step in delete.index.tolist():
            print("\nDeleting files from "+step+"...")
            
            jobID = delete.loc[step,"jobID"]
            cluster = delete.loc[step,"cluster"]
            logs = list(set(getLogs(logsDir,sess,step,jobID,cluster)))
            outs = getOutputFiles(cursor,sess,step,location,optional=True)
            print("rm -f "+" ".join(logs+outs)+"\n")
            if not test_mode:
                os.system("rm -f "+" ".join(logs+outs))
                
            print("done")
    
    updateDB(cnx,cursor)
    cursor.close()
    cnx.close()
    
    return True

def getSummaryTable(cursor,pipeline="",subjects_cols=[],sessions_cols=[],demographic_vars=[],behavioral_vars=[],results_meas=[],include={},exclude={},sort_cols=[]):
    if pipeline!="":
        df = pd.DataFrame({"sess":getSessList(cursor,pipeline)})
    else:
        df = pd.DataFrame({"sess":getColumnFromTable(cursor,"sess","sessions")})
    df.set_index('sess',inplace=True)
    
    df["sbjID"] = [getSbjID(cursor,sess) for sess in df.index.values.tolist()]
    
    for col in subjects_cols:
        df[col] = [getValueFromTable(cursor,col,"subjects",["sbjID"],[sbjID]) for sbjID in df["sbjID"].values.tolist()]
        
    for col in sessions_cols:
        df[col] = [getValueFromTable(cursor,col,"sessions",["sess"],[sess]) for sess in df.index.values.tolist()]
        
    for var in demographic_vars:
        df[var] = [getValueFromTable(cursor,"value","demographicData",["sbjID","measure"],[sbjID,var]) for sbjID in df["sbjID"].values.tolist()]
        
    for var in behavioral_vars:
        df[var] = [getValueFromTable(cursor,"value","behavioralData",["sess","measure"],[sess,var]) for sess in df.index.values.tolist()]
        
    for result in results_meas:
        df[result] = [getValueFromTable(cursor,"value","results",["sess","result"],[sess,result]) for sess in df.index.values.tolist()]
        
    for col,values in include.items():
        df = df[df[col].isin(values)]
        
    for col,values in exclude.items():
        df = df[not df[col].isin(values)]
    
    if len(sort_cols)>0:
        df.sort_values(by=sort_cols,inplace=True)
        
    return df

def printError(msg):
    print("\nXXXX ERROR: "+str(msg)+" XXXX")
    
def printWarning(msg):
    print("\nXX "+str(msg)+" XX")
    
def getUniqueValsDFcol(df,col_name):
    return np.unique(df.loc[:,col_name].values.tolist())

#############################################################################
def deleteEntryFromTable(cursor,table,columns_in,values_in,test_mode=False):
    n = deleteEntriesFromTable(cursor,table,columns_in,values_in,test_mode=test_mode)
    if n>1:
        printError(str(n)+" records where updated. rolling back.")
        return -1
    return n

def deleteEntriesFromTable(cursor,table,columns_in,values_in,test_mode=False):
    if len(columns_in)!=len(values_in) or len(columns_in)==0:
        printError("size of columns_edit and values_edit dont match or are empty")
        return -1
    
    # Make sure the table exists
    table_info = getTableInfo(cursor,table)
    if len(table_info)==0:
        printError(table+" table doesnt exist")
        return -1
    
    # Make sure all values in columns_in are a column in the table
    if not set(columns_in)<=set(getColsFromTable(table_info)):
        printError("one or more values in ["+",".join(columns_in)+"] are not valid fields for table "+table)
        return -1
    
    cmd = "delete from "+table+whereStatement(columns_in,values_in)
    
    # Print if in test mode
    if test_mode:
        print(cmd)
        return 0
    
    # Execute if not in test mode
    try:
        cursor.execute(cmd)
        return cursor.rowcount
    except Error as err:
        printError(err)
        return -1

def deleteEntriesFromDF(cnx,cursor,table,df,test_mode=False):
    # Make sure the table exists
    table_info = getTableInfo(cursor,table)
    if len(table_info)==0:
        printError(table+" table doesnt exist")
        return -1
    
    # Set the index of the DF as the primary key(s) of the table in the DB
    # A primary key is needed to make sure we delete only the correct entries
    prim_keys = getPrimaryKeys(table_info)
    df_cols = df.columns.values.tolist()
    if len(prim_keys)==0:
        printError(table+" is missing a primary key")
        return -1
    if not set(prim_keys)<=set(df.columns.values.tolist()):
        printError(table+" has primary keys ["+",".join(prim_keys)+"] but one or more of them is missing in DF columns: ["+",".join(df_cols)+"]")
        return -1
    df.set_index(prim_keys,inplace=True)
    
    # Make sure all the columns in the DF are present in the table
    db_cols = getColsFromTable(table_info)
    del_cols = [col for col in df_cols if col not in db_cols]
    if len(del_cols)>0:
        df.drop(del_cols,axis=1,inplace=True)
        printWarning("The following columns where removed from the DF because they dont exist in "+table+": "+",".join(del_cols))
    df_cols = df.columns.values.tolist()
    if len(df_cols)==0:
        printError("none of the columns of the df belong to table "+table)
        return -1
    
    n_deleted = 0
    for index,line in df.iterrows():
        values = np.array(index).tolist() if len(prim_keys)>1 else [index]
        n = deleteEntryFromTable(cursor,table,prim_keys,values,test_mode=test_mode)
        if n>=0:
            n_deleted+=n
        else:
            printWarning(f"entry with values {values} for keys {prim_keys} was not deleted")
    return n_deleted
        
def updateEntriesFromDF(cnx,cursor,table,df,update_all=False,test_mode=False):
    # Make sure the table exists
    table_info = getTableInfo(cursor,table)
    if len(table_info)==0:
        printError(table+" table doesnt exist")
        return [-1,-1]
    
    # Set the index of the DF as the primary key(s) of the table in the DB
    # A primary key is needed to make sure we edit only the correct entries
    prim_keys = getPrimaryKeys(table_info)
    df_cols = df.columns.values.tolist()
    if len(prim_keys)==0:
        printError(table+" is missing a primary key")
        return [-1,-1]
    if not set(prim_keys)<=set(df.columns.values.tolist()):
        printError(table+" has primary keys ["+",".join(prim_keys)+"] but one or more of them is missing in DF columns: ["+",".join(df_cols)+"]")
        return [-1,-1]
    df.set_index(prim_keys,inplace=True)
    
    # Make sure all the columns in the DF are present in the table
    db_cols = getColsFromTable(table_info)
    del_cols = [col for col in df_cols if col not in db_cols]
    if len(del_cols)>0:
        df.drop(del_cols,axis=1,inplace=True)
        printWarning("The following columns where removed from the DF because they dont exist in "+table+": "+",".join(del_cols))
    df_cols = df.columns.values.tolist()
    if len(df_cols)==0:
        printError("none of the columns of the df belong to table "+table)
        return [-1,-1]

    n_added = 0
    n_updated = 0
    for index,line in df.iterrows():
        # Get database values for each column
        idx_array = np.array(index).tolist() if len(prim_keys)>1 else [index]
        db_values = getEntryFromTable(cursor,df_cols,table,prim_keys,idx_array)
        
        keys = []
        values = []
        for key,value in line.items():
            keys+=[key]
            values+=[value]
        
        # The entry does not exist in the DB
        if len([val for val in db_values if val!=None])==0:
            if insertEntryToTable(cursor,table,prim_keys+keys,idx_array+values,test_mode)==-1:
                return [-1,-1]
            updateDB(cnx,cursor)
            n_added+=1
            continue
        
        # Update each line
        if update_all:
            if updateEntryFromTable(cursor,table,keys,values,prim_keys,idx_array,test_mode)==-1:
                cnx.rollback()
                return [-1,-1]
        
        # Check each column and update if desired
        else:
            for i in range(len(df_cols)):
                col = df_cols[i]
                db_value = str(db_values[i])
                if db_value=="None":
                    db_value = ""
                df_value = str(line[col]) if "date" not in str(table_info.loc[col,"Type"]) else str(line[col]).split(" ")[0]
                if df_value in empty_values:
                    df_value = ""
                
                if db_value==df_value:
                    continue
                
                if df_value=="" and table_info.loc[col,"Null"]=="NO":
                    printWarning("Cant update "+col+" for "+",".join(idx_array)+": value cant be None")
                    continue
                
                printWarning("inconsistency between DF and DB")
                print("column: "+col)
                print(",".join(prim_keys)+": "+",".join(idx_array))
                print("DF value: "+df_value)
                print("DB value: "+db_value)
                        
                if input("Update DB value? [Y]: ")=='Y' and updateEntryFromTable(cursor,table,[col],[df_value],prim_keys,idx_array,test_mode)==-1:
                    cnx.rollback()
                    return [-1,-1]
                
        updateDB(cnx,cursor)
        n_updated+=1
        
    return n_added,n_updated

def recoverTable(cnx,cursor,table,cols_update="",location="petrov",col_filter="",val_filter="",add_new_lines=True,delete_missing_lines=False,test_mode=False):
    cursor.execute("select database()")
    database = cursor.fetchone()[0]
    db_dir = getValueFromTable(cursor,"db_dir","generalInfo.locations",["location"],[location])
    recov_file = f"{db_dir}/recovery/{database}.{table}.txt"
    
    return recoverTable2(cnx,cursor,recov_file,cols_update,col_filter,val_filter,add_new_lines,delete_missing_lines,test_mode)

# add_new_lines: 
### add to the DB lines that appear in the recovery file but not currently in the DB (using the primary keys of the table)
### it would be True if I want to add any lines currently missing in the DB but present in the recovery file
### it would be False if I onbly want to update the current entries in the DB but not add anything new
# delete_missing_lines: 
### delete from the DB lines that do not appear in the recovery file but are currently in the DB
# If cols_update is empty, update all columns in the DF
def recoverTable2(cnx,cursor,recov_file,cols_update="",col_filter="",val_filter="",add_new_lines=True,delete_missing_lines=False,test_mode=False):
    if not os.path.isfile(recov_file):
        return [-1,-1,-1]
    
    # Create the DF to update the DB
    df_file = fileslib.fileToDataFrame(recov_file,True,"\t")
    table = os.path.basename(recov_file).split(".")[1]
    prim_keys = getPrimaryKeys(getTableInfo(cursor,table))
    cols_df = prim_keys+[col for col in cols_update if col not in prim_keys]
    if cols_update!="":
        cols_df = prim_keys+[col for col in cols_update if col not in prim_keys]
    else:
        cols_df = prim_keys+[col for col in df_file.columns if col not in prim_keys]
    df_file = df_file[cols_df]
    
    # Keep only the desired lines based on filter if any
    if col_filter!="" and val_filter!="":
        df_file = df_file[df_file[col_filter]==val_filter]
        
    # Get the DF with the current status of the table in the DB
    df_table = getEntriesFromTable(cursor,cols_df,table)
    
    # Create a new column with all the index columns merged into one
    cols_drop = []
    if len(prim_keys)>1:
        new_index = "_".join(prim_keys)
        df_file[new_index] = df_file[prim_keys].agg("_".join,axis=1)
        df_table[new_index] = df_table[prim_keys].agg("_".join,axis=1)
        cols_drop+=[new_index]
    else:
        new_index = prim_keys[0]
        
    # Remove from the DB any lines that are not in the recovery file
    if delete_missing_lines:
        df_table["delete"] = list(map(lambda index: index not in df_file[new_index].values.tolist(),df_table[new_index].values.tolist()))
        df_delete = df_table[df_table["delete"]]
        df_delete = df_delete.drop(columns=["delete"], axis=1)
        n_deleted = deleteEntriesFromDF(cnx,cursor,table,df_delete,test_mode=test_mode)
    else:
        n_deleted = 0
    
    if not add_new_lines:
        # Only include those lines of df_file for which the index is also present in df_table
        df_file["include"] = list(map(lambda index: index in df_table[new_index].values.tolist(),df_file[new_index].values.tolist()))
    else:
        # Update all entries in df_file
        df_file["include"] = [True for i in range(len(df_file))]
    cols_drop+=["include"]
    
    df_file = df_file[df_file["include"]]
    df_file.drop(columns=cols_drop, axis=1, inplace=True)
            
    [n_added,n_updated] = updateEntriesFromDF(cnx,cursor,table,df_file,True,test_mode)
    
    print(f"{n_added} entries added, {n_updated} entries updated, {n_deleted} entries deleted")
    return [n_added,n_updated,n_deleted]

# Takes back the table to the state it was when the recovery file was created
def recoverTables(username,password,hostdir,db_name,location="petrov"):
    print("\n### Recovering database from backup ###")
          
    cnx =  connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    
    # Wont be able to add any tables/columns if this script doesn't exist
    db_dir = getValueFromTable(cursor,"db_dir","generalInfo.locations",["location"],[location])
    recov_dir = f"{db_dir}/recovery/"
    recreate_script = recov_dir+"recreateTables.sql"
    if not os.path.isfile(recreate_script):
        return
        
    print("Select the DB to recover:")
    print("\n".join(getDatabases(cursor)))
    print("[Enter] to create a new database from backup")
    db_recover = input(">> ")
    while db_recover=="":
        db_recover = input("Name of the new DB: ")
    
    # Create new DB if it doesnt exist
    cursor.execute("create database if not exists "+db_recover)
    updateDB(cnx,cursor)
    
    # Get the list of tables that already exist in the DB
    tables_db = getTablesDB(cursor,db_recover)
        
    # Add any missing tables
    all_tables = []
    for recov_file in glob.glob(recov_dir+db_recover+".*.txt"):
        new_table = os.path.basename(recov_file).replace(db_recover+".","").replace(".txt","")
        
        if new_table in tables_db:
            all_tables+=[new_table]
        elif input("Create "+new_table+"? [Y]")=='Y':
            cmd = systemOut("grep "+new_table+" "+recreate_script+" | head -n1",False).replace(";","")
            if cmd!="":
                print(cmd)
                cursor.execute(cmd)
                all_tables+=[new_table]
                
    # Wont be able to create any tables/columns if there's no recovery files for this new DB
    if len(all_tables)==0:
        return
    
    print("\nSelect tables to recover (divided by comma):")
    print("\n".join(all_tables))
    print("[Enter] to include all tables")
    tables_recover = input(">> ")
    if tables_recover=="":
        tables_recover = all_tables
    else:
        tables_recover = tables_recover.split(",")
        
    for table_recover in tables_recover:
        print(f'\nSelect columns from {db_recover}.{table_recover} to recover (divided by comma):')
        all_cols = getColsFromTable(getTableInfo(cursor,table_recover))
        print("\n".join(all_cols))
        print("[Enter] to update all columns")
        cols_recover = input(">> ")
        cols_recover = all_cols if cols_recover=="" else cols_recover.split(",")
        
        [n_added,n_edited] = recoverTable(cnx,cursor,table_recover,cols_recover,location)
        print(f'\n{n_added} entries were added')
        print(f'{n_edited} entries were updated')
        
    cursor.close()
    cnx.close()

# This function doesn't add the sessions of each subject! Call the corresponding function after this.
def addSubjects(username,password,hostdir,db_name,project):
    cnx = connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    sites = colvalsTable(cursor,"sites","site")
    subjects = colvalsTable(cursor,"subjects","sbjID")
    
    for site in sites:
        print("\nsbj_list can have either one column with just sbjIDs or two columns sbjID,grp")
        sbj_list = input(f"Path to sbj_list for {site} site ([Enter] if not adding subjects for this site): ")
        if not os.path.isfile(sbj_list):
            continue
        
        fin = open(sbj_list,'r')
        for line in fin:
            info = line.replace("\n","").split(",")
            sbj = info[0]
            
            if sbj not in subjects:
                if len(info)==1:
                    subjects+=[sbj]
                    insertEntryToTable(cursor,"subjects",["sbjID","site",project],[sbj,site,1])
                elif len(info)==2:
                    subjects+=[sbj]
                    insertEntryToTable(cursor,"subjects",["sbjID","site","grp",project],[sbj,site,info[1],1])
                else:
                    printError(f"Line not added (wrong format): {line}")
            else:
                updateEntriesFromTable(cursor,"subjects",["site",project],[site,"1"],["sbjID"],[sbj])
                if len(info)==2:
                    updateEntriesFromTable(cursor,"subjects",["grp"],[info[1]],["sbjID"],[sbj])
                    
            updateDB(cnx,cursor)
        fin.close()
    
    cursor.close() 
    cnx.close()
    return subjects

# The corresponding sbjs must have been added before this!
def addSessions(username,password,hostdir,db_name,sess_list):
    cnx = connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    all_sessions = colvalsTable(cursor,"sessions","sess")
    subjects = colvalsTable(cursor,"subjects","sbjID")
    
    print("Select an option:\n"
          "1. sess==sessID\n"
          "2. sess==sbjID_sessID (default)")
    opt = input(">> ")
        
    fin = open(sess_list,'r')
    for line in fin:
        info = line.replace("\n","").split(":")
        if len(info)!=2:
            printError(f"Line not added (wrong format): {line}")
            continue
            
        sbj = info[0]
        if sbj not in subjects:
            printError(f"Subject {sbj} not found. Line not added: {line}")
            continue
            
        sessions = info[1].split(",")
        for sessID in sessions:
            sess = sessID if opt=="1" else sbj+"_"+sessID
            if sess in all_sessions:
                printWarning(f"Session already exists: {sess}")
            else:
                insertEntryToTable(cursor,"sessions",["sess","sbjID","sessID"],[sess,sbj,sessID])
                all_sessions+=[sess]
                    
        updateDB(cnx,cursor)
    
    fin.close()
    cursor.close()
    cnx.close()
    
    return all_sessions

def addPipes(cnx,cursor,database,project,new_pipes,pipe_type="ECP"):
    pipelines = getPipelines(cursor,project)
    sess_cols = getColsFromTable(getTableInfo(cursor,"sessions"))
    
    for new_pipe in [pipe for pipe in new_pipes if pipe not in pipelines]:
        print(f"\nAdding {new_pipe} pipelines...")
        steps = getDefaultSteps(pipe_type)
        if len(steps)>0:
            print(f"Suggested 1st step: {steps[0]}")
        first_step = input(f"1st step in {new_pipe}: ")
        cmd1 = "pipename,project,first_step"
        cmd2 = "'"+new_pipe+"','"+project+"','"+first_step+"'"
                
        data_dir = input("data_dir (remove project_dir i.e.'Data_Release'): ")
        if data_dir!="":
            cmd1+=",data_dir"
            cmd2+=",'"+data_dir+"'"
            
        sess_dir = input("sess_dir (remove data_dir i.e. '$sbj/$sess')")
        if sess_dir!="":
            cmd1+=",sess_dir"
            cmd2+=",'"+sess_dir+"'"
            
        scripts_dir = input("scripts_dir (remove project_dir i.e. 'Scripts'): ")
        if scripts_dir!="":
            cmd1+=",scripts_dir"
            cmd2+=",'"+scripts_dir+"'"
        
        insertEntryToTable(cursor,"pipelines",[cmd1],[cmd2])
        pipelines+=[new_pipe]
        
        if new_pipe not in sess_cols:
            cursor.execute(f"alter table sessions add column {new_pipe} boolean")
            
        if input(f"Does ALL sessions belong to {new_pipe} pipeline? [Y]: ")=='Y':
            updateEntriesFromTable(cursor,"sessions",[new_pipe],["1"])
        elif input(f"Add sessions to {new_pipe} from file? [Y]: ")=='Y':
            sess_file = input(f"List of sess to add to {new_pipe}: ")
            sess_array = "("+",".join(["'"+sess+"'" for sess in fileslib.file2array(sess_file)])+")"
            n = updateEntriesFromTable(cursor,"sessions",[new_pipe],["1"],["sess"],[sess_array])
            if n==-1:
                printError("No sessions added to {new_pipe}")
            else:
                print(f"{n} sessions added to {new_pipe}")
            
        updateDB(cnx,cursor)
    
    return pipelines

def addImages(cnx,cursor,db,pipe_type):
    if pipe_type=="ECP":
        print("\nSuggested imageTypes: AP_1,AP_2,PA_1,PA_2")
    elif pipe_type=="DTI":
        print("\nSuggested imageTypes: DTI")
    elif pipe_type=="RS":
        print("\nSuggested imageTypes: RS")
        
    # Add one column per img to the sessions table to mark if each session has the corresponding img
    imgTypes = input("Add image types divided by comma: ").split(",")
    cols = getColsFromTable(getTableInfo(cursor,"sessions"))
    for img in imgTypes:
        if img not in cols:
            cursor.execute(f"alter table sessions add column {img} boolean")
            updateDB(cnx,cursor)
            
    return imgTypes

# Adds one or more steps to the pipeline when the pipeline IS BEING CREATED
# If the pipeline has already been created, better use function addStep()
# Serves also to modify column values for any step after the pipeline has been created
def addSteps(cnx,cursor,location,project,pipe_type):
    dic_desc = getDescDic(pipe_type)
    dic_input = getInputDic(pipe_type)
    dic_output = getOutputDic(pipe_type)
    
    for pipe, line in getEntriesFromTable(cursor,["pipename","scripts_dir","data_dir","first_step"],"pipelines",["project"],[project],index="pipename").iterrows():
        if input(f"\nDo you which to add/modify steps in {pipe} pipeline? [N]: ")!='N':
            proj_dir = getProjDir(cursor,pipe,location)
            scripts_dir = f"{proj_dir}/{line['scripts_dir']}/"
            data_dir = f"{proj_dir}/{line['data_dir']}/" 
                    
            step = line["first_step"] if line["first_step"]!=None else ""
            extra_steps = []
            while step!="":
                print("\n---------------------------\n"
                      "# "+step)
                # Add the new step to the DB
                if getNumberEntries(cursor,"pipesteps",["step"],[step])==0:
                    insertEntryToTable(cursor,"pipesteps",["step","pipeline"],[step,pipe])
                    updateDB(cnx,cursor)
                
                # Change move_files if desired
                if input(f'Change move_files to False for {step}? [Y]: ')=='Y':
                    updateEntriesFromTable(cursor,"pipesteps",["move_files"],["False"],["step"],[step])
                    updateDB(cnx,cursor)
                
                [description,curr_script,curr_local,curr_indir,curr_outdir,ns] = getEntryFromTable(cursor,["description","script","local_script","input_folder","output_folder","next_step"],"pipesteps",["step"],[step])
                
                # Add the description
                if description!=None:
                    print('\nCurrent description: {description}')
                if step in dic_desc.keys():
                    print(f'\nSuggested description: {dic_desc[step]}')
                desc = input(f'\nNew description for {step}: ')
                if desc!="":
                    updateEntriesFromTable(cursor,"pipesteps",["description"],[desc],["step"],[step])
                    updateDB(cnx,cursor)
                
                # Add rcc script
                if curr_script!=None:
                    print(f'\nCurrent **SUBMISSION** script: {curr_script}')
                script = input(f'\nNew **SUBMISSION** script for {step}: {scripts_dir}')
                if script!="":
                    print("Write the following lines to declare sbj and sess:\n"
                          "sbj=sbjID\n"
                          "sess=sess\n"
                          "proj=proj\n"
                          "step=step\n"
                          "sessID=sessID\n"
                          f'data_dir={data_dir}\n'
                          f'vi {scripts_dir}{script}')
                    input("[Enter]")
                    if os.path.isfile(scripts_dir+script):
                        updateEntriesFromTable(cursor,"pipesteps",["script"],[script],["step"],[step])
                        updateDB(cnx,cursor)
                
                # Add local script
                if curr_local!=None:
                    print(f'\nCurrent **LOCAL** script: {curr_local.split(" ")[0]}')
                local_script = input(f'\nNew **LOCAL** script for {step}: {scripts_dir}')
                if local_script!="":
                    print(f'vi {scripts_dir}{local_script}\n'
                          "** remember to update permissions **")
                    input("[Enter]")
                    if os.path.isfile(scripts_dir+local_script):
                        params = input(f'local_script parameters separated by SPACE (use $sbj or $sess if necessary) for {step}: ')
                        local_cmd = local_script if params=="" else local_script+" "+params
                        updateEntriesFromTable(cursor,"pipesteps",["local_script"],[local_cmd],["step"],[step])
                        updateDB(cnx,cursor)
                 
                # Add input folder and files
                if curr_indir!=None:
                    print(f'\nCurrent input folder: {curr_indir}')
                print("Use $sbj or $sess if needed")
                print("Use $sess_dir for the session directory inside the data folder")
                print("Use $scratch_dir for the scratch directory")
                input_folder = input(f'\nNew input folder for {step}: {data_dir}')
                if input_folder!="":
                    updateEntriesFromTable(cursor,"pipesteps",["input_folder"],[input_folder],["step"],[step])
                    
                    if step in dic_input.keys():
                        print(f'Suggested inputs: {dic_input[step]}')
                    input_files = input("New input files (separated by comma, use $sbj or $sess if needed): ")
                    if input_files!="":
                        updateEntriesFromTable(cursor,"pipesteps",["input_files"],[input_files],["step"],[step])
                            
                    updateDB(cnx,cursor)
                
                # Add output folder and files
                if curr_outdir!=None:
                    print(f'\nCurrent output folder: {curr_outdir}')
                print("Use $sbj or $sess if needed")
                print("Use $sess_dir for the session directory inside the data folder")
                print("Use $scratch_dir for the scratch directory")
                output_folder = input(f'\nNew output folder for {step}: {data_dir}')
                if output_folder!="":
                    updateEntriesFromTable(cursor,"pipesteps",["output_folder"],[output_folder],["step"],[step])
                    
                    if step in dic_output.keys():
                        print(f'Suggested outputs: {dic_output[step]}')
                    output_files = input(f'New output files (separated by comma, use $sbj or $sess if needed) for {step}: ')
                    if output_files!="":
                        updateEntriesFromTable(cursor,"pipesteps",["output_files"],[output_files],["step"],[step])
                        print("done")
                            
                    updateDB(cnx,cursor)
                    
                # Get the next step
                array = getDefaultSteps(pipe_type)
                if len(array)>0:
                    i = array.index(step)
                    if i<=len(array)-2:
                        print(f"\nSuggested next_step: {array[i+1]}")
                    else:
                        print("\nSuggestion: this could be the last step")
                    
                next_step = input(f"\nNext step to {step} (if different than '{ns}'): ")
                if next_step!="":
                    updateEntriesFromTable(cursor,"pipesteps",["next_step"],[next_step],["step"],[step])
                    updateDB(cnx,cursor)
                else:
                    next_step = ns if ns!=None else ""
                
                # There's no extra steps to add to the DB (from a previous step that had more than one next_step) before adding next_step
                ns_array = next_step.split(",")        
                if len(extra_steps)==0:
                    # Current step has only one next_step
                    if len(ns_array)<2:
                        step = next_step
                        
                    # Current step has more than one next_step
                    else:
                        # Add the first of those next_steps in the next while loop
                        step = ns_array[0]
                        del ns_array[0]
                        # Add the other next_steps to extra_steps
                        for item in ns_array:
                            if item not in extra_steps:
                                extra_steps+=[item]
                
                # There's other steps that need to be added to the DB before adding next_step (from a previous step that had more than one next_step)           
                else:
                    # Add the first of those extra_steps waiting to be added in the next while loop
                    step = extra_steps[0]
                    del extra_steps[0]
                    # Add all next_steps to extra_steps
                    for item in ns_array:
                        if item!="" and (item not in extra_steps):
                            extra_steps+=[item]
                
                print("---------------------------")

def startPipelineAll(cnx,cursor,pipeline,imgTypes):
    first_step = getValueFromTable(cursor,"first_step","pipelines",["pipename"],[pipeline])
    if first_step==None:
        return
    
    for sess in getSessList(cursor,pipeline):
        exist = 1
        for img in imgTypes:
            if not getValueFromTable(cursor,img,"sessions",["sess"],[sess]):
                exist = 0;
                break
        if exist==1:
            insertEntryToTable(cursor,"procs",["sess","step"],[sess,first_step])
            updateDB(cnx,cursor)

# Creates a new project from zero or modifies an existing project
def createProject(username,password,hostdir,location):
    ### Open a new connection with the DB ###
    cnx = connect(username,password,hostdir)
    
    ### Create the cursor ###
    cursor = cnx.cursor()
    # Set the number of seconds the server waits for activity on a connection before closing it to 86400 (24 hours)
    cursor.execute("SET GLOBAL wait_timeout=86400")
    # Set the number of seconds to wait for a lock timeout to 1 minute
    cursor.execute("SET GLOBAL innodb_lock_wait_timeout=60")
    
    ### Create the DB if not exist ###
    print("---------------------------\n# Current databases:")
    cursor.execute("show databases")
    databases = [item[0] for item in cursor.fetchall()]
    print("\n".join(databases)+"\n---------------------------")
    
    db = input("Database name: ")
    while db=="":
        db = input("Database name: ")
    if db not in databases:
        cursor.execute("create database "+db)
    cursor.execute("use "+db)
    updateDB(cnx,cursor)
    
    ### Create the basic tables ###
    print("\nCreating tables if not exist...")
    for cmd in commandsNewTables(db):
        print("\n"+cmd+"\n")
        cursor.execute(cmd)
    insertEntryToTable(cursor,"string2num",["variable","string_val","numerical_val","from_table"],['grp','CON',0,'subjects'])
    insertEntryToTable(cursor,"string2num",["variable","string_val","numerical_val","from_table"],['grp','PAT',1,'subjects'])
    updateDB(cnx,cursor)
    print("done")
    
    ### Populate projects (add new project) ###
    current_projects = colvalsTable(cursor,"projects","title",True)
    
    print("** The folder in $scratch must have the same name as the new project **")
    print("** Be sure to mount $scratch from the rcc and NOT just create it locally **")
    new_project = input("Project name ([Enter] for the first of current projects): ")
    if new_project=="":
        new_project = current_projects[0]
    
    if new_project not in current_projects:
        print("** Use $scratch or $home_dir if needed **")
        print("Select an existing folder:")
        insertEntryToTable(cursor,"projects",["title"],[new_project])
        updateDB(cnx,cursor)
    
    # Add one column per project to the subjects table to mark if that sbj belongs to the project    
    if new_project not in getColsFromTable(getTableInfo(cursor,"subjects")):
        cursor.execute("alter table subjects add column "+new_project+" boolean")
    
    ### Populate subjects ###
    subjects = addSubjects(username,password,hostdir,db,new_project,location) if input(f"\nAdding subjects to {new_project}? [Y]: ")=='Y' else colvalsTable(cursor,"subjects","sbjID")
    if len(subjects)==0:
        printError("At least one subject must be added")
        return
    
    ### Populate sessions ###
    print("\nEach line of sessions file: sbjID:sessID1,sessID2,...")
    sess_list = input("Path to sess_list ([Enter] if not adding any sessions): ")
    all_sessions = addSessions(username,password,hostdir,db,location,sess_list) if os.path.isfile(sess_list) else colvalsTable(cursor,"sessions","sess")
    if len(all_sessions)==0:
        print("At least one session must be added")
        return
    
    ### Populate pipelines ###
    print(f"Enter new pipeline names for project {new_project} divided by comma")
    print("[Enter] if not adding any pipeline")
    new_pipes = input(">> ")
    new_pipes = new_pipes.split(",") if new_pipes!="" else []
    pipe_type = getDefaultPipeline()
    pipelines = addPipes(cnx,cursor,db,new_project,new_pipes,pipe_type) if len(new_pipes)>0 else getPipelines(cursor,new_project)
    if len(pipelines)==0:
        print("You must add at least one pipelines")
        return
    
    ### Populate pipesteps from first step ###
    addSteps(cnx,cursor,location,new_project,pipe_type)
    
    ### Copy sbjs files ###
    imgTypes = input("list of images used divided by comma: ").split(",")
    for img in imgTypes:
        if img not in getColsFromTable(getTableInfo(cursor,"sessions")):
            print(f"Adding {img} column to sessions table...")
            cursor.execute("alter table sessions add column "+img+" tinyint(1)")
            print("done")
    
    if input("\nCopying original files to sbj folder? [Y]: ")=='Y':
        copyFiles(cnx,cursor,location,new_project,imgTypes)
                            
    ### Start running pipeline in sessions that are ready ###
    for pipe in pipelines:
        if input(f'\nStart {pipe} pipeline on all sessions with all files ({",".join(imgTypes)})? [Y]: ')=='Y':
            startPipelineAll(cnx,cursor,pipe,imgTypes)
            
    ### Create recovery files for the new DB ###
    if input("\nWas recreateTables.sql and reloadTables.sql already updated? [N]: ")=='N':
        createRecoveryFiles(cursor,location)
        
    ### Create the skeletton of the python script ###
    if input(f"\nCreate skeletton python script for {new_project}? [Y]: ")=='Y':
        for pipe in pipelines:
            createSkltPython(cursor,db,new_project,pipe)
    
    ### Close the connection ###
    print("###### DONE ######")
    cursor.close()
    cnx.close()

# This is not working well    
def delStep(username,password,hostdir,db_name,step):
    cnx = connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    
    # Unlink the step
    [next_step,pipeline] = getEntryFromTable(cursor,["next_step","pipeline"],"pipesteps",["step"],[step])
    if next_step!=None and getStringFromTable(cursor,"first_step","pipelines",["pipename"],[str(pipeline)])==step:
        updateEntriesFromTable(cursor,"pipelines",["first_step"],[next_step],["pipename"],[pipeline])
    else:
        # Get the previous steps
        prev_steps = getPreviousSteps(cursor,step)

        # Fix the links    
        for prev_step in prev_steps:
            array = getNextSteps(cursor,prev_step)
            array.remove(step)
            
            # The step we're removing had next step(s)
            # Add those steps to the next_steps of the previous step
            if next_step!=None:
                array+=[next_step]
                updateEntriesFromTable(cursor,"pipesteps",["next_step"],[",".join(array)],["step"],[prev_step])
                
            # The step we're removing didn't have a next step
            else:
                if len(array)==0:
                    updateEntriesFromTable(cursor,"pipesteps",["next_step"],[""],["step"],[prev_step])
                else:
                    updateEntriesFromTable(cursor,"pipesteps",["next_step"],[",".join(array)],["step"],[prev_step])
            
    # Remove processes
    deleteEntriesFromTable(cursor,"procs",["step"],[step])
    
    # Delete step from DB
    deleteEntryFromTable(cursor,"pipesteps",["step"],[step])
    
    # Should also delete the logs and output files if desired
    
    updateDB(cnx,cursor)
    cursor.close()
    cnx.close()

def insertResultFromDF(cursor,df):
    for sess,line in df.iterrows():
        for result,value in line.items():
            value = str(value)
            if value in empty_values:
                continue
            insertResult(cursor,sess,result,value)

def insertDemographicDataFromDF(cursor,df):
    for sbjID,line in df.iterrows():
        for measure,value in line.items():
            value = str(value)
            if value in empty_values:
                continue
            insertDemographicData(cursor,sbjID,measure,value)
        break
         
def copyFiles(cnx,cursor,location,project,imgTypes):
    for pipe in getPipelines(cursor,project):
        data_dir = getDataDir(cursor,pipe,location)
        if not os.path.isdir(data_dir):
            printError(data_dir+" not found")
            continue
        
        for img in imgTypes:
            print("Add *, ?, $sbj, $sess or $sessID when needed (or $data_dir for the destination folder)")
            img_path = input(f"{img} original path WITH NO EXTENSION: ")
            if img_path=="":
                continue
            ext_paths = input(f"{img} extensions in the original folder divided by comma (including the dot): ").split(",")
            dest_path = input(f"{img} destination folder: ")
            if not dest_path.endswith("/"):
                dest_path = dest_path+"/"
                
            for sess in getSessList(cursor,pipe):
                if isExcluded(cursor,sess):
                    continue
                print("Copying "+sess+" files...")
                        
                sessDir = getSessDir(cursor,pipe,location,sess)
                print("\n### sessDir: "+sessDir)
                if not os.path.isdir(sessDir):
                    os.makedirs(sessDir)
                if not os.path.isdir(sessDir):
                    printError("Could not create {sessDir}")
                    continue
                
                OK = True
                for ext in ext_paths:
                    print("\n## Copying "+ext[1:]+"...")
                    sbjID = getSbjID(cursor,sess)
                    sessID = getStringFromTable(cursor,"sessID","sessions",["sess"],[sess])
                    
                    orig = img_path.replace("$sessID",sessID).replace("$sess",sess).replace("$sbj",sbjID)
                    print("# Original file:")
                    print(orig)
                    if len(glob.glob(orig))==0:
                        printError(sess+" ("+img+"): file not found")
                        OK = False
                        break
                    
                    dest = dest_path.replace("$data_dir",data_dir).replace("$sessID",sessID).replace("$sess",sess).replace("$sbj",sbjID)+img+ext
                    print("# Destination file:")
                    print(dest)
                    if os.path.isfile(dest):
                        continue
                            
                    os.system(f"cp {orig} {dest}")
                    print("done")
                
                if OK:
                    updateEntriesFromTable(cursor,"sessions",[img],[1],["sess"],[sess])
                    updateDB(cnx,cursor)

def createSkltPython(cursor,db,project,pipeline):
    python_script = "/home/mkeith/Insync/DB/"+project+".py"
    if os.path.isfile(python_script):
        "/home/mkeith/Insync/DB/"+project+"_"+pipeline+".py"
    fout = open(python_script,'w')
    
    # Write the heather
    fout.write("#!/usr/bin/env python3\n"
               "from mysql.connector import Error\n"
               "from mysql.connector import errorcode\n"
               "import datetime\n"
               "import time\n"
               "import mylib\n\n"
    
               "# Global variables\n"
               f"pipeline = '{pipeline}'")
    
    # Write each done function
    steps = getStepsPipeline(cursor,pipeline)
    for step in steps:
        fout.write(f"def done{step.capitalize()}(cnx,location):\n"
                   "\tcursor = cnx.cursor()\n"
                   f"\tdatadir = mylib.getDataDir(cursor,'{pipeline}',location)\n"
                   f"\tlogdir = mylib.getScriptsDir(cursor,'{pipeline}',location)+\"/logs\"\n\n"
                   
                   "\t# Get the list of sessions that need to be checked\n"
                   f"\tarray = mylib.doneSbjs(cursor,'{step}')\n"
                   "\tif len(array)>0:\n"
                   f"\t\tmylib.printDescription(cursor,'{pipeline}')\n\n"
                   
                   "\tn = 0\n"
                   "\tfor sess in array:\n"
                   "\t\tprint(\"Checking \"+sess+\"...\")\n"
                   "\t\tn+=1\n\n"
                   
                   "\t\t# Commit all changes\n"
                   "\t\tmylib.updateDB(cnx,cursor)\n"
                   "\t\tif n>15:\n"
                   "\t\t\tbreak\n\n"
                   
                   "\tcursor.close()\n"
                   "\treturn len(array)-n+1\n\n")
            
    # Create main checkDone function
    fout.write("def checkDone(cnx,location,auto,screen):\n"
               "\tprint(\"\\n### DONE PROCS ###\")\n\n")
        
    fout.write("\twhile True:\n")
    i = 1
    line = "n1"
    for step in steps:
        fout.write("\t\tn"+str(i)+" = done"+step.capitalize()+"(cnx,location)\n")
        line+="+n"+str(i)
        i+=1
    fout.write("\t\tif "+line+"==0:\n")
    fout.write("\t\t\tbreak\n\n")
        
    # Create the function that runs the pipeline
    fout.write("def runPipeline(location,screen):\n")
    fout.write("\tprint(67 * '-')\n")
    fout.write("\tprint('### "+pipeline.upper()+" ###')\n")
    fout.write("\tprint(67 * '-')\n\n")
            
    fout.write("\topt = mylib.menu({}) if location!='hpc' else mylib.hpcMenu({})\n")
    fout.write("\twhile opt!='1':\n")
    fout.write(f"\t\tmylib.checkReadyProcs('{db}','{pipeline}',opt=='0',location,screen)\n")
    fout.write(f"\t\tmylib.RunningProcs('{db}','{pipeline}',opt=='0',location,test_mode=test_mode)\n")
    fout.write("\t\tcheckDone(cnx,location,opt=='0',screen)\n\n")
            
    fout.write("\t\tif opt!='0':\n")
    fout.write("\t\t\tbreak\n\n")
            
    fout.write("\t\t# Wait one hour between automatic loops\n")
    fout.write("\t\tprint(\"Went to sleep at \"+str(datetime.datetime.today().strftime(\"%H:%M:%S\")))\n")
    fout.write("\t\tprint(\"ZZZZZZ...\")\n")
    fout.write("\t\ttime.sleep(3600)\n")
    fout.write("\t\tprint(\"Wake up!\")\n\n")
            
    # Create main
    fout.write("def main():\n")
    fout.write("\ttry:\n")
    fout.write("\t\tlocation = input(\"Location? [hpc/petrov/petrov_damaged/oneDrive*]: \")\n")
    fout.write("\t\tif location=='':\n")
    fout.write("\t\t\tlocation = \"oneDrive\"\n")
    fout.write("\t\tscreen = False if input(\"Screen available? [N]: \")=='N' else True\n")
    fout.write("\t\trunPipeline(location,screen)\n\n")
        
    fout.write("\texcept Error as err:\n")
    fout.write("\t\tif err.errno == errorcode.ER_ACCESS_DENIED_ERROR:\n")
    fout.write("\t\t\tprintError(\"Wrong username or password\")\n")
    fout.write("\t\telse:\n")
    fout.write("\t\t\tprintError(err)\n\n")
            
    # Create executable
    fout.write("if __name__ == \"__main__\":\n")
    fout.write("\tmain()\n\n")
            
    fout.close()
    
def checkErrors(cnx,screen,pipeline,location):
    print(f'### ERRORS {pipeline} ###')
    cursor = cnx.cursor()
    
    # Get the list of steps for the pipeline
    pipesteps = getStepsPipeline(cursor,pipeline)

    # Get the information of each error
    for sess,line in getEntriesFromTable(cursor,["sess","step","notes"],"procs",["status"],["error"],["step"],"sess").iterrows():
        if line["step"] not in pipesteps:
            continue

        # Show the session and step
        print(sess+' '+line["step"]+" note: "+str(line["notes"]))

        # Ask new status and update DB
        valid_opts = ["ready","hold","done","checked"]
        newstatus = input("New status [ "+" / ".join(valid_opts)+" ]: ")
        if newstatus in valid_opts:
            if newstatus=="checked":
                checked(cursor,sess,line["step"],False,'',location)
            elif newstatus=="hold":
                note = input("Hold note: ")
                updateProcStatus(cursor,sess,line["step"],newstatus,location,note)
            else:
                updateProcStatus(cursor,sess,line["step"],newstatus,location)

    # Save changes
    updateDB(cnx,cursor)
    cursor.close()
    
def CheckOnHold(username,password,hostname,db_name,pipeline,location):
    print(f'\n### PROCS ON HOLD {pipeline} ###')
          
    # Open connection & cursor
    cnx = connect(username,password,hostname,db_name)
    cursor = cnx.cursor()

    # Get the list of steps for the pipeline
    pipesteps = getStepsPipeline(cursor,pipeline)

    # Get the list of sessions and steps that need to be cheched
    for sess,line in getEntriesFromTable(cursor,["sess","step","notes"],"procs",["status"],["hold"],["step"],"sess").iterrows():
        if (line["step"] not in pipesteps) or line["notes"]==None:
            continue
            
        # Step was missing input files in the sbj directory
        if line["notes"]=="missing input_files" and checkInputFiles(cursor,sess,line["step"],location):
            updateProcStatus(cursor,sess,line["step"],'ready',location)
            
        # Step was missing input files in the sbj scratch
        elif line["notes"]=="missing input_files (in scratch)" and checkInputFiles(cursor,sess,line["step"],location):
            updateProcStatus(cursor,sess,line["step"],'ready',location)
            
        # sess folder was being used by another process
        elif line["notes"].endswith("is currently running on this session"):
            if getStringFromTable(cursor,"status","procs",["sess","step"],[sess,line["notes"].split(" ")[0]])!="running":
                updateProcStatus(cursor,sess,line["step"],'ready',location)
                 
        # Step was waiting for some previous step(s) to finish
        elif line["notes"]=="previous step not checked yet" and checkPrevStepsCompleted(cursor,sess,line["step"]):
            updateProcStatus(cursor,sess,line["step"],'ready',location)
    
        updateDB(cnx,cursor)
        
    # Close cursor & connection
    cursor.close()
    cnx.close()
       
def progressReport(username,password,hostdir,db_name,project,pipeline,location):
    print(f"\n### PROGRESS REPORT {project} ({pipeline}) ###")
          
    cnx = connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    outDir = getProjDir(cursor,pipeline,location)+"/docs/"
    os.system(f"mkdir -p {outDir}")
    csv = outDir+pipeline.lower()+"_exclusion_report.csv"
    pdf = outDir+pipeline.lower()+"_progress_report.pdf"
    exclrep = open(csv,'w')
    
    # Group the subjects that are not excluded in lists of 10 subjects
    # To show one group per page
    included = []
    grp = []
    n_excl = 0
    for sess in getSessList(cursor,pipeline):
        if pipeline=="TBSS":
            DIFFsbj = getValueFromTable(cursor,"DIFFsbj","subjects",["sbjID"],[getSbjID(cursor,sess)])
            if DIFFsbj==None or DIFFsbj==0:
                continue
        if not isExcluded(cursor,sess):
            if len(grp)<10:
                grp+=[sess]
            else:
                included+=[grp]
                grp = []
                grp+=[sess]
        else:
            n_excl+=1
            exclrep.write(sess+','+getStringFromTable(cursor,"criteria","excluded",["sess"],[sess])+'\n')
    if len(grp)>0:
        included+=[grp]
    exclrep.close()
    
    # Get the list of steps for the pipeline
    steps = getStepsPipeline(cursor,pipeline)
    if len(steps)==0:
        printWarning("Cannot generate progress report. Probably the firt step is missing in the DB.")
        return
    
    # Generate the progress map for each subject and add it to the report
    merger = PdfFileMerger()  
    n = 1
    for grp in included:
        outfile = outDir+pipeline.lower()+"_progress_report_"+str(n)
        statslib.grpimg(cursor,grp,steps,outfile+".png")
        fileslib.convertToPDF(outfile+".png")
        merger.append(PdfFileReader(outfile+".pdf"))
        os.remove(outfile+".png")
        os.remove(outfile+".pdf")
        n+=1
    merger.write(pdf)
    merger.close()
    
    print(f"Progress report saved in: {pdf}")
    print(f"Exclusion report saved in: {csv}")
  
    cursor.close()
    cnx.close()
    
# Order a dictionary by its values in decreasing order (higher to lower)
# INPUT: dictionary
# OUTPUT: list of keys order by their value in dic
def orderDicByVals_dec(dic):
    sorted_keys = []
    sorted_values = []
    
    for key,value in dic.items():
        for i in range(len(sorted_values)+1):
            if i>len(sorted_values)-1 or value>sorted_values[i]:
                break
        sorted_keys = sorted_keys[0:i]+[key]+sorted_keys[i:len(sorted_keys)]
        sorted_values = sorted_values[0:i]+[value]+sorted_values[i:len(sorted_values)]
        
    return sorted_keys
    
# Order a dictionary by its values in increasing order (lower to higher)
# INPUT: dictionary
# OUTPUT: list of keys order by their value in dic
def orderDicByVals_inc(dic):
    sorted_keys = []
    sorted_values = []
    
    for key,value in dic.items():
        for i in range(len(sorted_values)+1):
            if i>len(sorted_values)-1 or value<sorted_values[i]:
                break
        sorted_keys = sorted_keys[0:i]+[key]+sorted_keys[i:len(sorted_keys)]
        sorted_values = sorted_values[0:i]+[value]+sorted_values[i:len(sorted_values)]
        
    return sorted_keys

def printStepInfo(username,password,hostdir,pipeline,database):
    cnx = connect(username,password,hostdir,database)
    cursor = cnx.cursor()
    
    print(f"Available steps: {getStepsString(cursor,pipeline)}")
    step = input("Select step: ")
    info = getColsFromTable(getTableInfo(cursor,"pipesteps"))
    print(f'Available columns: {",".join(info)}')
    cols = input("Select columns (divided by comma): ").split(",")
    
    for index,line in getEntriesFromTable(cursor,cols,"pipesteps",["step"],[step]).iterrows():
        for col,val in line.items():
            print(f"\n{step} {col}: {str(val)}")
            if input(f"Modify value of {col}? [Y]: ")=='Y':
                new_val = input("New value: ")
                updateEntriesFromTable(cursor,"pipesteps",[col],[new_val],["step"],[step])
                updateDB(cnx,cursor)
            
    cursor.close()
    cnx.close()
    
def cleanLogs(username,password,hostdir,db_name,pipeline,location):
    print("\n### Cleaning log files ###")
    cnx = connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    
    logDir = getScriptsDir(cursor,pipeline,location)+"/logs/"
    for f in glob.glob(logDir+'*'):
        fname = os.path.basename(f)
        parts = fname.split(".")
        info = parts[0].split("_")
        
        ## The log is NOT from TBSS
        if (not fname.startswith("randomise")) and (not fname.startswith("preStats_part1")) and (not fname.startswith("preStats_part2")):
            # The log is from a job that ran locally: step.sess
            
            if len(parts)==2 and getValueFromTable(cursor,"status","procs",["step","sess"],[parts[0],parts[1]]) in [None,"ready"]:
                    print("rm "+f)
                    os.system("rm "+f)
            
            # The log is from a job that ran in the cluster: step_sess.sh.[eo]jobID
            elif len(parts)!=2:
                if len(info)<2:
                    print("rm "+f)
                    os.system("rm "+f)
                if info[-1]=="1":
                    sess = "_".join(info[-2:])
                elif info[-1]=="A":
                    sess = "_".join(info[-3:])
                else:
                    sess = info[-1]
                step = parts[0].replace("_"+sess,"")
                jobID = parts[-1].replace("o","").replace("e","")
                if getValueFromTable(cursor,"status","procs",["step","sess","jobID"],[step,sess,jobID]) in [None,"ready"]:
                    print("rm "+f)
                    os.system("rm "+f)
                
        # The log is from TBSS
        else:
            # The log is from a job that run locally: [randomise,preStats_1,preStats_2]_tbssname.sh.log
            if fname.endswith(".log"):
                if getNumberEntries(cursor,"tbss",["name"],[info[-1]])==0:
                    print("rm "+f)
                    os.system("rm "+f)
            
            # The log is from a job that ran in the cluster: [randomise,preStats_1,preStats_2]_tbssname.sh.[eo]jobID
            elif fname.startswith("randomise"):
                tbssname = parts[0].replace("randomise_","")
                if getNumberEntries(cursor,"tbss",["name"],[tbssname])==0:
                    print("rm "+f)
                    os.system("rm "+f)
                elif getNumberEntries(cursor,"tbss",["name","randomise"],[tbssname,"%"+parts[-1].replace("e","").replace("o","")+"%"])==0:
                    print("rm "+f)
                    os.system("rm "+f)
                        
            elif fname.startswith("preStats_part1"):
                tbssname = parts[0].replace("preStats_part1_","")
                if getNumberEntries(cursor,"tbss",["name"],[tbssname])==0:
                    print("rm "+f)
                    os.system("rm "+f)
                elif getNumberEntries(cursor,"tbss",["name","preStats_part1"],[tbssname,"%"+parts[-1].replace("e","").replace("o","")+"%"])==0:
                    print("rm "+f)
                    os.system("rm "+f)
    
    cursor.close()
    cnx.close()
    print("done")
    
def insertProc(cursor,sess,nextstep,location,hold=False,note='',test_mode=False):
    # Are all previous steps are checked? (some steps are next_step of more than one step)
    # It is not enough that all the input files be there, previous steps should be checked before marking it as ready
    all_checked = checkPrevStepsCompleted(cursor,sess,nextstep)
    
    # Check if the step is already added for that session, and in that case what is the status and notes
    [status,notes] = getProcStatus(cursor,sess,nextstep)
    
    # Cut note if it's too long (if it's empty it won't do anything)
    type_notes = str(getTableInfo(cursor,"procs").loc["notes","Type"])
    length_notes = int(type_notes[type_notes.index("(")+1:type_notes.index(")")])
    if len(note)>length_notes:
        note = note[0:length_notes-1]
    
    # Step has not been added for that session
    marked_ready = False
    if status=='':
        # Put nextstep on hold 
        if hold:
            insertEntryToTable(cursor,"procs",["sess","step","status"],[sess,nextstep,'hold'],test_mode=test_mode)
            
        # All previous steps are checked, mark as ready
        elif all_checked:
            insertEntryToTable(cursor,"procs",["sess","step","status"],[sess,nextstep,"ready"],test_mode=test_mode)
            marked_ready = True
            
        # Some previous steps are not checked, put on hold
        else:
            insertEntryToTable(cursor,"procs",["sess","step","status"],[sess,nextstep,'hold'],test_mode=test_mode)
            note = 'previous step not checked yet'
    
    # Step was already added for that session, but has to be put on hold
    elif hold:
        updateProcStatus(cursor,sess,nextstep,'hold',location,test_mode=test_mode)
        
    # Step was already added for that session, but not all previous steps are checked
    elif not all_checked:
        updateProcStatus(cursor,sess,nextstep,'hold',location,test_mode=test_mode)
        note = "previous step not checked yet"
        
    # Step was already added for that session and all previous steps are checked
    else:
        updateProcStatus(cursor,sess,nextstep,'ready',location,test_mode=test_mode)
        marked_ready = True
        
    if note!="":
        updateNotesProc(cursor,sess,nextstep,note,test_mode=test_mode)
            
    return marked_ready

def addStep(username,password,hostdir,db_name,pipeline,location):
    cnx = connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    
    steps = getStepsPipeline(cursor,pipeline)
    step = input("new step: ")
    if step not in steps:
        # Get previous steps
        print(f'\nCurrent steps: {",".join(steps)}')
        prev_steps = input(f"Previous steps for {step} divided by comma ([Enter] if {step} will be the 1st step): ")
        while prev_steps!="" and not (set(prev_steps.split(",")) <= set(steps)):
            print(f"\nOne or more elements in {prev_steps} dont exist")
            print(f'Current steps: {",".join(steps)}')
            prev_steps = input(f"\nPrevious steps for {step} divided by comma ([Enter] if {step} will be the 1st step): ")
        # If previous steps is empty, this is the new first step
        if prev_steps=="":
            updateEntriesFromTable(cursor,"pipelines",["first_step"],[step],["pipename"],[pipeline])
        
        # Next steps of the new step
        next_steps = ""
        # There are no previous steps
        # step is the new first step and its next_step should be the old first step
        # if len(steps)==0, there was no first step, the new step is the first & only step & prev_next_steps & next_steps stay empty
        if prev_steps=="" and len(steps)>0:
            next_steps = steps[0]
        # There are previous steps, ns is the next step of that previous step
        elif prev_steps!="":
            for prev_step in prev_steps.split(","):
                ns = getValueFromTable(cursor,"next_step","pipesteps",["step"],[prev_step])
        
                # If the previous step didn't have a next_step
                if ns==None:
                    updateEntryFromTable(cursor,"pipesteps",["next_step"],[step],["step"],[prev_step])
                # The previous step had a next_step
                else:
                    ns1 = ns+","+step
                    print("What should be the next step for "+prev_step+"?")
                    print("1. "+ns1)
                    print("*2. "+step)
                    opt = input(">> ")
                    
                    # Option 1 was selected
                    # Concatenate the new step to the previous next_step for previous step
                    # In this case, next_step for the new step can be left empty (for now)
                    if opt=="1":
                        updateEntryFromTable(cursor,"pipesteps",["next_step"],[ns1],["step"],[prev_step])
                    
                    # Option 2 was selected
                    # Put the new step between the previous step & it's current next step
                    # In this case next_step for the new step needs to be the old next_step for the previous step
                    else:
                        updateEntryFromTable(cursor,"pipesteps",["next_step"],[step],["step"],[prev_step])
                        if next_steps=="":
                            next_steps = ns
                        elif ns not in next_steps.split(","):
                            next_steps+=","+ns
                    
        help_txt = "(use $sess/$sbj/$sess_dir as needed)"        
        desc = input(f"{step} description: ")
        local_script = input(f"{step} LOCAL script {help_txt}: ")
        scripts_dir = getScriptsDir(cursor,pipeline,location)+"/"
        if os.path.isdir(scripts_dir):
            if local_script!="" and not os.path.isfile(scripts_dir+local_script.split(" ")[0]):
                os.system("vi "+local_script)
            script = input(f"{step} CLUSTER script: ")
            if script!="" and not os.path.isfile(scripts_dir+script.split(" ")[0]):
                os.system("vi "+script)
            
        input_folder = input(f"{step} input folder {help_txt}: ")
        if input_folder!="":
            input_files = input(f"{step} mandatory input files {help_txt}: ")
            opt_input_files = input(f"{step} optional input files {help_txt}: ")
                
        output_folder = input(f"{step} output folder {help_txt}: ")
        if output_folder!="":
            output_files = input(f"{step} mandatory output files {help_txt}: ")
            opt_output_files = input(f"{step} optional output files {help_txt}: ")
            
        if next_steps=="":
            print(f'\nCurrent seps: {",".join(steps)}')
            next_steps = input(f"\nNext step for {step} ([Enter] if it's still unknown or hasn't been added yet): ")
        while next_steps!="" and not (set(next_steps.split(",")) <= set(steps)):
            print(f"\none or more steps in {next_steps} dont exist")
            print(f'Current seps: {",".join(steps)}')
            next_steps = input(f"\nNext step for {step} ([Enter] if next_step is still unknown or hasn't been added yet): ")
        
        # Create new step
        cols = ["step","pipeline","description","local_script","script","input_folder","input_files","input_files_opt","output_folder","output_files","output_files_opt","next_step"]
        vals = [step,pipeline,desc,local_script,script,input_folder,input_files,opt_input_files,output_folder,output_files,opt_output_files,next_steps]
        insertEntryToTable(cursor,"pipesteps",cols,vals)
    else:
        print(step+" already exists for "+pipeline)
        
    # We're not adding a first step
    mark_hold = []
    if prev_steps!="":
        mark_ready = input(f'\nList of sbjs divided by comma to mark {step} as ready ([Enter] to mark on all that have {prev_steps} checked): ')
        mark_ready = [] if mark_ready=="" else mark_ready.split(",")
        
        add_all = len(mark_ready)==0
        checked = []
        for prev_step in prev_steps.split(","):
            for sess in getColumnFromTable(cursor,"sess","procs",["step","status"],[prev_step,"checked"]):
                if sess not in checked:
                    checked+=[sess]
                if (sess in mark_ready) or (sess in mark_hold):
                    continue
                elif add_all:
                    mark_ready+=[sess]
                else:
                    mark_hold+=[sess]
        print(f'{len(checked)} sessions had at least one previous step marked as checked')
        if len(mark_ready)>0:
            print(f'Will attempt to mark the following sessions as ready: {",".join(mark_ready)}')
        if len(mark_hold)>0:
            print(f'Will attempt to mark the following sessions as hold: {",".join(mark_hold)}')
    
    # We're adding a first step
    else:
        mark_ready = getSessList(cursor,pipeline)
    
    # Marking sessions as hold        
    for sess in mark_hold:
        insertProc(cursor,sess,step,location,True,"other sbjs running first")
    print(f'{len(mark_hold)} sessions marked as hold')
    if len(mark_hold)>0:
        print(f'\nNow on hold: {",".join(mark_hold)}')
    
    # Marking sessions as hold
    new_hold = []
    for sess in mark_ready:
        if not insertProc(cursor,sess,step,location):
            mark_ready.remove(sess)
            new_hold+=[sess]
    print(f'{len(mark_ready)} sessions marked as ready')
    if len(mark_ready)>0:
        print(f'\nNow ready: {",".join(mark_ready)}')
    print(f'{len(new_hold)} sessions marked as hold because not all previous steps are checked')
    if len(new_hold)>0:
        print(f'\nNew on hold: {",".join(new_hold)}')
    
    # Is better to update at the end in case anything goes wrong then nothing should be saved
    updateDB(cnx,cursor)
    cursor.close()
    cnx.close()

def error(cursor,notes,sess,step,location,scratchdir="",ofile="",efile="",replace_notes=True,test_mode=False):
    updateProcStatus(cursor,sess,step,'error',location,notes,replace=replace_notes,test_mode=test_mode)
    
    # mv the logs and the session scratch dir to the failed folder if running in the cluster
    # ofile and efile are '' when calling from runningLocal
    if ofile!="" and efile!="" and os.path.isdir(scratchdir):
        if os.path.isfile(ofile):
            print(f"Moving {ofile} to failed location: {scratchdir}/failed")
            if not test_mode:
                os.system(f"mv {ofile} {scratchdir}/failed")
        if os.path.isfile(efile):
            print(f"Moving {efile} to failed location: {scratchdir}/failed")
            if not test_mode:
                os.system(f"mv {efile} {scratchdir}/failed")
        
        sess_scratch = scratchdir+'/'+sess
        if mvFiles(cursor,step) and os.path.isdir(sess_scratch):
            print(f"Moving {sess_scratch} scratch folder to failed location:{scratchdir}/failed/{sess}_{step}...")
            
            if not test_mode:
                os.system(f"mkdir -p {scratchdir}/failed/{sess}_{step}")
                os.system(f"mv {sess_scratch}/* {scratchdir}/failed/{sess}_{step}")
                os.system(f"rm -r {sess_scratch}")
    
    printError(notes)

def checked(cursor,sess,step,hold,note,location,test_mode=False):
    updateProcStatus(cursor,sess,step,'checked',location,test_mode=test_mode)
    
    pipeline = getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    next_steps = getNextSteps(cursor,step)
    print(f'next steps: {",".join(next_steps)}')
    if test_mode:
        print("Previous steps will not be checked since its in test mode")
    if len(next_steps)>0:
        move_files = mvFiles(cursor,step)
        
        for next_step in next_steps:
            if move_files and (not checkInputFiles(cursor,sess,next_step,location)):
                printWarning(f"MISSING INPUT FILES FOR NEXT STEP: {next_step}")
                insertProc(cursor,sess,next_step,location,True,"missing input_files",test_mode=test_mode)
            elif (not move_files) and (not checkInputFiles(cursor,sess,next_step,location)):
                printWarning(f"MISSING INPUT FILES FOR NEXT STEP: {next_step}")
                insertProc(cursor,sess,next_step,location,True,"missing input_files (in scratch)",test_mode=test_mode)
            elif not hold:
                insertProc(cursor,sess,next_step,location,test_mode=test_mode)
            else:
                insertProc(cursor,sess,next_step,location,True,note,test_mode=test_mode)
    
    saveExecTime(cursor,sess,step,getScriptsDir(cursor,pipeline,location)+"/logs",test_mode=test_mode)
    print("OK")

def done(cursor,sess,step,scratchdir,ofile,efile,location,move_files=True,test_mode=False):
    OK = True
    if os.path.isdir(scratchdir)and move_files:
        OK = moveFromScratch(cursor,sess,step,scratchdir,ofile,efile,location,test_mode=test_mode)
    
    if OK:
        updateProcStatus(cursor,sess,step,'done',location,test_mode=test_mode)
        print("OK")

def checkPrevStepsCompleted(cursor,sess,step):
    for prev in getPreviousSteps(cursor,step):
        if getStringFromTable(cursor,"status","procs",["step","sess"],[prev,sess])!="checked":
            return False
    return True

def updateNotesProc(cursor,sess,step,new_notes,replacenote=True,test_mode=False):
    new_notes = new_notes.replace('\n','').replace('\'','').replace('\t','').strip()
    old_notes = getNotesProc(cursor,sess,step).strip()
    if replacenote or old_notes=="" or old_notes=="None":
        notes = new_notes
    elif new_notes=="":
        notes = old_notes
    else:
        notes = old_notes+". "+new_notes if not old_notes.endswith(".") else old_notes+" "+new_notes
    
    type_notes = str(getTableInfo(cursor,"procs").loc["notes","Type"])
    length_notes = int(type_notes[type_notes.index("(")+1:type_notes.index(")")])
    if len(notes)>length_notes:
        notes = notes[0:length_notes-1]
    updateEntryFromTable(cursor,"procs",["notes"],[notes],["step","sess"],[step,sess],test_mode)
        
def updateNotesSess(cursor,sess,new_notes,replace=True):
    new_notes = new_notes.replace('\n','').replace('\'','').replace('\t','').strip()
    old_notes = getNotesSess(cursor,sess).strip()
    if replace or old_notes=="" or old_notes=="None":
        notes = new_notes
    else:
        notes = old_notes+". "+new_notes if not old_notes.endswith(".") else old_notes+" "+new_notes
    
    type_notes = str(getTableInfo(cursor,"sessions").loc["notes","Type"])
    length_notes = int(type_notes[type_notes.index("(")+1:type_notes.index(")")])
    if len(notes)>length_notes:
        notes = notes[0:length_notes-1]
    updateEntryFromTable(cursor,"sessions",["notes"],[notes],["sess"],[sess])
        
def updateNotesSbj(cursor,sbjID,new_notes,replace=True):
    new_notes = new_notes.replace('\n','').replace('\'','').replace('\t','').strip()
    old_notes = getNotesSbj(cursor,sbjID).strip()
    if replace or old_notes=="" or old_notes=="None":
        notes = new_notes
    else:
        notes = old_notes+". "+new_notes if not old_notes.endswith(".") else old_notes+" "+new_notes
    
    type_notes = str(getTableInfo(cursor,"subjects").loc["notes","Type"])
    length_notes = int(type_notes[type_notes.index("(")+1:type_notes.index(")")])
    if len(notes)>length_notes:
        notes = notes[0:length_notes-1]
    updateEntryFromTable(cursor,"subjects",["notes"],[notes],["sbjID"],[sbjID])

def updateNotesTBSS(cursor,tbss_name,new_notes,replace=True):
    new_notes = new_notes.replace('\n','').replace('\'','').replace('\t','').strip()
    old_notes = getNotesTBSS(cursor,tbss_name).strip()
    if replace or old_notes=="" or old_notes=="None":
        notes = new_notes
    else:
        notes = old_notes+". "+new_notes if not old_notes.endswith(".") else old_notes+" "+new_notes
    
    type_notes = str(getTableInfo(cursor,"tbss").loc["notes","Type"])
    length_notes = int(type_notes[type_notes.index("(")+1:type_notes.index(")")])
    if len(notes)>length_notes:
        notes = notes[0:length_notes-1]
    updateEntryFromTable(cursor,"tbss",["notes"],[notes],["name"],[tbss_name])
        
def updateProcStatus(cursor,sess,step,status,location,notes="",replace=False,test_mode=False):
    # An invalid status was entered
    if status not in ["hold","ready","running","done","checked","error"]:
        printError(f"unrecognized status: {status}")
        return False
    
    # The process does not exist
    if getNumberEntries(cursor,"procs",["sess","step"],[sess,step])==0:
        printError(f"couldnt update proc, doesnt exist ({sess},{step})")
        return False
    
    # Update the status
    updateEntryFromTable(cursor,"procs",["status"],[status],["sess","step"],[sess,step],test_mode)
    
    # Unless the process is now on hold or error, any error message should be deleted
    old_notes = getNotesProc(cursor,sess,step).lower()
    error_msgs = ["output file missing","missing input_files","previous step not checked yet","running on this session","in queue for","error","bad output"]
    if status not in ["error","hold"]:
        for msg in error_msgs:
            if msg in old_notes:
                replace = True
                break
            
    # Remove previous logs if they exist if status is now ready
    if status=="ready":
        pipeline = getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
        localLog = getScriptsDir(cursor,pipeline,location)+'/logs/'+step+'.'+sess
        scratchdir = getScratchDirPipe(cursor,location,pipeline)
        ofile = scratchdir+'/'+step+'_'+sess+".sh.o*"
        efile = scratchdir+'/'+step+'_'+sess+".sh.e*"
        os.system("rm -f "+localLog+" "+ofile+" "+efile)
    
    # Update notes for the process
    updateNotesProc(cursor,sess,step,notes,replace,test_mode)
    
    return True

def saveExecTime(cursor,sess,step,logDir,test_mode=False):
    [jobID,cluster] = getEntryFromTable(cursor,["jobID","cluster"],"procs",["sess","step"],[sess,step])
    [ofile,efile] = getLogs(logDir,sess,step,jobID,cluster)
        
    if os.path.isfile(ofile) and os.path.isfile(efile):
        # Get execution time
        try:
            fo = open(ofile,'r')
            foundline = False
            for line in fo:
                if re.match(r"Total execution time was [0-9][0-9] hrs [0-9][0-9] mins [0-9][0-9] secs",line):
                    foundline = True
                    break         
        except ValueError as err:
            printError(f'Could not save execution time. Error reading log file: {err}')
            
        # Save execution time if found in the file               
        else:
            if foundline:
                array = line.split(' ')
                if len(array)==10:
                    print(f'Saving execution time on {step} {sess}...')               
                    # Transform it to seconds
                    exectimesec = int(array[8])+int(array[6])*60+int(array[4])*3600              
                    # Save it to the DB
                    updateEntryFromTable(cursor,"procs",["exectimesec"],[exectimesec],["sess","step"],[sess,step],test_mode=test_mode)
                else:
                    printWarning("Could not save execution time. Information not found in log file.")
            else:
                printWarning("Could not save execution time. Information not found in log file.")
        finally:
            fo.close()
    else:
        printWarning(f"Could not save execution time. Log files not found: {ofile} {efile}")
  
def getDefaultPipeline():
    print("Select one option:\n"
          "1. ECP pipeline\n"
          "2. Basic DTI pipeline\n"
          "3. Basic RS pipeline\n"
          "4. Other*")
    opt = input(">> ")
    
    if opt=="1":
        return "ECP"
    elif opt=="2":
        return "DTI"
    elif opt=="3":
        return "RS"
    else:
        return "other"

def getStepsInfo(username,password,hostdir,db_name,pipeline="",automatic=False,step=""):
    if pipeline=="" and step=="":
        return
    elif pipeline!="":
        print(f"\n### Steps info for {pipeline} ###")
    elif step!="":
        print(f"\n### Step info for {step} ###")
              
    cnx = connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    
    if pipeline!="":
        steps = getStepsPipeline(cursor,pipeline)
    else:
        steps = [step]
        
    for step in steps:
        print(f"\n## {step}")
        for status in getColumnFromTable(cursor,"status","procs",["step"],[step],["status"],distinct=True):
            print(status+": "+str(getNumberEntries(cursor,"procs",["step","status"],[step,status])))
        if not automatic:
            input("[Enter]")
        
    cursor.close()
    cnx.close()

def getProcStatus(cursor,sess,step):
    if not inProcs(cursor,sess,step):
        return ['','']
    return [x if x!=None else "" for x in  getEntryFromTable(cursor,["status","notes"],"procs",["sess","step"],[sess,step])]

def getProjDir(cursor,pipeline,location):
    home_dir = getStringFromTable(cursor,"home_dir","generalInfo.locations",["location"],[location])
    scratch_dir = getStringFromTable(cursor,"scratch","generalInfo.locations",["location"],[location])
    
    project_name = getStringFromTable(cursor,"project","pipelines",["pipename"],[pipeline])
    project_dir = getStringFromTable(cursor,"project_dir","projects",["title"],[project_name]).replace("$home_dir",home_dir).replace("$scratch",scratch_dir)
    
    return project_dir

def getDataDir(cursor,pipeline,location):
    project_dir = getProjDir(cursor,pipeline,location)
    datadir_name = getValueFromTable(cursor,"data_dir","pipelines",["pipename"],[pipeline])
    return f"{project_dir}/{datadir_name}"

def getSessDir(cursor,pipeline,location,sess):
    datadir = getDataDir(cursor,pipeline,location)
    sessdir_path = getStringFromTable(cursor,"sess_dir","pipelines",["pipename"],[pipeline]).replace("$sess",sess).replace("$sbj",getSbjID(cursor,sess))
    return f"{datadir}/{sessdir_path}"

def getScratchDirPipe(cursor,location,pipeline):
    scratch_dir = getValueFromTable(cursor,"scratch","generalInfo.locations",["location"],[location])
    project_name = getValueFromTable(cursor,"project","pipelines",["pipename"],[pipeline])
    return f"{scratch_dir}/{project_name}"

def getScriptsDir(cursor,pipeline,location):
    project_dir = getProjDir(cursor,pipeline,location)
    fname = getValueFromTable(cursor,"scripts_dir","pipelines",["pipename"],[pipeline])
    return f"{project_dir}/{fname}"

def getNextSteps(cursor,step):
    next_steps = getStringFromTable(cursor,"next_step","pipesteps",["step"],[step])
    return next_steps.split(",") if next_steps!="None" else []

def doneSess(cursor,step):
    return getColumnFromTable(cursor,"sess","procs",["step","status"],[step,"done"])

def getPreviousSteps(cursor,step):
    prev_steps = []
    pipeline = getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    for pipeline_step in getStepsPipeline(cursor,pipeline):
        if step in getStringFromTable(cursor,"next_step","pipesteps",["step"],[pipeline_step]).split(","):
            prev_steps+=[pipeline_step]
    return prev_steps

def getStepsPipeline(cursor,pipeline="",joint=False,step_start=""):
    ## Get the first step for the pipeline
    if step_start=="":
        first_step = getStringFromTable(cursor,"first_step","pipelines",["pipename"],[pipeline])
    elif pipeline!="":
        # Do this to make sure step_start is part of the pipeline
        first_step = getStringFromTable(cursor,"step","pipesteps",["step","pipeline"],[step_start,pipeline])
    else:
        first_step = step_start
    if first_step=="None":
        return []
    
    ## Get the next steps of the first step
    next_steps = getNextSteps(cursor,first_step)
    
    ## Add each level of the tree as a different element of the array steps
    steps = [first_step]
    while len(next_steps)>0:
        # Add next tree level to steps
        if joint:
            steps+=[",".join(next_steps)]
        else:
            steps+=next_steps
        
        # Get the array of steps for the next level
        array = []
        for next_step in next_steps:
            array+=getNextSteps(cursor,next_step)
            
        # Next level to add
        next_steps = []
        for step in array:
            if (not step in steps) and (step not in next_steps):
                # A step could have a next step that goes into a different pipeline
                if pipeline=="" or getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])==pipeline:
                    next_steps+=[step]
    
    return steps

def getStepsString(cursor,pipeline):
    return " > ".join(getStepsPipeline(cursor,pipeline,joint=True))

def checkRunningLocal(step,sess,cursor,automatic,location,test_mode=False):
    localLog = getScriptsDir(cursor,getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step]),location)+'/logs/'+step+'.'+sess
    print(f"Local log: {localLog}")
    
    # If the local log file exists
    if os.path.isfile(localLog):
        # Check if there's an obvious error in the log
        try:
            errorMsg = fileslib.errorInLog(localLog)
        except ValueError as err:
            error(cursor,"exception reading log: "+str(err),sess,step,location,test_mode=test_mode)
        else:
            if errorMsg!='':
                error(cursor,errorMsg,sess,step,location,test_mode=test_mode)
            else:
                # Check if script seem to have reach the end with no errors
                try:
                    isDone = fileslib.lineInFile(localLog,"DONE "+step)
                except ValueError as err:
                    error(cursor,"exception reading log: "+str(err),sess,step,location,test_mode=test_mode)
                else:
                    if isDone:
                        done(cursor,sess,step,'','','',location,test_mode=test_mode)
                    else:
                        print(f"{sess} ({step}) appears to still be running")
            
    # If there is no local log file, it's running in the command line
    elif not automatic:
        print("It's supposed to be running locally but no local log file was found. Is it running in the comand line?")
        print(f"1. {sess} Still running [Enter]")
        print(f"2. {sess} Done running")
        print(f"3. {sess} Failed")
        opt = input(">> ")
        if opt=='2':
            done(cursor,sess,step,'','','',location,test_mode=test_mode)
        elif opt=='3':
            error(cursor,'failed locally with no log',sess,step,location,test_mode=test_mode)

def checkDone(username,password,hostdir,db_name,pipeline,location):
    print(f"\n### DONE PROCS {pipeline} ###")

    cnx = connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    pipesteps = getStepsPipeline(cursor,pipeline)
    step_sessions = {}
    for key,line in getEntriesFromTable(cursor,["step"],"procs",["status"],["done"],["step"],distinct=True).iterrows():
        step = line["step"]
        if step not in pipesteps:
            continue
        step_sessions[step] = doneSess(cursor,step)
    cursor.close()
    cnx.close()
    
    for step,sess_list in step_sessions.items():
        if step=="start":
            ECP.doneStart(sess_list,username,password,hostdir,db_name,location)
        elif step=="camino_wlf":
            ECP.doneCamino(sess_list,username,password,hostdir,db_name,True,location,False)
        elif step=="finaldtifit":
            ECP.doneDTIFIT(sess_list,username,password,hostdir,db_name,True,False,location)
    
    print(f"\n### FINISHED CHECKING DONE PROCS FOR {pipeline} ###\n")

def RunningProcs(username,password,hostdir,db_name,pipeline,automatic,location,test_mode=False):
    print(f"\n### RUNNING PROCS {pipeline} ###")
    
    # Open connection & cursor
    cnx = connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    
    scratchdir = getScratchDirPipe(cursor,location,pipeline)
    print(f"## Scratch dir: {scratchdir}")
    logDir =  getScriptsDir(cursor,pipeline,location)+"/logs"
    print(f"## Logs dir: {logDir}")

    # Get the list of steps for the pipeline
    pipesteps = getStepsPipeline(cursor,pipeline)
    
    for key,line in getEntriesFromTable(cursor,["sess","step","jobID","qsub_time","cluster"],"procs",["status"],["running"],["step"]).iterrows():
        if line["step"] not in pipesteps:
            continue
        
        ## Check if files should be rm from scratch after proc runs ##
        move_files = mvFiles(cursor,line["step"])
        
        ## GET THE LOGS PATH ##
        [ofile,efile] = getLogs(scratchdir,line['sess'],line["step"],line["jobID"],line["cluster"]) if line["jobID"]!=None else getLogs(logDir,line['sess'],line["step"],line["jobID"],line["cluster"])
        
        print(f"\nChecking {line['step']} on {line['sess']}...")
        print(f"Ran in {line['cluster']} cluster") if line["cluster"]!=None else print("Ran locally")
        print(f"Output log: {ofile}")
        print(f"Error log: {efile}")
        
        ## Check local job ##
        if str(line["jobID"]) in empty_values:
            checkRunningLocal(line["step"],line["sess"],cursor,automatic,location,test_mode=test_mode)
            updateDB(cnx,cursor)
            
        ## Check if job is still in queque or didn't run (log file doesn't exist) ##
        elif (not os.path.isfile(ofile)) or (not os.path.isfile(efile)):
            # If there is information about the submission date/time, check how long it's being quequed
            if str(line["qsub_time"]) not in empty_values:
                # Get current date and time
                now_date = str(datetime.date.today().strftime("%Y_%m_%d"))
                now_time = str(datetime.datetime.today().strftime("%H:%M:%S"))
                now_dt = now_date+'_'+now_time
                
                # Calculate the difference between current date/time and submission date/time
                tdelta = str(datetime.datetime.strptime(now_dt,"%Y_%m_%d_%H:%M:%S") - datetime.datetime.strptime(line["qsub_time"],"%Y_%m_%d_%H:%M:%S"))
                tdelta_array = tdelta.split(',')
                if len(tdelta_array)==1:
                    diff_days = 0
                    diff_hours = int(tdelta_array[0].split(':')[0])
                else:
                    diff_days = int(tdelta_array[0].split(' ')[0])
                    diff_hours = int(tdelta_array[1].split(':')[0])
                print(f"in queque for {tdelta}")
                
                # Add message if it has been in queque for more than 10 hours
                if diff_days>0 or diff_hours>12:
                    updateNotesProc(cursor,line['sess'],line["step"],f"in queue for {tdelta}",replacenote=True,test_mode=test_mode)
                    updateDB(cnx,cursor)
            else:
                print("in queue for unknown time")
        
        ## Job is running or ran in the cluster ##
        else:
            # Check if job finished running
            try:
                finished = procFinishedRunning(line["step"],ofile,location,str(line["jobID"]))
            except ValueError as err:
                error(cursor,str(err),line['sess'],line["step"],location,scratchdir,ofile,efile,test_mode=test_mode)
                updateDB(cnx,cursor)
                continue
            
            # Job finished running
            if finished:
                print("finished")
                # Check if script seem to have reach the end with no errors
                try:
                    isDone = fileslib.lineInFile(ofile,"DONE "+line["step"])
                except ValueError as err:
                    error(cursor,str(err),line['sess'],line["step"],location,scratchdir,ofile,efile,test_mode=test_mode)
                    updateDB(cnx,cursor)
                    continue
                
                if isDone:
                    done(cursor,line['sess'],line["step"],scratchdir,ofile,efile,location,move_files,test_mode=test_mode)
                else:
                    error(cursor,'check log files',line['sess'],line["step"],location,scratchdir,ofile,efile,test_mode=test_mode)
                updateDB(cnx,cursor)
            
            # Job is still running or never finished because it encountered an error
            else:
                print("not finished")
                # Check if job produced an error
                try:
                    errorMsg1 = fileslib.errorInLog(efile)
                    errorMsg2 = fileslib.errorInLog(ofile)
                except ValueError as err:
                    error(cursor,str(err),line['sess'],line["step"],location,scratchdir,ofile,efile,test_mode=test_mode)
                    updateDB(cnx,cursor)
                    continue
                
                # There is an error in the elog file
                if errorMsg1!='':
                    error(cursor,errorMsg1,line['sess'],line["step"],location,scratchdir,ofile,efile,test_mode=test_mode)
                    updateDB(cnx,cursor)
                        
                # There is an error in the olog file
                elif errorMsg2!='':
                    error(cursor,errorMsg2,line['sess'],line["step"],location,scratchdir,ofile,efile,test_mode=test_mode)
                    updateDB(cnx,cursor)
                    
                # Job is still running
                else:
                    print("running")

    # Close the cursor & connection
    cursor.close()
    cnx.close()
    print(f"\n### FINISHED CHECKING RUNNING PROCS FOR {pipeline} ###\n")

def getScratchOutputFiles(cursor,sess,step,scratchdir):
    [output_files,opt_output_files] = getEntryFromTable(cursor,["output_files","output_files_opt"],"pipesteps",["step"],[step])
    output_files = [] if output_files==None else output_files.replace("$sess",sess).replace("$sbj",getSbjID(cursor,sess)).split(',')
    opt_output_files = [] if opt_output_files==None else opt_output_files.replace("$sess",sess).replace("$sbj",getSbjID(cursor,sess)).split(',')
    scratchdir = scratchdir[:-1] if scratchdir.endswith("/") else scratchdir
    
    return [f"{scratchdir}/{sess}/{os.path.basename(outfile)}" for outfile in output_files+opt_output_files]

def getOutputDir(cursor,sess,step,location):
    # it's ok if it doesn't exist because when I'm moving files to the output folder it'll be created
    [output_folder,pipe] = getEntryFromTable(cursor,["output_folder","pipeline"],"pipesteps",["step"],[step])
    if mvFiles(cursor,step):
        output_folder = str(output_folder).replace("$sess_dir",getSessDir(cursor,pipe,location,sess))
    else:
        output_folder = str(output_folder).replace("$scratch_dir",getScratchDirPipe(cursor,location,pipe))
    
    return output_folder.replace("$sess",sess).replace("$sbj",getSbjID(cursor,sess))

def getOutputFiles(cursor,sess,step,location,optional=False):
    [output_files,opt_output_files] = getEntryFromTable(cursor,["output_files","output_files_opt"],"pipesteps",["step"],[step])
    output_folder = getOutputDir(cursor,sess,step,location)
    sbjID = getSbjID(cursor,sess)
    output_files = [] if output_files==None else output_files.replace("$sess",sess).replace("$sbj",sbjID).split(',')
    opt_output_files = [] if opt_output_files==None or not optional else opt_output_files.replace("$sess",sess).replace("$sbj",sbjID).split(',')
    
    return [f"{output_folder}/{outfile}" for outfile in output_files+opt_output_files]

def checkOutputFiles(cursor,sess,step,location):
    return checkFilesExist(getOutputFiles(cursor,sess,step,location))

def getScratchInputFiles(cursor,sess,step,scratchdir):
    [input_files,opt_input_files] = getEntryFromTable(cursor,["input_files","input_files_opt"],"pipesteps",["step"],[step])
    input_files = [] if input_files==None else input_files.replace("$sess",sess).replace("$sbj",getSbjID(cursor,sess)).split(',')
    opt_input_files = [] if opt_input_files==None else opt_input_files.replace("$sess",sess).replace("$sbj",getSbjID(cursor,sess)).split(',')
    scratchdir = scratchdir[:-1] if scratchdir.endswith("/") else scratchdir
    
    return [f"{scratchdir}/{sess}/{os.path.basename(infile)}" for infile in input_files+opt_input_files]

def getInputDir(cursor,sess,step,location):
    [indir,pipename] = getEntryFromTable(cursor,["input_folder","pipeline"],"pipesteps",["step"],[step])
    sess_dir = getSessDir(cursor,pipename,location,sess)
    scratch = getScratchDirPipe(cursor,location,pipename)
    
    return str(indir).replace("$sess_dir",sess_dir).replace("$scratch_dir",scratch).replace("$sess",sess).replace("$sbj",getSbjID(cursor,sess))
    
def getInputFiles(cursor,sess,step,location,optional=False):
    [input_files,opt_input_files] = getEntryFromTable(cursor,["input_files","input_files_opt"],"pipesteps",["step"],[step])
    sbjID = getSbjID(cursor,sess)
    input_files = [] if input_files==None else input_files.replace("$sess",sess).replace("$sbj",sbjID).split(',')
    opt_input_files = [] if opt_input_files==None or not optional else opt_input_files.replace("$sess",sess).replace("$sbj",sbjID).split(',')
    input_folder = getInputDir(cursor,sess,step,location)
    
    return [f"{input_folder}/{infile}" for infile in input_files+opt_input_files]
    
def checkInputFiles(cursor,sess,step,location):
    for infile in getInputFiles(cursor,sess,step,location):
        if (not os.path.isfile(infile)) and (not os.path.isdir(infile)):
            printWarning(infile)
            return False
    
    return True

def checkFilesExist(list_files):
    for f in list_files:
        if not os.path.isfile(f) and not os.path.isdir(f) and not os.path.islink(f):
            printWarning(f)
            return False
    return True

# Returns which step is currently running in the cluster for that sess
# no more than one step at a time would be running since a step uses the scratch folder of the sess        
def runningInCluster(cursor,sess):
    return getValueFromTable(cursor,"step","procs",["sess","status","jobID"],[sess,"running","*"])

# Number of procs currently running in the cluster
def numRunningInCluster(cursor):
    return getNumberEntries(cursor,"procs",["status","jobID"],["running","*"])

# Number of threads currently running locally
def numRunningLocally(cursor):
    return getNumberEntries(cursor,"procs",["status","jobID"],["running",""])

def checkReadyProcs(username,password,hostdir,db_name,pipeline,automatic,location,screen=True,emergency=False,priority=[],absolute=[],test_mode=False,force_manual=False):
    print(f"\n### READY PROCS {pipeline} ###")
    
    # Calculate the number of jobs that can be submitted at this time
    # The limit must be a multiple of 3 because I am submitting max limit/3 jobs per step
    now_hour = int(datetime.datetime.today().strftime("%H"))
    if now_hour>=9 and now_hour<=17:
        jobs_limit = 36 if not emergency else 63
        print("\nIt's "+str(now_hour)+"hr (high traffic)")
    elif (now_hour>17 and now_hour<22) or (now_hour>5 and now_hour<9):
        jobs_limit = 39 if not emergency else 66
        print("\nIt's "+str(now_hour)+"hr (medium traffic)")
    else:
        jobs_limit = 45 if not emergency else 69
        print("\nIt's "+str(now_hour)+"hr (low traffic)")
          
    # Calculate the number of threads that can run locally
    # The limit must be a multiple of 3 because I am running max limit/3 threads per step
    threads_limit = 6 if not emergency else 9
    
    # Open a new connection and cursor
    cnx = connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    
    # Get the number of jobs available to run in the cluster
    n_cluster = numRunningInCluster(cursor)
    print(f'## Allowing {jobs_limit} jobs (max {int(jobs_limit/3)} per step)')
    print(f'## Currently running {n_cluster} jobs in the cluster')
    jobs_available = jobs_limit-n_cluster
    print(f'## Number of submissions available: {jobs_available}\n')
          
    # Get the number of threads available to run locally
    print(f'## Allowing {threads_limit} threads (max {int(threads_limit/3)} per step)')
    n_local = numRunningLocally(cursor)
    print(f"## Currently running {n_local} threads locally")
    threads_available = threads_limit-n_local
    print(f"## Number of threads available: {threads_available}\n")
    
    scriptsDir = getScriptsDir(cursor,pipeline,location)
    print(f'## Scripts dir: {scriptsDir}')
    if not os.path.isdir(scriptsDir):
        printError("scriptsDir does not exist")
        return

    # Get the list of steps for the pipeline
    pipesteps = getStepsPipeline(cursor,pipeline)
    
    # Get the list of sessions and steps that can run, and organize by step
    dic_steps = {}
    total_cluster_jobs = 0
    total_local_threads = 0
    for key,line in getEntriesFromTable(cursor,["sess","step"],"procs",["status"],["ready"],["step"]).iterrows():
        sess = line["sess"]
        step = line["step"]
        
        if len(absolute)>0 and (sess not in absolute):
            continue
        
        if step in pipesteps:
            # initialize the array of sessions for that step if it's not initizalized
            array = dic_steps[step] if step in dic_steps.keys() else []
            
            # Do not submit more than the number of available jobs
            # Limit to (jobs_limit/3) sess per step so the other steps don't stay waiting forever in queue
            if getValueFromTable(cursor,"script","pipesteps",["step"],[step])!=None:
                if total_cluster_jobs>=jobs_available or len(array)>=int(jobs_limit/3):
                    continue
                total_cluster_jobs+=1
                
            # Do not run more than the number of available threads
            # Limit to (threads_limit/3) sess per step so the other steps dont stay waiting forever
            elif getValueFromTable(cursor,"local_script","pipesteps",["step"],[step])!=None:
                if total_local_threads>=threads_available or len(array)>=int(threads_limit/3):
                    continue
                total_local_threads+=1
            
            # add session to the step array
            array+=[sess]
            dic_steps[step] = array
            
    # For each step that has subjects to run, submit the job or run locally
    for step,array in dic_steps.items():
        # If there are any subjects that should be prioritized, put them at the beginning of the array
        for sess in priority:
            if sess in array:
                array.remove(sess)
                array.insert(0,sess)
        
        print(f'\n## Attempting to run {len(array)} proc(s) from {step}: {",".join(array)}')
        
        # There is a cluster script
        # Give the option to submit manually if there is a cluster script and we're not in automatic mode
        # Submit automatically if we are in the cluster running automatic mode
        # If we're running locally on automatic mode, do not submit cluster job (cluster scripts are not submitted in automatic mode unless we are in the cluster)
        
        # If there's local and cluster scripts:
        # If it's not in automatic mode, it will try first to submit, then to run locally.
        # If it's in automatic mode, it will only run the cluster job (if the program is running in the cluster)
        # If its in automatic mode, running locally, and with local and cluster scripts, it wont run. Otherwise, it could run locally while also being submitted in the cluster.
        # If we're running locally and not on automatic mode and there's also a local script, give the option to run the local script instead of submitting
        
        [cluster_script,local_script] = getEntryFromTable(cursor,["script","local_script"],"pipesteps",["step"],[step])
        if cluster_script!=None and (location=="hpc" or (location!="hpc" and (not automatic or force_manual))):
            submit_manual = force_manual or input(f"Submit {step} [Y]? ")=='Y' if not automatic else False
            
            # Submit job, either manually or automatic
            if automatic or submit_manual:
                # Copy subjects files to scratch dir
                array_ok = []
                for sess in array:
                    if copyFilesToScratch(cnx,step,sess,location,test_mode=test_mode):
                        array_ok+=[sess]
                        
                # Submit job for sbjs that copied with no error
                if len(array_ok)>0:
                    scriptpath = scriptsDir+'/'+cluster_script
                    submitIndividual(cnx,array_ok,scriptpath,step,location,test_mode=test_mode)
                    
            # We're not in automatic mode and decided to not submit manually and there's also a local script
            elif local_script!=None and input(f"Run {step} local [N]? ")!='N':
                runLocally(cnx,step,array,local_script,location,"ran locally instead of cluster",test_mode=test_mode)
                
        # There is no cluster script, but there is a local script
        # If we're in automatic mode, only run locally (otherwise they could run both places at the same time)
        # If we're not in automatic mode, give the option to run eitherway
        elif cluster_script==None and local_script!=None and ((automatic and location!="hpc") or (not automatic and input(f"Run {step} local [N]? ")!='N')):
            runLocally(cnx,step,array,local_script,location,test_mode=test_mode)
            
        else:
            print("not ran, it's either in automatic mode in the cluster, or user chose not to run.")
            
        if test_mode and input("Stop the test? [Y]: ")=='Y':
            break
    
    # Close cursor and connection
    cursor.close()
    cnx.close()
    print(f"\n### FINISHED CHECKING READY PROCS FOR {pipeline} ###\n")

def copyFilesToScratch(cnx,step,sess,location,test_mode=False):
    cursor = cnx.cursor()
    sbjScratch = getScratchDirPipe(cursor,location,getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step]))+'/'+sess
    prev_steps = getPreviousSteps(cursor,step) # previous steps
    running_step = runningInCluster(cursor,sess) # steps that are running for this session
    infiles_all = getInputFiles(cursor,sess,step,location,optional=True) # list of all input files
    infiles_must = getInputFiles(cursor,sess,step,location) # list of only mandatory input files
    cursor.close()
    
    if len(infiles_all)==0:
        printError("No input files")
        return False
    
    # If the previous steps (given that there's any) has move_files as false, there's nothing to do here
    print(f'prev step: {prev_steps}')
    if len(prev_steps)>0: 
        for prev_step in prev_steps:
            cursor = cnx.cursor()
            mv = mvFiles(cursor,prev_step)
            cursor.close()
            if not mv:
                return True
    
    # Remove sbj scratch dir if exists and no other process is running on it
    if os.path.isdir(sbjScratch):
        if running_step!=None:
            cursor = cnx.cursor()
            updateProcStatus(cursor,sess,step,"hold",location,running_step+" is currently running on this session",True,test_mode=test_mode)
            updateDB(cnx,cursor)
            cursor.close()
            return False
        elif not test_mode:
            os.system("rm -r "+sbjScratch)
        else:
            print("rm -r "+sbjScratch)
            
    # Re-create the folder
    print(f"Scratch: {sbjScratch}")
    if not test_mode:
        os.system(f"mkdir -p {sbjScratch}")
    else:
        print(f"mkdir -p {sbjScratch}")
    
    # Copy input files to sess scratch dir
    print(f'Copying {len(infiles_all)} {sess} files to scratch for {step}...')
    for infile in infiles_all:
        if os.path.exists(infile):
            if not test_mode:
                os.system("cp -r "+infile+' '+sbjScratch)
            else:
                print("cp -r "+infile+' '+sbjScratch)
        elif infile in infiles_must:
            printWarning(infile)
            cursor = cnx.cursor()
            error(cursor,"input file missing",sess,step,location,test_mode=test_mode)
            updateDB(cnx,cursor)
            cursor.close()
            return False
    print("done")
    
    return True

def runLocally(cnx,step,array,local_script,location,notes="",test_mode=False):
    # If the script doesn't exist, it cannot run for any sess in array
    cursor = cnx.cursor()
    scriptsDir = getScriptsDir(cursor,getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step]),location)+"/"
    cursor.close()
    if not os.path.isfile(scriptsDir+'/'+local_script.split(" ")[0]):
        printError(f'{scriptsDir}/{local_script.split(" ")[0]} does not exist')
        return False
    
    remove = []
    cursor = cnx.cursor()
    for sess in array:
        # Check that the input files are present
        if not checkInputFiles(cursor,sess,step,location):
            updateProcStatus(cursor,sess,step,'hold',location,'missing input_files',True,test_mode=test_mode)
            remove+=[sess]
            continue
        
        # Set all sessions as running so they don't run at the same time elsewhere
        updateProcStatus(cursor,sess,step,'running',location,notes,replace=False,test_mode=test_mode)
        updateEntriesFromTable(cursor,"procs",["jobID","exectimesec","qsub_time","cluster"],["","","",""],["sess","step"],[sess,step],test_mode=test_mode)
        updateDB(cnx,cursor)
        
    for sess in [sess for sess in array if sess not in remove]:
        # Create the temporary dir where the script will run and change the current dir to that path
        cmd = f"rm -rf tmp_{step}_{sess} && "\
        f"mkdir tmp_{step}_{sess} && "\
        f"cd tmp_{step}_{sess} && "\
        "nohup "+scriptsDir+local_script.replace("$sess",sess).replace("$sbj",getSbjID(cursor,sess))+" && "\
        f"[ -f nohup.out ] && mv nohup.out {scriptsDir}logs/{step}.{sess} && "\
        f"cd .. && "\
        f"rm -rf tmp_{step}_{sess}"
        print(cmd)
        if not test_mode:
            subprocess.Popen(cmd,shell=True)
        else:
            continue
    
    # Returns true even if there was some error and for some reason nohup was not created
    cursor.close()
    return True

def createIndividualJob(cursor,scriptpath,sess,step,location,test_mode=False):
    print(f'Creating individual job for {step} {sess}...')
    
    if not os.path.isfile(scriptpath):
        printError(f'{scriptpath} does not exist')
        return ""
    
    pipe = getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    proj = getStringFromTable(cursor,"project","pipelines",["pipename"],[pipe])
    out = getScratchDirPipe(cursor,location,pipe)+'/'+step+'_'+sess+".sh"
    if test_mode:
        return out
    
    [sbjID,sess] = getEntryFromTable(cursor,["sbjID","sess"],"sessions",["sess"],[sess])
    fin = open(scriptpath,'r')
    fout = open(out,'w')    
    for line in fin:
        line = line.replace('\n','')
        if line=="sbj=${SUBJECTS[PBS_ARRAYID-1]}":
            fout.write("sbj="+sbjID+'\n')
        elif line=="sbj=sbjID":
            fout.write("sbj="+sbjID+'\n')
        elif line=="sess=sess":
            fout.write("sess="+sess+'\n')
        elif line=="proj=proj":
            fout.write("proj="+proj+'\n')
        elif line=="step=step":
            fout.write("step="+step+'\n')
        elif (not line.startswith("#SBATCH -t")) and ("cat sbj_list.txt" not in line):
            fout.write(line+'\n')
    
    fout.close()
    fin.close()  
    print("done")
    
    return out
    
def submitIndividual(cnx,array_ok,scriptpath,step,location,notes="",test_mode=False):
    cursor = cnx.cursor()
    
    n = 0
    for sess in array_ok:
        # Re-write the script for each subject in the scratch directory
        out = createIndividualJob(cursor,scriptpath,sess,step,location,test_mode=test_mode)
        if not os.path.isfile(out) and not test_mode:
            printError(f"job not submitted for {sess}({step})")
            continue
        info = out.split('/')
        fname = info[-1]
        
        # Submit the job
        # Do not check if submission_scratch exists because if I'm not in the hpc it's not going to exist
        submission_scratch = getScratchDirPipe(cursor,"hpc",getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step]))
        print(f"cd {submission_scratch}")
        print(f"qsub {fname}")
        
        jobID = ""
        if location=="hpc":
            os.chdir(submission_scratch)
            try:
                if not test_mode:
                    jobID = systemOut(f"qsub {fname}")
                    print(jobID)
            except ValueError as err:
                printError(f"exception while submitting {fname}: {err}")
                continue
        else:
            jobID = input("jobID: ")
        if jobID=="" and not test_mode:
            printError(f"job not submitted for {sess}({step})")
            continue
        
        # Get current date and time
        now_date = str(datetime.date.today().strftime("%Y_%m_%d"))
        now_time = str(datetime.datetime.today().strftime("%H:%M:%S"))
        now = now_date+'_'+now_time
        
        # Set the status as running
        updateProcStatus(cursor,sess,step,'running',location,notes,test_mode=test_mode)
        n+=updateEntriesFromTable(cursor,"procs",["jobID","qsub_time","cluster"],[jobID,now,"hpc"],["sess","step"],[sess,step],test_mode=test_mode)
        updateDB(cnx,cursor)
    
    print(f"{n} procs were marked as running")
    cursor.close()

def menu(projopts,showAuto=True):
    print("\nChoose an option:")
    if showAuto:
        print("0. Run all (automatic mode)")
    print("1. Exit\n"
          "2. Check running procs\n"
          "3. Check done procs\n"
          "4. Check ready procs\n"
          "6. Check steps status\n"
          "7. Backup DB\n" 
          "8. Recover tables from backup\n"
          "9. Clean logs\n"
          "10. Check on hold\n"
          "11. Create new step / add sbjs to step\n"
          "12. Delete step\n"
          "13. Print/modify step information")
    #print("5. Progress report")
    
    for key,val in projopts.items():
        print(key+". "+val)
    
    return input(">> ")

def hpcMenu(projopts,showAuto=True):
    print("Choose an option:")
    if showAuto:
        print("0. Run all (automatic mode)")
    print("1. Exit\n"
          "2. Check running procs\n"
          "3. Check done procs\n"
          "4. check ready procs\n"
          "6. Check steps status\n"
          "10. Check on hold")
    
    for key,val in projopts.items():
        print(key+". "+val)
    
    return input(">> ")

def backupDB(username,password,hostdir,database,location="petrov",automatic=False):
    print("\n### BACKUP DB ###")
          
    # Binder is just to create the connection but it backups all the databases
    cnx = connect(username,password,hostdir,database)
    cursor = cnx.cursor()
    
    ### Get the path to the backup directory ###
    dbDir = getValueFromTable(cursor,"db_dir","generalInfo.locations",["location"],[location])
    recovDir = f"{dbDir}/recovery/"
    if not os.path.isdir(recovDir):
        return False
    
    ### Create the backup directory and move the current backup files there ###
    now = str(datetime.date.today()).replace('-','_')
    outdir = recovDir+now
    if automatic and (os.path.isdir(outdir) or location=="hpc"):
        return
    while os.path.isdir(outdir):
        outdir+='+'
    
    os.system(f"mkdir -p {outdir}")
    for txt in glob.glob(recovDir+"*.txt"):
        os.system("mv "+txt+" "+outdir)
    for qsl in glob.glob(recovDir+"*.sql"):
        os.system("mv "+qsl+" "+outdir)
    
    ### Get the list of columns per database.table ###
    tables = getColumnsDB(cursor)
    
    ### Create new backup files ####
    for table,cols in tables.items():
        print("Backing up "+table+"...")
        fout = open(recovDir+table+".txt",'w')
        fout.write("\t".join(cols)+"\n")
        cursor.execute("select "+",".join(cols)+" from "+table)
        for record in cursor.fetchall():
            line = "\\N" if record[0] is None else "\""+str(record[0])+"\""
            for i in range(1,len(record)):
                value = "\\N" if record[i] is None else "\""+str(record[i])+"\""
                line+="\t"+value
            fout.write(line+"\n")
        fout.close()
    
    ### Update recreateTables.sql and reloadTables.sql ###
    fout1 = open(recovDir+"recreateTables.sql",'w')
    fout2 = open(recovDir+"reloadTables.sql",'w')
    for table,cols in tables.items():
        create = ""
        cursor.execute(f"show create table {table}")
        for(table_raw,cmd) in cursor.fetchall():
            create = cmd.replace("\n","").replace("`"+table_raw+"`","`"+table+"`")
        if create!="":
            fout1.write(create.replace("`","")+";\n\n")
            
        load = "load data local infile '"+recovDir+table+".txt' into table "+table+" fields optionally enclosed by '\"' ("+cols[0]
        for i in range(1,len(cols)):
            load+=","+cols[i]
        load+=");"
        fout2.write(load+"\n\n")
    fout2.close()
    fout1.close()
    
    cursor.close()
    cnx.close()

def moveFromScratch(cursor,sess,step,scratchdir,ofile,efile,location,test_mode=False):
    print(f"Moving {sess} files from scratch after {step}...")
    
    ## Remove individual sess/step script
    if not test_mode:
        os.system(f"rm -f {scratchdir}/{step}_{sess}.sh")
    else:
        print(f"rm -f {scratchdir}/{step}_{sess}.sh")
    
    ## Move logs to project log folder
    pipeline = getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    outlogs = getScriptsDir(cursor,pipeline,location)+"/logs"
    if os.path.isdir(outlogs):
        if os.path.isfile(ofile):
            print(f"Moving {ofile} to {outlogs}")
            if not test_mode:
                os.system("mv "+ofile+' '+outlogs)
        else:
            printError(f"output log doesnt exist: {ofile}")
            
        if os.path.isfile(efile):
            print(f"Moving {efile} to {outlogs}")
            if not test_mode:
                os.system(f"mv {efile} {outlogs}")
        else:
            printError(f"error log doesnt exist: {efile}")
    else:
        printError("logs folder doesn't exist, leaving logs in scratch")
        
    ## Remove input files from scratch directory
    scratch_outputs = getScratchOutputFiles(cursor,sess,step,scratchdir)
    for scrin in getScratchInputFiles(cursor,sess,step,scratchdir):
        if scrin not in scratch_outputs:
            if not test_mode:
                os.system(f"rm -rf {scrin}") 
            else:
                print(f"rm -rf {scrin}") 
        
    ## Move output files to subject directory (the reminding files after removing the inputs)
    output_folder = getOutputDir(cursor,sess,step,location)
    if not test_mode:
        os.system("mkdir -p "+output_folder)
    else:
        print("mkdir -p "+output_folder)
        
    sbjscratch = scratchdir+'/'+sess
    for f in glob.glob(sbjscratch+"/*"):
        if not test_mode:
            os.system(f"mv {f} {output_folder}")
        else:
            print(f"mv {f} {output_folder}")
        
    ## For preTopupEddy a few more things need to be done
    if step=="preTopupEddy":
        print("Moving preTopupEddy files...")
        
        # Create output folders
        topupdir = output_folder+"/topup/"
        eddydir = output_folder+"/eddy/"
        preddydir = output_folder+"/preEddy/"
        if not test_mode:
            os.system(f"mkdir -p {topupdir} {eddydir} {preddydir}")
        else:
            print(f"mkdir -p {topupdir} {eddydir} {preddydir}")
        
        # Move files into topup, eddy and preEddy dirs
        dic = {topupdir:["Pos_b0.nii.gz","Neg_b0.nii.gz","Pos_Neg_b0.nii.gz","acqparams.txt"],
               eddydir:["Pos.nii.gz","Pos.bval","Pos.bvec","Neg.nii.gz","Neg.bval","Neg.bvec","Pos_Neg.nii.gz","Pos_Neg.bval","Pos_Neg.bvec","index.txt"],
               preddydir:["Pos_b0_0000.nii.gz","Pos_b0_0001.nii.gz","Neg_b0_0000.nii.gz","Neg_b0_0001.nii.gz"]}
        for target,array in dic.items():
            for f in array:
                ff = output_folder+'/'+f
                if os.path.isfile(ff) and not test_mode:
                    os.system(f"mv {ff} {target}")
                elif os.path.isfile(ff):
                    print(f"mv {ff} {target}")
        
        # Copy files into eddy
        if os.path.isfile(topupdir+"acqparams.txt") and not test_mode:
            os.system(f"ln -s {topupdir}acqparams.txt {eddydir}")
        elif os.path.isfile(topupdir+"acqparams.txt"):
            print(f"ln -s {topupdir}acqparams.txt {eddydir}")
                
    ## Remove subject scratch
    if not test_mode:
        os.system(f"rm -r {sbjscratch}")
    else:
        print(f"rm -r {sbjscratch}")
    print("Done moving files, sbj scratch directory deleted")
    return True

def mergeColsDF(df,cols,rm=True):
    new_col = "_".join(cols)
    df = pd.concat([pd.DataFrame({new_col:[]}),df])
    for index,line in df.iterrows():
        vals = []
        for col in cols:
            vals+=[df.loc[index,col]]
        df.loc[index,new_col] = "_".join(vals)
        
    if rm:
        df.drop(columns=cols,axis=1,inplace=True)
    
    return df
