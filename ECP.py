#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Monica Keith"
__status__ = "Production"

import mylib
import os
import mrilib
import fileslib
import numpy as np
import glob
import statistics

empty_values = ["nan","NaT","None","\\N","NULL","NaN","none","nat",""]

def doneStart(sess_list,username,password,hostdir,db_name,location,absolute=[]):
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        mylib.checked(cursor,sess,"start",False,'',location)
        mylib.updateDB(cnx,cursor)
        
    cursor.close()
    cnx.close()
        
def doneCamino(sess_list,username,password,hostdir,db_name,step,screen,location,automatic,absolute=[],priority=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    soft = "freeview" if location!="petrov" else "fsleyes"
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        # Check output files
        if not mylib.checkOutputFiles(cursor,sess,step,location):
            mylib.error(cursor,"output file missing",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        if automatic or not screen or location=="hpc":
            continue
        
        # Get fa, md and mask from output files
        indir = mylib.getInputDir(cursor,sess,step,location)+"/"
        print("input dir: "+indir)
        mask = indir+"nodif_brain_mask.nii.gz"
        
        outdir = mylib.getOutputDir(cursor,sess,step,location)
        print("output dir: "+outdir)
        fa = outdir+"/wlf/fa.nii.gz"
        os.system("fslmaths "+mask+" -mul "+fa+" "+fa)
        
        opt = mrilib.faMenu(fa,mask,screen,location) if location=="petrov" else mrilib.faMenu(fa,mask,screen,location,soft="freeview")
        
        # Check FA
        if opt=="3":
            mylib.error(cursor,"Bad fa",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        mylib.checked(cursor,sess,step,False,'',location)
        mylib.updateDB(cnx,cursor)
        
        # Check MD
        md = outdir+"/wlf/md.nii.gz"
        os.system("fslmaths "+mask+" -mul "+md+" "+md)
        os.system("fslmaths "+md+" -thr 0 "+md)
            
        mrilib.fsleyes(mrilib.fsleyescmd_brainN([md],software=soft))
            
        if input("Error? [Y]: ")=='Y':
            mylib.updateNotesProc(cursor,sess,step,"Bad md")
            mylib.updateDB(cnx,cursor)
            continue
            
    cursor.close()
    cnx.close()

def doneDTIFIT(sess_list,username,password,hostdir,db_name,step,screen,automatic,location,other_imgs=[],priority=[],absolute=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    software = "freeview" if location!="petrov" else "fsleyes"
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        if not mylib.checkOutputFiles(cursor,sess,step,location):
            mylib.error(cursor,"output file missing",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        if automatic or not screen or location=="hpc":
            continue
        
        # Check FA
        mylib.printWarning("REMEMBER TO ERODE FA IF NECESSARY")
        outDir = mylib.getOutputDir(cursor,sess,step,location)+'/'
        mask = mylib.getInputDir(cursor,sess,step,location)+"/nodif_brain_mask.nii.gz"
        FA = outDir+"dti_FA.nii.gz"
        V1 = outDir+"dti_V1.nii.gz"
        
        opt = mrilib.vecmenu(FA,V1,mask,screen,location,soft=software)
        if opt=="3":
            mylib.error(cursor,'bad fa',sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Check vectors
        mrilib.fsleyes(mrilib.fsleyescmd_brainN([V1],software=software))
                
        if input("Error? [Y]: ")=='Y':
            mylib.error(cursor,'bad vectors',sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
            
        if location=="petrov":
            mrilib.fsleyes(mrilib.fsleyescmd_vecs(V1))
                
            if input("Error? [Y]: ")=='Y':
                mylib.error(cursor,'bad vectors',sess,step,location)
                mylib.updateDB(cnx,cursor)
                continue
        
        # Check other images
        for img in other_imgs:
            nifti = f"{outDir}dti_{img}.nii.gz"
            os.system("fslmaths "+nifti+" -mul "+mask+' '+nifti)
            os.system("fslmaths "+nifti+" -thr 0 "+nifti)
            
            mrilib.fsleyes(mrilib.fsleyescmd_brainN([nifti],software=software))
            
            if input("Error [Y]? ")=='Y':
                mylib.error(cursor,'bad '+img,sess,step,location)
                mylib.updateDB(cnx,cursor)
                continue
        
        mylib.checked(cursor,sess,step,False,'',location)
        mylib.updateDB(cnx,cursor)
        
    cursor.close()
    cnx.close()
    
def diffdirsInfo(cursor,sess,pipe,location):
    df = mylib.getEntriesFromTable(cursor,["75_AP","75_PA","76_AP","76_PA"],"sessions",["sess"],[sess]).swapaxes('columns', 'rows')
    df.loc[:,0] = df.loc[:,0]==1
    sbjDir = mylib.getSessDir(cursor,pipe,location,sess)+"/raw/"
    sbjID = mylib.getSbjID(cursor,sess)
    df.loc[:,1] = [sbjDir+sbjID+"_3T_DWI_dir"+idx for idx in df.index.tolist()]
    df.columns = ["included","prefix"]    
    return df
    
def doneVolReg(sess_list,username,password,hostdir,db_name,step,location,priority=[],absolute=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    pipe = mylib.getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{len(array)})...')
        
        ok = True
        df = diffdirsInfo(cursor,sess,pipe,location)
        for diffdir,line in df[df["included"]].iterrows():
            for ext in [".1D","_enorm.1D"]:
                if not os.path.isfile(line["prefix"]+ext):
                    ok = False
                    mylib.error(cursor,"output file missing for "+diffdir,sess,step,location)
                    mylib.updateDB(cnx,cursor)
                    break
            if not ok:
                break
            
            nvols = mrilib.getDims(line["prefix"]+".nii.gz")[3]
            nlines = fileslib.linesCols(line["prefix"]+"_enorm.1D")[0]
            if nvols!=nlines:
                ok = False
                mylib.error(cursor,f"{diffdir}: {nvols} vols, {nlines} lines",sess,step,location)
                mylib.updateDB(cnx,cursor)
                break
            
            array = [float(num) for num in fileslib.file2array(line["prefix"]+"_enorm.1D")]
            mylib.insertResult(cursor,sess,diffdir+"_avg_enorm_motion",np.mean(array))
            mylib.updateDB(cnx,cursor)
            
        if ok:
            mylib.checked(cursor,sess,step,False,'',location)
            mylib.updateDB(cnx,cursor)
    
    cursor.close()
    cnx.close()
    
def doneCaminoSeries(sess_list,username,password,hostdir,db_name,step,location,priority=[],absolute=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    pipe = mylib.getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{len(array)})...')
        
        ok = True
        df = diffdirsInfo(cursor,sess,pipe,location)
        for diffdir,line in df[df["included"]].iterrows():
            for ext in [".scheme",".Bfloat","_snr.txt"]:
                if not os.path.isfile(line["prefix"]+ext):
                    ok = False
                    mylib.error(cursor,"output file missing for "+diffdir,sess,step,location)
                    mylib.updateDB(cnx,cursor)
                    break
                
            if not ok:
                break
            
            try:
                txt = open(line["prefix"]+"_snr.txt",'r')
                snr = "NA"
                for line in txt:
                    if line.startswith("SNR mult:"):
                        snr = float(line.replace("\n","").split("\t")[-1])
                        break
                txt.close()
                
                if snr=="NA":
                    ok = False
                    mylib.error(cursor,diffdir+" snr file potentially empty",sess,step,location)
                    mylib.updateDB(cnx,cursor)
                    break
                
                mylib.insertResult(cursor,sess,diffdir+"_SNR",snr)
            except ValueError as err:
                ok = False
                mylib.error(cursor,str(err),sess,step,location)
                mylib.updateDB(cnx,cursor)
                break
            
        if ok:
            mylib.checked(cursor,sess,step,False,'',location)
            mylib.updateDB(cnx,cursor)
        
    cursor.close()
    cnx.close()
    
def generate4d(mask,nvols,tr):
    print("Generating 4D file...")
    res = mask.replace("_mask.nii.gz","_4dmask.nii.gz")
    
    cmd="fslmerge -tr "+res
    for i in range(int(nvols)):
        cmd=cmd+' '+mask
    cmd=f"{cmd} {tr}"
    os.system(cmd)
    
    return mrilib.getDims(res)[3]==nvols and mrilib.getPixDims(res)[3]==tr

def done3Dmask(sess_list,username,password,hostdir,db_name,step,screen,automatic,location,priority=[],absolute=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    pipe = mylib.getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    soft = "fsleyes" if location=="petrov" else "freeview"
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        # Check output files
        if not mylib.checkOutputFiles(cursor,sess,step,location):
            mylib.error(cursor,"output file missing",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        if automatic or not screen or location=="hpc":
            continue
        
        ok = True
        df = diffdirsInfo(cursor,sess,pipe,location)
        for diffdir,line in df[df["included"]].iterrows():
            orig = line["prefix"]+".nii.gz"
            mask = line["prefix"]+"_bet_mask.nii.gz"
            if not os.path.isfile(mask):
                mylib.error(cursor,"output file missing for "+diffdir,sess,step,location)
                mylib.updateDB(cnx,cursor)
                ok = False
                break
            
            if mrilib.maskmenu(orig,mask,location,screen=True,opt='7',softwareP=soft,step=step,sess=sess)[0]=='3':
                mylib.error(cursor,'bad mask for '+diffdir,sess,step,location)
                mylib.updateDB(cnx,cursor)
                ok = False
                break
            
            print("Masking file...")
            masked = line["prefix"]+"_masked.nii.gz"
            os.system(f"fslmaths {orig} -mul {mask} {masked}")
            nvols = mrilib.getDims(orig)[3]
            tr = mrilib.getPixDims(orig)[3]
            if not generate4d(mask,nvols,tr):
                mylib.error(cursor,'error with '+diffdir+' 4d mask',sess,step,location)
                mylib.updateDB(cnx,cursor)
                ok = False
                break
            
        if ok:
            mylib.checked(cursor,sess,step,False,'',location)
            mylib.updateDB(cnx,cursor)
                                    
    cursor.close()
    cnx.close()
    
def doneDS(sess_list,username,password,hostdir,db_name,step,location,priority=[],absolute=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    pipe = mylib.getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        # Check output files
        if not mylib.checkOutputFiles(cursor,sess,step,location):
            mylib.error(cursor,"output file missing",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        ok = True
        df = diffdirsInfo(cursor,sess,pipe,location)
        for diffdir,line in df[df["included"]].iterrows():
            outfile = line["prefix"]+"_masked_ds.nii.gz"
            if not os.path.isfile(outfile):
                mylib.error(cursor,"output file missing for "+diffdir,sess,step,location)
                mylib.updateDB(cnx,cursor)
                ok = False
                continue
            
            if mrilib.hasNegVals(outfile):
                mylib.error(cursor,diffdir+' output has negative values',sess,step,location)
                mylib.updateDB(cnx,cursor)
                ok = False
                break
            
        if ok:
            mylib.checked(cursor,sess,step,False,'',location)
            mylib.updateDB(cnx,cursor)
            
    cursor.close()
    cnx.close()

def fixAPPA(df):
    prefixes = [df.loc[idx,"prefix"] for idx in df.index.tolist()]
    rawDir = os.path.dirname(prefixes[0])+"/"
    tmpDir = rawDir+"/tmp/"
    os.system("rm -rf "+tmpDir)
    os.system("mkdir "+tmpDir)
    
    # Move all files to a temporary dir
    for prefix in prefixes:
        for f in glob.glob(prefix+"*"):
            os.system("mv "+f+" "+tmpDir)
            
    # Move back to raw dir chaging the name
    for orig,new in {"75_AP":"75_PA","76_AP":"76_PA","75_PA":"75_AP","76_PA":"76_AP"}.items():
        for f in glob.glob(tmpDir+"*"+orig+"*"):
            os.system("mv "+f+" "+rawDir+os.path.basename(f).replace(orig,new))
            
    os.system("rm -r "+tmpDir)
    
def readBvals(bvals,delim):
    df = fileslib.fileToDataFrame(bvals,False,delim)
    if len(df)==0:
        return []
    
    try:
        return [int(x) for x in df.iloc[0,:].values.tolist() if str(x)!="nan"]
    except:
        return []

def getDelim(bvals,nifti):
    nii_vols = int(mrilib.getDims(nifti)[3])
    
    for delim in [" ","\t"]:
        if len(readBvals(bvals,delim))==nii_vols:
            return delim
        
    return None

def checkVectors(nifti,bvecs,bvals):
    # Check that the 3 files exist
    if not os.path.isfile(nifti) or not os.path.isfile(bvecs) or not os.path.isfile(bvals):
        return False
    
    # Check that there are no double spaces
    fileslib.changeDelim(bvecs,{"  ":" "},True)
    fileslib.changeDelim(bvals,{"  ":" "},True)
    
    # Make sure bvecs has the same delimiter
    fileslib.changeDelim(bvecs,{" ":"\t"},True)
    fileslib.changeDelim(bvals,{" ":"\t"},True)
    
    # Get the correct delimiter from bvals. It should give the same number of columns as volumes in the nifti
    delim = getDelim(bvals,nifti)
    
    # check that the number of lines in the vectors is correct
    return delim!=None and fileslib.linesCols(bvecs,delim)[0]==3 and fileslib.linesCols(bvals,delim)[0]==1

def n_zeros(bvals,delim):
    if not os.path.isfile(bvals):
        return 0
    
    return len([bval for bval in readBvals(bvals,delim) if bval==0])
    
def doneECC(sess_list,username,password,hostdir,db_name,step,screen,automatic,location,priority=[],absolute=[],test_mode=False):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    pipe = mylib.getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    soft = "fsleyes" if location=="petrov" else "freeview"
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        # Check output files
        if not mylib.checkOutputFiles(cursor,sess,step,location):
            mylib.error(cursor,"output file missing",sess,step,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue

        if automatic or not screen or location=="hpc":
            continue
        
        # Check AP_PA
        df = diffdirsInfo(cursor,sess,pipe,location)
        df = df[df["included"]]
        
        idx_75 = [idx for idx in df.index.tolist() if idx.startswith("75")]
        check_75 = [df.loc[idx,"prefix"]+"_masked_ds_ecc.nii.gz" for idx in idx_75]
        if len(check_75)>0:
            input("\nOnly check AP/PA [Enter]")
            if location!="oneDrive":
                mrilib.fsleyes(mrilib.fsleyescmd_brainN(check_75,screen,software=soft))
            else:
                mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[check_75],[]))
                input("[Enter]")
                mrilib.fsleyesOneDrive_delete(step,sess)
                
        if input("Error in AP/PA? [Y]: ")=='Y':
            fixAPPA(df.loc[idx_75])
            if location!="oneDrive":
                mrilib.fsleyes(mrilib.fsleyescmd_brainN(check_75,screen,software=soft))
            else:
                mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[check_75],[]))
                input("[Enter]")
                mrilib.fsleyesOneDrive_delete(step,sess)
                
            if input("Error? [Y]: ")=='Y':
                mylib.error(cursor,'could not fix 75 AP/PA',sess,step,location,test_mode=test_mode)
                mylib.updateDB(cnx,cursor)
                continue
        
        idx_76 = [idx for idx in df.index.tolist() if idx.startswith("76")]
        check_76 = [df.loc[idx,"prefix"]+"_masked_ds_ecc.nii.gz" for idx in idx_76]
        if len(check_76)>0:
            input("\nOnly check AP/PA [Enter]")
            if location!="oneDrive":
                mrilib.fsleyes(mrilib.fsleyescmd_brainN(check_76,screen,software=soft))
            else:
                mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[check_76],[]))
                input("[Enter]")
                mrilib.fsleyesOneDrive_delete(step,sess)
            
        if input("Error? [Y]: ")=='Y':
            fixAPPA(df.loc[idx_76])
            if location!="oneDrive":
                mrilib.fsleyes(mrilib.fsleyescmd_brainN(check_76,screen,software=soft))
            else:
                mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[check_76],[]))
                input("[Enter]")
                mrilib.fsleyesOneDrive_delete(step,sess)
            
            if input("Error? [Y]: ")=='Y':
                mylib.error(cursor,'could not fix 76 AP/PA',sess,step,location,test_mode=test_mode)
                mylib.updateDB(cnx,cursor)
                continue
        
        error = False
        for diffdir,line in df.iterrows():
            # Check output files
            nifti = line["prefix"]+"_masked_ds_ecc.nii.gz"
            bvals = line["prefix"]+".bval"
            bvecs = line["prefix"]+".bvec" 
            if not os.path.isfile(bvecs) or not os.path.isfile(bvals) or not os.path.isfile(nifti):
                mylib.error(cursor,'output file missing for '+diffdir,sess,step,location,test_mode=test_mode)
                mylib.printWarning(f"{bvecs}")
                mylib.printWarning(f"{bvals}")
                mylib.printWarning(f"{nifti}")
                #error = True
                break
            
            # Check that the vectors are correct
            if not checkVectors(nifti,bvecs,bvals):
                mylib.error(cursor,'bad vectors for '+diffdir,sess,step,location,test_mode=test_mode)
                error = True
                break
            
            # Check the number of b0s
            input("\nCount number of b0 images [Enter]")
            if location!="oneDrive":
                mrilib.fsleyes(mrilib.fsleyescmd_brainN([nifti],screen,software=soft))
            else:
                mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[nifti],[]))
                input("[Enter]")
                mrilib.fsleyesOneDrive_delete(step,sess)
                
            b0s = input("# b0s: ")
            while b0s=="":
                b0s = input("# b0s: ")
            b0s = int(b0s)
            
            if b0s!=n_zeros(bvals,getDelim(bvals,nifti)):
                mylib.error(cursor,'wrong number of b0s in '+diffdir,sess,step,location,test_mode=test_mode)
                error = True
                break
            mylib.insertResult(cursor,sess,diffdir+"_b0bvals",b0s)
            
            # Get & save the original number of slices and volumes
            dims = mrilib.getDims(nifti)
            mylib.insertResult(cursor,sess,diffdir+"_nvols",dims[3],test_mode=test_mode)
            
            nslices = dims[2]
            mylib.insertResult(cursor,sess,diffdir+"_nslices",nslices,test_mode=test_mode)
            df.loc[diffdir,"nslices"] = nslices
        
        mylib.updateDB(cnx,cursor)
        if error:
            continue
        
        # The minimum number of slices must be an even number and divisible by 3 (in order to create the slspec file)
        min_slices = min(df["nslices"].astype(int).values.tolist())
        while min_slices>0 and ((min_slices % 2)!=0 or (min_slices % 3)!=0):
            min_slices-=1
        if min_slices==0:
            mylib.error(cursor,'could not fix the number of slices '+diffdir,sess,step,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Make sure all existing series have the required number of slices
        notes = ""
        for diffdir,line in df[df["nslices"].astype(int)>min_slices].iterrows():
            if notes=="":
                notes = f"#slices changed: {diffdir}: {line['nslices']} to {min_slices}"
            else:
                notes+=f",{diffdir}: {line['nslices']} to {min_slices}"
                
            print(f"\nReducing the number of slices for {diffdir} from {line['nslices']} to {min_slices}...")
            
            nifti = line["prefix"]+"_masked_ds_ecc.nii.gz"
            mask = line["prefix"]+"_bet_mask.nii.gz"
                
            orig = nifti.replace(".nii.gz","_orig.nii.gz")
            orig_mask = mask.replace(".nii.gz","_orig.nii.gz")
                
            if not test_mode:
                if not os.path.isfile(orig):
                    os.system(f"cp {nifti} {orig}")
                if not os.path.isfile(orig_mask):
                    os.system(f"cp {mask} {orig_mask}")
                        
                os.system(f"chmod 775 {nifti}")
                os.system(f"chmod 775 {mask}")
                os.system(f"fslroi {nifti} {nifti} 0 -1 0 -1 0 {min_slices}")
                os.system(f"fslroi {mask} {mask} 0 -1 0 -1 0 {min_slices}")
            else:
                if not os.path.isfile(orig):
                    print(f"cp {nifti} {orig}")
                if not os.path.isfile(orig_mask):
                    print(f"cp {mask} {orig_mask}")
                        
                print(f"chmod 775 {nifti}")
                print(f"chmod 775 {mask}")
                print(f"fslroi {nifti} {nifti} 0 -1 0 -1 0 {min_slices}")
                print(f"fslroi {mask} {mask} 0 -1 0 -1 0 {min_slices}")
                
            print("done")
                
            if mrilib.getDims(nifti)[2]!=min_slices and not test_mode:
                mylib.error(cursor,'could not fix the number of slices on '+diffdir,sess,step,location,test_mode=test_mode)
                mylib.updateDB(cnx,cursor)
                error = True
                break
            
        if error:
            continue
        
        # Save the final number of slices & series
        mylib.insertResult(cursor,sess,'final_slices_series',min_slices,test_mode=test_mode)
        mylib.updateDB(cnx,cursor)
        
        mylib.checked(cursor,sess,step,False,notes,location,test_mode=test_mode)
        mylib.updateDB(cnx,cursor)

    cursor.close()
    cnx.close()
    
def save_vec(lines,outfile,delim):
    fout = open(outfile,'w')
    for i in range(len(lines)):
        if i>0:
            fout.write("\n")
        line = lines[i]
        for j in range(len(line)):
            if j>0:
                fout.write(delim)
            fout.write(str(line[j]))
    fout.close()
    
def invert_sign(vec_line):
    return [-x if x!=0.0 else 0.0 for x in vec_line]

def readBvecs(bvecs,delim):
    df = fileslib.fileToDataFrame(bvecs,False,delim)
    line1 = [float(x) for x in df.iloc[0,:].values.tolist() if not np.isnan(float(x))]
    line2 = [float(x) for x in df.iloc[1,:].values.tolist() if not np.isnan(float(x))]
    line3 = [float(x) for x in df.iloc[2,:].values.tolist() if not np.isnan(float(x))]
    return [line1,line2,line3]

# Invert X axis (red)
def bvecs_invX(bvecs,delim):
    oldbvecs = bvecs+"_old"
    os.system("mv "+bvecs+" "+oldbvecs)
    
    [line1,line2,line3] = readBvecs(oldbvecs,delim)
    line1 = invert_sign(line1)
    save_vec([line1,line2,line3],bvecs,delim)
    
# Swap X axis (red) with Z axis (blue)
def bvecs_swapXZ(bvecs,delim):
    oldbvecs = bvecs+"_old"
    os.system("mv "+bvecs+" "+oldbvecs)
    
    [line1,line2,line3] = readBvecs(oldbvecs,delim)
    save_vec([line3,line2,line1],bvecs,delim)

def doneDtifitPreTopUp(sess_list,username,password,hostdir,db_name,step,screen,automatic,location,priority=[],absolute=[],test_mode=False):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    pipe = mylib.getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    soft = "fsleyes" if location=="petrov" else "freeview"
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        # Check output files
        if not mylib.checkOutputFiles(cursor,sess,step,location):
            mylib.error(cursor,"output file missing",sess,step,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue

        if automatic or not screen or location=="hpc":
            continue
        
        ok = True
        new_notes = ""
        df = diffdirsInfo(cursor,sess,pipe,location)
        for diffdir,line in df[df["included"]].iterrows():
            rawDir = os.path.dirname(line["prefix"])+"/"
            FA = rawDir+"dtifit_ecc/dti_"+diffdir+"_FA.nii.gz"
            V1 = rawDir+"dtifit_ecc/dti_"+diffdir+"_V1.nii.gz"
            if not os.path.isfile(FA) or not os.path.isfile(V1):
                mylib.error(cursor,'output file missing for '+diffdir,sess,step,location,replace_notes=False,test_mode=test_mode)
                mylib.updateDB(cnx,cursor)
                ok = False
                break
            
            # Visually inspect the V1 map
            # first inspect vectors because if these are wrong the FA will be wrong
            if location=="petrov":
                mrilib.fsleyes(mrilib.fsleyescmd_colors(V1))
                mrilib.fsleyes(mrilib.fsleyescmd_vecs(V1))
            elif location!="oneDrive":
                input("Choose the option to display as vectors, Ctrl+Lmouse [Enter]")
                mrilib.fsleyes(mrilib.fsleyescmd_brainN([V1],software=soft))
            else:
                mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[V1],[]))
                input("[Enter]")
                mrilib.fsleyesOneDrive_delete(step,sess)
            
            if input("Invert X (if directions are wrong independent of the colors) [Y]? ")=='Y':
                new_notes = ", ".join([new_notes,f"inverted x on {diffdir}"])
                bvecs_invX(line["prefix"]+".bvec","\t")
                
            if input("Swap X and Y (if colors are wrong independent of the directions) [Y]? ")=='Y':
                new_notes = ", ".join([new_notes,f"swapped x,y on {diffdir}"])
                bvecs_swapXZ(line["prefix"]+".bvec","\t")
                
            mask = line["prefix"]+"_bet_mask.nii.gz"
            opt = mrilib.vecmenu(FA,V1,mask,screen,location,soft=soft) if location!="petrov" else mrilib.vecmenu(FA,V1,mask,screen,location)
            if opt=='3':
                new_notes = ", ".join([new_notes,f"bad outputs for {diffdir}"])
                ok = False
            if new_notes.startswith(", "):
                new_notes = new_notes[2:]
        
        if new_notes!="":
            old_notes = mylib.getNotesProc(cursor,sess,step).strip()
            old_notes = old_notes[:-1] if old_notes.endswith(".") else old_notes
            if ("inverted" in old_notes) or ("swapped" in old_notes):
                new_notes = ". ".join([old_notes,"re-run and was not fixed"])
                ok = False
            elif not old_notes in empty_values:
                new_notes = ". ".join([old_notes,new_notes])
        
        if not ok:
            mylib.error(cursor,new_notes,sess,step,location,replace_notes=True,test_mode=test_mode)
            
        elif new_notes!="":
            mylib.updateProcStatus(cursor,sess,step,'ready',location,notes=new_notes,replace=True,test_mode=test_mode)
            
        else:
            mylib.checked(cursor,sess,step,False,"",location,test_mode=test_mode)
            
        mylib.updateDB(cnx,cursor)

    cursor.close()
    cnx.close()
    
def doneIN(sess_list,username,password,hostdir,db_name,step,location,priority=[],absolute=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    pipe = mylib.getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    dic = {"75_AP":"AP_1","75_PA":"PA_1","76_AP":"AP_2","76_PA":"PA_2"}
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        error = False
        df = diffdirsInfo(cursor,sess,pipe,location)
        for diffdir,line in df[df["included"]].iterrows():
            outfiles = os.path.dirname(line["prefix"]).replace("/raw","/preEddy")+"/"
            nifti = outfiles+dic[diffdir]+".nii.gz"
            bvecs = outfiles+dic[diffdir]+".bvec"
            bvals = outfiles+dic[diffdir]+".bval"
            
            # Check output files
            if not os.path.isfile(nifti) or not os.path.isfile(bvecs) or not os.path.isfile(bvals):
                mylib.error(cursor,dic[diffdir]+' output not generated',sess,step,location)
                mylib.updateDB(cnx,cursor)
                error = True
                break
            
            # check bvecs/bvals
            if not checkVectors(nifti,bvecs,bvals):
                mylib.error(cursor,'wrong vectors for '+dic[diffdir],sess,step,location)
                mylib.updateDB(cnx,cursor)
                error = True
                break
            
            # Check that the output has no negative values
            if mrilib.hasNegVals(nifti):
                mylib.error(cursor,dic[diffdir]+' has negative values',sess,step,location)
                mylib.updateDB(cnx,cursor)
                error = True
                break
        
        if not error:
            mylib.checked(cursor,sess,step,False,'',location)
            mylib.updateDB(cnx,cursor)
        
    cursor.close()
    cnx.close()
    
def checkVol(vol,expected_nvols):
    # Check that the number of expected volumes is correct
    nvols = mrilib.getDims(vol+".nii.gz")[3]
    if expected_nvols!=int(nvols):
        mylib.printError(f'{expected_nvols} vs {nvols}')
        return False
    
    # Check that the vectors are correct
    return checkVectors(vol+".nii.gz",vol+".bvec",vol+".bval")
    
def checkIndexFile(nvols_series,index):
    series = 1
    # Check that the number of lines per series corresponds to the number of volumes in that series
    for n_vols in nvols_series.values():
        if n_vols==0:
            continue
        if n_vols!=fileslib.nlinesval(index,series):
            return False
        series+=1
     
    # Check that the number of unique values in the index file corresponds with the number of series present
    n_uniq_lines = mylib.systemOut("cat "+index+" | sort -nu | wc -l")
    return True if series-1==int(n_uniq_lines) else False
 
def donePreTopupEddy(sess_list,username,password,hostdir,db_name,step,location,priority=[],absolute=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    pipe = mylib.getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{len(array)})...')
        
        # Check output files
        if not mylib.checkOutputFiles(cursor,sess,step,location):
           mylib.error(cursor,"output file missing",sess,step,location)
           mylib.updateDB(cnx,cursor)
           continue
        
        # Get number of volumes per each series
        nvols_series = {}
        total = 0
        for fin in ["PA_1","PA_2","AP_1","AP_2"]:
            fin_nvols = mrilib.getDims(mylib.getSessDir(cursor,pipe,location,sess)+"/preEddy/"+fin+".nii.gz")[3]
            if fin_nvols>0:
                nvols_series[fin] = fin_nvols
                total+=fin_nvols
        
        # Check the index file
        eddyDir = mylib.getSessDir(cursor,pipe,location,sess)+"/eddy/"
        if not checkIndexFile(nvols_series,eddyDir+"index.txt"):
            mylib.error(cursor,"Bad index file",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Check acquisition paremeters file
        n_series = int(mylib.getStringFromTable(cursor,"75_AP","sessions",["sess"],[sess])=="1")+int(mylib.getStringFromTable(cursor,"75_PA","sessions",["sess"],[sess])=="1")+int(mylib.getStringFromTable(cursor,"76_AP","sessions",["sess"],[sess])=="1")+int(mylib.getStringFromTable(cursor,"76_PA","sessions",["sess"],[sess])=="1")
        n_lines = fileslib.linesCols(eddyDir+"acqparams.txt",' ')[0]
        if n_lines!=n_series:
            mylib.error(cursor,f"Bad acqparam file: {n_lines} lines, {n_series} series",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        total_pvols = 0
        total_nvols = 0
        for i in ["1","2"]:
            if "PA_"+i in nvols_series.keys():
                total_pvols+=nvols_series["PA_"+i]
            if "AP_"+i in nvols_series.keys():
                total_nvols+=nvols_series["AP_"+i]
                
        # Check that the positive is correct
        if not checkVol(eddyDir+"Pos",total_pvols):
            mylib.error(cursor,"positive output series incorrect",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Check that the negative is correct
        if not checkVol(eddyDir+"Neg",total_nvols):
            mylib.error(cursor,"negative output series incorrect",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Check that the pos_neg file is correct
        if not checkVol(eddyDir+"Pos_Neg",total):
            mylib.error(cursor,"combined output series incorrect",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        mylib.checked(cursor,sess,step,False,'',location)
        mylib.updateDB(cnx,cursor)
          
    cursor.close()
    cnx.close()
    
# Create the slspec file for eddy
def createSpectFile(cursor,sess,file_path):
    n1 = mylib.getValueFromTable(cursor,"value","results",["result","sess"],['final_slices_series',sess])
    if None:
        return
    
    i2 = int(int(n1)/3)
    i3 = i2*2
    col1 = list(range(0,i2))
    col2 = list(range(i2,i3))
    col3 = list(range(i3,int(n1)))
    
    specFile = open(file_path,'w')
    for i in col1:
        specFile.write(str(col1[i])+' '+str(col2[i])+' '+str(col3[i])+'\n')
    specFile.close()

def doneTopUp(sess_list,username,password,hostdir,db_name,step,screen,automatic,location,priority=[],absolute=[],test_mode=False):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    pipe = mylib.getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        # Check output files
        if not mylib.checkOutputFiles(cursor,sess,step,location):
           mylib.error(cursor,"output file missing",sess,step,location,test_mode=test_mode)
           mylib.updateDB(cnx,cursor)
           continue
       
        # Check movpar file (#lines==#series)
        nseries = sum(mylib.getEntryFromTable(cursor,["75_AP","75_PA","76_AP","76_PA"],"sessions",["sess"],[sess]))
        sbjDir = mylib.getSessDir(cursor,pipe,location,sess)+"/"
        nlines = fileslib.linesCols(sbjDir+"topup/topup_Pos_Neg_b0_movpar.txt",' ')[0]
        if nlines!=nseries:
            mylib.error(cursor,'error in movpar output: '+str(nlines)+" lines, "+str(nseries)+" series",sess,step,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue
       
        # Check img (#vols==total volumes)
        imgvols = mrilib.getDims(sbjDir+"eddy/Pos_Neg.nii.gz")[3]
        n1 = int(mylib.systemOut("cat "+sbjDir+"eddy/index.txt | grep 1 | wc -l"))
        n2 = int(mylib.systemOut("cat "+sbjDir+"eddy/index.txt | grep 2 | wc -l"))
        n3 = int(mylib.systemOut("cat "+sbjDir+"eddy/index.txt | grep 3 | wc -l"))
        n4 = int(mylib.systemOut("cat "+sbjDir+"eddy/index.txt | grep 4 | wc -l"))
        nvols = n1+n2+n3+n4
        if imgvols!=nvols:
            mylib.error(cursor,'wrong # of vols in Pos_Neg',sess,step,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Check bvecs (#cols=total volumes, #lines==3)
        [nlines1, ncols1] = fileslib.linesCols(sbjDir+"eddy/Pos_Neg.bvec",' ')
        ncols2 = fileslib.linesCols(sbjDir+"eddy/Pos_Neg.bvec",'\t')[1]
        if nlines1!=3 or (ncols1!=nvols and ncols2!=nvols):
            mylib.error(cursor,f"bad Pos_Neg bvec nlines:{nlines},ncols:{ncols1}/{ncols2},nvols:{nvols}",sess,step,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Check bvals (#cols=total volumes, #lines==1)
        [nlines1, ncols1] = fileslib.linesCols(sbjDir+"eddy/Pos_Neg.bval",' ')
        ncols2 = fileslib.linesCols(sbjDir+"eddy/Pos_Neg.bval",'\t')[1]
        if nlines1!=1 or (ncols1!=nvols and ncols2!=nvols):
            mylib.error(cursor,f"bad Pos_Neg bval nlines:{nlines},ncols:{ncols1}/{ncols2},nvols:{nvols}",sess,step,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue
            
        if automatic or not screen or location=="hpc":
            continue
                        
        # Check mask and generate spec file
        if mrilib.maskmenu(sbjDir+"topup/hifib0.nii.gz",sbjDir+"topup/nodif_brain_mask.nii.gz",location)[0]=='3':
            mylib.error(cursor,'nodif_brain_mask error',sess,step,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue
        
        createSpectFile(cursor,sess,sbjDir+"eddy/slspecFile.txt")
        if not os.path.isfile(sbjDir+"eddy/slspecFile.txt"):
            mylib.error(cursor,'could not generate spec file',sess,step,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue
        
        mylib.checked(cursor,sess,step,False,'',location,test_mode=test_mode)
        mylib.updateDB(cnx,cursor)
    
    cursor.close()
    cnx.close()
    
def doneEddy(sess_list,username,password,hostdir,db_name,step,screen,automatic,location,priority=[],absolute=[],test_mode=False):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    software = "freeview" if location!="petrov" else "fsleyes"
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        # Check output files
        if not mylib.checkOutputFiles(cursor,sess,step,location):
           mylib.error(cursor,"output file missing",sess,step,location,test_mode=test_mode)
           mylib.updateDB(cnx,cursor)
           continue
        
        [AP_75,PA_75,AP_76,PA_76] = mylib.getEntryFromTable(cursor,["75_AP","75_PA","76_AP","76_PA"],"sessions",["sess"],[sess])
        dim2_exp = 0
        if AP_75==1:
            nvols_AP75 = mylib.getValueFromTable(cursor,"value","results",["sess","result"],[sess,'75_AP_nvols'])
            if nvols_AP75!=None:
                dim2_exp+=int(nvols_AP75)
            else:
                mylib.error(cursor,'cannot check output because 75_AP_nvols result is missing in the DB',sess,step,location,test_mode=test_mode)
                mylib.updateDB(cnx,cursor)
                continue
        if PA_75==1:
            nvols_PA75 = mylib.getValueFromTable(cursor,"value","results",["sess","result"],[sess,'75_PA_nvols'])
            if nvols_PA75!=None:
                dim2_exp+=int(nvols_PA75)
            else:
                mylib.error(cursor,'cannot check output because 75_PA_nvols result is missing in the DB',sess,step,location,test_mode=test_mode)
                mylib.updateDB(cnx,cursor)
                continue
        if AP_76==1:
            nvols_AP76 = mylib.getValueFromTable(cursor,"value","results",["sess","result"],[sess,'76_AP_nvols'])
            if nvols_AP76!=None:
                dim2_exp+=int(nvols_AP76)
            else:
                mylib.error(cursor,'cannot check output because 76_AP_nvols result is missing in the DB',sess,step,location,test_mode=test_mode)
                mylib.updateDB(cnx,cursor)
                continue
        if PA_76==1:
            nvols_PA76 = mylib.getValueFromTable(cursor,"value","results",["sess","result"],[sess,'76_PA_nvols'])
            if nvols_PA76!=None:
                dim2_exp+=int(nvols_PA76)
            else:
                mylib.error(cursor,'cannot check output because 76_PA_nvols result is missing in the DB',sess,step,location,test_mode=test_mode)
                mylib.updateDB(cnx,cursor)
                continue
        
        sbjDir = mylib.getOutputDir(cursor,sess,step,location)+'/'    
        nifti_final = sbjDir+"eddy_unwarped_images.nii.gz"
        dim1 = mrilib.getDims(nifti_final)[3]
        nseries = AP_75+PA_75+AP_76+PA_76
        dim1_exp = dim2_exp if (nseries==1 or nseries==3) else dim2_exp/2    
        
        nifti_outlierfree = sbjDir+"eddy_unwarped_images.eddy_outlier_free_data.nii.gz"
        dim2 = mrilib.getDims(nifti_outlierfree)[3]
        if dim1!=dim1_exp or dim2!=dim2_exp:
            print(f'dim2: {dim2}')
            print(f'dim2_exp: {dim2_exp}')
            
            print(f'dim1: {dim1}')
            print(f'dim1_exp: {dim1_exp}')
            
            mylib.error(cursor,'bad output dimension',sess,step,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue
            
        if automatic or not screen or location=="hpc":
            continue
        
        if location!="oneDrive":
            mrilib.fsleyes(mrilib.fsleyescmd_brainN([nifti_final],software=software))
        else:
            mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[nifti_final],[]))
            input("[Enter]")
            mrilib.fsleyesOneDrive_delete(step,sess)
        
        if input("Error? [Y]: ")=='Y':
            mylib.error(cursor,'bad output',sess,step,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Mean RMS movement respect to first volume
        rms = fileslib.fileToDataFrame(sbjDir+"eddy_unwarped_images.eddy_movement_rms",False,' ').iloc[:,0].values.tolist()[1:]
        M = str(statistics.mean(rms))
        mylib.insertResult(cursor,sess,'mean_rms_mov',M,test_mode=test_mode)
        
        # Total percentage of outliers
        outliers = sbjDir+"eddy_unwarped_images.eddy_outlier_report"
        noutliers = fileslib.linesCols(outliers)[0]
        dims = mrilib.getDims(nifti_final)
        total = dims[2]*dims[3]
        perc_outliers = 0 if noutliers==0 else noutliers*100/total
        mylib.insertResult(cursor,sess,'perc_outliers',perc_outliers,test_mode=test_mode)
        
        mylib.checked(cursor,sess,step,False,'',location,test_mode=test_mode)
        
        # Save bad outliers info
        dic = {}
        fin = open(outliers,'r')
        for line in fin:
            array = line.split(' ')
            scan = array[4]
            if scan not in dic.keys():
                dic[scan] = 1
            else:
                dic[scan]+=1
        fin.close()
        
        final_slices = mylib.getValueFromTable(cursor,"value","results",["sess","result"],[sess,'final_slices_series'])
        if final_slices==None:
            continue
        max_ok = 10*int(final_slices)/100
        bad_vols = []
        for scan,slc in dic.items():
            if slc>max_ok:
                bad_vols+=[int(scan)]
        if len(bad_vols)>0:
            bad_vols.sort()
            perc_vols = len(bad_vols)*100/dim1
            val = str(int(perc_vols))+"% of vols had more than 10% outlier slices: "+str(bad_vols[0])
            for i in range(1,len(bad_vols)):
                val = val+","+str(bad_vols[i])
            mylib.insertResult(cursor,sess,'eddy_outliers',val,test_mode=test_mode)
        
        mylib.updateDB(cnx,cursor)
        
    cursor.close()
    cnx.close()
    
def donePostEddy(sess_list,username,password,hostdir,db_name,step,screen,automatic,location,priority=[],absolute=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    nn = 0
    pipe = mylib.getStringFromTable(cursor,"pipeline","pipesteps",["step"],[step])
    soft = "freeview" if location!="petrov" else "fsleyes"
    for sess in array:
        nn+=1
        print(f'\nChecking {step} on {sess} ({nn}/{n_max})...')
        
        outdir = mylib.getOutputDir(cursor,sess,step,location)+'/'
        preData = mylib.getSessDir(cursor,pipe,location,sess)+"/eddy/Pos_Neg.nii.gz"
        data = outdir+"data.nii.gz"
        mask = outdir+"nodif_brain_mask.nii.gz"
        bvecs = outdir+"bvecs"
        bvals = outdir+"bvals"
        
        if not mylib.checkOutputFiles(cursor,sess,step,location):
            mylib.error(cursor,'One or more output files missing',sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
            
        if not checkVectors(data,bvecs,bvals):
            mylib.error(cursor,'Wrong number of cols or lines in vectors',sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        df = diffdirsInfo(cursor,sess,pipe,location)
        n = len(df[df["included"]])
        nvols = mrilib.getDims(data)[3]
        rawvols = mrilib.getDims(preData)[3]/2
        if not (n!=4 or (nvols==rawvols and n==4)):
            mylib.error(cursor,'Wrong number of volumes in data',sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
       
        if automatic or not screen or location=="hpc":
            continue
                    
        if location!="oneDrive":
            mrilib.fsleyes(mrilib.fsleyescmd_maskN(data,[mask],software=soft))
        else:
            mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[data],[mask]))
            input("[Enter]")
            mrilib.fsleyesOneDrive_delete(step,sess)
        
        if input("Error? [Y]: ")=='Y':
            mylib.error(cursor,'Bad output',sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
            
        mylib.insertResult(cursor,sess,"final_slices",mrilib.getDims(data)[2])
        mylib.insertResult(cursor,sess,"final_volumes",mrilib.getDims(data)[3])
        mylib.checked(cursor,sess,step,False,'',location)          
        mylib.updateDB(cnx,cursor)

    cursor.close()
    cnx.close()
    
def doneTransfSklt(sess_list,username,password,hostdir,db_name,step,screen,automatic,location,dtifit_step,priority=[],absolute=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    soft = "freeview" if location!="petrov" else "fsleyes"
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        # Check output files
        if not mylib.checkOutputFiles(cursor,sess,step,location):
            mylib.error(cursor,"output file missing",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Check the output dimensions
        fa = mylib.getOutputDir(cursor,sess,dtifit_step,location)+"/wlf/fa.nii.gz"
        sklt = mylib.getOutputDir(cursor,sess,step,location)+"/FMRIB58_FA-skeleton_1mm_bin.nii.gz"
        if not mrilib.check3Ddimentions(fa,[sklt]):
            mylib.error(cursor,"bad output dimensions",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        if automatic or not screen or location=="hpc":
            continue
        
        if location!="oneDrive":
            mrilib.fsleyes(mrilib.fsleyescmd_maskN(fa,[sklt],software=soft))
        else:
            mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[fa],[sklt]))
            input("[Enter]")
            mrilib.fsleyesOneDrive_delete(step,sess)
        
        if input("Error? [Y]: ")=='Y':
            mylib.error(cursor,"Bad registration",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Get the mean FA
        os.system("fslmaths "+sklt+" -mul "+fa+" "+sess+"_tmp.nii.gz")
        M = mrilib.getMean(sess+"_tmp.nii.gz")
        os.system("rm "+sess+"_tmp.nii.gz")
        mylib.insertResult(cursor,sess,"mean_fa_sklt",M)
        
        mylib.checked(cursor,sess,step,False,'',location)
        mylib.updateDB(cnx,cursor)
            
    cursor.close()
    cnx.close()
    
def doneWithPDForPNG(sess_list,username,password,hostdir,db_name,screen,antstep,automatic,location,priority=[],absolute=[],test_mode=False):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,antstep)
    print(f'Checking {antstep} on: {",".join(array)}')
    
    n = 0
    for sess in array:
        n+=1
        print(f'\nChecking {antstep} on {sess} ({n}/{n_max})...')
        
        # Check output files
        if not mylib.checkOutputFiles(cursor,sess,antstep,location):
            mylib.error(cursor,"output file missing",sess,antstep,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue
        
        outscheck = []
        for outfile in mylib.getOutputFiles(cursor,sess,antstep,location,optional=True):
            if outfile.endswith(".pdf") or outfile.endswith(".png"):
                outscheck+=[outfile]
            
        if len(outscheck)>0 and (automatic or not screen):
            continue
        
        error = False
        for outfile in outscheck:
            print("gio open "+outfile)
            os.system("gio open "+outfile)
            if input("Error? [Y]: ")=='Y':
                error = True
                break
                
        if error:
            mylib.error(cursor,"Bad registration",sess,antstep,location,test_mode=test_mode)
            mylib.updateDB(cnx,cursor)
            continue
            
        mylib.checked(cursor,sess,antstep,False,'',location,test_mode=test_mode)
        mylib.updateDB(cnx,cursor)
    
    cursor.close()
    cnx.close()
    
def doneBedpostEnd(sess_list,username,password,hostdir,db_name,step,screen,automatic,location,priority=[],absolute=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    soft = "freeview" if location!="petrov" else "fsleyes"
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        # Check output files
        out = mylib.getOutputDir(cursor,sess,step,location)+"/"
        outputs = [out+"bvals",out+"bvecs",out+"commands.txt",out+"logs",out+"dyads2_thr0.05_modf2.nii.gz",out+"dyads2_thr0.05.nii.gz",out+"dyads3_thr0.05_modf3.nii.gz",out+"dyads3_thr0.05.nii.gz",out+"mean_dsamples.nii.gz",out+"mean_d_stdsamples.nii.gz",out+"mean_fsumsamples.nii.gz",out+"mean_S0samples.nii.gz",out+"monitor",out+"nodif_brain_mask.nii.gz",out+"xfms"]
        for i in range(1,4):
            outputs+=[out+"dyads"+str(i)+"_dispersion.nii.gz"]
            outputs+=[out+"dyads"+str(i)+".nii.gz"]
            outputs+=[out+"mean_f"+str(i)+"samples.nii.gz"]
            outputs+=[out+"mean_ph"+str(i)+"samples.nii.gz"]
            outputs+=[out+"mean_th"+str(i)+"samples.nii.gz"]
            outputs+=[out+"merged_f"+str(i)+"samples.nii.gz"]
            outputs+=[out+"merged_ph"+str(i)+"samples.nii.gz"]
            outputs+=[out+"merged_th"+str(i)+"samples.nii.gz"]
        OK = True
        for output in outputs:
            if not os.path.isfile(output) and not os.path.isdir(output):
                OK = False
                mylib.printWarning(output)
                break
        if not OK:
            mylib.error(cursor,"output file missing",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Check size of outputs
        if not mrilib.check3Ddimentions(out+"nodif_brain_mask.nii.gz",[output for output in outputs if output.endswith(".nii.gz")]):
            mylib.error(cursor,"output files have wrong dimentions",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        if automatic or not screen or location=="hpc":
            continue
        
        # Check dyads
        if location!="oneDrive":
            mrilib.fsleyes(mrilib.fsleyescmd_brainN([out+"dyads1.nii.gz"],software=soft))
        else:
            mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[out+"dyads1.nii.gz"],[]))
            input("[Enter]")
            mrilib.fsleyesOneDrive_delete(step,sess)
        
        if input("Error? [Y]: ")=='Y':
            mylib.error(cursor,"bad output",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        mylib.checked(cursor,sess,step,False,'',location)        
        mylib.updateDB(cnx,cursor)
        
    cursor.close()
    cnx.close()

def doneCaminoSNR(sess_list,username,password,hostdir,db_name,step,location,priority=[],absolute=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        if not mylib.checkOutputFiles(cursor,sess,step,location):
            mylib.error(cursor,"One or more output files missing",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        outdir = mylib.getOutputDir(cursor,sess,step,location)+'/'
        SNR = mylib.systemOut("cat "+outdir+"estimatesnr.txt | grep \"SNR mult\" | awk '{print $3}'")
        sigma = mylib.systemOut("cat "+outdir+"estimatesnr.txt | grep \"sigma mult\" | awk '{print $3}'")
        b0s = mylib.systemOut("cat "+outdir+"estimatesnr.txt | grep \"Number of b=0 images\" | awk '{print $5}'")
        if SNR=="" or sigma=="" or b0s=="":
            mylib.error(cursor,"empty estimatesnr file",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Save results in DB
        mylib.insertResult(cursor,sess,"SNR_new",SNR) 
        mylib.insertResult(cursor,sess,"sigma_new",sigma) 
        mylib.insertResult(cursor,sess,"final_b0s",b0s) 
        mylib.checked(cursor,sess,step,False,'',location)
        mylib.updateDB(cnx,cursor)
        
        # If SNR<10 it should not run bedpost
        if float(SNR)<10:
            mylib.updateNotesProc(cursor,sess,step,'SNR less than 10')
            next_step = mylib.getStringFromTable(cursor,"next_step","pipesteps",["step"],[step])
            mylib.deleteEntryFromTable(cursor,"procs",["sess","step"],[sess,next_step])
            mylib.updateDB(cnx,cursor)
            continue
    
    cursor.close()
    cnx.close()
    
def doneProbtrackX(sess_list,username,password,hostdir,db_name,screen,automatic,location,step,ants_step,brain_name,forcepdf=False,pdf_outdir="",priority=[],absolute=[]):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    soft = "freeview" if location!="petrov" else "fsleyes"
    for sess in array:
        n+=1
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        # Check output files
        if not mylib.checkOutputFiles(cursor,sess,step,location):
            mylib.error(cursor,"output file missing",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Check that waytotal is greater than zero
        out = mylib.getOutputDir(cursor,sess,step,location)+"/"
        fin = open(out+"waytotal",'r')
        waytotal = fin.readline()
        fin.close()
        if float(waytotal)==0:
            mylib.error(cursor,"waytotal zero",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        # Check that fdt_paths is not empty
        paths = out+"fdt_paths_norm.nii.gz"
        max_val = mrilib.getRange(paths)[1]
        if float(max_val)==0:
            os.system("rm "+paths)
            mylib.error(cursor,"Empty fdt_paths",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        if automatic or not screen or location=="hpc" or (location=="oneDrive" and forcepdf):
            continue            
        
        # Get the list of masks
        masks_tracto = []
        for infile in mylib.getInputFiles(cursor,sess,step,location):
            if infile.endswith(".nii.gz"):
                masks_tracto+=[infile]
        
        # Threshold the paths
        thr_paths = out+"fdt_paths_norm_thr001.nii.gz"
        if not os.path.isfile(thr_paths):
            os.system("fslmaths "+paths+" -thr 0.001 "+thr_paths)
        
        # Get path to the underlay brain
        brain = mylib.getOutputDir(cursor,sess,ants_step,location)+"/"+brain_name
        
        # Check paths
        if not forcepdf:
            input("paths are already thresholded. Remember to visualize in different views if necessary. [Enter]")
            if location!="oneDrive":
                mrilib.fsleyes(mrilib.fsleyescmd_maskN(brain,masks_tracto+[thr_paths],force_pdf=forcepdf,software=soft))
            else:
                mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[brain],masks_tracto+[thr_paths]))
                input("[Enter]")
                mrilib.fsleyesOneDrive_delete(step,sess)
            
            if input("Error? [Y]: ")=='Y':
                mylib.error(cursor,"bad output",sess,step,location)
                mylib.updateDB(cnx,cursor)
                continue
        else:
            if pdf_outdir=="":
                pdf_outdir = mylib.getOutputDir(cursor,sess,step,location)
            
            input("Axial view ('normal') [Enter]")
            output = f"{pdf_outdir}/axial_{step}_{sess}.pdf"
            mrilib.fsleyesPDF_probtrackX(brain,masks_tracto,thr_paths,output)
            input("[Enter]")
            
            input("Sagittal view (profile) [Enter]")
            output = f"{pdf_outdir}/sagittal_{step}_{sess}.pdf"
            mrilib.fsleyesPDF_probtrackX_2(brain,masks_tracto,output,thr_paths)
            input("[Enter]")
            
            input("Coronal view (LGN drawing view) [Enter]")
            output = f"{pdf_outdir}/coronal_{step}_{sess}.pdf"
            mrilib.fsleyesPDF_probtrackX_3(brain,masks_tracto,output,thr_paths)
            input("[Enter]")
            
            if input("Error? [Y]: ")=='Y':
                mylib.error(cursor,"bad output",sess,step,location)
                mylib.updateDB(cnx,cursor)
                continue
        
        mylib.checked(cursor,sess,step,False,'',location)        
        mylib.updateDB(cnx,cursor)
        
    cursor.close()
    cnx.close()
    
# paths_std: file name (WITHOUT extension) of the output file for the paths in std space

# masks_std_step =step that has masks_std as output files
# masks_std: file names (WITHOUT extension) of the output mask files in std space (output file of masks_std_step)

# ants_step: i.e. ants_wlf. step where ants_brain is as an output file (in std space)
# ants_brain: i.e. fa (WITHOUT extension) in std space (output file of ants_step)
    
# paths_diff: file name (WITHOUT extension) of the output file for the paths in diff space

# masks_diff_step = step that has masks_diff as output files
# masks_diff: file names (WITHOUT extension) of the output mask files in diff space (output file of masks_diffstep)

# diff_step: i.e. postEddy step where diff_brain is an output file (in diff space)
# diff_brain: i.e. data_b0 (WITHOUT extension) in diff space (output file of diff_step)
def doneProbmap(sess_list,username,password,hostdir,db_name,step,screen,automatic,location,priority=[],absolute=[],forcepdf=False,masks_std_step="",masks_std=[],paths_std="",ants_step='',ants_brain='',paths_diff="",masks_diff_step="",masks_diff=[],diff_step="",diff_brain=""):
    array = []
    for sess in sess_list:
        if len(absolute)>0 and (sess not in absolute):
            continue
        array+=[sess]
    if len(array)==0:
        return
        
    # Put the priority sessions at the beginning of the array
    for sess in priority:
        if sess in array:
            array.remove(sess)
            array.insert(0,sess)
    
    n_max = len(array) if len(array)<5 else 5
    array = array[:n_max]
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')
    
    n = 0
    soft = "freeview" if location!="petrov" else "fsleyes"
    for sess in array:
        n+=1
        if n>n_max:
            break
        print(f'\nChecking {step} on {sess} ({n}/{n_max})...')
        
        # Check output files
        if not mylib.checkOutputFiles(cursor,sess,step,location):
            mylib.error(cursor,"output file missing",sess,step,location)
            mylib.updateDB(cnx,cursor)
            continue
        
        if not screen or automatic or location=="hpc":
            continue
        
        # Threhold paths to 1% for visualization
        stdPaths = ""
        diffPaths = ""
        for fout in mylib.getOutputFiles(cursor,sess,step,location):
            if fout.endswith(".nii.gz"):
                fout_name = os.path.basename(fout).replace(".nii.gz","")
                if fout_name==paths_std:
                    stdPaths = fout
                elif fout_name==paths_diff:
                    diffPaths = fout
                
        if stdPaths=="" and diffPaths=="":
            mylib.printError("could not find paths in output files")
            continue
        
        print("thresholding paths...")
        if stdPaths!="":
            print(f"STD paths: {stdPaths}")
            os.system("fslmaths "+stdPaths+" -thr 0.01 "+stdPaths.replace(".nii.gz","_thr.nii.gz"))
            stdPaths = stdPaths.replace(".nii.gz","_thr.nii.gz")
        if diffPaths!="":
            print(f"DIFF paths: {diffPaths}")
            os.system("fslmaths "+diffPaths+" -thr 0.01 "+diffPaths.replace(".nii.gz","_thr.nii.gz"))
            diffPaths = diffPaths.replace(".nii.gz","_thr.nii.gz")
        print("done")
            
        # Check paths in FSL standard space space
        if len(masks_std)>0 and ants_step!="" and ants_brain!="":
            print("\nPaths in STANDARD SPACE:")
            
            brain = ""
            for fout in mylib.getOutputFiles(cursor,sess,ants_step,location):
                if fout.endswith(".nii.gz") and os.path.basename(fout).replace(".nii.gz","")==ants_brain:
                    brain = fout
                    break
            if brain=="":
                mylib.printError(f"could not find {ants_brain} in {ants_step} outputs")
                continue
            print(f"underlay brain: {brain}")
            
            stdMasks = []
            for fout in mylib.getOutputFiles(cursor,sess,masks_std_step,location):
                if fout.endswith(".nii.gz"):
                    fout_name = os.path.basename(fout).replace(".nii.gz","")
                    if fout_name in masks_std:
                        stdMasks+=[fout]
            print(f"masks: {','.join(stdMasks)}")
            
            input("paths are already thresholded [Enter]")
            if location!="oneDrive":
                mrilib.fsleyes(mrilib.fsleyescmd_maskN(brain,stdMasks+[stdPaths],force_pdf=forcepdf,software=soft),delete=False)
            else:
                mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[brain],stdMasks+[stdPaths]))
                input("[Enter]")
                mrilib.fsleyesOneDrive_delete(step,sess)
                
            if input("Error? [Y]: ")=='Y':
                mylib.error(cursor,"bad output",sess,step,location)
                mylib.updateDB(cnx,cursor)
                continue
        
        # Check paths in diffusion space
        if len(masks_diff)>0 and diff_step!="" and diff_brain!="":
            print("\nPaths in DIFFUSION SPACE:")
            
            brain = ""
            for fout in mylib.getOutputFiles(cursor,sess,diff_step,location):
                if fout.endswith(".nii.gz") and os.path.basename(fout).replace(".nii.gz","")==diff_brain:
                    brain = fout
                    break
            if brain=="":
                mylib.printError(f"could not find {diff_brain} in {diff_step} outputs")
                continue
            print(f"underlay brain: {brain}")
            
            diffMasks = []
            for fout in mylib.getOutputFiles(cursor,sess,masks_diff_step,location):
                if fout.endswith(".nii.gz"):
                    fout_name = os.path.basename(fout).replace(".nii.gz","")
                    if fout_name in masks_diff:
                        diffMasks+=[fout]
            print(f"masks: {','.join(diffMasks)}")
            
            input("\npaths are already thresholded [Enter]")
            if location!="oneDrive":
                mrilib.fsleyes(mrilib.fsleyescmd_maskN(brain,diffMasks+[diffPaths],force_pdf=forcepdf,software=soft),delete=False)
            else:
                mrilib.fsleyes(mrilib.fsleyesOneDrive_maskN(step,sess,[brain],diffMasks+[diffPaths]))
                input("[Enter]")
                mrilib.fsleyesOneDrive_delete(step,sess)
            
            if input("Error? [Y]: ")=='Y':
                mylib.error(cursor,"bad output",sess,step,location)
                mylib.updateDB(cnx,cursor)
                continue
        
        mylib.checked(cursor,sess,step,False,'',location)
        mylib.updateDB(cnx,cursor)
        
    cursor.close()
    cnx.close()
