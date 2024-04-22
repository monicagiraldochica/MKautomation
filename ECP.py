#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import mylib
import os
import mrilib
import fileslib

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
        
def doneCamino(sess_list,username,password,hostdir,db_name,screen,location,automatic,absolute=[],priority=[]):
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
    step = "camino_wlf"
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

def doneDTIFIT(sess_list,username,password,hostdir,db_name,screen,automatic,location,other_imgs=[],priority=[],absolute=[]):
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
    step = "finaldtifit"
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
def doneProbmap(sess_list,username,password,hostdir,db_name,screen,automatic,location,step,priority=[],absolute=[],forcepdf=False,masks_std_step="",masks_std=[],paths_std="",ants_step='',ants_brain='',paths_diff="",masks_diff_step="",masks_diff=[],diff_step="",diff_brain=""):
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
    
def doneCaminoSNR(sess_list,username,password,hostdir,db_name,location,priority=[],absolute=[]):
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
    step = "camino_snr"
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
    
def doneTransfSklt(sess_list,username,password,hostdir,db_name,screen,automatic,location,priority=[],absolute=[]):
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
    step = "transform_sklt"
    dtifit_step="camino_wlf"
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

def diffdirsInfo(cursor,sess,pipe,location):
    df = mylib.getEntriesFromTable(cursor,["75_AP","75_PA","76_AP","76_PA"],"sessions",["sess"],[sess]).swapaxes('columns', 'rows')
    df.loc[:,0] = df.loc[:,0]==1
    sbjDir = mylib.getSessDir(cursor,pipe,location,sess)+"/raw/"
    sbjID = mylib.getSbjID(cursor,sess)
    df.loc[:,1] = [sbjDir+sbjID+"_3T_DWI_dir"+idx for idx in df.index.tolist()]
    df.columns = ["included","prefix"]    
    return df
    
def donePostEddy(sess_list,username,password,hostdir,db_name,screen,automatic,location,priority=[],absolute=[]):
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
    step = "postEddy"
    cnx = mylib.connect(username,password,hostdir,db_name)
    cursor = cnx.cursor()
    mylib.printDescription(cursor,step)
    print(f'Checking {step} on: {",".join(array)}')

    nn=0
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
    
def doneCaminoSeries(sess_list,username,password,hostdir,db_name,location,priority=[],absolute=[]):
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
    step = "camino_series"
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
    
def done3Dmask(sess_list,username,password,hostdir,db_name,screen,automatic,location,priority=[],absolute=[]):
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
    step = "3dmask"
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