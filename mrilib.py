#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import mylib
import fileslib
import statslib
import numpy as np
from PyPDF2 import PdfFileMerger, PdfFileReader
import pdfkit

colors = ["red","blue","green","yellow","random","copper","cool","blue-lightblue","red-yellow","brain_colours_1hot","brain_colours_1hot_iso","brain_colours_2winter","brain_colours_2winter_iso","brain_colours_3warm","brain_colours_3warm_iso","brain_colours_4cool","brain_colours_4cool_iso","brain_colours_5redyell","brain_colours_5redyell_iso","brain_colours_6bluegrn","brain_colours_6bluegrn_iso","brain_colours_actc","brain_colours_actc_iso","brain_colours_blackbdy","brain_colours_blackbdy_iso","brain_colours_bluegray","brain_colours_bluegray_iso","brain_colours_bone","brain_colours_bone_iso","brain_colours_cardiac","brain_colours_cardiac_iso","brain_colours_cortex","brain_colours_cortex_iso","brain_colours_diverging_bgy","brain_colours_diverging_bgy_iso","brain_colours_diverging_bwr","brain_colours_diverging_bwr_iso","brain_colours_flow","brain_colours_flow_iso","brain_colours_french","brain_colours_french_iso","brain_colours_ge_color","brain_colours_ge_color_iso","brain_colours_gold","brain_colours_gold_iso","brain_colours_gooch","brain_colours_gooch_iso","brain_colours_greengray","brain_colours_greengray_iso","brain_colours_hotiron","brain_colours_hotiron_iso","brain_colours_nih","brain_colours_nih_fire","brain_colours_nih_fire_iso","brain_colours_nih_ice","brain_colours_nih_ice_iso","brain_colours_nih_iso","brain_colours_nih_new","brain_colours_nih_new_iso", "brain_colours_pink","brain_colours_pink_iso","brain_colours_rainramp","brain_colours_rainramp_iso","brain_colours_redgray","brain_colours_redgray_iso","brain_colours_spectrum","brain_colours_spectrum_iso","brain_colours_surface","brain_colours_surface_iso","brain_colours_x_hot","brain_colours_x_hot_iso","brain_colours_x_rain","brain_colours_x_rain_iso","cortical","greyscale","hot","hsv","render1","render1t","render2","render2t","render3","retino","subcortical","pink"]

#############################################################################
############################### Visualization ###############################
#############################################################################  
def fsleyes(cmd,delete=True):
    print(cmd)
    if cmd.startswith("fsleyes"):
        os.system(f"{cmd} >/dev/null 2>&1")
    else:
        os.system(cmd)
    if cmd.startswith("gio open"):
        input("[Enter]")
        if delete:
            os.system("rm "+cmd.replace("gio open",""))

def fsleyescmd_colors(v1):
    return "fsleyes" if not os.path.isfile(v1) else "fsleyes -xh -yh -hc "+v1+" -ot rgbvector"

def fsleyescmd_vecs(v1):
    return "fsleyes" if not os.path.isfile(v1) else "fsleyes -xh -yh -hc "+v1+" -ot linevector"

# Compare an array of registered brains with the target brain to see if they align correctly
def checkRegistration(target_brain,registered_brains,output_dir):
    if not output_dir.endswith("/"):
        output_dir+="/"
    
    cmd = "slicesdir -p "+target_brain+" -S"
    outputs = []
    for brain in registered_brains:
        cmd+=" "+brain
        outputs+=[output_dir+brain.replace("/","_").replace(".nii.gz",".png")]
        
    print("Getting screenshots of registration...")
    print(cmd)
    os.system(cmd)
    print("done")
    
    os.system(f"mv slicesdir/*.png {output_dir}")
    os.system("rm -r slicesdir")
    os.system(f"rm {output_dir}grota.png")
    
    return outputs

# All images should be in the same space
def compare_imgs(images,output_prefix,std=False):
    # Get the number of slices
    n_slices = getDims(images[0])[2]
    
    # Get the png for each slice in each img
    slice_pngs_dic = {}
    for slice_num in range(0,n_slices):
        pngs = []
        for img in images:
            print(f'Rendering slice {slice_num+1} of {n_slices}...')
            
            com = getCOM(img,True) if not std else [89, 105, 83]
            outslice = img.replace(".nii.gz",f"slice_{slice_num}.png")
            outslice = slicePNG(outslice,True,True,False,com[0],com[1],slice_num,img)
            if outslice!="":
                pngs+=[outslice]
            
        slice_pngs_dic[slice_num] = pngs
    print("DONE rendering slices\n")
    
    # Merge horizontally the pngs of all imgs for each slice   
    for slice_num in range(0,n_slices):
        print(f'Merging slices horizontally: {slice_num+1} of {n_slices}...')
        fileslib.mergeHorizontalImgs(slice_pngs_dic[slice_num],output_prefix+"_slice"+str(slice_num)+".png",True)
    print("DONE merging slices\n")
    
    # Merge vertically all slices / create pdf
    print("Creating pdf...")
    merger = PdfFileMerger()
    for slice_num in range(0,n_slices):
        print(f'Merging page {slice_num+1} of {n_slices}...')
        
        slice_path = output_prefix+"_slice"+str(slice_num)+".png"
        while os.path.isfile(slice_path):
            slice_path = slice_path.replace(".png","+.png")
        
        pdfout = fileslib.convertToPDF(slice_path,True)
        merger.append(pdfout)
        os.remove(pdfout)
    
    output = output_prefix+".pdf"
    while os.path.isfile(output):
        output = output.replace(".pdf","+.pdf")
    merger.write(output)
    merger.close()
    print(f"\npdf created: {output}")
    
    return output

def slicePNG(outfile,hidex,hidey,hidez,x,y,z,nifti):
    if not os.path.isfile(nifti) or (hidex and hidey and hidez):
        return ""
    
    if not outfile.endswith(".png"):
        outfile+=".png"
    
    os.system(f"rm -f {outfile}")
    
    [minval,maxval] = getRange(nifti,True)
    cmd = f"fsleyes render --outfile {outfile}"
    if hidex:
        cmd+=" -xh"
    if hidey:
        cmd+=" -yh"
    if hidez:
        cmd+=" -zh"
    cmd+=f" -hc -vl {x} {y} {z} {nifti} -or {minval} {maxval} >/dev/null 2>&1"
    
    os.system(cmd)
    if os.path.isfile(outfile):
        return outfile
    return ""

def fsleyesPDF_brain(brain,std=True,select=[],output=""):
    # number of slices
    dim3 = getDims(brain)[2]
    # center of mass
    com = getCOM(brain,True) if not std else [89, 105, 83]
    merger = PdfFileMerger()
    
    for slice_num in range(0,dim3):
        if len(select)>0 and slice_num not in select:
            continue
        
        print(f'Rendering slice {slice_num+1} of {dim3}...')
        
        outslice = brain.replace(".nii.gz",f"slice_{slice_num}.png") if not std else f"slice_{slice_num}.png"
        outslice = slicePNG(outslice,True,True,False,com[0],com[1],slice_num,brain)
        if outslice!="":
            print("done")
        else:
            mylib.printError("failed")
            break
        
        pdfout = fileslib.convertToPDF(outslice,True)
        merger.append(pdfout)
        os.remove(pdfout)
    
    if output=="":    
        output = brain.replace(".nii.gz",".pdf")
    while os.path.isfile(output):
        output = output.replace(".pdf","+.pdf")
    merger.write(output)
    merger.close()
    
    return output

# All imgs should be in the same space
# on_top will only work when there's two imgs or more
def fsleyescmd_brainN(brains,screen=True,standard=False,software="fsleyes",force_pdf=False,on_top=False,pdfout=""):
    dic_range = {}
    for brain in brains:
        if not os.path.isfile(brain):
            mylib.printError(f"FILE NOT FOUND {brain}")
            return software
        dic_range[brain] = getRange(brain)
    
    if not screen or force_pdf:
        print("Generating pdf...")
        
        # Generate the PDF for one brain
        if len(brains)==1:
            output = fsleyesPDF_brain(brain,std=standard,output=pdfout)
        else:
            prefix = brains[0].replace(".nii.gz","").replace(".nii","")+"_tmp" if pdfout=="" else pdfout.replace(".pdf","")
            
            # Generate the PDF for more than one brain, showing each brain next to each other
            if not on_top:
                output = compare_imgs(brains,prefix)
            
            # Generate the PDF for more than one brain, overloading the outline of the first one to the others
            else:
                cmd = "slicesdir -p "+brains[0]+" -S "+" ".join(brains[1:])
                print(cmd)
                os.system(cmd)
                output = prefix+".pdf"
                pdfkit.from_file('slicesdir/index.html',output)
                os.system("rm -r slicesdir")
            
        print(f"done: {output}")
        return "gio open "+output
    
    # Visualize brain(s) using fsleyes
    if software=="fsleyes":
        cmd = "fsleyes -xh -yh -hc"
        for brain in brains:
            [minval,maxval] = dic_range[brain]
            cmd+=f" {brain} -or {minval} {maxval}"
        return cmd
    
    # If it's only one brain, it can be opened with any software (i.e. freeview, afni, mricron) for my visualization purposes
    if len(brains)==1:
        cmd = software+" "+brains[0]
        return cmd
            
    # If there's more than one brain & I can't use fsleyes, force to open with freeview
    cmd = "freeview "+" ".join(brains)
    return cmd

# Works with software: fsleyes, freeview
# All imgs should be in the same space
def fsleyescmd_maskN(brain,masks,screen=True,standard=False,transp=50,software="fsleyes",force_pdf=False,pdfout=""):
    for img in [brain]+masks:
        if not os.path.isfile(img):
            mylib.printError(f"FILE NOT FOUND {img}")
            return software
        
    if not screen or force_pdf:
        print("Generating overlay pdf...")
        if pdfout=="":
            pdfout = "tmp.pdf"
        output = fsleyesPDF_maskN(brain,masks,pdfout,std=standard,trans=transp)
        print("done: "+output)
        return "gio open "+output
    
    elif software=="fsleyes":
        cmd = "fsleyes -xh -yh -hc "+brain
        j = 0
        for i in range(len(masks)):
            cmd+=" "+masks[i]+" -cm "+colors[j]+" -a "+str(transp)
            j = j+1 if j<len(colors)-2 else 0
        return cmd
    
    else:
        cmd = "freeview "+brain+" "+" ".join(masks)
        return cmd

# Remember to delete files after calling this function by calling the function below
def fsleyesOneDrive_maskN(step,sess,brains,masks,transp=50,leaveAs4D=False):
    # Create the output dir in oneDrive
    oneDir = f"fsleyes/{step}_{sess}/"
    os.system(f"rclone mkdir {oneDir}")
    
    # Extract the first volume for faster transmission of the brain files
    rm_files = []
    final_brains = []
    for brain in brains:
        if getDims(brain)[3]>1 and not leaveAs4D:
            orig_brain = brain
            brain = brain.replace(".nii.gz","_v0_tmp.nii.gz")
            os.system(f"fslroi {orig_brain} {brain} 0 1")
            rm_files+=[brain]
        final_brains+=[brain]
    
    # Copy files to oneDrive
    print(f"Copying files to oneDrive ({oneDir})...")
    copied = []
    for img in final_brains+masks:
        if not os.path.isfile(img):
            mylib.printError(f"FILE NOT FOUND {img}")
            os.system(f"rclone purge oneDrive:{oneDir}")
            return ""
        
        img_name = os.path.basename(img)
        print(img_name)
        os.system(f"rclone copyto --copy-links {img} oneDrive:{oneDir}{img_name}")
        copied+=[img_name]
    print("done")
    
    # Remove the extracted first volume
    if len(rm_files)>0:
        os.system(f"rm {' '.join(rm_files)}")
    
    # Create cmd to visualize in Mac
    basedir = f"/Users/monica/Documents/OneDrive - mcw.edu/{oneDir}"
    cmd = f"/usr/local/fsl/bin/fsleyes -xh -yh -hc"
    
    for brain in final_brains:
        cmd+=" '{basedir}{os.path.basename(brain)}'"
    
    j = 0
    for i in range(len(masks)):
        cmd+=f" '{basedir}{os.path.basename(masks[i])}' -cm {colors[j]} -a {transp}"
        j = j+1 if j<len(colors)-2 else 0
        
    return cmd

def fsleyesOneDrive_delete(step,sess):
    os.system(f"rclone purge oneDrive:fsleyes/{step}_{sess}/")

# Axial view
def fsleyesPDF_maskN(brain,masks,output,std=True,trans=50,silence=False,colors_opt=[]):
    for img in [brain]+masks:
        if not os.path.isfile(img):
            mylib.printError(f"FILE NOT FOUND: {img}")
            return ""
    
    # Slices
    if not silence:
        print("Getting slices range...")
    box_min = []
    box_max = []
    for mask in masks:
        box = getROIbox(mask)[2]
        box_min+=[box[0]]
        box_max+=[box[1]]
    dim3_min = int(min(box_min))
    dim3_max = int(max(box_max))
    if dim3_max==dim3_min:
        dim3_max+=1
    if not silence:
        print(f'Slices range: {dim3_min}-{dim3_max}')   
    
    # center of mass
    com = getCOM(brain,silence) if not std else [89, 105, 83]
    if not silence:
        print(f"Center of mass: {com[0]},{com[1]},{com[2]}")
    
    # Range values
    if not silence:
        print("Getting range values...")
    minvals = []
    maxvals = []
    for mask in masks:
        [minval,maxval] = getRange(mask,silence)
        if minval==0 and maxval==0:
            mylib.printError(f"EMPTY MASK: {mask}")
            return ""
        minvals+=[minval]
        maxvals+=[maxval]
    if not silence:
        print("done")
    
    # Name of the final output
    if not output.endswith(".pdf"):
        output = output+".pdf"
    while os.path.isfile(output):
        output = output.replace(".pdf","+.pdf")
    
    # Create directory to store temporary files    
    tmp_dir = output.replace(".pdf","")
    while os.path.isdir(tmp_dir):
        tmp_dir+="+"
    os.system("mkdir "+tmp_dir)
    
    # Replace colors array if necessary
    if len(masks)==len(colors_opt):
        colors_use = colors_opt
    else:
        colors_use = colors
    
    # Create PDF
    merger = PdfFileMerger()
    for slice in range(dim3_min,dim3_max):
        if not silence:
            print(f'Rendering slice {slice+1} of {dim3_max}...')
        
        outslice = tmp_dir+"/slice_"+str(slice)+".png"
        while os.path.isfile(outslice):
            outslice = outslice.replace(".png","+.png")
        
        cmd = f"fsleyes render --outfile {outslice} -xh -yh -hc -vl {com[0]} {com[1]} {slice} {brain} "
        i = 0
        j = 0
        for mask in masks:
            cmd+=f"{mask} -or {minvals[i]} {maxvals[i]} -cm {colors_use[j]} -a {trans} "
            i+=1
            j+=1
            if j>len(colors_use):
                j=0
        cmd+=" >/dev/null 2>&1"
        os.system(cmd)
        
        if os.path.isfile(outslice) and not silence:
            print("done: "+outslice)
        elif not os.path.isfile(outslice):
            print(cmd)
            mylib.printError(f"failed: {outslice}")
            merger.close()
            os.system("rm -f "+tmp_dir+"/*")
            os.system("rm -rf "+tmp_dir)
            return ""
        
        pdfout = fileslib.convertToPDF(outslice,True)
        merger.append(pdfout)
        os.remove(pdfout)
    
    if not silence:
        print(f"PDF generated: {output}")    
    merger.write(output)
    merger.close()
    os.system(f"rm -f {tmp_dir}/*")
    os.system(f"rm -rf {tmp_dir}")
    
    return output

# Coronal view
def fsleyesPDF_probtrackX_3(brain,masks,output,fdt_paths="",silence=False):
    allfiles = [brain]+masks if fdt_paths=="" else [brain]+masks+[fdt_paths]
    for img in allfiles:
        if not os.path.isfile(img):
            mylib.printError("FILE NOT FOUND: "+img)
            return ""
    
    # Slices    
    box_min = []
    box_max = []
    if not silence:
        print("Getting min and max slices...")
    all_masks = masks if fdt_paths=="" else masks+[fdt_paths]
    for mask in all_masks:
        box = getROIbox(mask)[1]
        box_min+=[box[0]]
        box_max+=[box[1]]
    dim2_min = int(min(box_min))
    dim2_max = int(max(box_max))
    if not silence:
        print(f'done: {dim2_min},{dim2_max}')
    
    # center of mass
    if not silence:
        print(f"\nGetting {brain} COM...")
    com = getCOM(brain,silence)
    if not silence:
        print(f"done: {com[0]},{com[1]},{com[2]}")
    
    # Range values
    minvals = {}
    maxvals = {}
    if not silence:
        print("\nGetting range values for each mask...")
    for mask in masks:
        [minval,maxval] = getRange(mask,silence)
        if minval==0 and maxval==0:
            mylib.printError(f"EMPTY MASK: {mask}")
            return ""
        minvals[mask] = minval
        maxvals[mask] = maxval
    [min_fdt,max_fdt] = getRange(fdt_paths) if fdt_paths!="" else [0.0,0.0]
    if not silence:
        print("done")
    
    # Name of the final output
    if not output.endswith(".pdf"):
        output = output+".pdf"
    while os.path.isfile(output):
        output = output.replace(".pdf","+.pdf")
        
    # Create directory to store temporary files
    tmp_dir = output.replace(".pdf","")
    while os.path.isdir(tmp_dir):
        tmp_dir+="+"
    os.system("mkdir "+tmp_dir)
    
    # Create PDF
    merger = PdfFileMerger()
    for slice in range(dim2_min,dim2_max):
        if not silence:
            print(f'\nRendering slice {slice+1} of {dim2_max}...')
        
        outslice = tmp_dir+"/slice_"+str(slice)+".png"
        while os.path.isfile(outslice):
            outslice = outslice.replace(".png","+.png")
        
        cmd = f"fsleyes render --outfile {outslice} -xh -zh -hc -vl {com[0]} {slice} {com[2]} {brain} "
        i = 0
        for mask in masks:
            cmd+=f"{mask} -or {minvals[mask]} {maxvals[mask]} -cm {colors[i]} -a 50 "
            i+=1
            if i>len(colors):
                i=0
        if fdt_paths!="":
            cmd+=f"{fdt_paths} -or {min_fdt} {max_fdt} -cm red-yellow -a 80 >/dev/null 2>&1"
        else:
            cmd+=">/dev/null 2>&1"
        os.system(cmd)
        
        if os.path.isfile(outslice) and not silence:
                print("done: "+outslice)
        elif not os.path.isfile(outslice):
            print(cmd)
            mylib.printError(f"failed: {outslice}")
            merger.close()
            os.system(f"rm -f {tmp_dir}/*")
            os.system(f"rm -rf {tmp_dir}")
            return ""
        
        pdfout = fileslib.convertToPDF(outslice,True)
        merger.append(pdfout)
        os.remove(pdfout)
    
    if not silence:
        print(f"PDF generated: {output}")    
    merger.write(output)
    merger.close()
    os.system(f"rm -f {tmp_dir}/*")
    os.system(f"rm -rf {tmp_dir}")
    
    return output

# Sagittal view
def fsleyesPDF_probtrackX_2(brain,masks,output,fdt_paths="",silence=False):
    allfiles = [brain]+masks if fdt_paths=="" else [brain]+masks+[fdt_paths]
    for img in allfiles:
        if not os.path.isfile(img):
            mylib.printError("FILE NOT FOUND: "+img)
            return ""
    
    # Slices    
    box_min = []
    box_max = []
    if not silence:
        print("Getting min and max slices...")
    all_masks = masks if fdt_paths=="" else masks+[fdt_paths]
    for mask in all_masks:
        box = getROIbox(mask)[0]
        box_min+=[box[0]]
        box_max+=[box[1]]
    dim1_min = int(min(box_min))
    dim1_max = int(max(box_max))
    if not silence:
        print(f'done: {dim1_min},{dim1_max}')
    
    # center of mass
    if not silence:
        print(f"\nGetting {brain} COM...")
    com = getCOM(brain,silence)
    if not silence:
        print(f"done: {com[0]},{com[1]},{com[2]}")
    
    # Range values
    minvals = {}
    maxvals = {}
    if not silence:
        print("\nGetting range values for each mask...")
    for mask in masks:
        [minval,maxval] = getRange(mask,silence)
        if minval==0 and maxval==0:
            mylib.printError(f"EMPTY MASK: {mask}")
            return ""
        minvals[mask] = minval
        maxvals[mask] = maxval
    [min_fdt,max_fdt] = getRange(fdt_paths) if fdt_paths!="" else [0.0,0.0]
    if not silence:
        print("done")
    
    # Name of the final output
    if not output.endswith(".pdf"):
        output = output+".pdf"
    while os.path.isfile(output):
        output = output.replace(".pdf","+.pdf")
        
    # Create directory to store temporary files
    tmp_dir = output.replace(".pdf","")
    while os.path.isdir(tmp_dir):
        tmp_dir+="+"
    os.system("mkdir "+tmp_dir)
    
    # Create PDF
    merger = PdfFileMerger()
    for slice in range(dim1_min,dim1_max):
        if not silence:
            print(f'\nRendering slice {slice+1} of {dim1_max}...')
        
        outslice = tmp_dir+"/slice_"+str(slice)+".png"
        while os.path.isfile(outslice):
            outslice = outslice.replace(".png","+.png")
        
        cmd = f"fsleyes render --outfile {outslice} -yh -zh -hc -vl {slice} {com[1]} {com[2]} {brain} "
        i = 0
        for mask in masks:
            cmd+=f"{mask} -or {minvals[mask]} {maxvals[mask]} -cm {colors[i]} -a 50 "
            i+=1
            if i>len(colors):
                i=0
        if fdt_paths!="":
            cmd+=f"{fdt_paths} -or {min_fdt} {max_fdt} -cm red-yellow -a 80 >/dev/null 2>&1"
        else:
            cmd+=">/dev/null 2>&1"
        os.system(cmd)
        
        if os.path.isfile(outslice) and not silence:
                print("done: "+outslice)
        elif not os.path.isfile(outslice):
            print(cmd)
            mylib.printError(f"failed: {outslice}")
            merger.close()
            os.system(f"rm -f {tmp_dir}/*")
            os.system(f"rm -rf {tmp_dir}")
            return ""
        
        pdfout = fileslib.convertToPDF(outslice,True)
        merger.append(pdfout)
        os.remove(pdfout)
    
    if not silence:
        print(f"PDF generated: {output}")    
    merger.write(output)
    merger.close()
    os.system(f"rm -f {tmp_dir}/*")
    os.system(f"rm -rf {tmp_dir}")
    
    return output

# Axial view        
def fsleyesPDF_probtrackX(brain,masks,fdt_paths,output):
    for img in [brain]+masks+[fdt_paths]:
        if not os.path.isfile(img):
            mylib.printError("FILE NOT FOUND")
            return ""
    
    # Slices    
    box_min = []
    box_max = []
    for mask in masks+[fdt_paths]:
        print("getting box for "+mask+"...")
        box = getROIbox(mask)[2]
        box_min+=[box[0]]
        box_max+=[box[1]]
    dim3_min = int(min(box_min))
    dim3_max = int(max(box_max))
    
    # center of mass
    com = getCOM(brain,True)
    
    # Range values
    minvals = []
    maxvals = []
    for mask in masks:
        [minval,maxval] = getRange(mask)
        if minval==0 and maxval==0:
            mylib.printError(f"EMPTY MASK: {mask}")
            return ""
        minvals+=[minval]
        maxvals+=[maxval]
    [min_fdt,max_fdt] = getRange(fdt_paths)
    
    # Name of the final output
    if not output.endswith(".pdf"):
        output = output+".pdf"
    os.system(f"rm -f {output}")
    
    # Create directory to store temporary files
    tmp_dir = output.replace(".pdf","")
    os.system(f"rm -rf {tmp_dir}")
    os.system("mkdir "+tmp_dir)
    
    # Create PDF
    merger = PdfFileMerger()
    for slice in range(dim3_min,dim3_max):
        print(f'Rendering slice {slice+1} of {dim3_max}...')
        
        outslice = f"{tmp_dir}/slice_{slice}.png"
        while os.path.isfile(outslice):
            outslice = outslice.replace(".png","+.png")
        
        cmd = f"fsleyes render --outfile {outslice} -xh -yh -hc -vl {com[0]} {com[1]} {slice} {brain} "
        i = 0
        j = 0
        for mask in masks:
            cmd+=f"{mask} -or {minvals[i]} {maxvals[i]} -cm {colors[j]} -a 50 "
            i+=1
            j+=1
            if j>len(colors):
                j=0
        cmd+=f"{fdt_paths} -or {min_fdt} {max_fdt} -cm red-yellow -a 80 >/dev/null 2>&1"
        os.system(cmd)
        
        if os.path.isfile(outslice):
            print("done: "+outslice)
        else:
            print(cmd)
            mylib.printError(f"failed: {outslice}")
            merger.close()
            os.system(f"rm -f {tmp_dir}/*")
            os.system(f"rm -rf {tmp_dir}")
            return ""
        
        pdfout = fileslib.convertToPDF(outslice,True)
        merger.append(pdfout)
        os.remove(pdfout)
        
    print(f"PDF generated: {output}")    
    merger.write(output)
    merger.close()
    os.system(f"rm -f {tmp_dir}/*")
    os.system(f"rm -rf {tmp_dir}")
    
    return output

# sess & step parameters are only used when location=oneDrive
def vecmenu(FA,V1,mask,screen,location,soft="fsleyes",sess="",step=""):
    opt=''
    while opt!='2' and opt!='3':
        if opt=='1':
            os.system(f"3dmask_tool -overwrite -input {mask} -prefix {mask} -dilate_input -1 >/dev/null 2>&1")
        elif opt=='0':
            os.system(f"3dmask_tool -overwrite -input {mask} -prefix {mask} -dilate_input 1 >/dev/null 2>&1")
                
        os.system(f"fslmaths {FA} -mul {mask} {FA}")
        if os.path.isfile(V1):
            os.system(f"fslmaths {V1} -mul {mask} {V1}")
        
        if location!="oneDrive":     
            fsleyes(fsleyescmd_brainN([FA],screen,software=soft))
        else:
            fsleyesOneDrive_maskN(step,sess,[FA],[])
            input("[Enter]")
            fsleyesOneDrive_delete(step,sess)
       
        print("0. expand")
        print("1. erode")
        print("2. done*")
        print("3. bad")
        opt = input(">> ")
        if opt=="":
            opt = "2"
            
    return opt
    
# sess & step parameters are only used when location=oneDrive
def faMenu(FA,mask,screen,location,soft="fsleyes",sess="",step=""):
    if not os.path.isfile(FA):
        mylib.printError(f"{FA} file not found")
        return '3'
    if not os.path.isfile(mask):
        mylib.printError(f"{mask} file not found")
        return '3'
    
    opt=''
    while opt!='2' and opt!='3':
        if opt=='1':
            os.system(f"3dmask_tool -overwrite -input {mask} -prefix {mask} -dilate_input -1 >/dev/null 2>&1")
        elif opt=='0':
            os.system(f"3dmask_tool -overwrite -input {mask} -prefix {mask} -dilate_input 1 >/dev/null 2>&1")
            
        os.system(f"fslmaths {FA} -mul {mask} {FA}")
        
        if location!="oneDrive":     
            fsleyes(fsleyescmd_brainN([FA],screen,software=soft))
        else:
            fsleyesOneDrive_maskN(step,sess,[FA],[])
            input("[Enter]")
            fsleyesOneDrive_delete(step,sess)
       
        print("0. expand")
        print("1. erode")
        print("2. done*")
        print("3. bad")
        opt = input(">> ")
        if opt=="":
            opt = '2'
        
    return opt

# step & sess only used if location==oneDrive
def maskmenu(brain,mask,location,screen=True,opt='',step="",sess=""):
    dil_area = mask.replace(".nii.gz","_dil_area.nii.gz")
    erode_area = mask.replace(".nii.gz","_erode_area.nii.gz")
    edited = False
    
    while opt!='2' and opt!='3':
        if opt=='0':
            os.system(f"3dmask_tool -overwrite -input {mask} -prefix {mask} -dilate_input 1 >/dev/null 2>&1")
            edited = True
        elif opt=='1':
            os.system(f"3dmask_tool -overwrite -input {mask} -prefix {mask} -dilate_input -1 >/dev/null 2>&1")
            edited = True
        elif opt=='4':
            os.system(f"3dmask_tool -overwrite -input {mask} -prefix {mask} -fill_holes -fill_dirs xy >/dev/null 2>&1")
            edited = True
        elif opt=='5':
            if not os.path.isfile(dil_area):
                mylib.printError(f"{dil_area} mask area does not exist")
            
            else:
                # Expand original mask
                expanded = mask.replace(".nii.gz","_expanded.nii.gz")
                os.system(f"3dmask_tool -overwrite -input {mask} -prefix {expanded} -dilate_input 1 >/dev/null 2>&1")
                
                # Calculate the difference
                diff = mask.replace(".nii.gz","_dif.nii.gz")
                os.system(f"fslmaths {expanded} -sub {mask} {diff}")
                
                # Calculate the area to be added to the mask
                toAdd = mask.replace(".nii.gz","_toAdd.nii.gz")
                os.system(f"fslmaths {dil_area} -mul {diff} {toAdd}")
                
                # Add selected dilated area to mask & binarize
                os.system(f"fslmaths {mask} -add {toAdd} -bin {mask}")
                os.system(f"rm {expanded} {diff} {toAdd}")
                
                edited = True
        elif opt=='6':
            if not os.path.isfile(erode_area):
                mylib.printError(f"{erode_area} mask area does not exist")
            
            else:
                # Dilate original mask
                eroded = mask.replace(".nii.gz","_eroded.nii.gz")
                os.system(f"3dmask_tool -overwrite -input {mask} -prefix {eroded} -dilate_input -1 >/dev/null 2>&1")
                
                # Calculate the difference
                diff = mask.replace(".nii.gz","_dif.nii.gz")
                os.system(f"fslmaths {mask} -sub {eroded} {diff}")
                
                # Calculate the area to be subtracted from the mask
                toDel = mask.replace(".nii.gz","_toDel.nii.gz")
                os.system(f"fslmaths {erode_area} -mul {diff} {toDel}")
                
                # Delete selected eroded area from mask & binarize
                os.system(f"fslmaths {mask} -sub {toDel} -bin {mask}")
                os.system(f"rm {eroded} {diff} {toDel}")
                
                edited = True
        elif opt=='7':
            for i in range(4):
                # Fill any wholes
                os.system(f"3dmask_tool -overwrite -input {mask} -prefix {mask} -fill_holes -fill_dirs xy >/dev/null 2>&1")
                # Expand mask
                os.system(f"3dmask_tool -overwrite -input {mask} -prefix {mask} -dilate_input 1 >/dev/null 2>&1")
            for i in range(4):
                # Dilatate mask
                os.system(f"3dmask_tool -overwrite -input {mask} -prefix {mask} -dilate_input -1 >/dev/null 2>&1")
            edited = True
        elif opt=='8':
            os.system(f"fslmaths {mask} -bin {mask}")
            os.system(f"3dcalc -overwrite -a {mask} -expr '1-a' -prefix {mask}")
        elif opt=='9':
            os.system(f"3dmask_tool -overwrite -input {mask} -prefix {mask} -fill_holes >/dev/null 2>&1")
        
        if location!="oneDrive":
            softwareP = "fsleyes" if location=="petrov" else "freeview"
            fsleyes(fsleyescmd_maskN(brain,[mask],screen,software=softwareP))
        else:
            fsleyes(fsleyesOneDrive_maskN(step,sess,[brain],[mask]))
            input("[Enter]")
            fsleyesOneDrive_delete(step,sess)
       
        print("0. Expand")
        print("1. Erode")
        print("2. Good*")
        print("3. Bad")
        print("4. Fill holes per slice")
        print("5. Expland in _dil_area")
        print("6. Erode in _erode_area (not tested)")
        print("7. Fix coliflower look")
        print("8. Invert mask")
        print("9. Fill holes in volume")
        opt = input(">> ")
        if opt=="":
            opt = '2'
    
    os.system("rm -f "+dil_area+" "+erode_area)    
    return [opt,edited]

#############################################################################
############################ Get info from MRI ##############################
############################################################################# 
def getROIbox(roiPath):
    box = np.zeros((3,2))
    if not os.path.isfile(roiPath):
        return box
    
    # Get the coordinates for all non-zero voxels
    ext = fileslib.getFileExt(roiPath)
    coords = roiPath.replace(ext,".txt")
    os.system(f"rm -f {coords}")
    os.system(f"3dmaskdump -nozero -quiet -xyz -o {coords} {roiPath}")
    if fileslib.linesCols(coords,' ')[0]==0:
        os.remove(coords)
        return box
    
    # Get the min and max X coordinates
    box[0,0] = int(mylib.systemOut("cat "+coords+" | awk '{print $1}' | sort -nu | head -n1"))
    box[0,1] = int(mylib.systemOut("cat "+coords+" | awk '{print $1}' | sort -nu | tail -n1"))
    
    # Get the min and max Y coordinates
    box[1,0] = int(mylib.systemOut("cat "+coords+" | awk '{print $2}' | sort -nu | head -n1"))
    box[1,1] = int(mylib.systemOut("cat "+coords+" | awk '{print $2}' | sort -nu | tail -n1"))
    
    # Get the min and max Z coordinates
    box[2,0] = int(mylib.systemOut("cat "+coords+" | awk '{print $3}' | sort -nu | head -n1"))
    box[2,1] = int(mylib.systemOut("cat "+coords+" | awk '{print $3}' | sort -nu | tail -n1"))
    
    os.remove(coords)
    return box

def getnvoxels(filepath):
    return int(mylib.systemOut("fslstats "+filepath+" -V | awk '{print $1}'")) if os.path.isfile(filepath) else 0

def getVolmm(filepath):
    return float(mylib.systemOut("fslstats "+filepath+" -V | awk '{print $2}'")) if os.path.isfile(filepath) else 0.0

def getPixDims(filepath):
    if os.path.isfile(filepath):
        dim1 = float(mylib.systemOut(f"fslval {filepath} pixdim1"))
        dim2 = float(mylib.systemOut(f"fslval {filepath} pixdim2"))
        dim3 = float(mylib.systemOut(f"fslval {filepath} pixdim3"))
        # dim4 is TR
        dim4 = float(mylib.systemOut(f"fslval {filepath} pixdim4"))
        return [dim1,dim2,dim3,dim4]
    else:
        return [0.0,0.0,0.0,0.0]

def getDims(filepath):
    if os.path.isfile(filepath):
        dim1 = int(mylib.systemOut(f"fslval {filepath} dim1"))
        dim2 = int(mylib.systemOut(f"fslval {filepath} dim2"))
        # dim3 is nslices
        dim3 = int(mylib.systemOut(f"fslval {filepath} dim3"))
        # dim4 is nvols
        dim4 = int(mylib.systemOut(f"fslval {filepath} dim4"))
        return [dim1,dim2,dim3,dim4]
    else:
        return [0,0,0,0]
    
def copyDims(source,destination):
    os.system("fslcpgeom "+source+" "+destination)

def getMean(filepath):
    return float(mylib.systemOut(f"fslstats {filepath} -M")) if os.path.isfile(filepath) else 0.0

def getMaskedMean(niftipath,maskpath):
    if not os.path.isfile(niftipath) or not os.path.isfile(maskpath):
        return 0.0
    
    filepath = niftipath.replace(".nii.gz","_tmp.nii.gz")
    os.system(f"fslmats {niftipath} -mul {maskpath} {filepath}")
    res = getMean(filepath)
    os.system(f"rm {filepath}")
    return res

def getSTD(filepath):
    return float(mylib.systemOut(f"fslstats {filepath} -S")) if os.path.isfile else 0.0

def getMaskedSTD(niftipath,maskpath):
    if not os.path.isfile(niftipath) or not os.path.isfile(maskpath):
        return 0.0
    
    filepath = niftipath.replace(".nii.gz","_tmp.nii.gz")
    os.system(f"fslmats {niftipath} -mul {maskpath} {filepath}")
    res = getSTD(filepath)
    os.system(f"rm {filepath}")
    return res

def getCOM(filepath,silence=False):
    if os.path.isfile(filepath):
        x = round(float(mylib.systemOut("fslstats "+filepath+" -C | awk '{print $1}'")))
        y = round(float(mylib.systemOut("fslstats "+filepath+" -C | awk '{print $2}'")))
        z = round(float(mylib.systemOut("fslstats "+filepath+" -C | awk '{print $3}'")))
        if not silence:
            print(f"COM {filepath}: {x},{y},{z}")
        return [x,y,z]
    else:
        mylib.printError(f"{filepath} not found")
        return [0.0,0.0,0.0]

def getRange(filepath,quiet=False):
    if os.path.isfile(filepath):
        if not quiet:
            print("Getting file range, might take a while if file is big...")
        minval = float(mylib.systemOut("fslstats "+filepath+" -R | awk '{print $1}'"))
        maxval = float(mylib.systemOut("fslstats "+filepath+" -R | awk '{print $2}'"))
        return [minval,maxval]
    else:
        mylib.printError("File not found: "+filepath)
        return [0.0,0.0]

def getClusterLocations(cursor,cluster,location,pipeline):
    if not os.path.isfile(cluster):
        return
    print("cluster:\n"+cluster)
    
    # Save the location information for each atlas in a temporary file
    print("Checking atlases (might take a while)...")
    atlases = open(mylib.getScratchDirPipe(cursor,location,pipeline)+"/atlases.txt",'r')
    os.system("rm -f tmp.txt")
    for atlas in [atlas.replace("\n","") for atlas in atlases]:
        os.system("echo \"#"+atlas+"\" >> tmp.txt")
        os.system("atlasquery -a \""+atlas+"\" -m "+cluster+" >> tmp.txt")
    atlases.close()
    print("done")
    
    # Create a dictionary for each non-empty atlas
    # atlases key: atlas name, value: atlas dictionary
    # atlases dictionary in atlases: key: labelled area, value: probability (in percentage) of the cluster being a member of that area
    atlases = {}
    atlas = ""
    dic = {}
    # dictionary with the total percentage of each atlas (to order)
    # key: atlas name, value: total percentage
    tot_perc = {}
    total = 0
    tmp = open("tmp.txt",'r')
    for line in [line.replace("\n","") for line in tmp]:
        if line.startswith("#"):
            if len(dic)>0:
                atlases[atlas] = dic
                tot_perc[atlas] = total
            atlas = line.replace("#","")
            dic = {}
            total = 0
        else:
            # perc: probability (in percentage) of the cluster being a member of array[0]
            val = line.split(":")[-1]
            perc = float(val)
            label = line.replace(":"+str(val),"")
            if perc>=1 and ("Unclassified" not in label) and len(label.replace("*","").replace(".","").replace("0",""))>0:
                dic[label] = perc
                total+=perc
    tmp.close()
    os.system("rm tmp.txt")
    
    # Sort the dictionaries by the total percentage of labelled areas
    sorted_atlases = mylib.orderDicByVals_dec(tot_perc)
    
    # Create the html file
    html = cluster.replace(".nii.gz",".html")
    fout = open(html,'w')
    for atlas in sorted_atlases:
        fout.write("<h1>"+atlas+"</h1>\n")
        dic = atlases[atlas]
        fout.write("<p>\n")
        for label in mylib.orderDicByVals_dec(dic):
            perc = dic[label]
            fout.write(str(perc)+"% "+label+"<br>\n")
        fout.write("</p>\n\n")
    fout.close()
    
    print(f"atlas summary:\n{html}")
    return html

#############################################################################
############################ Check info from MRI ##############################
############################################################################# 
def checkBrainOutside(brain,mask):
    if not os.path.isfile(brain) or not os.path.isfile(mask):
        return False
    else:
        # Multiply brain by mask
        os.system(f"fslmaths {brain} -mul {mask} {brain}")
        # Invert the mask
        os.system(f"fslmaths {mask} -sub 1 -mul -1 rm1.nii.gz")
        # Check that the brain has no voxels in the inverted mask region
        os.system(f"fslmaths rm1.nii.gz -mul {brain} rm2.nii.gz")
        nvox = getnvoxels("rm2.nii.gz")
        # Remove temporary files
        os.system("rm rm1.nii.gz rm2.nii.gz")
        return True if nvox==0 else False

# 1. Check that the # of lines in each time series is the same as number of volumes in the nifti file
# 2. Generate png files and temporary pdf file to check the actual time series   
def checkTimeSeries(niftipath,tsarray,screen,automatic=False):
    pdfout = niftipath.replace('.nii.gz','.pdf')
    while os.path.isfile(pdfout):
        pdfout = pdfout.replace(".pdf","+.pdf")
    
    if not os.path.isfile(pdfout):
        print("Generating time series graphs...")
        pngs = []
        for ts in tsarray:
            if fileslib.linesCols(ts,' ')[0]!=getDims(niftipath)[3]:
                return ["","Wrong number of lines in "+ts]
            statslib.plotTS(ts)
            if os.path.isfile(ts.replace('.txt','.png')):
                pngs+=[ts.replace('.txt','.png')]
            else:
                return ["","Could not generate graph for "+ts]
    
        print("Generating QC pdf...")
        merger = PdfFileMerger()
        for img in pngs:
            out = fileslib.convertToPDF(img)
            merger.append(PdfFileReader(out))
            os.remove(out)
        
        merger.write(pdfout)
        merger.close()
    
    if not automatic:    
        fileslib.openPDFfile(pdfout,screen)
        if input("Error [Y]? ")!='Y':
            msg = ''
            os.remove(niftipath.replace('.nii.gz','.pdf'))
        else:
            msg = "Error in some time series. Check "+pdfout
    
        return [pdfout,msg]
    
    else:
        return [pdfout,''] if os.path.isfile(pdfout) else ['','']
    
def check3Ddimentions(img,array):
    expected_dims = getDims(img)
    del expected_dims[-1]
    expected_pixDims = getPixDims(img)
    del expected_pixDims[-1]
    
    for outfile in array:
        res_dims = getDims(outfile)
        del res_dims[-1]
        res_pixDims = getPixDims(outfile)
        del res_pixDims[-1]
        
        if res_dims!=expected_dims or res_pixDims!=expected_pixDims:
            mylib.printError(f"expected dimensions: {img}\n{expected_dims}\n{expected_pixDims}")
            mylib.printError(f"output dimensions: {outfile}\n{res_dims}\n{res_pixDims}")
            return False
        
    return True

def check4Ddimentions(img,array):
    expected_dims = getDims(img)
    expected_pixDims = getPixDims(img)
    
    for outfile in array:
        res_dims = getDims(outfile)
        res_pixDims = getPixDims(outfile)
        
        if res_dims!=expected_dims or res_pixDims!=expected_pixDims:
            mylib.printError(f"expected dimensions: {img}\n{expected_dims}\n{expected_pixDims}")
            mylib.printError(f"output dimensions: {outfile}\n{res_dims}\n{res_pixDims}")
            return False
        
    return True