#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Monica Keith"
__status__ = "Production"

import pathlib
from PIL import Image, ImageDraw, ImageFont
import os
import PyPDF2 
import mylib
import pandas as pd
import glob
import cv2

#############################################################################
########################## Manipulating text files ##########################
#############################################################################
def fileToDataFrame(infile,header,separator):
    return pd.read_csv(infile, sep=separator) if header else pd.read_csv(infile, sep=separator, header=None)

def lineInFile(infile,line):
    if not os.path.isfile(infile):
        return False
    
    try:
        ans = False
        fo = open(infile,'r')
        for lastline in fo:
            if line in lastline:
                ans = True
                break
    except ValueError as err:
        error = ("an exception occured while reading "+os.path.basename(infile)+": "+str(err)).replace('\n','').replace('\'','')
        raise ValueError(error)
    finally:
        fo.close()
        return ans
    
def getLinesFromFile(infile,text):
    lines_found = []
    if not os.path.isfile(infile):
        return lines_found
    
    try:
        fo = open(infile,'r')
        for line in fo:
            if text in line:
                lines_found+=[line.replace("\n","")]
    except ValueError as err:
        error = ("an exception occured while reading "+infile+": "+str(err)).replace('\n','').replace('\'','')
        raise ValueError(error)
    finally:
        fo.close()
        return lines_found
    
def warningInLog(log):
    if not os.path.isfile(log):
        return ''
    
    try:
        ans = ''
        fo = open(log,'r')
        for lastline in fo:
            if "warning" in lastline.lower():
                ans = lastline
                break
        warning = ans.replace("\n",'')
    except ValueError as err:
        warning = ("an exception occured while reading "+os.path.basename(log)+": "+str(err)).replace('\n','').replace('\'','')
    finally:
        fo.close()
        if len(warning)>250:
            warning = warning[0:249]
        return warning

def errorInLog(efile):
    if not os.path.isfile(efile):
        return "log file not found: "+efile
    try:
        ans = ''
        fo = open(efile,'r')
        errorList = ["exception","missing operand","offending address","address not mapped","segmentation violation","segmentation fault","error","cannot open","command not found","no such file or directory","core dumped","job killed: walltime","failed"]
        for lastline in fo:
            for error in errorList:
                if error in lastline.lower():
                    ans = lastline
                    break
        error = ans.replace("\n",'')
    except ValueError as err:
        error = ("an exception occured while reading "+os.path.basename(efile)+": "+str(err)).replace('\n','').replace('\'','')
    finally:
        fo.close()
        if len(error)>250:
            error = error[0:249]
        return error

def file2array(infile,separator=' '):
    if not os.path.isfile(infile):
        return []
    
    fin = open(infile,'r')
    array = []
    for line in fin:
        for item in list(filter(None,line.replace("\n","").split(separator))):
            array+=[item]
    fin.close()
    
    return array

def array2file(array,txtfile,sep,append=False,addFinalEnter=False):
    if len(array)==0:
        return False
    
    if (not append) or (append and not os.path.isfile(txtfile)):
        fout = open(txtfile,'w')
    else:
        fout = open(txtfile,'a')
        
    fout.write(array[0])
    for i in range(1,len(array)):
        fout.write(sep+array[i])
    if addFinalEnter:
        fout.write("\n")
        
    fout.close()
    return True
        
def nlinesval(textFile,val):
    return 0 if not os.path.isfile(textFile) else int(mylib.systemOut("cat "+textFile+" | grep \""+str(val)+"\" | wc -l"))

def printFile(introline,filepath):
    if not os.path.isfile(filepath):
        return
    
    print(introline)
    try:
        fin = open(filepath,'r')
        for line in fin:
            print(line)
    except ValueError as err:
        mylib.printError(f'Error reading file: {err}')
    finally:
        fin.close()

def changeDelim(txt,old_new,rm_orig=False):
    os.system(f"mv {txt} {txt}_orig")
    
    fin = open(txt+"_orig",'r')
    fout = open(txt,'w')
    for line in fin:
        for old_delim,new_delim in old_new.items():
            line = line.replace(old_delim,new_delim)
        fout.write(line)
    fout.close()
    fin.close()
    
    if rm_orig:
        os.system(f"rm {txt}_orig")

def linesCols(textFile,delim=","):
    if not os.path.isfile(textFile):
        return [0, 0]
    
    ncols = int(mylib.systemOut("awk '{print NF;exit}' "+textFile))
    nlines = int(mylib.systemOut("cat "+textFile+" | wc -l"))
    return [nlines,ncols]

# Returns True if the two files have the same content, False otherwise
def compareFiles(file1,file2):
    if not os.path.isfile(file1) or not os.path.isfile(file2):
        return False
    
    with open(file1) as f1, open(file2) as f2:
        for x,y in zip(f1,f2):
            if x.replace("\n","")!=y.replace("\n",""):
                return False
    return True

# Retruns True if concatFile is the horizontal concatenation of each line in file1 and file2
def isConcat(concatFile,file1,file2,delim):
    if not os.path.isfile(file1) or not os.path.isfile(file2):
        return False
    
    n = 1
    ok = True
    with open(file1) as f1, open(file2) as f2, open(concatFile) as f3:
        for ln1,ln2,ln in zip(f1,f2,f3):
            ln1 = ln1.replace("\n","")
            if ln1.endswith(delim):
                ln1 = ln1[:-1]
                
            ln2 = ln2.replace("\n","")
            if ln2.endswith(delim):
                ln2 = ln2[:-1]
        
            ln = ln.replace("\n","")
            if ln.endswith(delim):
                ln = ln[:-1]
                
            A = ln.split(delim)
            B = ln1.split(delim)+ln2.split(delim)
            if (A!=B):
                mylib.printWarning(f'line{n}')
                for i in range(max(len(A),len(B))):
                    a = "" if i>=len(A) else A[i]
                    b = "" if i>=len(B) else B[i]
                    if a==b:
                        print(f'{i}. {a}=={b}')  
                    else: 
                        print(f'{i}. {a}!={b}')
                        break
                ok = False
                
            n+=1
            if not ok:
                break
    
    f1.close()
    f2.close()
    f3.close()        
    return ok

def concat(concatFile,file1,file2,delim):
    fout = open(concatFile,'w')
    
    with open(file1) as f1, open(file2) as f2:
        for ln1,ln2 in zip(f1,f2):
            ln1 = ln1.replace("\n","")
            if ln1.endswith(delim):
                ln1 = ln1[:-1]
                
            ln2 = ln2.replace("\n","")
            if ln2.endswith(delim):
                ln2 = ln2[:-1]
                
            fout.write(ln1+delim+ln2+"\n")
            
    f1.close()
    f2.close()
    fout.close()

#############################################################################
############################# Manipulating PDFs #############################
#############################################################################
def convertToPDF(img,rmorigs=False):
    if not img.endswith(".png"):
        mylib.printError(f"Wrong input format for img {img}")
        return ""
    
    im = Image.open(img)
    im = im.convert("RGB")
    new = img.replace(".png",".pdf")
    im.save(new)
    im.close()
    
    if rmorigs:
        os.remove(img)
    return new

def getNumPages(pdf):
    return PyPDF2.PdfFileReader(pdf).numPages

def mergeVerticalImgsPDF(imgsList,output,rmorigs):
    merger = PyPDF2.PdfFileMerger()
    
    for img in imgsList:
        if not img.endswith(".pdf"):
            img = convertToPDF(img,rmorigs)
            merger.append(img)
        else:
            merger.append(img)
            
        if rmorigs:
            os.remove(img)
               
    merger.write(output)
    merger.close()
    
def mergeVerticalImgs(imgsList,output,rmorigs):
    if output.endswith(".pdf"):
        mergeVerticalImgsPDF(imgsList,output,rmorigs)
    else:
        mergeVerticalImgsPNG(imgsList,output,rmorigs)

def openPDFfile(pathfile,screen):
    if os.path.isfile(pathfile):
        if screen:
            os.system(f"gio open {pathfile}")
        else:
            print(pathfile)
            input("[Enter]")
    
# dic[page] = title
# Pages start at zero!!
def addSimpleOutline(inpdf,outpdf,dic,rmorig):
    writer = PyPDF2.PdfFileWriter()
    reader = PyPDF2.PdfFileReader(open(inpdf,"rb"))
    npages = reader.getNumPages()
    
    for page in range(npages):
        writer.addPage(reader.getPage(page))
        if page in dic.keys():
            writer.addBookmark(dic[page],page)
    
    output = open(outpdf,"wb")
    writer.write(output)
    output.close()
    
    if rmorig:
        os.system("rm "+inpdf)
            
#############################################################################
############################# Manipulating PNGs #############################
#############################################################################    
def mergeVerticalImgsPNG(imgsList,output,rmorigs,space=0):
    imgs = [Image.open(i) for i in imgsList]
    min_img_width = min(i.width for i in imgs)
    
    # Re-size images if necessary
    total_height = 0
    for i, img in enumerate(imgs):
        # If the image is larger than the minimum width, resize it
        if img.width > min_img_width:
            imgs[i] = img.resize((min_img_width, int(img.height / img.width * min_img_width)), Image.ANTIALIAS)
        total_height += imgs[i].height
    total_height+=space*len(imgsList)
    img_merge = Image.new(imgs[0].mode, (min_img_width, total_height))
    
    # Concatenate vertically
    y = 0
    for img in imgs:
        img_merge.paste(img, (0, y))
        y += img.height+space
        img.close()
    
    img_merge.save(output)
    if rmorigs:
        for img in imgsList:
            os.remove(img)

def mergeHorizontalImgs(imgsList,output,rmorigs,space=0):
    imgs1 = [Image.open(i) for i in imgsList]
    min_img_height1 = min(i.height for i in imgs1)
    
    # Re-size images if necessary
    total_width1 = 0
    for i, img in enumerate(imgs1):
        # If the image is larger than the minimum height, resize it
        if img.height > min_img_height1:
            imgs1[i] = img.resize((min_img_height1, int(img.height / img.width * min_img_height1)), Image.ANTIALIAS)
        total_width1 += imgs1[i].width
    total_width1+=space*len(imgsList)
    img_merge1 = Image.new(imgs1[0].mode, (total_width1, min_img_height1))
    
    # Concatenate horizontally
    x = 0
    for img in imgs1:
        img_merge1.paste(img, (x, 0))
        x += img.width+space
        img.close()
    
    img_merge1.save(output)
    if rmorigs:
        for img in imgsList:
            os.remove(img)

def resizeImg(img,out,xdim,ydim):
    image = Image.open(img)
    new_image = image.resize((xdim,ydim))
    new_image.save(out)
    
def addTitle(inpng,outpng,text,rmorigs,fontname="BodoniFLF-Bold.ttf",fontsize=20):
    img = Image.open(inpng)
    I1 = ImageDraw.Draw(img)
    font = ImageFont.truetype("/home/mkeith/Insync/DB/fonts/"+fontname,size=fontsize)
    I1.text((10,10),text,fill=(255,0,0),font=font)
    img.save(outpng)
    if rmorigs:
        os.system("rm "+inpng)
        
def addTitleTop(inpng,outpng,text,rmorigs,fontname="BodoniFLF-Bold.ttf",fontsize=20):
    old_img = Image.open(inpng)
    new_img = Image.new(old_img.mode,(old_img.width,old_img.height+fontsize+15))
    new_img.paste(old_img,(0,fontsize+15))
    new_img.save(outpng)
    addTitle(outpng,outpng,text,rmorigs,fontname,fontsize)
        
def addBorder(img,imgout,rmorigs):
    virat_img = cv2.imread(img)
    borderoutput = cv2.copyMakeBorder(virat_img, 10, 10, 10, 10, cv2.BORDER_CONSTANT, value=[0, 0, 255])
    cv2.imwrite(imgout, borderoutput)
    if rmorigs:
        os.system("rm "+img)

#############################################################################
######################### Manipulating other files ##########################
############################################################################# 
def getFileExt(filePath):
    return ".nii.gz" if filePath.endswith(".nii.gz") else pathlib.Path(filePath).suffix

def mkDir(dirpath):
    if dirpath.endswith("/"):
        dirpath = dirpath[:-1]
    
    if os.path.isdir(dirpath) and len(glob.glob(dirpath+"/*"))==0:
        return
    
    if os.path.isdir(dirpath):
        dest = dirpath+"+"
        while os.path.isdir(dest):
            dest+="+"
        os.system("mv "+dirpath+" "+dest)
        
    os.system("mkdir "+dirpath)
