#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Monica Keith"
__status__ = "Production"

import os
import mylib

def doMetricMath(array,expr,out):
    mvars = " ".join(["-var x"+str(i)+" "+array[i] for i in range(len(array))])
    cmd = "wb_command -metric-math "+expr+" "+out+" "+mvars
    print("\n"+cmd)
    os.system(cmd)
    print("done\n")

def doMath(array,expr,out):
    mvars = " ".join(["-var x"+str(i)+" "+array[i] for i in range(len(array))])
    cmd = "wb_command -cifti-math "+expr+" "+out+" "+mvars
    print("\n"+cmd)
    os.system(cmd)
    print("done\n")
    
def avgSurfs(array,avg_out):
    nsbj = len(array)
    avg_expr = "'("+"+".join(["x"+str(i) for i in range(nsbj)])+")/"+str(nsbj)+"'"
    doMath(array,avg_expr,avg_out)
    
def stdSurfs(array,std_out):
    nsbj = len(array)
    std_expr = "'sqrt(("+"+".join(["x"+str(i) for i in range(nsbj)])+")/"+str(nsbj)+")'"
    doMath(array,std_expr,std_out)
    
def getMaxValSurface(surface):
    MAX = mylib.systemOut("wb_command -cifti-stats "+surface+" -reduce MAX")
    return float(MAX)

def smooth(surface,output,left_mid,right_mid):
    print("\nSmoothing "+surface)
    os.system("wb_command -cifti-smoothing "+surface+" 5 5 COLUMN "+output+" -fwhm -left-surface "+left_mid+" -right-surface "+right_mid)
    print("done\n")

def project(vol_in,surf_out,surf_mid,surf_white,surf_pial):
    os.system("rm -rf "+surf_out)
    print("\nMapping "+vol_in+" into "+surf_out+"...")
    os.system("wb_command -volume-to-surface-mapping "+vol_in+" "+surf_mid+" "+surf_out+" -ribbon-constrained "+surf_white+" "+surf_pial)
    print("done\n")
    
def generateCifti(cifti,L_metric,R_metric):
    print("Generating "+cifti+"...")
    os.system("wb_command -cifti-create-dense-scalar "+cifti+" -left-metric "+L_metric+" -right-metric "+R_metric)
    print("done\n")
