#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Monica Keith"
__status__ = "Production"

import tkinter as tk
from tkinter import ttk, messagebox
import mylib
from tkinter import filedialog
import pandas as pd
import os
import numpy as np
import math

DB_USER = ""
DB_PASS = ""
DB_HOST = ""
popup = None
user_field = None
passwd_field = None
host_field = None

test_mode = False
emergency = False

manual_textArea = None
auto_textArea = None
auto_button = None

menu_loc = None
radio_loc = None
loc_list = []
menu_db = None
radio_db = None
db_list = []
menu_proj = None
radio_proj = None
proj_list = []
menu_pipe = None
radio_pipe = None
pipe_list = []

excel_menu = None
grp_vars = []
menu_grp = None
grp_list = []
sel_grps = []
site_vars = []
menu_site = None
site_list = []
sel_site = []
res_vars = []
menu_res_array = []
sel_res = []
demo_vars = []
menu_demo_array = []
sel_demo = []

def connect():
    global DB_USER, DB_PASS, DB_HOST
    
    DB_USER = user_field.get()
    DB_PASS = passwd_field.get()
    DB_HOST = host_field.get()
    
    try:
        cnx = mylib.connect(DB_USER, DB_PASS, DB_HOST)
        connected = cnx.is_connected()
        cnx.close()
    except:
        connected = False
        messagebox.showerror("Error","Failed authentication")
    
    if connected:
        popup.destroy()
    else:
        DB_USER = ""
        DB_PASS = ""
        DB_HOST = ""

def toggleAuto():
    if auto_button!=None:
        auto_button.config(text="Stop") if auto_button.cget('text')=="Start" else auto_button.config(text="Start")

def manualTA_append(text):
    if manual_textArea!=None:
        manual_textArea.insert(tk.END, f"{text}")
        
def manualTA_insert(text):
    if manual_textArea!=None:
        manual_textArea.delete("1.0", tk.END)
        manual_textArea.insert(tk.END, f"{text}")
        
def autoTA_append(text):
    if auto_textArea!=None:
        auto_textArea.insert(tk.END, f"{text}")
        
def autoTA_insert(text):
    if auto_textArea!=None:
        auto_textArea.delete("1.0", tk.END)
        auto_textArea.insert(tk.END, f"{text}")
        
def initString():
    if radio_db!=None and len(db_list)>0:
        database = db_list[int(radio_db.get())-1]
    else:
        database = ""
        
    if radio_proj!=None and len(proj_list)>0:
        project = proj_list[int(radio_proj.get())-1]
    else:
        project = ""
        
    if radio_pipe!=None and len(pipe_list)>0:
        pipeline = pipe_list[int(radio_pipe.get())-1]
    else:
        pipeline = ""
        
    if radio_loc!=None and len(loc_list)>0:
        location = loc_list[int(radio_loc.get())-1]
    else:
        location = ""
        
    return f"Current database: {database}"\
           f"\nCurrent project: {project}"\
           f"\nCurrent pipeline: {pipeline}"\
           f"\nCurrent location: {location}"\
           f"\nTest mode: {test_mode}"\
           f"\nEmergency status: {emergency}"

def resetProject():
    global menu_proj, proj_list, radio_proj
    
    if menu_proj==None:
        return
    
    for proj in proj_list:
        menu_proj.delete(proj)
            
    if radio_db!=None and len(db_list)>0:
        database = db_list[int(radio_db.get())-1] 
        cnx = mylib.connect(DB_USER, DB_PASS, DB_HOST)
        cursor = cnx.cursor()
        proj_list = mylib.getColumnFromTable(cursor, "title", database+".projects")
        cursor.close()
        cnx.close()
               
        radio_proj = tk.StringVar()
        i = 1
        for proj in proj_list:
            menu_proj.add_radiobutton(label=proj, variable=radio_proj, value=i, command=changeProject)
            i+=1
        if len(proj_list)>0:
            radio_proj.set(1)
            
def resetPipeline():
    global menu_pipe, pipe_list, radio_pipe
        
    if menu_pipe==None:
        return
    
    for pipe in pipe_list:
        menu_pipe.delete(pipe)
            
    if radio_db!=None and len(db_list)>0:
        database = db_list[int(radio_db.get())-1]
        cnx = mylib.connect(DB_USER, DB_PASS, DB_HOST)
        cursor = cnx.cursor()
        pipe_list = mylib.getColumnFromTable(cursor, "pipename", database+".pipelines")
        cursor.close()
        cnx.close()
        
        radio_pipe = tk.StringVar()
        i = 1
        for pipe in pipe_list:
            menu_pipe.add_radiobutton(label=pipe, variable=radio_pipe, value=i, command=changePipeline)
            i+=1
        if len(pipe_list)>0:
            radio_pipe.set(1)
    
def resetGroups():
    global menu_grp, grp_vars, grp_list
    
    if menu_grp==None:
        return
    
    for grp in grp_list:
        menu_grp.delete(grp)
    grp_vars = []
    
    if radio_db!=None and len(db_list)>0:
        database = db_list[int(radio_db.get())-1]
        cnx = mylib.connect(DB_USER, DB_PASS, DB_HOST)
        cursor = cnx.cursor()
        grp_list = mylib.getColumnFromTable(cursor, "grp", database+".subjects", distinct=True)
        grp_list = [grp for grp in grp_list if grp!=None]
        cursor.close()
        cnx.close()
        
        for i in range(len(grp_list)):
            grp_vars+=[tk.StringVar()]
            menu_grp.add_checkbutton(label=grp_list[i], variable=grp_vars[i], onvalue=f"{grp_list[i]}1", offvalue=f"{grp_list[i]}0", command=changeGroups)
            
def resetSites():
    global menu_site, site_vars, site_list
    
    if menu_site==None:
        return
    
    for site in site_list:
        menu_site.delete(site)
    site_vars = []
    
    if radio_db!=None and len(db_list)>0:
        database = db_list[int(radio_db.get())-1]
        cnx = mylib.connect(DB_USER, DB_PASS, DB_HOST)
        cursor = cnx.cursor()
        site_list = mylib.getColumnFromTable(cursor, "site", database+".subjects", distinct=True)
        site_list = [site for site in site_list if site!=None]
        cursor.close()
        cnx.close()
        
        for i in range(len(site_list)):
            site_vars+=[tk.StringVar()]
            menu_site.add_checkbutton(label=site_list[i], variable=site_vars[i], onvalue=f"{site_list[i]}1", offvalue=f"{site_list[i]}0", command=changeSites)
            
def resetResults():
    global res_vars, menu_res_array
    
    for i in range(len(menu_res_array)):
        excel_menu.delete(f"Results{i+1}")
        
    res_vars = []
    menu_res_array = []
    
    if radio_pipe!=None and len(pipe_list)>0:
        pipeline = pipe_list[int(radio_pipe.get())-1]
    else:
        pipeline = ""
    
    if radio_db!=None and len(db_list)>0:
        database = db_list[int(radio_db.get())-1]
        cnx = mylib.connect(DB_USER, DB_PASS, DB_HOST, database)
        cursor = cnx.cursor()
        sess_list = mylib.getSessList(cursor,pipeline)
        if len(sess_list)>0:
            string_array = "("+','.join([f"'{sess}'" for sess in sess_list])+")"
            array = mylib.getColumnFromTable(cursor, "result", "results", ["sess"], [string_array], ["result"], True)
        else:
            array = []
        cursor.close()
        cnx.close()
        
        if len(array)>0:
            n_arrays = math.ceil(len(array)/10)
            arrays = np.array_split(array, n_arrays)
            
            j = 0
            for i in range(len(arrays)):
                menu_res = tk.Menu()
                menu_res_array+=[menu_res]
                excel_menu.add_cascade(menu=menu_res, label=f"Results{i+1}")
                for item in list(arrays[i]):
                    res_vars+=[tk.StringVar()]
                    menu_res.add_checkbutton(label=item, variable=res_vars[j], onvalue=f"{item}1", offvalue=f"{item}0", command=changeResults)
                    j+=1
                    
def resetDemographic():
    global demo_vars, menu_demo_array
    
    for i in range(len(menu_demo_array)):
        excel_menu.delete(f"Demographic{i+1}")
        
    demo_vars = []
    menu_demo_array = []
    
    if radio_pipe!=None and len(pipe_list)>0:
        pipeline = pipe_list[int(radio_pipe.get())-1]
    else:
        pipeline = ""
    
    if radio_db!=None and len(db_list)>0:
        database = db_list[int(radio_db.get())-1]
        cnx = mylib.connect(DB_USER, DB_PASS, DB_HOST, database)
        cursor = cnx.cursor()
        sbj_list = []
        for sess in mylib.getSessList(cursor,pipeline):
            sbj_list+=[mylib.getSbjID(cursor,sess)]
        if len(sbj_list)>0:
            string_array = "("+','.join([f"'{sbj}'" for sbj in sbj_list])+")"
            array = mylib.getColumnFromTable(cursor, "measure", "demographicData", ["sbjID"], [string_array], ["measure"], True)
        else:
            array = []
        cursor.close()
        cnx.close()
        
        if len(array)>0:
            n_arrays = math.ceil(len(array)/10)
            arrays = np.array_split(array, n_arrays)
            
            j = 0
            for i in range(len(arrays)):
                menu_demo = tk.Menu()
                menu_demo_array+=[menu_demo]
                excel_menu.add_cascade(menu=menu_demo, label=f"Demographic{i+1}")
                for item in list(arrays[i]):
                    demo_vars+=[tk.StringVar()]
                    menu_demo.add_checkbutton(label=item, variable=demo_vars[j], onvalue=f"{item}1", offvalue=f"{item}0", command=changeDemo)
                    j+=1
            
def resetLocation():
    global menu_loc, loc_list, radio_loc, DB_USER, DB_PASS
    
    if menu_loc!=None:
        for loc in loc_list:
            menu_loc.delete(loc)
    
    radio_loc = tk.StringVar()
    
    cnx = mylib.connect(DB_USER, DB_PASS)
    cursor = cnx.cursor()
    loc_list = mylib.getColumnFromTable(cursor,"location","generalInfo.locations")
    cursor.close()
    cnx.close()
    
    i = 1
    for loc in loc_list:
        menu_loc.add_radiobutton(label=loc, variable=radio_loc, value=i, command=changeLocation)
        i+=1
    if len(loc_list)>0:
        radio_loc.set(1)
        
def resetDB():
    global menu_db, db_list, radio_db, DB_USER, DB_PASS
    
    if menu_db!=None:
        for db in db_list:
            menu_db.delete(db)
            
    radio_db = tk.StringVar()
    
    cnx = mylib.connect(DB_USER, DB_PASS)
    cursor = cnx.cursor()
    db_list = mylib.getDatabases(cursor)
    cursor.close()
    cnx.close()
    
    i = 1
    for db in db_list:
        menu_db.add_radiobutton(label=db, variable=radio_db, value=i, command=changeDB)
        i+=1
    if len(db_list)>0:
        radio_db.set(1)

def changeDB():
    resetPipeline()
    resetProject()
    resetGroups()
    resetSites()
    resetResults()
    manualTA_insert(initString())
    autoTA_insert(initString())

def changeLocation():
    if radio_loc!=None and len(loc_list)>0:
        location = loc_list[int(radio_loc.get())-1]
        manualTA_append(f"\n\nCurrent location: {location}")
        autoTA_append(f"\n\nCurrent location: {location}")
        
def changeProject():
    if radio_proj!=None and len(proj_list)>0:
        project = proj_list[int(radio_proj.get())-1]
        manualTA_append(f"\n\nCurrent project: {project}")
        autoTA_append(f"\n\nCurrent project: {project}")
        
def changePipeline():
    if radio_pipe!=None and len(pipe_list)>0:
        resetResults()
        pipeline = pipe_list[int(radio_pipe.get())-1]
        manualTA_append(f"\n\nCurrent pipeline: {pipeline}")
        autoTA_append(f"\n\nCurrent pipeline: {pipeline}")

def changeGroups():
    global grp_vars, sel_grps
    
    sel_grps = []
    for var in grp_vars:
        val = var.get()
        if str(val).endswith("1"):
            sel_grps+=[val[:-1]]
            
def changeSites():
    global site_vars, sel_site
    
    sel_site = []
    for var in site_vars:
        val = var.get()
        if str(val).endswith("1"):
            sel_site+=[val[:-1]]

def changeResults():
    global sel_res
    
    sel_res = []
    for var in res_vars:
        val = var.get()
        if str(val).endswith("1"):
            sel_res+=[val[:-1]]
            
def changeDemo():
    global sel_demo
    
    sel_demo = []
    for var in demo_vars:
        val = var.get()
        if str(val).endswith("1"):
            sel_demo+=[val[:-1]]
            
def changeTestMode(status):
    global test_mode
    
    test_mode = True if int(status)==1 else False
    manualTA_append(f"\n\nTest mode: {test_mode}")
    
def changeEmergency(status):
    global emergency
    
    emergency = True if int(status)==1 else False
    manualTA_append(f"\n\nEmergency: {emergency}")

def checkReady():
    db_name = db_list[int(radio_db.get())-1] if radio_db!=None and len(db_list)>0 else ""
    pipeline = pipe_list[int(radio_pipe.get())-1]
    location = loc_list[int(radio_loc.get())-1]
    
    mylib.checkReadyProcs(DB_USER, DB_PASS, DB_HOST, db_name, pipeline, True, location, emergency=emergency, priority=[], absolute=[], test_mode=test_mode, force_manual=True)

def checkRunning():
    db_name = db_list[int(radio_db.get())-1] if radio_db!=None and len(db_list)>0 else ""
    pipeline = pipe_list[int(radio_pipe.get())-1]
    location = loc_list[int(radio_loc.get())-1]
    
    mylib.RunningProcs(DB_USER, DB_PASS, DB_HOST, db_name, pipeline, True, location, test_mode=test_mode)
    
def checkDone():
    db_name = db_list[int(radio_db.get())-1] if radio_db!=None and len(db_list)>0 else ""
    pipeline = pipe_list[int(radio_pipe.get())-1]
    location = loc_list[int(radio_loc.get())-1]
    
    mylib.checkDone(DB_USER, DB_PASS, DB_HOST, db_name, pipeline, location)
    
def checkSteps():
    db_name = db_list[int(radio_db.get())-1] if radio_db!=None and len(db_list)>0 else ""
    pipeline = pipe_list[int(radio_pipe.get())-1]
    
    mylib.getStepsInfo(DB_USER, DB_PASS, DB_HOST, db_name, pipeline, True)
    
def generateExcel():
    if radio_pipe==None or len(pipe_list)==0:
        return
    
    pipe = pipe_list[int(radio_pipe.get())-1]
    db_name = db_list[int(radio_db.get())-1] if radio_db!=None and len(db_list)>0 else ""
    sbjcols = []
    incl = {}
    if len(sel_grps)>0:
        sbjcols+=["grp"]
        incl["grp"] = sel_grps
    elif len(grp_list)>0:
        sbjcols+=["grp"]
    if len(sel_site)>0:
        sbjcols+=["site"]
        incl["site"] = sel_site
    elif len(site_list)>0:
        sbjcols+=["site"]
    
    cnx = mylib.connect(DB_USER, DB_PASS, DB_HOST, db_name)
    cursor = cnx.cursor()
    df = mylib.getSummaryTable(cursor, pipe, sbjcols, results_meas=sel_res, demographic_vars=sel_demo, include=incl)
    cursor.close()
    cnx.close()
    
    file_name = filedialog.asksaveasfilename(initialfile='database.xlsx',defaultextension=".xlsx",filetypes=[("Excel","*.xlsx")])
    if str(file_name)=="()" or str(file_name)=="":
        return
    if not os.path.isfile(file_name):
        writer = pd.ExcelWriter(file_name, engine='xlsxwriter')
        df.to_excel(writer, sheet_name=pipe)
        writer.close()
    else:
        with pd.ExcelWriter(
            file_name,
            mode="a",
            engine="openpyxl",
            if_sheet_exists="new",
        ) as writer:df.to_excel(writer, sheet_name=pipe)
    
def main():
    global manual_textArea, auto_textArea, auto_button, menu_db, menu_proj, menu_pipe, menu_loc, popup, user_field, passwd_field, host_field, DB_USER, DB_PASS, DB_HOST, menu_grp, menu_site, excel_menu
    
    # Popup to enter username & password for the DB
    popup = tk.Tk()
    popup.title("Login")
    popup.geometry("230x140")
    
    popup_line1 = ttk.Frame(popup)
    ttk.Label(popup_line1,text="Username:", font=('Helvetica',10), width=8).pack(side="left", pady=10, padx=5)
    user_field = ttk.Entry(popup_line1, width=20)
    user_field.pack(side="left", padx=5)
    popup_line1.pack(fill=tk.BOTH)
    
    popup_line2 = ttk.Frame(popup)
    ttk.Label(popup_line2,text="Password:", font=('Helvetica',10), width=8).pack(side="left", pady=10, padx=5)
    passwd_field = ttk.Entry(popup_line2, show="*", width=20)
    passwd_field.pack(side="left", padx=5)
    popup_line2.pack(fill=tk.BOTH)
    
    popup_line3 = ttk.Frame(popup)
    ttk.Label(popup_line3, text="Host:", font=('Helvetica',10), width=8).pack(side="left", pady=10, padx=5)
    host_field = ttk.Entry(popup_line3, width=20)
    host_field.pack(side="left", padx=5)
    host_field.insert(0,"localhost")
    popup_line3.pack(fill=tk.BOTH)
    
    popup_line4 = ttk.Frame(popup)
    ttk.Button(popup_line4, text="Connect", command=connect).pack(side="left", padx=5)
    popup_line4.pack(fill=tk.BOTH)
    
    popup.mainloop()
    
    if DB_USER=="" or DB_PASS=="" or DB_HOST=="":
        messagebox.showerror("Error","Failed authentication")
        exit()
    
    # Menu
    root = tk.Tk()
    root.title("MKautomation")
    menubar = tk.Menu(root)
    
    # General menu
    general_menu = tk.Menu(menubar)    
    menubar.add_cascade(menu=general_menu, label="General")
    
    menu_db = tk.Menu(general_menu)
    general_menu.add_cascade(menu=menu_db, label="Database")
    resetDB()
    
    menu_loc = tk.Menu(general_menu)
    general_menu.add_cascade(menu=menu_loc, label="Location")
    resetLocation()
    
    check1 = tk.StringVar()
    general_menu.add_checkbutton(label="Test mode", variable=check1, onvalue=1, offvalue=0, command=lambda: changeTestMode(check1.get()))
    
    check2 = tk.StringVar()
    general_menu.add_checkbutton(label="Emergency", variable=check2, onvalue=1, offvalue=0, command=lambda: changeEmergency(check2.get()))
    
    # Database menu
    dabase_menu = tk.Menu(menubar)
    menubar.add_cascade(menu=dabase_menu, label="Database")
    
    menu_proj = tk.Menu()
    dabase_menu.add_cascade(menu=menu_proj, label="Project")
    resetProject()
    
    menu_pipe = tk.Menu(dabase_menu)
    dabase_menu.add_cascade(menu=menu_pipe, label="Pipeline")
    resetPipeline()
    
    # Excel menu
    excel_menu = tk.Menu(menubar)
    menubar.add_cascade(menu=excel_menu, label="Excel")
    
    menu_grp = tk.Menu()
    excel_menu.add_cascade(menu=menu_grp, label="Groups")
    resetGroups()
    
    menu_site = tk.Menu()
    excel_menu.add_cascade(menu=menu_site, label="Sites")
    resetSites()
    
    resetResults()
    
    resetDemographic()
    
    # Manual execution
    nb = ttk.Notebook(root)
    manual_frame = ttk.Frame(nb)
        
    manual_topFrame = ttk.Frame(manual_frame)
    v = tk.Scrollbar(manual_topFrame, orient='vertical')
    v.pack(side=tk.RIGHT, fill='y')
    manual_textArea = tk.Text(manual_topFrame, height=55, width=95, yscrollcommand=v.set)
    v.config(command=manual_textArea.yview)
    manual_textArea.pack()
    manualTA_append(initString())
    manual_topFrame.pack(fill=tk.BOTH, expand=True)
        
    manual_bottomFrame = ttk.Frame(manual_frame)
    ttk.Button(manual_bottomFrame, text="Check ready", command=checkReady).pack(side="left", padx=5)
    ttk.Button(manual_bottomFrame, text="Check running", command=checkRunning).pack(side="left", padx=5)
    ttk.Button(manual_bottomFrame, text="Check done", command=checkDone).pack(side="left", padx=5)
    ttk.Button(manual_bottomFrame, text="Steps status", command=checkSteps).pack(side="left", padx=5)
    ttk.Button(manual_bottomFrame, text="Generate excel", command=generateExcel).pack(side="left", padx=5)
    manual_bottomFrame.pack(side="bottom", expand=True, padx=5, pady=5)
        
    manual_frame.pack(fill=tk.BOTH, expand=True)
    nb.add(manual_frame, text="Manual")
    
    # Automatic execution
    auto_frame = ttk.Frame(nb)
        
    auto_topFrame = ttk.Frame(auto_frame)
    auto_textArea = tk.Text(auto_topFrame, height=55, width=95)
    auto_textArea.pack()
    autoTA_append(initString())
    auto_topFrame.pack(fill=tk.BOTH, expand=True)
        
    auto_bottomFrame = ttk.Frame(auto_frame)
    auto_button = ttk.Button(auto_bottomFrame, text="Start", command=toggleAuto)
    auto_button.pack(side="left", padx=5)
    auto_bottomFrame.pack(side="bottom", expand=True, padx=5, pady=5)
        
    auto_frame.pack(fill=tk.BOTH, expand=True)
    nb.add(auto_frame, text="Automatic")
    
    nb.pack(padx=5, pady=5, expand=True)
    root.config(menu=menubar)
    root.mainloop()
    
if __name__ == "__main__":
    main()

