#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  1 10:58:00 2022

@author: ding
"""

from tkinter import *
from tkinter import ttk

# =============================================================================

root = Tk()
root.title("Setting")
root.geometry("500x720")
root.config(bg="skyblue")

# =============================================================================

class GUI:
    
    def __init__(self, ws):
        
        self.ws = ws
        self.font = ('Courier New','16')
        x_ini = 60
        x_des = 190
        y_ini = 60
        y_des = 0

# -----------------------------------------------------------------------------
        
        RL_Label = Label(ws, text="Agent",font=self.font)
        RL_Label.config(bg="skyblue")
        RL_Label.place(x=x_ini,y=y_ini+y_des,anchor="nw")
        
        self.RL_combo = ttk.Combobox(ws,width=15,font=self.font)
        self.RL_combo["values"] = ["Q-Learning","Deep Q Network"]
        self.RL_combo.current(0)
        self.RL_combo.place(x=x_ini+x_des,y=y_ini+y_des,anchor="nw")
        
# -----------------------------------------------------------------------------

        epoch_Label = Label(ws,text="Epoch",font=self.font)
        epoch_Label.config(bg="skyblue")
        epoch_Label.place(x=x_ini,y=y_ini+y_des+40,anchor="nw")
        
        self.epoch_entry = Entry(ws,font=self.font,width=16)
        self.epoch_entry.place(x=x_ini+x_des,y=y_ini+y_des+40,anchor="nw")

# -----------------------------------------------------------------------------

        Method_Label = Label(ws,text="Method",font=self.font)
        Method_Label.config(bg="skyblue")
        Method_Label.place(x=x_ini,y=y_ini+y_des+80,anchor="nw")
        
        self.Method_combo = ttk.Combobox(ws,width=15,font=self.font)
        if ws == CHAN_frame:
            self.Method_combo["values"] = ["ORI","JCP","ATM","MUSCL","THINCEM"]
        if ws == HAR_frame:
            self.Method_combo["values"] = ["Hybrid","ROE","AUSMD"]
        self.Method_combo.current(0)
        self.Method_combo.place(x=x_ini+x_des,y=y_ini+y_des+80,anchor="nw")

# -----------------------------------------------------------------------------

        ISAVE_Label = Label(ws,text="Isave",font=self.font)
        ISAVE_Label.config(bg="skyblue")
        ISAVE_Label.place(x=x_ini,y=y_ini+y_des+120,anchor="nw")
        
        self.ISAVE_entry = Entry(ws,font=self.font,width=16)
        self.ISAVE_entry.place(x=x_ini+x_des,y=y_ini+y_des+120,anchor="nw")

# -----------------------------------------------------------------------------

        ITMAX_Label = Label(ws,text="ITmax",font=self.font)
        ITMAX_Label.config(bg="skyblue")
        ITMAX_Label.place(x=x_ini,y=y_ini+y_des+160,anchor="nw")
        
        self.ITMAX_entry = Entry(ws,font=self.font,width=16)
        self.ITMAX_entry.place(x=x_ini+x_des,y=y_ini+y_des+160,anchor="nw")

# -----------------------------------------------------------------------------

        if ws == CHAN_frame:
            ngrd_Label = Label(ws, text="Mesh",font=self.font)
            ngrd_Label.config(bg="skyblue")
            ngrd_Label.place(x=x_ini,y=y_ini+y_des+200,anchor="nw")
            
            self.ngrd_combo = ttk.Combobox(ws,width=15,font=self.font)
            self.ngrd_combo["values"] = ["400x80","800x160"]
            self.ngrd_combo.current(0)
            self.ngrd_combo.place(x=x_ini+x_des,y=y_ini+y_des+200,anchor="nw")
            
# -----------------------------------------------------------------------------

        if ws == HAR_frame:
            cfl_Label = Label(ws, text="CFL",font=self.font)
            cfl_Label.config(bg="skyblue")
            cfl_Label.place(x=x_ini,y=y_ini+y_des+200,anchor="nw")
            
            self.cfl_entry = Entry(ws,font=self.font,width=16)
            self.cfl_entry.place(x=x_ini+x_des,y=y_ini+y_des+200,anchor="nw")

# -----------------------------------------------------------------------------

        alpha_ini_Label = Label(ws,text="Initial Alpha",font=self.font)
        alpha_ini_Label.config(bg="skyblue")
        alpha_ini_Label.place(x=x_ini,y=y_ini+y_des+240+20,anchor="nw")
        
        self.alpha_ini_entry = Entry(ws,font=self.font,width=16)
        self.alpha_ini_entry.place(x=x_ini+x_des,y=y_ini+y_des+240+20,anchor="nw")
        
        alpha_fin_Label = Label(ws,text="Final Alpha",font=self.font)
        alpha_fin_Label.config(bg="skyblue")
        alpha_fin_Label.place(x=x_ini,y=y_ini+y_des+280+20,anchor="nw")
        
        self.alpha_fin_entry = Entry(ws,font=self.font,width=16)
        self.alpha_fin_entry.place(x=x_ini+x_des,y=y_ini+y_des+280+20,anchor="nw")
        
        alpha_ran_Label = Label(ws,text="Alpha Range",font=self.font)
        alpha_ran_Label.config(bg="skyblue")
        alpha_ran_Label.place(x=x_ini,y=y_ini+y_des+320+20,anchor="nw")
        
        self.alpha_ran_entry = Entry(ws,font=self.font,width=16)
        self.alpha_ran_entry.place(x=x_ini+x_des,y=y_ini+y_des+320+20,anchor="nw")
        
# -----------------------------------------------------------------------------

        beta_ini_Label = Label(ws,text="Initial Beta",font=self.font)
        beta_ini_Label.config(bg="skyblue")
        beta_ini_Label.place(x=x_ini,y=y_ini+y_des+360+40,anchor="nw")
        
        self.beta_ini_entry = Entry(ws,font=self.font,width=16)
        self.beta_ini_entry.place(x=x_ini+x_des,y=y_ini+y_des+360+40,anchor="nw")
        
        beta_fin_Label = Label(ws,text="Final Beta",font=self.font)
        beta_fin_Label.config(bg="skyblue")
        beta_fin_Label.place(x=x_ini,y=y_ini+y_des+400+40,anchor="nw")
        
        self.beta_fin_entry = Entry(ws,font=self.font,width=16)
        self.beta_fin_entry.place(x=x_ini+x_des,y=y_ini+y_des+400+40,anchor="nw")
        
        beta_ran_Label = Label(ws,text="Beta Range",font=self.font)
        beta_ran_Label.config(bg="skyblue")
        beta_ran_Label.place(x=x_ini,y=y_ini+y_des+440+40,anchor="nw")
        
        self.beta_ran_entry = Entry(ws,font=self.font,width=16)
        self.beta_ran_entry.place(x=x_ini+x_des,y=y_ini+y_des+440+40,anchor="nw")

# -----------------------------------------------------------------------------
        
        Save_button = Button(ws,text="Save",command=self.save)
        Save_button.place(relx=0.25, y=y_ini+y_des+480+80, anchor="n")

# -----------------------------------------------------------------------------
        
        Default_button = Button(ws,text="Default",command=self.default)
        Default_button.place(relx=0.5, y=y_ini+y_des+480+80, anchor="n")

# -----------------------------------------------------------------------------
        
        Clear_button = Button(ws,text="Clear",command=self.clean)
        Clear_button.place(relx=0.75, y=y_ini+y_des+480+80, anchor="n")
        
# -----------------------------------------------------------------------------

    def save(self):
        
        save_data = []
        
        save_data.append(str(self.RL_combo.current()+1)+"\n")
        save_data.append(str(self.epoch_entry.get())+"\n")
        save_data.append(str(self.Method_combo.current()+1)+"\n")
        save_data.append(str(self.ISAVE_entry.get())+"\n")
        save_data.append(str(self.ITMAX_entry.get())+"\n")
        if self.ws == CHAN_frame:
            save_data.append(str(self.ngrd_combo.current()+1)+"\n")
            path = 'CHAN_setting.txt'
        if self.ws == HAR_frame:
            save_data.append(str(self.cfl_entry.get())+"\n")
            path = 'HAR_setting.txt'
        save_data.append(str(self.alpha_ini_entry.get())+"\n")
        save_data.append(str(self.alpha_fin_entry.get())+"\n")
        save_data.append(str(self.alpha_ran_entry.get())+"\n")
        save_data.append(str(self.beta_ini_entry.get())+"\n")
        save_data.append(str(self.beta_fin_entry.get())+"\n")
        save_data.append(str(self.beta_ran_entry.get())+"\n")
        
        f = open(path, 'w')
        f.writelines(save_data)
        f.close()
        
        self.message_box()
        
        return

# =============================================================================
    
    def message_box(self):
        
        message_font = ('Courier New','14')
        
        message = Toplevel()
        message.title("Message")
        message.geometry("400x120")
        message.config(bg="skyblue")
        
        message_Label = Label(message,text="New settings have been saved.",font=message_font)
        message_Label.config(bg="skyblue")
        message_Label.place(relx=0.1,rely=0.2,anchor="nw")
        
        close_message_button = Button(message, text="OK", command=message.destroy)
        close_message_button.place(relx=0.5,rely=0.6, anchor="n")
        
        return
    
# =============================================================================
    
    def default(self):
        
        self.clean()
        
        if self.ws == CHAN_frame:
        
            self.epoch_entry.insert(0, "500")
            self.ISAVE_entry.insert(0, "20")
            self.ITMAX_entry.insert(0, "8000")
            
            self.alpha_ini_entry.insert(0, "100")
            self.alpha_fin_entry.insert(0, "500")
            self.alpha_ran_entry.insert(0, "10")
            
            self.beta_ini_entry.insert(0, "1.2")
            self.beta_fin_entry.insert(0, "1.5")
            self.beta_ran_entry.insert(0, "0.01")
        
        if self.ws == HAR_frame:
            
            self.epoch_entry.insert(0, "100")
            self.ISAVE_entry.insert(0, "50")
            self.ITMAX_entry.insert(0, "9000")
            self.cfl_entry.insert(0, "1.9")
            
            self.alpha_ini_entry.insert(0, "0.9")
            self.alpha_fin_entry.insert(0, "0.99")
            self.alpha_ran_entry.insert(0, "0.01")
            
            self.beta_ini_entry.insert(0, "1.25")
            self.beta_fin_entry.insert(0, "1.60")
            self.beta_ran_entry.insert(0, "0.001")
            
        return
    
# =============================================================================
    
    def clean(self):
        
        if self.ws == CHAN_frame:
            
            self.epoch_entry.delete(0, END)
            self.ISAVE_entry.delete(0, END)
            self.ITMAX_entry.delete(0, END)
            
            self.alpha_ini_entry.delete(0, END)
            self.alpha_fin_entry.delete(0, END)
            self.alpha_ran_entry.delete(0, END)
            
            self.beta_ini_entry.delete(0, END)
            self.beta_fin_entry.delete(0, END)
            self.beta_ran_entry.delete(0, END)
            
        if self.ws == HAR_frame:
            
            self.epoch_entry.delete(0, END)
            self.ISAVE_entry.delete(0, END)
            self.ITMAX_entry.delete(0, END)
            self.cfl_entry.delete(0, END)
            
            self.alpha_ini_entry.delete(0, END)
            self.alpha_fin_entry.delete(0, END)
            self.alpha_ran_entry.delete(0, END)
            
            self.beta_ini_entry.delete(0, END)
            self.beta_fin_entry.delete(0, END)
            self.beta_ran_entry.delete(0, END)
        
        return
    
# =============================================================================

notebook = ttk.Notebook(root)
notebook.pack(pady=0)

CHAN_frame = Frame(notebook,takefocus=True,width=700,heigh=820,bg="skyblue")
HAR_frame = Frame(notebook,takefocus=True,width=700,heigh=820,bg="skyblue")

CHAN_frame.pack(fill="both",expand=True)
HAR_frame.pack(fill="both",expand=True)

notebook.add(CHAN_frame, text="2D detonation waves Case")
notebook.add(HAR_frame, text="Blunt-Body Flow Case")

CHAN = GUI(CHAN_frame)
HAR = GUI(HAR_frame)

# =============================================================================

root.mainloop()

