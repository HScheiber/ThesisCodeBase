# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 20:32:43 2020

@author: Hayden
"""
def loadbar_on():
    from tkinter import ttk, Tk, Button
    import time
    
    root = Tk()
    root.title('Loading trajectory...')
    root.geometry=("1000x1000")
    
    def step():
        #my_progress['value'] += 20
        my_progress.start(10)
    
    def stop():
        #my_progress['value'] += 20
        my_progress.stop()
    
    
    my_progress = ttk.Progressbar(root,orient='horizontal',
                                 length=300,mode='indeterminate')
    my_progress.pack(pady=40,padx=50)
    
    my_button = Button(root,text="Progress", command=step)
    my_button.pack(pady=20)
    
    my_button2 = Button(root,text="Stop", command=root.destroy)
    my_button2.pack(pady=20)
    
    
    root.mainloop()
    return root
    
def loadbar_off(object):
    object.destroy()
    
    
    
import PySimpleGUI as sg
pg.one_line_progress_meter('Loading...',
    0,
    100,
    key="OK for 1 meter",
    orientation="h",
    bar_color=(None, None),
    button_color=None,
    size=(20, 20),
    border_width=None,
    grab_anywhere=True,
    no_titlebar=False)

for i in range(1,100):
    sg.one_line_progress_meter('Calculating Order Parameters', i+1, 100, 'key',
                               'Calculating Order Parameters',no_titlebar=True,
                               orientation="h")
