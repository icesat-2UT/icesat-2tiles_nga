# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 21:37:08 2021

@author: malonzo
"""

# Import functions
from mt_utils import (handle_user_input,
                      handle_home_option, handle_exit_option,
                      handle_help_option, handle_list_option, 
                      handle_load_option, handle_process_option, 
                      handle_filter_option, handle_compute_option,
                      handle_plot_option, handle_export_option,
                      GridStruct)

# Show main menu
handle_home_option()

# Initialize parameters
userInpKey = 'none'
df_all = []
rasterData = GridStruct([],[],[],[])
pptkViewers = []

# Keep program running until 'exit' or 'quit' or 'stop' is entered
exitWords = ['exit','quit','stop']
while(userInpKey not in exitWords):
    
    # Wait for user input
    userInp = input('>> ')

    # Handle user input
    userInpSplit, userInpKey = handle_user_input(userInp)
        
    # Take action based on user input key
    if(userInpKey=='home'):
        
        # Show header and home screens
        handle_home_option()
        
    elif(userInpKey=='help'):
        
        handle_help_option(userInpSplit)
        
    elif(userInpKey=='list'):
        
        # Show .inp files
        handle_list_option(userInpSplit)
        
    elif(userInpKey=='load'):
        
        # Read .gzip file
        df_all = handle_load_option(userInpSplit, df_all)
     
    elif(userInpKey=='process'):
        
        # Process .inp file
        df_all = handle_process_option(userInpSplit, df_all)
        
    elif(userInpKey=='filter'):
        
        # Filter processed data
        df_all = handle_filter_option(userInpSplit, df_all)
    
    elif(userInpKey=='compute'):
        
        # Filter processed data
        handle_compute_option(userInpSplit, df_all)
        
    elif(userInpKey=='plot'):
            
        rasterData, pptkViewers = handle_plot_option(userInpSplit, df_all, rasterData, pptkViewers)
        
    elif(userInpKey=='export'):
                
        handle_export_option(userInpSplit, df_all, rasterData)
        
    elif(userInpKey in exitWords):
        
        handle_exit_option(pptkViewers)
        
    else:
        
        print('Unknown command.\n')
        
    # endIf
    
# endWhile
