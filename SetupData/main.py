# Import the parameters module
import sys, os
import pandas as pd
import numpy as np
# ------------- Self functions-----------------
import parameters
from src import setup_functions, setup_bspline, setup_layercake
#####################################################################################
def main():
    # Access the project directory from input
    if len(sys.argv[:]) < 1:
        print("Usage: python main.py {project_dir}")
        sys.exit(0)
    project_dir = str(sys.argv[1])
    print("Setup data for project directory: ",project_dir)
    # Now access model information and data directory
    pwd = os.getcwd()
    #
    main_dir = os.path.dirname(pwd)
    # main project directory
    full_project_dir = os.path.join(main_dir, project_dir)
    # print(full_project_dir)
    # station list file name
    stafile = os.path.join(full_project_dir, parameters.stafilename); 
    # print(stafile)
    # Data directory
    Vel_dir = os.path.join(full_project_dir,"Vel_mod"); 
    print(Vel_dir) 
    # data directory
    data_dir = os.path.join(full_project_dir , 'data'); 
    # script direcotry
    src_dir = os.path.join(main_dir, 'SetupData')
    src_dir = os.path.join(src_dir, 'src')
    #----------------------------------------------------------------------------------------------------
    # ----- Read the station list to build the model and data for each station ----------------
    sta_data = pd.read_csv(stafile,delim_whitespace=True, na_values=0,names=["name","long","lat","elv"])
    # sta_data['name'] = sta_data['name'].apply( lambda x: f"{x:05d}"); 
    # sta_data['name'] = sta_data['name'].astype(str)
   # force to string safely (handles numbers, NaN)
    sta_data["name"] = sta_data["name"].apply(lambda x: "" if pd.isna(x) else str(int(x)) if str(x).strip().replace(".0","").isdigit() else str(x).strip())
    # pad to 5 digits
    sta_data["name"] = sta_data["name"].astype(str).str.zfill(5)
    #
    # make the directory contained all check files
    query_dir = os.path.join(full_project_dir, 'query_data')
    
    #
    if not os.path.isdir(query_dir):
        os.mkdir(query_dir)
    # create the file for all excutable stations:
    if os.path.isfile(os.path.join(full_project_dir, 'runstations.lst')):
        os.remove(os.path.join(full_project_dir, 'runstations.lst'))
    # The MonteCarlo directory
    in_sub = os.path.join(full_project_dir, "MonteCarlo/indata")
    # Now loop for all station:
    for i,stanm in enumerate(sta_data["name"]):
      print("Processing the station: ",i,stanm)
      query_sub = os.path.join(query_dir, stanm)
      data_sub = os.path.join(data_dir, f'{stanm}_data')
      #
      if not os.path.isdir(query_sub):
        print("create: ",query_sub)
        os.mkdir(query_sub)
    # ----------- call function to settup the model --------------------------
      if (parameters.mod_type_flag == 1):
        print("> setup_layercake");
        setup_layercake.setup(sta_data.iloc[i],Vel_dir,src_dir,query_sub,data_sub,in_sub) 
      elif (parameters.mod_type_flag == 2):
        print("> setup_bspline"); 
        setup_bspline.setup(sta_data.iloc[i],Vel_dir,src_dir,query_sub,data_sub,in_sub)
      else:
          print("Wrong model type! Exit!")
          sys.exit(0)


        





    

if __name__ == "__main__":
    main()