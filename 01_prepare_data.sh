#!/bin/bash
echo " 1 .Are you certain you wish to proceed with execution? 

 2. Have you reviewed and made changes to the Setupdata/parameters.py file?

 3. Have you correct the project directory on this script?
"
echo ""
echo " If no -> Press Ctrl + C to abort"
echo ""
read -p "Press any key to continue... " -n1 -s
#
# Project directory: pass it as the first argument, e.g. "bash 01_prepare_data.sh CHT".
# Defaults to FM (the active project); pass CHT to go back to the old one.
project_dir="${1:-FM}"
#
cd SetupData/
python main.py $project_dir
cd ..
echo "Setup data finish!"

