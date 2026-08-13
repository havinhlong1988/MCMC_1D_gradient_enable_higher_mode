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
project_dir="CHT"
#
cd SetupData/
python main.py $project_dir
cd ..
echo "Setup data finish!"

