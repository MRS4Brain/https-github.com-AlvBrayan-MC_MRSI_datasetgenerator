#!/bin/bash 
if [ $# -ne 9 ]; then
    echo -e "Please give as argument:\n the name of the metabolite Brucker raw file ,\n the name of the Water Brucker raw file, \n the name of the matlab file containing the acquisition parameters (e.g.AcqParam.m ),\n the name given to the resulting dataset, \n the TGV regularization parameter value (1E-4 - 1E-3),\n the tolerence percentage for lipid suppression (goal of brain/skull signal ratio,range: 0.4-1, 0 to skip the lipid supp.),\n the number of components for water suppression (3T:16, 7T:32, 0 skip the water supp.),\n the number of components for the low-rank (2D: 20, 3D: 40) ,\n the simulated undersampling ratio (1 for none, 0.5 for x2 acceleration, 0.25 x4)."
else
current_dir=`dirname $0`
matlab_exec=matlab

TWIX=$1
WaterTWIX=$2
Param_File=$3
Data_name=$4
mu_TV=$5
Nbasis=$6
WSComp=$7
Comp=$8
Undersampling=$9

rm matlab_command.m
echo "addpath(genpath('${current_dir}/scripts'))" >> matlab_command.m
echo "MRSI_CSSENSELR_Recon('${TWIX}','${WaterTWIX}','${Param_File}','${Data_name}',${mu_TV},${Nbasis},${WSComp},${Comp},${Undersampling})" >> matlab_command.m
${matlab_exec} -nodisplay -nosplash < matlab_command.m  &> debug_${Data_name}.log
rm matlab_command.m

fi
