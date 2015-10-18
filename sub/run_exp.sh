#!/usr/bin/bash

# parameters for fms submission script
run_sub_folder=../../../sub
run_folder=../exp;     			# experiment folder
run_file=run_exp; 			# default run file name
in_file='list_exp';			# input list filename
out_file='list_exp_submitted';		# output list filename
wt=0; 					# waiting time between experiments

# loop over experiments in $in_file
while read line
do
	now=$(date +'%m-%d-%y-%T');		# set current time
	echo $now : $line >> $out_file; 	# write current instance of experiment and time to output file
	
	# create arrays with experiment namelist parameters and values for these parameters
	vars=(${line// / }); 			# get all experiment name and parameters from $line in $in_file
	exp_name=(${vars[0]}); 			# save name of experiment from $line
	vars=(${vars[@]:1}); 			# save only parameters and parameter values from $line in $vars

	# creat name for this particular instance of this experiment
	inst_name=$exp_name;
	for ((i=0; i<${#vars[@]}; i++))
	do
		inst_name=$inst_name'_'${vars[$i]};
	done
		
	# save namelist parameters and parameter values from $line in separate arrays
	for ((i=0; i<${#vars[@]}; i+=2))
        do
		pars[$(( $i / 2 ))]=${vars[$i]};
		vals[$(( $i / 2 ))]=${vars[$(( $i + 1 ))]};
	done
  
	# create a copy of the run file in the fms experiment folder with the standard parameters
	cp $run_folder/$exp_name/run/$run_file $run_folder/$exp_name/run/run_$inst_name;

	# replace the standard parameters in the run file copy by the modified parameters from $line
	for ((i=0; i<${#pars[@]}; i++))
	do
		sed "s/.*${pars[$i]}.*/${pars[$i]} = ${vals[$i]},/" $run_folder/$exp_name/run/run_$inst_name >> $run_folder/$exp_name/run/run_tmp;
	done
	mv $run_folder/$exp_name/run/run_tmp $run_folder/$exp_name/run/run_$inst_name;

	# replace default name in standard run file with name of this instance of the experiment $line
	sed "s/RUN_NAME/$inst_name/" $run_folder/$exp_name/run/run_$inst_name > $run_folder/$exp_name/run/run_tmp;
	mv $run_folder/$exp_name/run/run_tmp $run_folder/$exp_name/run/run_$inst_name;

	# submit modified run file
	cd $run_folder/$exp_name/run;
	bsub < run_$inst_name;
	cd $run_sub_folder;	

	# wait $wt seconds before repeat in order to avoid jamming the cluster when copying data (might lead to data loss)
	sleep $wt;
done < $in_file;
