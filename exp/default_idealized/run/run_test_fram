#!/bin/csh -f  
#PBS -l nodes=16  
#PBS -l walltime=04:00:00 
#PBS -N fms_default
#PBS -o out_err/out.$PBS_JOBID
#PBS -e out_err/err.$PBS_JOBID
#PBS -S /bin/csh    
#PBS -V 
#PBS -q default

# === Default test run script for idealized GCM ===

# See description at run_test.readme
 
# Ian Eisenman, Yohai Kaspi, Tim Merlis, November 2010
# Farid Ait Chaalal, Xavier Levine, Zhihong Tan, March 2012 
# Farid Ait Chaalal, September 2012
# Robb Wills, Ori Adam, May 2014

# change the working directory (default is home directory) 
cd $PBS_O_WORKDIR
echo Working directory is $cwd    
 

set model_type     = moist                          # if "moist", the moist model is run and if "dry, it is the dry. "moist_hydro" is for the bucket hydrology model. The namelists for the parameters are below (L212). 
set analysis_type  = 2d                             # choose type of analysis: 2d (zonally averaged) or 3d (zonally varying) outputs
set run_name       = default_test_moist_2d                   # label for run; output dir and working dir are run_name specific
set run_script     = $cwd/run_test_fram                  # path/name of this run script (for resubmit) 
set exp_home       = $cwd:h                         # directory containing run/$run_script and input/ 
set exp_name       = $exp_home:t                    # name of experiment (i.e., name of this model build) 
set fms_home       = $cwd:h:h:h/idealized           # directory containing model source code, etc, usually /home/$USER/fms/idealized

echo $fms_home 
 
set days            = 90                             # length of integration 
set runs_per_script = 4                             # number of runs within this script
set start_analysis  = 1                              # number of script run at which to start analysis (after spin-up) 
set num_script_runs = 1                              # how many times to resubmit script to queue
set days_per_segment = ${days}                       # days per segment of analysis (for seasonally-varying analysis)

@ num_segments       = ${days} / ${days_per_segment} # number of analysis segments
echo num_segments    = $num_segments

# find data directory location, 

set data_dir_base = $HOME 
set data_dir = ${data_dir_base}/fms_output/${exp_name}/${run_name}
 

set echo  
echo "*** Running ${run_script} on $HOSTNAME ***"
date

# Tell me which nodes it is run on; for sending messages to help-hpc           
echo " "
echo This jobs runs on the following processors: 
echo `cat $PBS_NODEFILE`
echo " "
#----------------------------------------------------------------------

# zonally averaged analysis 
set analysis_version = analysis_${analysis_type} 
set analysis_script = run_analysis_${model_type}_${analysis_type}
set diagtable   = $exp_home/input/diag_table_${model_type}_${analysis_type}     # path to diagnostics table

set analysis_dir = ${fms_home:h}/analysis/$analysis_version/run                 # location of analysis directory

#--------------------------------------------------------------------------------------------------------

source /etc/profile.d/modules.csh
module load intel/intel-12
module load intel/impi

echo "MPI Used:" `which mpirun`

set NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cpus

#--------------------------------------------------------------------------------------------------------

limit stacksize unlimited 

cd $exp_home 
 
# define variables
set tmpdir      = $home/fms_tmp/${exp_name}          # temporary directory for model workdir, output, etc
set run_dir     = $tmpdir/$run_name                  # tmp directory for current run
set workdir     = $run_dir/workdir                   # where model is run and model output is produced
set output_dir  = $run_dir/output                    # output directory will be created here
set platform    = ifc                                # a unique identifier for your platform

# note the following init_cond's are overwritten later if reload_commands exists
set init_cond   = ""

set pathnames   = $exp_home/input/path_names                       # path to file containing list of source paths


set namelist    = $exp_home/input/namelists_${model_type}          # path to namelist file
set fieldtable  = $exp_home/input/field_table_${model_type}        # path to field table (specifies tracers)
set execdir     = $tmpdir/exe.fms                                  # where code is compiled and executable is created
set run_analysis = $run_dir/analysis                               # where analysis is run
set analysis_out_err = $run_analysis/out_err                       # out and err for analysis
set mppnccombine = $tmpdir/mppnccombine.$platform                  # path to executable mppnccombine
set template    = $fms_home/bin/mkmf.template.${platform}_fram_mpi # path to template for your platform
set mkmf        = $fms_home/bin/mkmf                               # path to executable mkmf
set sourcedir   = $fms_home/src                                    # path to directory containing model source code
set time_stamp  = $fms_home/bin/time_stamp.csh                     # generates string date for file name labels

set ireload     = 1                                  # counter for resubmitting this run script
set irun        = 1                                  # counter for multiple model submissions within this script
#--------------------------------------------------------------------------------------------------------

# if exists, load reload file 

set reload_file = ${run_dir}/reload_commands 

if ( -d $run_dir )  then
  if ( -f $reload_file ) then
     # set irun, ireload, init_cond
     source $reload_file
  endif
endif

#--------------------------------------------------------------------------------------------------------

# setup directory structure
if ( ! -d $execdir ) mkdir -p $execdir
if ( ! -d $run_analysis ) mkdir -p $run_analysis
if ( ! -d $analysis_out_err ) mkdir -p $analysis_out_err

if ( ! -e $workdir ) then
  mkdir $workdir $workdir/INPUT $workdir/RESTART
else
  rm -rf $workdir
  mkdir $workdir $workdir/INPUT $workdir/RESTART
  echo "WARNING: Existing workdir $workdir removed."
endif

if ( ! -d $output_dir )  then
  mkdir -p $output_dir
  mkdir -p $output_dir/combine  
  mkdir -p $output_dir/logfiles
  mkdir -p $output_dir/restart
endif

#--------------------------------------------------------------------------------------------------------

# compile mppnccombine.c, needed only if $npes > 1
if ( ! -f $mppnccombine ) then
  gcc -O -o $mppnccombine -I$fms_home/bin/nc_inc -L$fms_home/bin/nc_lib $fms_home/postprocessing/mppnccombine.c -lnetcdf
endif

#--------------------------------------------------------------------------------------------------------

# compile the model code and create executable

# append fms_home (containing netcdf libraries and include files) to template
/bin/cp $template $workdir/tmp_template
echo "fms_home = $fms_home" >> $workdir/tmp_template


# Prepend fortran files in srcmods directory to pathnames. 
# Use 'find' to make list of srcmod/*.f90 files. mkmf uses only the first instance of any file name.
cd $sourcedir
find $exp_home/srcmods/ -maxdepth 1 -iname "*.f90" -o -iname "*.inc" -o -iname "*.c" -o -iname "*.h" > $workdir/tmp_pathnames
echo "Using the following sourcecode modifications:"
cat $workdir/tmp_pathnames
cat $pathnames >> $workdir/tmp_pathnames

cd $execdir
$mkmf -p fms.x -t $workdir/tmp_template -c "-Duse_libMPI -Duse_netCDF" -a $sourcedir $workdir/tmp_pathnames $sourcedir/shared/include $sourcedir/shared/mpp/include
make -f Makefile

cd $workdir/INPUT

#--------------------------------------------------------------------------------------------------------

# set initial conditions and move to executable directory

if ( $init_cond != "" ) then
  cp $init_cond $init_cond:t
  cpio -iv  < $init_cond:t
#  rm -f $init_cond:t
endif

# name of ocean mask file, will only be used if load_mask = .true. in atmosphere_nml
set ocean_mask = ocean_mask_T42.nc

# if ocean_mask exists, move it to workdir/INPUT folder O.A. May 2014
if (-e $exp_home/input/${ocean_mask}) then
   cp $exp_home/input/${ocean_mask} ocean_mask.nc
   cd $output_dir
   cp $exp_home/input/${ocean_mask} ocean_mask.nc
endif

#--------------------------------------------------------------------------------------------------------

#  --- begin loop over $irun ---                                     
while ($irun <= $runs_per_script)

    cd $workdir  

    # set run length and time step, get input data and executable
    if ( $ireload == 1 && $irun == 1 ) then
      cat > input.nml <<EOF
	&main_nml
         current_time = 0,
         override = .true.,
         days   = $days,
         dt_atmos = 600 /
EOF
    else
      cat > input.nml <<EOF
	&main_nml
         days   = $days,
         dt_atmos = 600 /
EOF
    endif

    if (${model_type} == dry) then 
    cat >> input.nml <<EOF

      &atmosphere_nml      
	two_stream           = .false.,
	turb                 = .true.,
	ldry_convection      = .true.,
	dry_model            = .true.,
	lwet_convection      = .false.,
	mixed_layer_bc       = .false.,
	do_virtual           = .false.,
	tapio_forcing        = .true.,
	hs                   = .false.,
	atmos_water_correction = .false.,
	roughness_mom        = 0.05, 
	roughness_heat       = 0.05,
	roughness_moist      = 0.05,
	bucket               = .false./

EOF
    else if (${model_type} == moist) then

    cat >> input.nml <<EOF    
  
      &atmosphere_nml
	two_stream           = .true.,
	turb                 = .true.,
	ldry_convection      = .false., 
        dry_model            = .false.,
	lwet_convection      = .true.,
	mixed_layer_bc       = .true.,
	do_virtual           = .true.,
	tapio_forcing        = .false., 
	hs                   = .false.,
	atmos_water_correction = .false.,
	roughness_mom        = 5e-03,
	roughness_heat       = 1e-05,
	roughness_moist      = 1e-05,
	bucket               = .false./

EOF

    else if (${model_type} == moist_hydro) then

    cat >> input.nml <<EOF    
  
      &atmosphere_nml
	two_stream           = .true.,
	turb                 = .true.,
	ldry_convection      = .false., 
        dry_model            = .false.,
	lwet_convection      = .true.,
	mixed_layer_bc       = .true.,
	do_virtual           = .true.,
	tapio_forcing        = .false., 
	hs                   = .false.,
	atmos_water_correction = .false.,
	roughness_mom        = 5e-03,
	roughness_heat       = 1e-05,
	roughness_moist      = 1e-05,
	bucket               = .true., 
        load_mask            = .true.,
        init_bucket_depth    = 1000., 
        init_bucket_depth_land = 1., 
        land_left            = 0.,
        land_right           = 360.,
        land_bottom          = 10.,
        land_top             = 30.,
        max_bucket_depth_land= 2., 
        robert_bucket        = 0.04,   
        raw_bucket           = 0.53,       
        damping_coeff_bucket = 200/

EOF
   endif

    cat >> input.nml <<EOF


      &grid_phys_list 
	tsfc_sp                  = 260.0, 
	delh                     = 90.,
	ka_days                  = 50.0,
	ks_days                  = 7.0,
	Cdrag                    = 0.0e-5,
	t_strat                  = 200.0,
	sigma_b                  = 0.85,
	scale_height_ratio       = 3.5,
	reference_sea_level_press = 100000.,
	phi0                     = 0.0/

      &dry_convection_nml
	gamma                    = 0.7,
	tau                      = 14400.0/

      &spectral_init_cond_nml
	initial_temperature  = 280.0 /

      &diag_manager_nml
        mix_snapshot_average_fields = .true. /

      &radiation_nml
        albedo_value                 = 0.38,  
        lw_linear_frac               = 0.2,
        perpetual_equinox            =.true.,
        annual_mean                  =.false.,
        fixed_day                    =.false.,
        fixed_day_value              = 90.0,
        solar_constant               = 1360,
        lw_tau_exponent              = 4.0, 
        sw_tau_exponent              = 2.0,
        lw_tau0_pole                 = 1.8,
        lw_tau0_eqtr                 = 7.2,
        del_sol                      = 1.2,
        atm_abs                      = 0.22,
        days_in_year                 = 360,
        orb_long_perh                = 0,
        orb_ecc                      = 0.0,
        orb_obl                      = 23.5 /

      
      &mixed_layer_nml
        depth              = 40.0,
        qflux_amp          = 0.0,
        qflux_width        = 16.0,
        ekman_layer        = .false.,
        load_qflux         = .false.,
        evaporation        = .true.,
	depth_land         = 1.0/


      &qe_moist_convection_nml
        tau_bm               = 7200.0,
        rhbm                 = 0.7,
        val_inc              = 0.01,
        Tmin                 = 50.,
        Tmax                 = 450. /

      &lscale_cond_nml
	do_evap              = .false./

      # requires topography_option = 'gaussian' in spectral_dynamics_nml
      &gaussian_topog_nml
        height		   = 0.0,
        olon		   = 90.0,
        olat		   = 35.0,
	rlon               = 0.0,
	rlat               = 2.5,
        wlon		   = 4.95,
        wlat   		   = 4.95 /

EOF
    endif

    cat $namelist >> input.nml
    cp $diagtable diag_table
    cp $fieldtable field_table
    cp $execdir/fms.x fms.x

    cp input.nml $run_dir

    
    #--------------------------------------------------------------------------------------------------------
  
    # run the model with mpirun
    set MX_RCACHE=2
    mpirun -np $NPROCS -machinefile $PBS_NODEFILE ${workdir}/fms.x
    
    #--------------------------------------------------------------------------------------------------------

    #   --- generate date for file names ---

    set date_name = `$time_stamp -eh`
    if ( $date_name == "" ) set date_name = tmp`date '+%j%H%M%S'`
    if ( -f time_stamp.out ) rm -f time_stamp.out

    #--------------------------------------------------------------------------------------------------------

    #   --- move output files to their own directories (don't combine) --- 

    mkdir $output_dir/combine/$date_name

    foreach ncfile ( `/bin/ls *.nc *.nc.????` )
	mv $ncfile $output_dir/combine/$date_name/$date_name.$ncfile
    end

    #   --- save ascii output files to local disk ---

    foreach out (`/bin/ls *.out`)
	mv $out $output_dir/logfiles/$date_name.$out
    end

    #   --- move restart files to output directory --- 

    cd $workdir/RESTART
    set resfiles = `/bin/ls *.res*`
    if ( $#resfiles > 0 ) then
	#     --- desired filename for cpio of output restart files ---	
	set restart_file = $output_dir/restart/$date_name.cpio
	if ( ! -d $restart_file:h ) mkdir -p $restart_file:h
	#     --- also save namelist and diag_table ---
	cp $workdir/{*.nml,diag_table} .
	set files = ( $resfiles input.nml diag_table )
	/bin/ls $files | cpio -ocv > $restart_file:t
	mv $restart_file:t $restart_file
	#     --- set up restart for next run ---
	if ( $irun < $runs_per_script ) then 
	    mv -f *.res*  ../INPUT
	endif
    endif

    cd $workdir

    #--------------------------------------------------------------------------------------------------------

    #   --- write new reload information ---
    # for comparison with $start_analysis,  run_number = (ireload-1)*runs_per_script + irun
    set run_number = `expr $ireload \* $runs_per_script - $runs_per_script + $irun`
    echo Completed run $irun of $runs_per_script in bsub $ireload.
    set irun_prev = $irun
    @ irun++

    # remove restart file (init_cond) that is no longer in {reload_file} or ${reload_file}_prev
    if ( -f $reload_file"_prev" ) then
       set irun_tmp = $irun
       set ireload_tmp = $ireload
       set init_cond_tmp = $init_cond
       source $reload_file"_prev"
       rm -r $init_cond
       set irun = $irun_tmp
       set ireload = $ireload_tmp
       set init_cond = $init_cond_tmp
    endif
    if ( -f $reload_file ) mv -f $reload_file $reload_file"_prev"
    
    if ( $irun <= $runs_per_script ) then
	echo "set irun         =  $irun"          >  $reload_file
    else
	@ ireload++
	echo "set irun         =  1"              >  $reload_file
    endif

    echo     "set init_cond    =  $restart_file"  >> $reload_file
    echo     "set ireload      =  $ireload"       >> $reload_file
 
    ############################# post processing ############################
    
    cd $run_analysis
    if (${run_number} >= ${start_analysis}) then # combine data and do analysis

        # need to be careful not to write on top of file for analysis job currently pending in queue.
        # put each job in separate directory.
        set postproc_dir = ${run_analysis}/${date_name} # directory for this analysis run
        
        if ( ! -e $postproc_dir ) then
          mkdir -p ${postproc_dir}
        else
          rm -rf ${postproc_dir}
          mkdir ${postproc_dir}
          echo "WARNING: Existing analysis directory ${postproc_dir} removed."
        endif
        cd ${postproc_dir}
    
	echo "set exp_name = $exp_name" > post_processing_info 
	echo "set data_dir = $data_dir" >> post_processing_info
	echo "set run_name = $run_name" >> post_processing_info
	echo "set date_name = $date_name" >> post_processing_info
	echo "set run_analysis = $run_analysis" >> post_processing_info
	echo "set fms_home = $fms_home" >> post_processing_info
	echo "set tmpdir = $tmpdir" >> post_processing_info
        # specify model resolution, which is set in spectral_dynamics_nml section of input/namelists
        echo "set `grep num_fourier $namelist | tr ',' ' ' `" >> post_processing_info
	echo "set irun = $irun_prev" >> post_processing_info 
	echo "set runs_per_script = $runs_per_script" >> post_processing_info
        # information for segmentation of analysis
        echo "set days_per_segment = $days_per_segment" >> post_processing_info
	echo "set num_segments = $num_segments" >> post_processing_info
	echo "set isegment = 1" >> post_processing_info
    
        cp $analysis_dir/$analysis_script ./

        # ssh to head node and submit analysis script
        ssh -t -l $USER fram "cd ${postproc_dir}; /opt/torque/bin/qsub $analysis_script"
    else
	rm -rf $output_dir/combine/$date_name
    endif

    # don't resubmit if model failed to build a restart file
    if ( ! -f $restart_file ) then
      echo "FATAL ERROR: model restart file not saved. Try moving ${reload_file}_prev to ${reload_file} and re-running."
      echo "ireload = $ireload, irun = $irun"
      set irun = `expr $runs_per_script + 1 `
      set ireload = `expr $num_script_runs + 1 `
    endif


end # --- loop over $irun ended ---

rm -rf $workdir

cd $exp_home/run

if ($ireload > $num_script_runs) then
  echo "Note: not resubmitting job."
else
  echo "Submitting run $ireload."
  qsub $run_script
endif

date
 
 
