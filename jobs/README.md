**jobs**
---------

Contains top level scripts for running openquake models on the NCI. The main script that you use is `run_oq_model.py`. This is executable and therefore you should set your PATH environment variable to point to this directory. E.g. add this line to our `.bashrc` file:
	 
    export PATH=$PATH:/<path_to_sandpit>/sandpit/<user_name>/NSHA2018/jobs	

**Usage**

The script uses a parameter file params.txt and is run as:
    run_oq_model.py -param_file params.txt

Based on the parameters defined in `params.txt`, an output directory for the job will be created capturing timestamp, model name and user information. All relevant input files (`.xml` and `.ini`) are copied to the output directory to ensure provenance is captured.

An NCI run script is created and the job submitted to the NCI queue.

The user should edit `params.txt` to defined the NCI resources required and the location of the input data. Note that as currently written, it assumes you have a sandpit folder specified by your username under the folder `sandpit_path`.

    ncpus = 256 
    mem = 256GB
    walltime = 02:00:00
    jobfs = 300GB
    sandpit_path = /short/w84/NSHA18/sandpit/
    model_rel_path = NSHA2018/source_models/faults/Adelaide/Adelaide_faults
    shared_rel_path = NSHA2018/shared
    gsim_lt_filename = NSHA18_Aus_GMPE_75thp_logic_tree_cal_power_p0.4.xml
    model_output_base = /short/w84/NSHA18/PSHA_modelling/
    job_file = job.ini
    
You should not need to change anything in run_oq_model.py, however additional NCI job parameters are defined towards the end of the file. All your source mode input data adn the job file (i.e. `source model.xml`, `source_model_logic_tree.xml` and `job.ini`) should be under the directory specified by `/<sandpit_path>/<username>/<model_rel_path>`. Ground motion model files (`gsim_logic_tree`) and sites files (if required; `stes.csv`) are generally shared bewteen many models and therefore these are stored in the path defined by `shared_rel_path`. These are standard OpenQuake input files.

Standard output is written to `parjob.log` in the output directory. The output directory will be  a timestamped diretory, in this case under `/short/w84/NSHA18/PSHA_modelling/Adelaide_faults`.

After running the model the NCI will report on the resources actually used in a file in the outut directory with a name like `oq512c512ht.o<job_number>`. This can be used to adjust the resources erquested for future runs to ensure efficient use of the NCI. The default parameters given above are reasonable for a high resolution national scale model with a reasonable complex logic tree. Other models (e.g. a simple smoothed seismicity model) will use much less memory and wall time.