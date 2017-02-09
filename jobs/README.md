**jobs**
---------

Contains top level scripts for running openquake models on the NCI. The main script that you will need to edit is `run_oq_model.py`. 

Based on the parameters defined in the top part of this script, an output directory for the job will be created capturing timestamp, model name and user information. All relevant input files (`.xml` and `.ini`) are copied to the output directory to ensure provenance is captured.

An NCI run script is created and the job submitted to the NCI queue.

**Usage**

The user should edit the top part of `run_oq_model.py` to defined the NCI resources required and the location of the input data. Note that as currently written, it assumes you have a sandpit folder specified by your username under the folder `<sandpit_path>`.

    # Parameters for building run_<model>.sh qsub job
    # No rmally you would only need to change these, although
    # other parameters are hard coded at the bottom of the file
    ncpus = 256 
    mem = '512GB'
    walltime = '03:00:00'
    jobfs='200GB' 

    # Paths to sandpit and input path, assuming you have
    # a sandpit named after your NCI username
    sandpit_path = '/short/w84/NSHA18/sandpit/'
    model_rel_path = 'NSHA2018/source_models/smoothed_seismicity/Cuthbertson2016'
    # All outputs will be under this directory
    model_output_base = '/short/w84/NSHA18/PSHA_modelling/'

Below these lines you should not need to change anything, however additional NCI job parameters are defined towards the end of the file. All your input data (i.e. `source model.xml`, `source_model_logic_tree.xml`, 'gsim_logic_tree' and `job.ini`) should be under the directory specified by `/<sandpit_path>/<username>/<model_rel_path>`. These are standard OpenQuake input files.

The model is then simply run as `python run_oq_model.py`

Standard output is written to `parjob.log` in the output directory. The output directory will be  a timestamped diretory, in this case under `/short/w84/NSHA18/PSHA_modelling/Cuthbertson2016`.

After running the model the NCI will report on the resources actually used in a file in the outut directory with a name like `oq512c512ht.o<job_number>`. This can be used to adjust the resources erquested for future runs to ensure efficient use of the NCI. The default parameters given above are reasonable for a high resolution national scale model with a reasonable complex logic tree. Other models (e.g. a simple smoothed seismicity model) will use much less memory and wall time.