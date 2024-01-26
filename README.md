# AstridData
Utilities for processing the data from the ASTRID simulation


## compress_snapshot
Turn full ASTRID snapshot into float16 formatting for most columns to save disk space. Should work for both PART and PIG files. Currently the algorithm guarantees a 2% loss in all fields, and almost zero loss in Position and Velocities (an other important fields for restarting the simulation).
- see `job_example_selcols.slurm` for compressing selected data columns
- see `job_example_allcols.slurm` for compressing full snapshots
- current ASTRID compression/recovery jobs are stored in `job_scipts`


## merger_catalog



## subcat_bigfile
