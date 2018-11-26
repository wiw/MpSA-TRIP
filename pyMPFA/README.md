For run TRIP analysis needs perform next steps:


I. Prepare folder structure:

    1. For each experiment make folder with name of experiment name

    2. Into experiment directory create three folders named as: expression, mapping and normalization

    3. Into each previous folder copy content from github: "https://github.com/wiw/pyMPFA/tree/trip_0.3/pyMPFA"

    4. Then in each folder of expression, mapping and normalization edit file param.py for our conditions


II. Run load preliminary data and prepare it:
    
    1. If data is expression and/or normalization - it's OK. Use this command to run in common view:
        > ./command_to_run.sh genome <PATH_TO_FASTQ/FASTQ.GZ> <PATH_TO_EXPERIMENT_RUN_FOLDER> <PATH_TO_EXPERIMENT_OUTPUT_FOLDER> <WHICH,DATA,IS,USE> rev
        EXAMPLE:
        ./command_to_run.sh genome /home/anton/backup/input/trip/RUN_2018-06-07/1_S1_L001_R1_001.fastq.gz /home/anton/backup/input/trip/RUN_2018-06-07/run/trip_6_2 /home/anton/backup/input/trip/RUN_2018-06-07/results/trip_6_2 expression,normalization rev
    
    2. If data is mapping and you have paired sequences, then use next command:
    > ./command_to_run.sh paired <PATH_TO_FWD_FASTQ,PATH_TO_REV_FASTQ> <PATH_TO_EXPERIMENT_RUN_FOLDER> <PATH_TO_EXPERIMENT_OUTPUT_FOLDER> mapping rev
    EXAMPLE:
    ./command_to_run.sh paired /home/anton/backup/input/trip/RUN_2018-06-07/1_S1_L001_R1_001.fastq.gz,/home/anton/backup/input/trip/RUN_2018-06-07/1_S1_L001_R2_001.fastq.gz /home/anton/backup/input/trip/RUN_2018-06-07/run/trip_6_2 /home/anton/backup/input/trip/RUN_2018-06-07/results/trip_6_2 mapping rev


III. Analysing of received data:

    1. Go to folder "/home/anton/data/TRIP/pyMPFA" and edit variable `CONFIG` in file "AnalysisMain_trip_fork.py" to your parameters.

    2. Then execute this file in python environment:
        > python AnalysisMain_trip_fork.py