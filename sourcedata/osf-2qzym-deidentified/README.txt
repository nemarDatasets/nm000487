READ ME:

This data set is associated with the following manuscript :

Title:
'Concurrent category-selective neural activity across the ventral occipito-temporal cortex supports a non-hierarchical view of human visual recognition'.

Authors: 
Corentin Jacques, Jacques Jonas, Sophie Colnat-Coulbois, Bruno Rossion.

doi:     https://doi.org/10.7554/eLife.109640.1


Corresponding Author for article and dataset:
Bruno Rossion


--- 1. Experiment 

The experiment is described in the published manuscript. Each visual stimulation sequence starts with 4 non-face object images before the first face image, so as to result in a sequence of  [O O O O F O O O O F ...] .



--- 2. Electrophysiological data ---

Each folder contains the SEEG data from one participant. The data consist in recordings of intracerebral EEG in epileptic patient viewing visual stimuli. Each channel corresponds to the recording (i.e. voltage over time) of one intracerebral contact. 


Data is spread over two files (.mat and .lw5) that can be open in Letswave 5, a free software for EEG data analysis (http://www.nocions.org/letswave/) running on Matlab; or simply using matlab (see below). 

From matlab, use the following syntax:
[header, data] = LW_load('filename.lw5'); 

THe 'LW_load.m' function is provided.

Letswave 6 can also be used simply by changing the file extension from '.lw5' to '.lw6'. 


Recorded sEEG data (recorded as TRC files using Micromed system) has been minimaly pre-processed in the following way:

1. Imported in Letswave 5, which has been used for the subsequent processing steps
2. Each recording was segmented as sequences from 2 seconds before the onset of the visual stimulation sequence (but see some exceptions below) until 2 seconds after the end of the FPVS sequence. 
3. Sequences (varying number across participants) were merged in the same file. 


Further data transform and analyses have not been performed.


Letswave files can be open in Matlab without 'Letswave'. They consist of one header file (the .lw5 file) and one data file (the .mat file). 
THe header file (.lw5) is saved as version anterior to v7.3 while the data file (the .mat file) is save as a '-v7.3' file. This is important if loading the files outside Matlab (e.g. Python). 
The dimensions of the 'data' variable are organized as follow:

-dim 1: sequences/epochs
-dim 2: channels
-dim 3: frequency (only used for time-frequency analyses)
-dim 4: unused
-dim 5: unused
-dim 6: time

The header file contains the following fields

    filetype: 
        name: 
        tags: 
     history: 
    datasize: size of data variable
      xstart: time (sec) relative to event onset (here the onset of the stimulation sequence
      ystart: lowest frequency (if time-frequency data)
      zstart: 
       xstep:  = 1/time_samplingRate
       ystep: = 1/frequency_samplingRate
       zstep: 
    chanlocs: channel names
      events: contains codes and latencies for events (for most patient, the event code is meaningless)

*The 'chanloc' field contains:
-Contact label ( contacts with a ' are in the left hemisphere). 
-3D talairach coordinate
-Anatomical labels for main region (e.g. OCC, PTL, ATL, MTL, TP) - WM refers to white matter.
-Anatomical labels for sub region
-Whether the contact is within the VOTC region.

*Sampling rate is defined as 1/header.xstep. It was 256Hz for 3 participants, 500Hz for 11 participants and 512Hz for 126 participants.

*header.xstart defines the time frame of the recording. For a majority of participant header.xstart = -2 (recording starts 2 seconds before the onset of the visual stimulation). For other participants, header.xstart = -1.8333 . This is because for some participants, due to technical error, the very first event of each sequence was missing (even though all stimuli were presented correctly). The correction of header.xstart accounts for that (i.e. offsetting the recording by 1/6 s given the 6 Hz stimulation). THe information about missing first trigger or not is presented in the '140_participants_seeg_info.txt' (or .xlsx) file.

*In the manuscript, whole recording sequences were segmented in short epoch centered on face onset based timing only (given that image have predictable/periodic time onset), not based on triggers/events. 


--- 3. Demographics
Demographics about participants is contained in '140_participants_demographics.txt' 









