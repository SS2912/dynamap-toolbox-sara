# dynamap-toolbox-sara
Basic functions created by SS for analysis of PSD, FC and file organization

**Tools for file organization / modification**
  - excel2mtg.m
  - extractmontageEImax.m
  - load_files.m          --> Loads only files with chosen string in the file name

**Functions for data (pre-)processing**
- Delphos
  - select_markers.m      --> Analyze DELPHOS results with sliding window (VLM, SS)
  - window_spikes.m       --> Original by VLM. SS: added counts_window + changes 'spikes' into 'events' so that calculation of HFOs events is possible
- PSD
  - plot_psd_diff         --> support function for PSD comparison in PSD pipeline RF-TC project
  - zplot.m               --> computes and plot the zscore values for the difference/ ratio matrices of spectral density for each channel.
- FC
  - FCmatrix_no_rep.m     --> calculates mean/median on FC matrix for each anat label
  - load_conn_matrix.m    --> organizes FC and other data into melt tables (mixed models-like). Made specifically for PS2
  - MMTable.m             --> same as above, made specifically for RF-TC project
  - reduce_h2matrix.m     --> Reduce h2 matrix to contain only the regions selected in "roi" (or other table column of choice). Support for analysis pipelines
- NNEI
  - regpeaks.m            --> finds NNEI response and plots it in same plot for each detected stim peak
  - NNEItrial.m           --> plots stim signal and SEEG signal for each NNEI trial of stim
