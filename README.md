# READ ME FOR SAMPLE CODE TO PLOT CELL TRACTIONS, STRESSES, AND VELOCITIES

*Written by Notbohm Research Group, University of Wisconsin-Madison.* https://notbohm.ep.wisc.edu 

This document gives information about representative code to plot results
of traction force microscopy, monolayer stress microscopy, and computed
cell velocities. Scripts to compute tractions, stresses, and velocities are
in the repository Cell-Traction-Stress ( https://github.com/jknotbohm/Cell-Traction-Stress ).

**IMPORTANT: These scripts give examples of how to plot outputs. You will have to modify them
or write your own scripts to analyze the quantity of interest in your experiment.** 

## List of Files to Run

**plot_cell_vel.m**: Plots cell velocities using both a color map showing magnitudes and quivers
showing relative magnitudes and directions.

**compute_cell_trajectories.m**: Interpolates image correlation data to stitch together approximate
trajectories traveled by different cells.

**plot_displ_tractions.m**: Plots substrate displacements computed by image correlation
and tractions computed by traction force microscopy.

**plot_stress.m**: Plots principal stresses, mean of prinipal stresses, max in-plane shear stress, and prinicpal
orientations computed by monolayer stress microscopy.

Comments in all files state required inputs.

## Additional Subfunctions

**make_fig.m**: Subfunction used to make figure window with desired size and settings.

**map_cold-hot.dat**: Text file containing colormap of cold colors -> black -> hot colors. Use
Colormaps like this one for data having positive and negative numbers, setting black to 0.

**map_cold-hot-wrap.dat**: Text file containing colormap of colors that wrap from blue to black
to red to white. Useful for plotting angles.

Some scripts require m files available from the Matlab
file repository. See comments of each script for more information.















