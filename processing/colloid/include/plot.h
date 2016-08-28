#ifndef PLOT_H
#define PLOT_H

#include "colloid_base.h"

void plot_config(colloid_base &, const int, bool=true);

void plot_config2(colloid_base &, const int, 
		const char *, bool=true);

void config_plot(const char *file, const char *str_frameindex, 
				 int framerate, bool Tracked);

// I have 2 different traj_plot in src/plot/
// One is use 'plot' to plot dots with color	
// now the one we use is 'splot' to plot images
// hence it should be very careful to choose twindow and tstep
// Also a good frameindex should be assigned
void traj_plot(const char *file, int twindow, int tstep,
				 const char *str_frameindex, 
				 int framerate, bool Tracked);

#endif /* PLOT_H */
