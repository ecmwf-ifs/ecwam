# ---------------------------------------------------------
import softwaremodules

softwaremodules.module("load ecmwf-toolbox/new")

import os
import metview as mv
import numpy as np
# ---------------------------------------------------------

# Variables to be plotted
vvars=['swh']

# Tests to be plotted
tests=['etopo1_oper_an_fc_O48_cy50r1.yml']

rundir='../ecwam_runs'             # base directory of the model runs
plotdir='../ecwam_plots'           # directory where the plots will be saved
datestr='20230101060000'           # forecast date to be plotted

# Ensure output directory exists before saving plots.
if not os.path.exists(plotdir):
    os.makedirs(plotdir, exist_ok=True)

# ---------------------------------------------------------
# metview plot settings

coast = mv.mcoast(
    map_coastline_land_shade        = "on",
    map_coastline_land_shade_colour = "grey",
    map_coastline_sea_shade         = "on",
    map_coastline_sea_shade_colour  = "RGB(0.8944,0.9086,0.933)",
    map_coastline_thickness         =  2,
    map_boundaries                  = "on",
    map_boundaries_colour           = "charcoal",
    map_grid_colour                 = "charcoal",
    map_grid_longitude_increment    = 10
    )

view = mv.geoview(map_area_definition = 'corners',coastlines = coast)

shade_swh = mv.mcont(
    legend                       = "on",
    contour_line_colour          = "navy",
    contour_highlight            = "off",
    contour_level_selection_type = "level_list",
    contour_level_list           = [0,1,2,3,4,5,6,8,10],
    contour_label                = "off",
    contour_shade                = "on",
    contour_shade_colour_method  = "palette",
    contour_shade_method         = "area_fill",
    contour_shade_palette_name   = "eccharts_rainbow_blue_red_7"
    )

# ---------------------------------------------------------

# extract experiment names from test file names
exps = [os.path.basename(test).split('.')[0] for test in tests]

# plot (loop over variables and experiments)
for var in vvars:
    for exp in exps:
        ddir=f'{rundir}/{exp}/output/'
        
        file=f'MPP{datestr}_000000000000'
    
        print(f"Processing {var} for experiment {exp}")
        
        g_fc     = mv.read(ddir + file)
        fc       = g_fc.select(shortName=var)
        dw       = mv.plot_superpage(pages = mv.mvl_regular_layout(view,1,1,1,1))
        fc_title = mv.mtext(text_lines=[f"{var} {exp}"], text_font_size=0.5)

        prob_shade = globals()[f"shade_{var}"]
        
        mv.setoutput(mv.png_output(output_name=f"{plotdir}/{var}_{exp}"))
        mv.plot(dw, fc, prob_shade, fc_title)


print(f'\nPlots created successfully and saved into plotdir: {plotdir}')