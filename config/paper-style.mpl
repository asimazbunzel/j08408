# Matplotlib style sheet
# For more info, see: https://matplotlib.org/users/customizing.html

#### LINES
lines.linewidth      : 1.5     ## line width in points
lines.marker         : None    ## the default marker
lines.markersize     : 4       ## markersize, in points

#### FONT
font.family          : STIXGeneral
font.style           : normal
font.size            : 10.0
font.serif           : DejaVu Serif
font.sans-serif      : DejaVu Sans
font.monospace       : DejaVu Sans Mono

#### TEXT
text.usetex          : False
mathtext.rm          : sans
mathtext.fontset     : stix

#### AXES
axes.labelsize       : medium   ## fontsize of the x any y labels

#### TICKS
## see http://matplotlib.org/api/axis_api.html#matplotlib.axis.Tick
xtick.top            : True   ## draw ticks on the top side
xtick.bottom         : True   ## draw ticks on the bottom side
xtick.labelbottom    : True   ## draw label on the bottom
xtick.labelsize      : medium ## fontsize of the tick labels
xtick.direction      : inout  ## direction: in, out, or inout
xtick.minor.visible  : True   ## visibility of minor ticks on x-axis

ytick.left           : True   ## draw ticks on the left side
ytick.right          : True   ## draw ticks on the right side
ytick.labelleft      : True   ## draw tick labels on the left side
ytick.labelsize      : medium ## fontsize of the tick labels
ytick.direction      : inout  ## direction: in, out, or inout
ytick.minor.visible  : True   ## visibility of minor ticks on y-axis

#### LEGEND
legend.frameon       : False     ## if True, draw the legend on a background patch
legend.fontsize      : 8.0
# legend.borderpad     : 0.4
# legend.labelspacing  : 0.4
# legend.handletextpad : 0.5
# legend.columnspacing : 0.5

#### FIGURE
## See http://matplotlib.org/api/figure_api.html#matplotlib.figure.Figure
figure.figsize       : 3.54399, 3   ## figure size in inches
figure.dpi           : 300

#### SAVING FIGURES
## the default savefig params can be different from the display params
## e.g., you may want a higher resolution, or to make the figure
## background white
savefig.dpi         : figure
savefig.format      : pdf
savefig.directory   : ~/plots
savefig.bbox        : tight
savefig.pad_inches  : 0.1 # Padding when bbox is tight