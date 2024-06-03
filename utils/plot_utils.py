import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size

# The following two functions are written by Nathaniel Linden, provided May 2024
def export_legend(legend, filename="legend.pdf"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)

def get_sized_fig_ax(width, height, hspace=1.0, vspace=0.5, fig_height=6, fig_width=6):
    fig = plt.figure(figsize=(fig_width, fig_height))
    # The first items are for padding and the second items are for the axes.
    # sizes are in inch.
    h = [Size.Fixed(hspace), Size.Fixed(width)]
    v = [Size.Fixed(vspace), Size.Fixed(height)]

    divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
    # The width and height of the rectangle are ignored.
    ax = fig.add_axes(divider.get_position(),
                    axes_locator=divider.new_locator(nx=1, ny=1))
    
    return fig, ax