import matplotlib.pyplot as plt


def ShowLegend(ax):

	# set some legend properties.  All the code below is optional.  The
	# defaults are usually sensible but if you need more control, this
	# shows you how
	leg = plt.legend(loc='best', shadow=True, fancybox=True)

	# the matplotlib.patches.Rectangle instance surrounding the legend
	frame = leg.get_frame()
	frame.set_facecolor('1.00')    # set the frame face color to light gray

	# matplotlib.text.Text instances
	for t in leg.get_texts():
		t.set_fontsize('x-large')    # the legend text fontsize

	# matplotlib.lines.Line2D instances
	for l in leg.get_lines():
		l.set_linewidth(4.0)  # the legend line width

