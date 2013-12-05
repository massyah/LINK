import pylab

# points=[(1,2),(2,3),(3,12),(0,0),(1,1)]
def plot_roc(points,show=True):
	if show:
		pylab.clf()

	points.sort()
	pylab.ylim((0,1))
	pylab.xlim((0,1))

	cax = pylab.gca()
	cax.set_aspect('equal')

	pylab.text(0.2,0.1,"test ROC",fontsize=8)
	pylab.plot([x[0] for x in points],[y[1] for y in points], 'r-',linewidth=2)
	if show:
		pylab.show()