import pylab

points=[(0,0,0),(1,1,0),(2,2,0),(3,2,1),(4,3,1),(5,4,1),(6,4,2)]
def plot_counts(points,xlim=None,ylim=None,show=True):
	if show:
		pylab.clf()

	points.sort()
	if xlim==None:
		xlim=points[-1][1]*2
		ylim=points[-1][1]*1.4

	pylab.ylim((0,ylim))
	pylab.xlim((0,xlim))

	cax = pylab.gca()
	cax.set_aspect('equal')

	pylab.text(0.2,0.1,"Count plot",fontsize=8)
	pylab.plot([x[0] for x in points],[y[1] for y in points], 'b-',linewidth=2)
	pylab.plot([x[0] for x in points],[y[2] for y in points], 'r-',linewidth=2)
	if show:
		pylab.show()