from matplotlib import pyplot, _pylab_helpers

def clearAllPlots():
    #figures = [manager.canvas.figure for manager in _pylab_helpers.Gcf.get_all_fig_managers()]
    #for f in figures:
    #    f.close()
    pyplot.close('all')

def configurePlots(number):
    figures = []
    axes = []
    nfigs = len(_pylab_helpers.Gcf.get_all_fig_managers())
    for i in range(number):
        figures.append(pyplot.figure(i+nfigs))
        figures[-1].clear()
        axes.append(figures[-1].add_axes([0.1, 0.1, 0.8, 0.8]))
        axes[-1].clear()

    return figures, axes
