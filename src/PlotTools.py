from matplotlib import pyplot

def configurePlots(number):
    figures = []
    axes = []
    for i in range(number):
        figures.append(pyplot.figure(i))
        figures[-1].clear()
        axes.append(figures[-1].add_axes([0.1, 0.1, 0.8, 0.8]))
        axes[-1].clear()

    return figures, axes
