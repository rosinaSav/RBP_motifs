'''
Author: Rosina Savisaar.
Module containing functions for plotting data.
'''

import matplotlib.pyplot as plt
import pylab

def add_labels(x_lab, y_lab, title):
    '''
    Add axis labels and title to plot.
    '''
    if x_lab:
        plt.xlabel(x_lab)
    if y_lab:
        plt.ylabel(y_lab)
    if title:
        plt.title(title)

def histogram(x, bins, x_lab = None, y_lab = None, title = None, colour = None):
    '''
    Plot histogram.
    '''
    if not colour:
        plt.hist(x, bins, alpha = 0.5)
    else:
        plt.hist(x, bins, alpha = 0.5, color = colour)
    add_labels(x_lab, y_lab, title)

def line_plot(x, y, x_lab = None, y_lab = None, title = None):
    '''
    Plot line plot.
    '''
    plt.plot(x, y, color = "black", linewidth = 1.5)
    add_labels(x_lab, y_lab, title)
    
def save_and_show(size, dpi, fig_name, show = True):
    '''
    Save a plot to file and display it on screen.
    '''
    fig = pylab.gcf()
    fig.set_size_inches(size[0], size[1])
    fig.savefig('{0}.png'.format(fig_name), dpi=dpi)
    if show:
        plt.show()
