import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from artists import SpineArtist
from style import Style
from variable import Variable
from utilities import mark_pot, mark_preliminary

class Ternary(SpineArtist):
    """
    A class to create a "ternary" plot, which is a type of barycentric
    plot that shows the relationship between three variables that sum
    to a constant value. This is a common plot type in the field of
    machine learning, where some types of classification problems
    produce scores expressing the network's confidence in each of the
    possible classes. The scores are normalized are often normalized
    to sum to one. Functionally, this means that one of the scores
    is dependent on the others.

    In the three-class case, a ternary plot can be used to visualize
    the separation between the classes that is achieved by the network
    scores. The plot is a triangle with each vertex representing one
    class. The points inside the triangle represent the scores for
    each of the classes. The closer a point is to a vertex, the more
    confident the network is in that class. 

    Attributes
    ----------
    _title : str
        The title of the artist. This will be placed at the top of the
        axis assigned to the artist.
    _variables : list
        A list of the Variable objects that represent the scores for
        the classes. The order of the variables in the list determines
        the order of the vertices in the ternary plot. The number of
        variables must be equal to three.
    _categories : dict
        A dictionary mapping the category key to the group name.
    _coords : np.ndarray
        An array of the coordinates of the data points in the ternary
        plot. This has shape (2, N), where N is the number of data
        points collected across the samples.
    """
    def __init__(self, var0, var1, var2, categories, title=None):
        """
        Parameters
        ----------
        var0 : Variable
            The variable that represents the score for the first class.
        var1 : Variable
            The variable that represents the score for the second class.
        var2 : Variable
            The variable that represents the score for the third class.
        categories : dict
            A dictionary mapping the category key to the group name.
        title : str
            The title of the artist. The default is None.

        Returns
        -------
        None.
        """
        super().__init__(title)
        self._variables = [var0, var1, var2]
        self._categories = categories
        self._coords = None

    def add_sample(self, sample, is_ordinate):
        """
        Add a sample to the ternary plot.

        Parameters
        ----------
        sample : Sample
            The sample to add to the ternary plot.
        is_ordinate : bool
            A boolean flag that determines if the sample is to be used
            as the ordinate. If False, the sample is used as the abscissa.

        Returns
        -------
        None.
        """
        super().add_sample(sample, is_ordinate)

        # Retrieve the data for this sample
        vars = [v._key for v in self._variables]
        data, _ = sample.get_data(vars)
        data = np.hstack([np.stack(data[x]) for x in data.keys() if x in self._categories.keys()])   

        # Normalize the data to a fixed sum of 1.
        score_sum = np.sum(data, axis=0)
        valid_mask = score_sum > 0
        data = data[:, valid_mask]

        # Calculate the coordinates of the data points in a triangular
        # coordinate system.
        x = 0.5 * (2 * data[2] + data[1]) / score_sum[valid_mask]
        y = (np.sqrt(3)/2) * data[1] / score_sum[valid_mask]

        # Extend the _coords array by these new coordinates
        if self._coords is None:
            self._coords = np.stack([x, y])
        else:
            self._coords = np.hstack([self._coords, np.stack([x, y])])
        
    def draw(self, ax, labels=None, axis_labels=None, cmap=None, contour=False, npts=10000, style=None):
        """
        Draw the ternary plot on the provided axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to draw the plot on.
        labels : list, optional
            A list of the labels to use for the vertices of the triangle.
            The default of None will use basic placeholder labels.
        axis_labels : list, optional
            A list of the labels to use for the axes of the triangle. The
            default of None will use basic placeholder labels.
        cmap : str, optional
            The colormap to use when drawing the plot. The default is None,
            which will use the default colormap.
        contour : bool, optional
            A boolean toggle that determines if the plot should be drawn
            as a contour plot. The default is False, which will draw the
            plot as a scatter plot.
        npts : int, optional
            The number of points to use when drawing the plot. The default
            is 10000.
        style : Style, optional
            The style to use when drawing the plot. The default is None.

        Returns
        -------
        None.
        """
        # Estimate the density of the data points, but use a random
        # sampling of the data to evaluate the density for the plot
        # (the contour drawing is very slow for large numbers of
        # points).
        gaus = gaussian_kde(self._coords)
        selection = np.random.choice(len(self._coords[0]), int(npts))
        if contour:
            ax.tricontourf(self._coords[0][selection],
                        self._coords[1][selection],
                        gaus(self._coords[:,selection]),
                        levels=100, cmap=cmap)
        else:
            ax.scatter(self._coords[0][selection], self._coords[1][selection], s=0.5, c=gaus(self._coords[:,selection]), cmap=cmap)

        # Draw the triangle
        labelsize = plt.rcParams['axes.labelsize']
        ax.plot([0, 0.5, 1, 0], [0, np.sqrt(3)/2, 0, 0], 'k-', lw=1.5)
        ax.text(-0.05, -0.05, 'Vtx. 0' if labels is None else labels[0], ha='center', fontsize=labelsize)
        ax.text(0.5, np.sqrt(3)/2 + 0.05, 'Vtx. 1' if labels is None else labels[1], ha='center', fontsize=labelsize)
        ax.text(1.05, -0.05, 'Vtx. 2' if labels is None else labels[2], ha='center', fontsize=labelsize)

        ax.text(0.5, -0.02, 'Axis A' if axis_labels is None else axis_labels[0], ha='center', va='top', fontsize=labelsize, rotation=0)
        ax.text(0.27, np.sqrt(3)/4, 'Axis B' if axis_labels is None else axis_labels[1], ha='right', va='center', fontsize=labelsize, rotation=+60)
        ax.text(0.7, np.sqrt(3)/4, 'Axis C' if axis_labels is None else axis_labels[2], ha='left', va='center', fontsize=labelsize, rotation=-60)

        # Final adjustments
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, np.sqrt(3)/2 + 0.1)
        ax.set_aspect('equal')
        ax.axis('off')
        Ternary.draw_grid(ax)
        if self._title is not None:
            ax.set_title(self._title)

        if style.mark_pot:
            mark_pot(ax, self._exposure, style.mark_pot_horizontal)
        if style.mark_preliminary is not None:
            mark_preliminary(ax, style.mark_preliminary)

    @staticmethod
    def draw_grid(ax, levels=5, color='black', linestyle='--', linewidth=1):
        """
        Draw a triangular grid on the provided axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to draw the grid on.
        levels : int, optional
            The number of levels to draw on the grid. The default is 5.
        color : str, optional
            The color to use when drawing the grid. The default is 'black'.
        linestyle : str, optional
            The linestyle to use when drawing the grid. The default is '--'.
        linewidth : float, optional
            The width of the lines in the grid. The default is 0.5.
        """
        for i in range(1, levels):
            fraction = i / levels

            # Lines parallel to base.
            ax.plot([fraction / 2, 1 - fraction / 2], [fraction * np.sqrt(3) / 2, fraction * np.sqrt(3) / 2],
                    color=color, linestyle=linestyle, linewidth=linewidth)
            
            # Lines parallel to left edge.
            ax.plot([1-fraction, 1 - fraction / 2], [0, fraction * np.sqrt(3) / 2],
                    color=color, linestyle=linestyle, linewidth=linewidth)

            # Lines parallel to right edge.
            ax.plot([fraction / 2, fraction], [fraction * np.sqrt(3) / 2, 0],
                    color=color, linestyle=linestyle, linewidth=linewidth)