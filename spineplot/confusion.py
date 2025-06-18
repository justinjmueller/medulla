import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix

from artists import SpineArtist
from style import Style
from variable import Variable
from utilities import mark_pot, mark_preliminary

class ConfusionMatrix(SpineArtist):
    """
    A class to create a confusion matrix given a pair of variables that
    represent predicted and true class labels. The confusion matrix
    provides a summary of the prediction results on a classification
    problem. It shows the number of correct and incorrect predictions
    with respect to the true class labels, and how they are distributed
    amongst the available classes.

    Attributes
    ----------
    _title : str
        The title of the artist. This will be placed at the top of the
        axis assigned to the artist.
    _categories : dict
        A dictionary mapping the category key to the category name.
    _true_labels : Variable
        The Variable that contains the true class labels.
    _predicted_labels : Variable
        The Variable that contains the predicted class labels.
    _xlabel : str
        The x-axis label for the confusion matrix.
    _ylabel : str
        The y-axis label for the confusion matrix.
    """
    def __init__(self, categories, true_labels, predicted_labels, **kwargs):
        """
        Parameters
        ----------
        categories : dict
            A dictionary mapping the category key to the category name.
        true_labels : Variable
            The Variable that contains the true class labels.
        predicted_labels : Variable
            The Variable that contains the predicted class labels.
        kwargs : dict, optional
            title : str, optional
                The title of the confusion matrix. If not provided, a
                default title will be used.
            xlabel : str, optional
                The x-axis label for the confusion matrix. If not
                provided, a default label will be used.
            ylabel : str, optional
                The y-axis label for the confusion matrix. If not
                provided, a default label will be used.
        """
        # Keyword arguments
        title = kwargs.get('title', None)
        self._xlabel = kwargs.get('xlabel', 'Class Prediction')
        self._ylabel = kwargs.get('ylabel', 'Class Label')

        super().__init__(title=title)
        self._categories = categories
        self._true_labels = true_labels
        self._predicted_labels = predicted_labels

        # Modify categories to include a placeholder for NaN values
        self._categories[-1] = 'Null'

        # Initialize the confusion matrix with zeros.
        self._confusion_matrix = np.zeros((len(categories), len(categories)), dtype=int)

    def add_sample(self, sample, is_ordinate):
        """
        Add a sample to the confusion matrix.

        Parameters
        ----------
        sample : Sample
            The sample containing the predicted and true class labels.
        is_ordinate : bool
            Indicates whether the sample is an ordinate (true label) or not.
        """
        super().add_sample(sample, is_ordinate)

        # Update the confusion matrix with the sample's true and
        # predicted labels.
        data, _ = sample.get_data([self._true_labels._key, self._predicted_labels._key])
        for category, values in data.items():
            if category not in self._categories.keys():
                continue
            true_labels = values[0]
            predicted_labels = values[1]

            # The predicted label may contain NaN values, which should
            # be replaced with a placeholder for the confusion matrix.
            predicted_labels = np.nan_to_num(predicted_labels, nan=-1)

            # Update the confusion matrix with the true and predicted
            # labels for the current sample.
            self._confusion_matrix += confusion_matrix(true_labels, predicted_labels, labels=list(self._categories.keys()))


    def draw(self, ax, style, **kwargs):
        """
        Draw the confusion matrix on the provided axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis on which to draw the confusion matrix.
        style : Style
            The style to use for drawing the confusion matrix.
        **kwargs : dict
            show_null_column : bool, optional
                If True, the null column (for NaN values) will be shown in
                the confusion matrix. The default is False.
        """
        # Keyword arguments
        show_null_column = kwargs.get('show_null_column', False)

        super().draw(ax, style)

        # Make a float copy of the matrix
        cm = self._confusion_matrix.astype(float)

        # Find the index of the NaN category (label -1)
        labels = list(self._categories.keys())
        nan_index = labels.index(-1)

        # Remove the row for true label -1 (NaN), but keep the column
        cm = np.delete(cm, nan_index, axis=0)

        # If show_null_column is False, remove the column for predicted
        # label -1 (NaN) as well.
        if not show_null_column:
            cm = np.delete(cm, nan_index, axis=1)
            labels.pop(nan_index)

        # Normalize rows (excluding the removed one)
        row_sums = cm.sum(axis=1, keepdims=True)
        cm = np.divide(cm, row_sums, where=row_sums != 0)

        # Build index and column labels
        row_labels = [self._categories[lbl] for i, lbl in enumerate(labels) if i != nan_index]
        col_labels = [self._categories[lbl] for lbl in labels]

        # Create the DataFrame
        cm_df = pd.DataFrame(cm, index=row_labels, columns=col_labels)
        cm_df = cm_df.T  # Transpose the matrix
        cm_df = cm_df.iloc[::-1]  # Reverse the order of the index (y-axis)

        # Plot the confusion matrix using a heatmap.
        im = ax.imshow(cm_df, interpolation='nearest', vmin=0, vmax=1)
        ax.set_title(self._title)
        ax.set_xlabel(self._xlabel)
        ax.set_ylabel(self._ylabel)

        # Add annotations.
        for i in range(len(cm_df)):
            for j in range(len(cm_df.columns)):
                ax.text(j, i, f'{cm_df.iloc[i, j]:.2%}', ha='center', va='center', color='white' if cm_df.iloc[i, j] > cm_df.max().max() / 2 else 'black')

        # Set the ticks and labels.
        ax.set_xticks(np.arange(len(cm_df.columns)))
        ax.set_yticks(np.arange(len(cm_df.index)))
        ax.set_xticklabels(cm_df.columns)
        ax.set_yticklabels(cm_df.index, rotation=90, va='center')
        
        # Disable minor ticks
        ax.tick_params(axis='both', which='minor', bottom=False, left=False)

        mark_pot(ax, self._exposure, style.mark_pot_horizontal)
        mark_preliminary(ax, style.mark_preliminary)