import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

from artists import SpineArtist
from style import Style
from variable import Variable
from utilities import mark_pot, mark_preliminary

class ROCCurve(SpineArtist):
    """
    A class to create a Receiver Operating Characteristic (ROC) curve
    given a set of variables that represent discriminant scores and a
    set of classes that represent the true class labels.

    The ROC curve is a graphical representation of the trade-off between
    the true positive rate and the false positive rate for every possible
    cut-off value. The area under the curve (AUC) is a measure of the
    performance of the model. The closer the AUC is to 1, the better the
    model's discriminatory power.

    Attributes
    ----------
    _title : str
        The title of the artist. This will be placed at the top of the
        axis assigned to the artist.
    _categories : dict
        A dictionary mapping the category key to the category name.
    _labels : str
        The name of the field that represents the truth labels.
    _pos_label : list
        The list of labels that correspond to the positive class. The
        correspondence is one-to-one with the number of groups.
    _discriminant_scores : dict
        A dictionary of variables that represent discriminant scores.
        The keys are the category enumeration values and the values are
        the Variable objects that implement the discriminant score for
        the corresponding category.
    _background : list
        A list of categories that must be included in the ROC curve to fully
        populate the FPR calculation.
    _all_variable_names : list
        A list of all the variable names used as discriminant scores in
        the full ROC curve set.
    _groups : list
        A list of the group names that comprise the set of ROC curves.
    _roc_data : dict
        A dictionary of DataFrames that store the ROC curve data for each
        group in the set of ROC curves. The keys are the group names and
        the values are the DataFrames that store the score (discriminant)
        and class (label) values for the ROC curve of the corresponding
        group.
    _samples : list
        A list of Samples that will be used to draw the ROC curves.
    """
    def __init__(self, categories, labels, pos_label, discriminant_scores,
                 background=None, title=None):
        """
        Parameters
        ----------
        categories : dict
            A dictionary mapping the category key to the category name.
        labels : str
            The name of the field that represents the truth labels.
        pos_label : list
            The list of labels that correspond to the positive class. The
            correspondence is one-to-one with the number of groups.
        discriminant_scores : dict
            A dictionary of Variables that represent discriminant scores.
            The keys are the group names and the values are the Variables
            that represent the discriminant scores for the corresponding
            group.
        background : list
            A list of categories that must be included in the ROC curve to fully
            populate the FPR calculation.
        title : str, optional
            The title of the artist. The default is None.

        Returns
        -------
        None.
        """
        super().__init__(title)
        self._categories = categories
        self._labels = labels
        self._pos_label = pos_label
        self._discriminant_scores = discriminant_scores
        self._background = background
        self._all_variable_names = [x._key for x in self._discriminant_scores.values()]
        self._groups = list(self._discriminant_scores.keys())
        self._roc_data = {g : pd.DataFrame({'score': [], 'class': []}) for g in self._groups}

    def add_sample(self, sample, is_ordinate):
        """
        Add a sample to the ROC curve collection. The ensemble of samples
        will be used to draw the set of ROC curves on the same axis.

        Parameters
        ----------
        sample : Sample
            The sample to add to the ROC curve collection.
        is_ordinate : bool
            A flag indicating whether the sample is the ordinate sample.
        
        Returns
        -------
        None.
        """
        super().add_sample(sample, is_ordinate)
        
        # Update the ROC curve data
        reqvars = self._all_variable_names + [self._labels]
        reqidx = {x: i for i, x in enumerate(reqvars)}
        data, _ = sample.get_data(reqvars)
        for group in self._groups:
            if self._background is None:
                selected = [x for x in self._categories.keys() if self._categories[x] == group]
            else:
                selected = self._background
            scores = pd.concat([data[x][reqidx[self._discriminant_scores[group]._key]] for x in selected])
            classes = pd.concat([data[x][reqidx[self._labels]] for x in selected])
            self._roc_data[group] = pd.concat([self._roc_data[group], pd.DataFrame({
                'score': scores.values,
                'class': classes.values.astype(int)
            })])

    def draw(self, ax, style=None):
        """
        Draw the set of ROC curve on the given axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to draw the ROC curve on.
        style : Style, optional
            The style to use when drawing the ROC curve. The default is None.
        
        Returns
        -------
        None.
        """
        for gi, (k, v) in enumerate(self._roc_data.items()):
            mask = (~np.isnan(v['class']) & ~np.isnan(v['score']))
            fpr, tpr, _ = roc_curve(v['class'][mask], v['score'][mask], pos_label=self._pos_label[gi])
            roc_auc = auc(fpr, tpr)
            ax.plot(fpr, tpr, label=f'{k} (AUC = {roc_auc:.2f})')
            ax.plot([0, 1], [0, 1], 'k--')
            ax.set_xlim([0.0, 1.0])
            ax.set_ylim([0.0, 1.0])
            ax.set_xlabel('False Positive Rate')
            ax.set_ylabel('True Positive Rate')
            ax.set_title(self._title)
            ax.legend(loc='lower right')
        
        if style.scilimits:
            ax.ticklabel_format(axis='y', scilimits=style.scilimits)
        if style.mark_pot:
            mark_pot(ax, self._exposure, style.mark_pot_horizontal)
        if style.mark_preliminary is not None:
            mark_preliminary(ax, style.mark_preliminary, hadj=0.035 if style.scilimits is not None else 0)