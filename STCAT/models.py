import os
import pickle
import numpy as np
from typing import Optional 
from scipy.special import expit
from sklearn import __version__ as skv
from . import logger

models_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

class Model():
    """
    Class that wraps the logistic Classifier and the StandardScaler.

    Parameters
    ----------
    clf
        A logistic Classifier incorporated in the loaded model.
    scaler
        A StandardScaler incorporated in the loaded model.
    description
        Description of the model as a dictionary.

    Attributes
    ----------
    classifier
        The logistic Classifier incorporated in the loaded model.
    scaler
        The StandardScaler incorporated in the loaded model.
    description
        Description of the loaded model.
    """
    def __init__(self, clf, scaler, description):
        self.classifier = clf
        self.scaler = scaler
        self.description = description

    @staticmethod
    def load(model: Optional[str] = None):
        """
        Load the desired model.

        Parameters
        ----------
        model
            Model name specifying the model you want to load. Default to `'Immune_All_Low.pkl'` if not provided.
            To see all available models and their descriptions, use :func:`~celltypist.models.models_description`.

        Returns
        ----------
        :class:`~celltypist.models.Model`
            A :class:`~celltypist.models.Model` object.
        """
        if model in get_all_models():
            model = get_model_path(model)
        with open(model, "rb") as fh:
            try:
                pkl_obj = pickle.load(fh)
                return Model(pkl_obj['Model'], pkl_obj['Scaler_'], pkl_obj['description'])
            except Exception as exception:
                raise Exception(
                        f"🛑 Invalid model: {model}. {exception}")

    @property
    def cell_types(self) -> np.ndarray:
        """Get cell types included in the model."""
        return self.classifier.classes_

    @property
    def features(self) -> np.ndarray:
        """Get genes included in the model."""
        return self.classifier.features

    def __repr__(self):
        base = f"CellTypist model with {len(self.cell_types)} cell types and {len(self.features)} features"
        for x in ['date', 'details', 'source', 'version']:
            if self.description[x] != '':
                base += f"\n    {x}: {self.description[x]}"
        if len(self.cell_types) == 2:
            base += f"\n    cell types: {self.cell_types[0]}, {self.cell_types[1]}\n    features: {self.features[0]}, {self.features[1]}, ..., {self.features[-1]}"
        elif len(self.cell_types) == 3:
            base += f"\n    cell types: {self.cell_types[0]}, {self.cell_types[1]}, {self.cell_types[2]}\n    features: {self.features[0]}, {self.features[1]}, ..., {self.features[-1]}"
        else:
            base += f"\n    cell types: {self.cell_types[0]}, {self.cell_types[1]}, ..., {self.cell_types[-1]}\n    features: {self.features[0]}, {self.features[1]}, ..., {self.features[-1]}"
        return base

    def predict_labels_and_prob(self, indata, mode: str = 'best match', p_thres: float = 0.5) -> tuple:
        """
        Get the decision matrix, probability matrix, and predicted cell types for the input data.

        Parameters
        ----------
        indata
            The input array-like object used as a query.
        mode
            The way cell prediction is performed.
            For each query cell, the default (`'best match'`) is to choose the cell type with the largest score/probability as the final prediction.
            Setting to `'prob match'` will enable a multi-label classification, which assigns 0 (i.e., unassigned), 1, or >=2 cell type labels to each query cell.
            (Default: `'best match'`)
        p_thres
            Probability threshold for the multi-label classification. Ignored if `mode` is `'best match'`.
            (Default: 0.5)

        Returns
        ----------
        tuple
            A tuple of decision score matrix, raw probability matrix, and predicted cell type labels.
        """
        if skv.split('.')[0] != '0' and isinstance(indata, np.matrix):
            scores = self.classifier.decision_function(np.asarray(indata))
        else:
            scores = self.classifier.decision_function(indata)
        if len(self.cell_types) == 2:
            scores = np.column_stack([-scores, scores])
        probs = expit(scores)
        if mode == 'best match':
            return scores, probs, self.classifier.classes_[scores.argmax(axis=1)]
        elif mode == 'prob match':
            flags = probs > p_thres
            labs = np.array(['|'.join(self.classifier.classes_[np.where(x)[0]]) for x in flags])
            labs[labs == ''] = 'Unassigned'
            return scores, probs, labs
        else:
            raise ValueError(
                    f"🛑 Unrecognized `mode` value, should be one of `'best match'` or `'prob match'`")

    def write(self, file: str) -> None:
        """Write out the model."""
        obj = dict(Model = self.classifier, Scaler_ = self.scaler, description = self.description)
        file = os.path.splitext(file)[0] + '.pkl'
        with open(file, 'wb') as output:
            pickle.dump(obj, output)

    def extract_top_markers(self, cell_type: str, top_n: int = 10, only_positive: bool = True) -> np.ndarray:
        """
        Extract the top driving genes for a given cell type.

        Parameters
        ----------
        cell_type
            The cell type to extract markers for.
        top_n
            Number of markers to extract for a given cell type.
            (Default: 10)
        only_positive
            Whether to extract positive markers only. Set to `False` to include negative markers as well.
            (Default: `True`)

        Returns
        ----------
        :class:`~numpy.ndarray`
            A list of marker genes for the query cell type.
        """
        if cell_type not in self.cell_types:
            raise ValueError(
                    f"🛑 '{cell_type}' is not found. Please provide a valid cell type name")
        if len(self.cell_types) == 2:
            coef_vector = self.classifier.coef_[0] if cell_type == self.cell_types[1] else -self.classifier.coef_[0]
        else:
            coef_vector = self.classifier.coef_[self.cell_types == cell_type][0]
        if not only_positive:
            coef_vector = np.abs(coef_vector)
        return self.features[np.argsort(-coef_vector)][:top_n]

def get_model_path(file: str) -> str:
    """
    Get the full path to a file in the `models` folder.

    Parameters
    ----------
    file
        File name as a string.
        To see all available models and their descriptions, use :func:`~celltypist.models.models_description`.

    Returns
    ----------
    str
        A string of the full path to the desired file.
    """
    return os.path.join(models_path, f"{file}")

def get_all_models() -> list:
    """
    Get a list of all the available models.

    Returns
    ----------
    list
        A list of available models.
    """
    # download_if_required()
    available_models = []
    for model_filename in os.listdir(models_path):
        if model_filename.endswith(".pkl"):
            model_name = os.path.basename(model_filename)
            available_models.append(model_name)
    return available_models
