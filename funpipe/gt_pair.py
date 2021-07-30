import os
from math import ceil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from funpipe.utils import run
import subprocess as sp


class gt_pair:
    """ genotype pair """
    def __init__(self, gt1, gt2, na_ignore=False):
        """
        Parameters
        ----------
        gt1, gt2: pd.Series
        
        na_ignore: bool
            whether to ignore nan in comparison, default = False
            
        Example
        -------
        >>> gt1 = pd.Series([0, 1, 2, 0, 1, 2, 0, 1, 2, np.nan])
        >>> gt2 = pd.Series([0, 1, 2, 1, 0, 1, np.nan, np.nan, np.nan, np.nan])
        >>> gt = gt_pair(gt1, gt2).get_n_unique()
        >>> print(gt.n_total, gt.n_unique, gt.n_share)
        7 5 2
        >>> gt = gt_pair(gt1, gt2, na_ignore=True).get_n_unique()
        >>> print(gt.n_total, gt.n_unique, gt.n_share)
        5 3 2

        """
        # helper method to convert gt datatype
        def _gt_type(gt):
            gt_type = type(gt).__name__
            if gt_type == 'Series':
                return gt
            elif gt_type in ['list', 'ndarray']:
                return pd.Series(gt_type)
            else:
                raise ValueError("Input gt vector should be either pandas series, list or numpy ndarray.")
        
        
        self.gt1 = _gt_type(gt1)
        self.gt2 = _gt_type(gt2)
        self.na_ignore = na_ignore
        self.n_total = None
        self.n_share = None
        self.n_unique = None
        self._not_both_ref = None

    def get_n_total(self):
        """ compute total number of non-monomorphic sites between two samples
        
        Returns
        -------
        int
            total number of non-monomorphic sites between two samples
            
        Example
        -------
        >>> gt1 = pd.Series([0, 1, 2, 0, 1, 2, 0, 1, 2, np.nan])
        >>> gt2 = pd.Series([0, 1, 2, 1, 0, 1, np.nan, np.nan, np.nan, np.nan])
        >>> gt_pair(gt1, gt2).get_n_total().n_total
        7
        >>> gt_pair(gt1, gt2).get_n_total().n_total
        5

        """
        if self.na_ignore:
            self.n_total = ((self.gt1+self.gt2).fillna(0)
                            .map(lambda x: 1 if x !=0 else 0).sum())
        else:
            self.n_total = ((self.gt1.fillna(0) + self.gt2.fillna(0))
                            .map(lambda x: 1 if x !=0 else 0).sum())
        return self.n_total

    def get_n_share(self):
        """ Compare genotypes between two columns within a VCF, and report shared
        variants between the two samples.

                           A B
        for example: site1 1 .
                     site2 . 1
                     site3 1 1

        The shared variant here will be 1 (site3).

        Returns
        -------
        int
            number of shared sites

        Example
        -------
        >>> gt1 = pd.Series([0, 1, 2, 0, 1, 2, 0, 1, 2, np.nan])
        >>> gt2 = pd.Series([0, 1, 2, 1, 0, 1, np.nan, np.nan, np.nan, np.nan])
        >>> gt_pair(gt1, gt2).get_n_share().n_share
        2

        Note
        ----

        This method is also cross-validated with GenotypeConcordance in GATK and
        bcftools stats.

        NaN will not be matched to any others.

        """
        # is a polymorphic site (non-reference sites)
        is_poly = (self.gt1 + self.gt2).map(lambda x: 1 if x != 0 else 0)
        # two sites are similar, include reference
        is_same = (self.gt1 - self.gt2).map(lambda x: 1 if x == 0 else 0)

        # number of shared alleles
        self.n_share = int((is_poly * is_same).sum())
        
        return self.n_share


    def get_n_unique(self):
        """
        Unique variants here mean a site that are private to either sample.

                               A B
        for example: site1 1 .
                     site2 . 1
                     site3 1 1

        The unique variants here will be 2 (site1 and site2). If ignore NA,
        the unique variants will be 0 (site1 and 2 will not be considered here).


        Returns
        -------
        int
            number of unique sites

        Example
        -------
        >>> gt1 = pd.Series([0, 1, 2, 0, 1, 2, 0, 1, 2, np.nan])
        >>> gt2 = pd.Series([0, 1, 2, 1, 0, 1, np.nan, np.nan, np.nan, np.nan])
        >>> gt_pair(gt1, gt2).get_n_unique()
        5
        >>> gt_unique(gt1, gt2)
        3

        """
        if self.n_total is None:
            self = self.get_n_total()
        if self.n_share is None:
            self = self.get_n_share()
        self.n_unique = self.n_total - self.n_share
        
        return self.n_unique
