#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 14:23:06 2024

@author: hmfausey
"""

import pytest
import numpy as np
import sys
sys.path.append('../')

import power_law as pl

def test_power_law():
    output = pl.fv_emit(np.array([0.5, 1, 2]), 1, 1, 1)
    expected = [0.5, 1, 2]
    np.testing.assert_allclose(output, expected)