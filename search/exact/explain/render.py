"""Shim: record loading and drawing for the explain/ scripts (from search/regularize/jstyle.py)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'regularize'))
import jstyle
from jstyle import load, draw
