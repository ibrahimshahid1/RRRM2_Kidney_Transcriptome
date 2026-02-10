#!/usr/bin/env python3
"""Wrapper — actual implementation moved to src/markers/discover_dct.py"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from src.markers.discover_dct import main
main()
