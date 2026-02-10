#!/usr/bin/env python3
"""Wrapper — actual implementation moved to src/data/build_id_map.py"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from src.data.build_id_map import main
main()
