"""Misc utils for recipes"""
import argparse

class TypeDefaultFormat(
        argparse.ArgumentDefaultsHelpFormatter,
        argparse.MetavarTypeHelpFormatter):
    pass
