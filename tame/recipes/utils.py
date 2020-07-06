"""Misc utils for recipes"""
import argparse

class TypeDefaultFormat(
        argparse.RawTextHelpFormatter,
        argparse.MetavarTypeHelpFormatter):
    pass
