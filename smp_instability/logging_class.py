# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------
Script Name: logging_classs.py

Description:
    This module provides logging configuration and a reusable
    error-handling decorator to support consistent and maintainable logging.

Author: Bastian Bergfeld
Email: bbergfeld@gmx.net
Date Created: Fri May 2 08:59:48 2025
-----------------------------------------------------------
"""

import logging
import functools


class LoggerConfig:
    """
    Configures logging for the application.

    The setup_logging method sets up console logging and optionally file logging.
    It ensures:
        - Consistent log formatting
        - Prevention of duplicate log handlers
        - Suppression of verbose logs from common third-party libraries
          (e.g., snowmicropyn, matplotlib, sklearn)
    """

    @staticmethod
    def setup_logging(log_to_file=False, log_filename="pipeline.log"):
        log_format = "%(asctime)s - %(levelname)s - %(message)s"

        # Clear existing handlers to prevent duplicate logs
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        logging.basicConfig(level=logging.INFO, format=log_format)

        if log_to_file:
            file_handler = logging.FileHandler(log_filename)
            file_handler.setFormatter(logging.Formatter(log_format))
            logging.getLogger().addHandler(file_handler)
            logging.info("File logging enabled: %s", log_filename)

        # Suppress logs from external libraries
        logging.getLogger("snowmicropyn").setLevel(logging.WARNING)
        logging.getLogger("matplotlib").setLevel(logging.WARNING)
        logging.getLogger("sklearn").setLevel(logging.WARNING)


def error_handling_decorator(func):
    """
    A decorator that wraps a function and automatically logs any exceptions
    that occur during its execution. This improves robustness and simplifies
    debugging by including full stack traces in the log output.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logging.error(f"Error in {func.__name__}: {e}", exc_info=True)
            return None
    return wrapper