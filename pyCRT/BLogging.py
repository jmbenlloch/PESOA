"""
Base class for logging. It defines a handeler (by default output stream)
and the format of the logging messages. 
"""
import logging
ch = logging.StreamHandler()
#formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)