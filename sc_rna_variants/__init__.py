import sys
import logging

# basic verifications
if sys.version_info.major != 3 or sys.version_info.minor < 7:
    sys.exit("Wrong python version: use python version >= 3.7 and < 4")


from importlib import reload
reload(logging)
logger = logging.getLogger(__name__)


def config_logging(log_file=None):
    logging_args = {
        "level": logging.INFO,
        "filemode": 'w',
        "format": '%(asctime)s - %(name)s - %(levelname)s: %(message)s',
        "datefmt": '%m/%d/%Y %H:%M:%S'
    }
    # If a log file is given, add it and changes logging level to debug
    if log_file is not None:
        logging_args["filename"] = log_file
        logging_args["level"] = logging.DEBUG

    print(logging_args)
    logging.basicConfig(**logging_args)

