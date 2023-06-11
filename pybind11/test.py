"""
##############################################################################
#
#   Test: interfacing from Python to C++
#   AUTHOR: Maciej_Bak
#   CONTACT: wsciekly.maciek@gmail.com
#
##############################################################################
"""

# imports
import time
import logging
import logging.handlers
from argparse import ArgumentParser, RawTextHelpFormatter
import numpy as np
import functions


def parse_arguments():
    """Parser of the command-line arguments."""
    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--verbosity",
        dest="verbosity",
        choices=("DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"),
        default="ERROR",
        help="Verbosity/Log level. Defaults to ERROR",
    )
    parser.add_argument(
        "-l", "--logfile", dest="logfile", help="Store log to this file."
    )
    return parser

##############################################################################

def main():
    """Main body of the script."""

    np.random.seed(0)

    x1 = functions.add2numbers(np.float64(10), y=-5.0)
    print(x1)

    x_array = np.random.rand(5)
    x2_array = functions.square_a_number(x_array)
    print(x2_array)

    A = np.random.rand(3, 3)
    Ainv = functions.inverse(A)
    print(Ainv)

    X = np.random.rand(5, 3)
    y = np.array([1, 2, 1, 2, 1])
    coefficients = functions.lm(X, y)
    print(coefficients)

    res, err = functions.gsl_deriv()
    print(res,err)

    boost_deriv = functions.boost_deriv(0.0)
    print(boost_deriv)

    mlpack = functions.mlpack_example()
    print(mlpack)

##############################################################################

if __name__ == "__main__":

    try:
        # parse the command-line arguments
        options = parse_arguments().parse_args()

        # set up logging during the execution
        formatter = logging.Formatter(
            fmt="[%(asctime)s] %(levelname)s - %(message)s",
            datefmt="%d-%b-%Y %H:%M:%S",
        )
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger = logging.getLogger("logger")
        logger.setLevel(logging.getLevelName(options.verbosity))
        logger.addHandler(console_handler)
        if options.logfile is not None:
            logfile_handler = logging.handlers.RotatingFileHandler(
                options.logfile, maxBytes=50000, backupCount=2
            )
            logfile_handler.setFormatter(formatter)
            logger.addHandler(logfile_handler)

        # execute the body of the script
        start_time = time.time()
        logger.info("Starting script")
        main()
        seconds = time.time() - start_time

        # log the execution time
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        logger.info(
            "Successfully finished in {hours}h:{minutes}m:{seconds}s",
            hours=int(hours),
            minutes=int(minutes),
            seconds=int(seconds) if seconds > 1.0 else 1,
        )
    # log the exception in case it happens
    except Exception as e:
        logger.exception(str(e))
        raise e
