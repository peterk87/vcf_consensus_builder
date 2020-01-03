import logging

LOG_FORMAT = '\033[2;36m{asctime}\033[0;1m {levelname:<6}\033[0m {message} \033[2m[{name}:{funcName}:{lineno}]\033[0m'


def init_console_logger(logging_verbosity=0):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    if logging_verbosity > (len(logging_levels) - 1):
        logging_verbosity = 3
    lvl = logging_levels[logging_verbosity]
    logging.basicConfig(format=LOG_FORMAT, level=lvl, style='{')
