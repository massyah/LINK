import logging
import datetime
__loggers = {}

def get_logger(name):
	global __loggers

	if __loggers.get(name):
		return __loggers.get(name)

	else:
		print "Building a new logger",name
		logger=logging.getLogger(name)
		logger.setLevel(logging.DEBUG)
		now = datetime.datetime.now()
		handler=logging.StreamHandler()

		formatter = logging.Formatter('%(asctime)s %(levelname)s - %(filename)s - %(funcName)s -  %(message)s')
		handler.setFormatter(formatter)
		logger.addHandler(handler)
		__loggers.update(dict(name=logger))

		return logger

# if "logger" not in globals():
# 	logger = logging.getLogger('link_logger')
# 	logger.setLevel(logging.DEBUG)

# 	# while len(logger.handlers()) > 0:
# 	#  	logger.pop()

# 	# create console handler and set level to debug
# 	ch.setLevel(logging.DEBUG)

# 	# create formatter
# 	formatter = logging.Formatter('%(asctime)s - %(filename)s - %(funcName)s - %(message)s',"%Y-%m-%d %H:%M:%S")
# 	# formatter = logging.Formatter('%(asctime)s - %(message)s')
# 	# add formatter to ch
# 	ch.setFormatter(formatter)

# 	# add ch to logger
# 	logger.addHandler(ch)


# logger = logging.getLogger('link_logger')

logger=get_logger("LINK")
