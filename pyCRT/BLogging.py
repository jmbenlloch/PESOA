"""
Base class for logging. It defines a handeler (by default output stream)
and the format of the logging messages. 
"""

import logging
ch = logging.StreamHandler()
#formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)

class Debug:
	"""
	A very simple class to control stepping
	"""
	def __init__(self,name,switch=1):
		self.name = name
		self.switch = switch

	def Wait(self):
		if self.switch == 1:
			raw_input("Press a key...")
	def SwitchOn(self):
		self.switch = 1
	def SwitchOff(self):
		self.switch = 0

