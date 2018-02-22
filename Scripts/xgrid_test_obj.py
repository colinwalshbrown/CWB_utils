#!/usr/bin/env python
"""
xgrid_test_obj.py
"""
import random
from time import sleep

class Test(object):
    
    def __init__(self,count):
        self.count = count
        self.done = None

    def hold(self):
        sleep(self.count)
        self.done = True

    def runXgrid(self):
        self.hold()
