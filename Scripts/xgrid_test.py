#!/usr/bin/env python
"""
xgrid_test.py
"""

import xgrid_tools
import random
from xgrid_test_obj import Test
from time import sleep

def _main():
    runner = xgrid_tools.xgrid_runner()
    jobs = [xgrid_tools.xgrid_object_job(Test(x)) for x in random.sample(xrange(20),20)]
    [runner.add_job(x) for x in jobs]
    runner.run_jobs()
    while(not runner.all_done):
        print runner.completed_jobs
        
        for x in runner.completed_jobs:
            print x

if __name__ == "__main__":
    _main()
