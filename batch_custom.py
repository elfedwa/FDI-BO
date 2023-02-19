"""
Running a batch optimization with a custom predictor.


Optimizing the 2D Rosenbrock function, which is a 2D
function with one objective to be minimized. There
are no Z descriptors so we use only the X coordinates
for learning.


We show two examples here:

1. Running a batch optimization with a builtin predictor.
2. Using your own custom predictor while still using
    batch optimization.

Change the USE_CUSTOM_PREDICTOR variable False
to use the builtin predictor.


See the documentation for more information on batch
optimization and how it runs.
"""
import time

import random

import numpy as np
from fireworks.core.firework import FireTaskBase, Firework, FWAction, Workflow
from fireworks.core.launchpad import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from fireworks.utilities.fw_utilities import explicit_serialize

from rocketsled.control import MissionControl
from rocketsled.task import OptTask
from rocketsled.utils import split_xz

import datetime
import matplotlib
matplotlib.use('WebAgg')
from fireworks.features.multi_launcher import launch_multiprocess

# Setting up the FireWorks LaunchPad
launchpad = LaunchPad(name="rsled")
opt_label = "opt_default"
db_info = {"launchpad": launchpad, "opt_label": opt_label}
x_dim = [(-5.0, 5.0), (-5.0, 5.0)]

@explicit_serialize
class BRANINTask(FireTaskBase):
    _fw_name = "BRANINTask"
    def run_task(self, fw_spec):

        x = fw_spec["_x"]
        y = (x[1]-(5.1/(4*np.pi**2))*x[0]**2+5*x[0]-6)**2+10*(1-1/(8*np.pi))*np.cos(x[0])+10

        return FWAction(update_spec={"_y": y})

@explicit_serialize
class RosenbrockTask(FireTaskBase):
    _fw_name = "RosenbrockTask"

    def run_task(self, fw_spec):
        x = fw_spec["_x"]
        y = (1 - x[0]) ** 2 + 100 * (x[1] - x[0] ** 2) ** 2
        return FWAction(update_spec={"_y": y})

def wf_creator_BRANIN(x):
    spec = {"_x": x}
    # ObjectiveFuncTask writes _y field to the spec internally.
    firework1 = Firework([BRANINTask(), OptTask(**db_info)], spec=spec)
    return Workflow([firework1])
def wf_creator_rosenbrock(x):
    spec = {"_x": x}
    # ObjectiveFuncTask writes _y field to the spec internally.
    firework1 = Firework([RosenbrockTask(), OptTask(**db_info)], spec=spec)
    return Workflow([firework1])


if __name__ == "__main__":
    mc = MissionControl(**db_info)
    todays_date=str(datetime.datetime.now().year)+'-'+str(datetime.datetime.now().month).zfill(2)+'-'+str(datetime.datetime.now().day).zfill(2)
    print("todays date",todays_date)
    launchpad.reset(password=todays_date, require_password=True)
    iteration_size=30
    #########################BATCH SIZE A ###################

    # batch_size_a = 1
    # mc.reset(hard=True)
    # mc.configure(
    #     wf_creator=wf_creator_rosenbrock,
    #     dimensions=x_dim,
    #     predictor="GaussianProcessRegressor",
    #     batch_size=batch_size_a,
    #     acq='ei',
    # )
    # for bs in range(batch_size_a):
    #     launchpad.add_wf(
    #         wf_creator_rosenbrock(
    #             [np.random.uniform(-5, 5), np.random.uniform(-5, 5)]
    #         )
    #     )
    # batch_of_a_initial = time.time()
    # rapidfire(launchpad, nlaunches=iteration_size, sleep_time=0)
    # batch_of_a_final = time.time()
    # time_a=batch_of_a_final - batch_of_a_initial
    # mc.results(batch_size_a,iteration_size,time_a,[1,1])
    #
    # plt = mc.plot()

    ######################### BATCH SIZE B ###################

    batch_size_b = 20
    launchpad.reset(password=todays_date, require_password=True)
    mc.reset(hard=True)
    mc.configure(
        wf_creator=wf_creator_rosenbrock,
        dimensions=x_dim,
        predictor="GaussianProcessRegressor",
        batch_size=batch_size_b,
        acq='ei',
    )

    for bs in range(batch_size_b):
        hh=np.random.uniform(-5, 5)
        jj=np.random.uniform(-5, 5)
        print("SEE IF IT MATCHES ",hh,jj)
        launchpad.add_wf(
            wf_creator_rosenbrock(
                [hh, jj]
            )
        )
    batch_of_b_initial = time.time()
    # rapidfire(launchpad, nlaunches=(batch_size_b*iteration_size), sleep_time=0)
    launch_multiprocess(launchpad,FWorker(),nlaunches=(numbatch*batch_size), sleep_time=0, loglvl="CRITICAL",num_jobs=batch_size)

    batch_of_b_final = time.time()
    time_b=batch_of_b_final - batch_of_b_initial

    mc.results(batch_size_b,iteration_size,time_b,[1,1])
    plt = mc.plot()

    #########################################################

    print("Time for Batch of {} vs {} is {} and {}".format(batch_size_a,batch_size_b, time_a,time_b))
    plt.show()
