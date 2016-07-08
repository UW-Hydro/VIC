import time

from vic.vic import ffi
from vic import lib as vic_lib


def test_vic_timers():

    timer = ffi.new('timer_struct *')
    sleeptime = 1.5
    delta = 0.1

    # init sets delta to zero
    vic_lib.timer_init(timer)

    # start gets current time
    vic_lib.timer_start(timer)

    # sleep
    time.sleep(sleeptime)

    # stop pauses the timer
    vic_lib.timer_stop(timer)
    assert timer[0].delta_wall >= sleeptime
    assert timer[0].delta_wall < sleeptime + delta
    assert timer[0].delta_cpu >= 0.

    # start the timer again
    vic_lib.timer_continue(timer)

    # sleep again
    time.sleep(sleeptime)

    # stop after the lap time sleep
    vic_lib.timer_stop(timer)
    assert timer[0].delta_wall >= 2 * sleeptime
    assert timer[0].delta_wall < 2 * (sleeptime + delta)
