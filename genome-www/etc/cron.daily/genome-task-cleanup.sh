#!/bin/sh

find /var/cache/genome/task_runner_output/ -ctime +7 -exec rm -f {} \;
