#!/bin/bash

sbatch  --job-name=plat -t 7-00:00:00 --cpus-per-task=16 --mem=50g --output=offinder.out \
./offinder.sh
