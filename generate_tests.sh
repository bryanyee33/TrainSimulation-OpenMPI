#!/bin/bash

# Number of Stations
mkdir -p testcases/stations
python3 gen_test.py 30 10 10 500 55 500 > testcases/stations/30.in
python3 gen_test.py 50 10 10 500 55 500 > testcases/stations/50.in
python3 gen_test.py 70 10 10 500 55 500 > testcases/stations/70.in
python3 gen_test.py 90 10 10 500 55 500 > testcases/stations/90.in
python3 gen_test.py 110 10 10 500 55 500 > testcases/stations/110.in

# Max line length
mkdir -p testcases/line_length
python3 gen_test.py 100 10 10 500 50 500 > testcases/line_length/50.in
python3 gen_test.py 100 10 10 500 60 500 > testcases/line_length/60.in
python3 gen_test.py 100 10 10 500 70 500 > testcases/line_length/70.in
python3 gen_test.py 100 10 10 500 80 500 > testcases/line_length/80.in
python3 gen_test.py 100 10 10 500 90 500 > testcases/line_length/90.in
python3 gen_test.py 100 10 10 500 100 500 > testcases/line_length/100.in

# Max num trains
mkdir -p testcases/num_trains
python3 gen_test.py 30 10 10 500 15 500 > testcases/num_trains/500.in
python3 gen_test.py 30 10 10 600 15 500 > testcases/num_trains/600.in
python3 gen_test.py 30 10 10 700 15 500 > testcases/num_trains/700.in
python3 gen_test.py 30 10 10 800 15 500 > testcases/num_trains/800.in
python3 gen_test.py 30 10 10 900 15 500 > testcases/num_trains/900.in
python3 gen_test.py 30 10 10 1000 15 500 > testcases/num_trains/1000.in

# Number of Ticks
mkdir -p testcases/num_ticks
python3 gen_test.py 30 10 10 1000 15 500 > testcases/num_ticks/500.in
python3 gen_test.py 30 10 10 1000 15 600 > testcases/num_ticks/600.in
python3 gen_test.py 30 10 10 1000 15 700 > testcases/num_ticks/700.in
python3 gen_test.py 30 10 10 1000 15 800 > testcases/num_ticks/800.in
python3 gen_test.py 30 10 10 1000 15 900 > testcases/num_ticks/900.in
python3 gen_test.py 30 10 10 1000 15 1000 > testcases/num_ticks/1000.in