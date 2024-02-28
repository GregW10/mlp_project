#!/bin/zsh

sleep 1

echo 'sleep 1 exited'

sleep 2

echo 'sleep 2 exited'

echo "Hallelujah, the pope is fine."

echo 'finishing...'

gen -vne --evolutions 12 --epochs 1000 --iterations 100000 -o first > log_first 2> errors_first

gen -vne --evolutions 12 --epochs 1000 --iterations 100000 -o second > log_second 2> errors_second

gen -vne --evolutions 12 --epochs 1000 --iterations 100000 -o third > log_third 2> errors_third

gen -vne --evolutions 12 --epochs 1000 --iterations 100000 -o fourth > log_fourth 2> errors_fourth

echo "done"
