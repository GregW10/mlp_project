#!/bin/zsh

if [ $# != 1 ]; then
	g++ -std=c++20 -O3 gen.cpp -o gen
else
	g++ -ferror-limit=$1 -std=c++20 -O3 gen.cpp -o gen
fi
