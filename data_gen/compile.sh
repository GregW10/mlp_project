#!/bin/zsh

if [ $# != 1 ]; then
	g++ -std=c++20 -O3 main.cpp -o main
else
	g++ -ferror-limit=$1 -std=c++20 -O3 main.cpp -o main
fi
