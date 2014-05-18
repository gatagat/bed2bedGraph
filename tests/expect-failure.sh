#!/bin/sh

name=$1
bin=$2
input=$3

if "$bin" < "$input" > /dev/null 2>&1; then
	echo "$(tput setaf 1)Test $name did not fail!$(tput sgr0)"
else
	echo "$(tput setaf 2)Test $name failed as required.$(tput sgr0)"
fi
