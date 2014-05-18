#!/bin/sh

name=$1
bin=$2
input=$3
expected=$4

"$bin" < "$input" > "${input}.test"
if diff -q "$expected" "${input}.test" > /dev/null; then
	echo "$(tput setaf 2)Test $name passed.$(tput sgr0)"
else
	echo "$(tput setaf 1)Test $name failed!$(tput sgr0)"
	echo "Differences:"
	diff "$expected" "${input}.test" | head
fi
rm -f "${input}.test" > /dev/null
