#!/bin/bash

for i in $@;do
	inkscape -z --export-png="$(basename $i .svg).png" $i
done
