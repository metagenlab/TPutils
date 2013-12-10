#!/bin/bash

for i in $@;do
	inkscape -z --export-pdf="$(basename $i .svg).pdf" $i
done
