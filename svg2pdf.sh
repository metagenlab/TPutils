#!/bin/bash

# convert svg 2 pdf using inkscape

# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 02.2015
# -----------------------------------------------------

for i in $@;do
	inkscape -z --export-pdf="$(basename $i .svg).pdf" $i
done
