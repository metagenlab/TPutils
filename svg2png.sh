#!/bin/bash

# convert svg 2 png using inkscape

# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 02.2015
# -----------------------------------------------------

for i in $@;do
	inkscape -z --export-png="$(basename $i .svg).png" $i
done
