#!/bin/bash

"""
convert (not too big) csv tables into pdf using csv2latex and xelatex
"""



csv2latex --separator t --colorrows 0.75 --reduce 5 $1 > ${1/.*/.temp}

sed -i 's/\\begin{document}/\\usepackage{pdflscape}\n&/' ${1/.*/.temp}
sed -i 's/\\begin{document}/\\pagenumbering{gobble}\n&/' ${1/.*/.temp}
sed -i 's/\\begin{document}/\\usepackage{pdflscape}\n&/' ${1/.*/.temp}


sed -i 's/\\begin{document}/\\usepackage[margin=0.85in]{geometry}\n&/' ${1/.*/.temp}
title=${1/.*/}
title=`echo $title | sed 's/_/\\_/g'`

echo title $title

sed -i 's/\\end{document}/\\end{landscape}\n&/' ${1/.*/.temp}
sed -i 's/\\begin{tabular}/\\begin{landscape}\n&/' ${1/.*/.temp}
#sed -i 's/\\begin{tabular}/\\verb|'$title'|\\newline\\newline\n&/' ${1/.*/.temp}
sed -i 's/\\begin{tabular}/'$title'\\newline\\newline\n&/' ${1/.*/.temp}
sed -i 's/\\begin{tabular}/\\hspace{2cm}\n&/' ${1/.*/.temp}
sed -i 's/_/ /g' ${1/.*/.temp}

xelatex ${1/.*/.temp}
#rm ${1/.*/.temp}
#rm ${1/.*/.aux}
#rm ${1/.*/.log}
