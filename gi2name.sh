#!/bin/bash


# basic script to grep gi from gi_taxid_prot.dmp db
# TODO should be replaced by conversion tool using eutils
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------

while getopts ":g:h" opt; do
  case $opt in
    g)
      echo "-g was triggered, Parameter: $OPTARG" >&2
      gi=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
    h)
      echo -e "utilisation:\ngi2name.sh -g [gi number]"
      exit 1
  esac
done

taxdb_dir="/home/trestan/Dropbox/taxdb"

taxid=$(grep "^$gi" $taxdb_dir/gi_taxid_prot.dmp | cut -f2)
name=$(grep -w "^$taxid" $taxdb_dir/taxdump/names.dmp | grep "scientific name" | cut -d "|" -f 2 | sed 's/\t//')
echo $name
