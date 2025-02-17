#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ROOT=$(dirname $SCRIPT_DIR)

OSTYPE=$(uname -s)

if [[ $OSTYPE == "Linux" ]]
then
	source ~/.bashrc
else
	source ~/.zshrc
fi

conda activate bagel

python3 $ROOT/src/BAGEL.py bf -i $ROOT/results/$1 -o $ROOT/results/$2 -c 1 -e $ROOT/data/CEGv2.txt -n $ROOT/data/NEGv1.txt

rm $ROOT/results/$1
