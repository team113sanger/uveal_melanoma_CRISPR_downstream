#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ROOT=$(dirname $SCRIPT_DIR)
cwd=$ROOT/results/$3/$4/

OSTYPE=$(uname -s)

if [[ $OSTYPE == "Linux" ]]
then
        source ~/.bashrc
else
        source ~/.zshrc
fi

conda activate mageckenv
mkdir -p $cwd
cd $cwd

mageck mle -k $ROOT/$1 -d $ROOT/$2 -n $3

rm $ROOT/$2
