#!/bin/bash
#PBS -A spectorl
#PBS -m abe
#PBS -q mesabi

# Format arguments
SCRIPT=`basename $0`
if [ "$#" -ne 1 ]; then
cat << EOF
Usage: qsub ${SCRIPT} -F \"<Command_File>\"
        The purpose of this script is to submit, as a task array, the commands found in the Command_File
        
EOF
exit 1
fi

if ! [ -e "$1" ]; then
    echo "Command file: ${1} not found" >&2
    exit 1
fi

CMD_LIST="$1"

CMD="$(sed "${PBS_ARRAYID}q;d" ${CMD_LIST})"

eval ${CMD} && sed "${PBS_ARRAYID}q;d" ${CMD_LIST}