#!/bin/bash

INPUTS=(
    "rootfiles_data-mc/gammaN-0nXn-25/draw_2025PbPb_BeamA-prompt.root,rootfiles_data-mc/gammaN-0nXn-25/dataset-cook_BeamA-prompt.root,rootfiles_data-mc/gammaN-0nXn-25/dataset-cook_BeamA-nonprompt.root gammaN-0nXn-25/hists__2025PbPb_BeamA-prompt_BeamA-nonprompt"
)

make read_datamc.exe || exit 

for input in "${INPUTS[@]}" ; do
    inputpars=($input)
    [[ ${1:-0} -eq 1 ]] && {
        set -x
        ./read_datamc.exe ${inputpars[0]} "rootfiles/"${inputpars[1]}
        set +x
    }
    
done
