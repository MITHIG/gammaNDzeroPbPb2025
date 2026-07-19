#!/bin/bash

INPUTS=(
    "rootfiles_data-mc/gammaN-0nXn-25/draw_2025PbPb_BeamA-prompt.root,rootfiles_data-mc/gammaN-0nXn-25/dataset-cook_BeamA-prompt.root,rootfiles_data-mc/gammaN-0nXn-25/dataset-cook_BeamA-nonprompt.root gammaN-0nXn-25/hists__2025PbPb_BeamA-prompt_BeamA-nonprompt"
)

VARS=(
    "Dip3D"
    "Dip3Dsig"
)

make read_datamc.exe fit_fprompt.exe draw_fprompt.exe || exit 

for input in "${INPUTS[@]}" ; do
    inputpars=($input)
    inputfiles=${inputpars[0]}

    for var in ${VARS[@]} ; do
        itag=${inputpars[1]}"_"$var
        
        [[ ${1:-0} -eq 1 ]] && {
            set -x
            ./read_datamc.exe $inputfiles "rootfiles/"$itag $var
            set +x
        }
        itag_fit=${itag/hists__/fits__}
        [[ ${2:-0} -eq 1 ]] && {
            ./fit_fprompt.exe "rootfiles/"$itag".root" $itag_fit
        }
        itag_draw=${itag/hists__/draw__}
        [[ ${3:-0} -eq 1 ]] && {
            ./draw_fprompt.exe "rootfiles/"$itag_fit".root" $itag_draw
        }
    done    
done
