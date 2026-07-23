#!/bin/bash

INPUTS=(
    "rootfiles_data-mc/gammaN-0nXn-25/draw_2025PbPb_BeamA-prompt.root,rootfiles_data-mc/gammaN-0nXn-25/dataset-cook_BeamA-prompt.root,rootfiles_data-mc/gammaN-0nXn-25/dataset-cook_BeamA-nonprompt.root gammaN-0nXn-25/hists__2025PbPb_BeamA-prompt_BeamA-nonprompt"
    "rootfiles_data-mc/twodirs-0nXn-25/draw_2025PbPb_BeamA-prompt.root,rootfiles_data-mc/twodirs-0nXn-25/dataset-cook_BeamA-prompt.root,rootfiles_data-mc/twodirs-0nXn-25/dataset-cook_BeamA-nonprompt.root twodirs-0nXn-25/hists__2025PbPb_BeamA-prompt_BeamA-nonprompt"
    "rootfiles_data-mc/gammaN-0nXn-25/draw_2025PbPb_BeamA-prompt_ycoarse.root,rootfiles_data-mc/gammaN-0nXn-25/dataset-cook_BeamA-prompt_ycoarse.root,rootfiles_data-mc/gammaN-0nXn-25/dataset-cook_BeamA-nonprompt_ycoarse.root gammaN-0nXn-25/hists__2025PbPb_BeamA-prompt_BeamA-nonprompt_ycoarse"
    "rootfiles_data-mc/twodirs-0nXn-25/draw_2025PbPb_BeamA-prompt_ycoarse.root,rootfiles_data-mc/twodirs-0nXn-25/dataset-cook_BeamA-prompt_ycoarse.root,rootfiles_data-mc/twodirs-0nXn-25/dataset-cook_BeamA-nonprompt_ycoarse.root twodirs-0nXn-25/hists__2025PbPb_BeamA-prompt_BeamA-nonprompt_ycoarse"
)

VARS=(
    "Dip3D"
    "Dip3Dsig"
)

make read_datamc.exe fit_fprompt.exe draw_fprompt.exe collect_fprompt.exe || exit 

for input in "${INPUTS[@]}" ; do
    inputpars=($input)
    inputfiles=${inputpars[0]}
    tag_input=${inputpars[1]}

    inputs_collect=
    output_collect=
    for var in ${VARS[@]} ; do
        itag=$tag_input"_"$var
        
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
        
        inputs_collect=$inputs_collect',rootfiles/'$itag_draw'.root'
        output_collect=$output_collect'-'$var
    done

    inputs_collect=${inputs_collect#*,}
    echo $inputs_collect
    output_collect=${tag_input/hists__/collect__}'_'${output_collect#*-}
    echo $output_collect

    [[ ${4:-0} -eq 1 ]] && {
        ./collect_fprompt.exe "$inputs_collect" $output_collect
    }
done
