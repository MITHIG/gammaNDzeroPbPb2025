#!/bin/bash

SAVE_PNG=1

TAG_TARGET_BINNING="b-default" ; TARGET_BINNING_Y='-2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.' ; TARGET_BINNING_PT='2., 5.' ;
# TAG_TARGET_BINNING="b-ptdiff" ; TARGET_BINNING_Y='-2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.' ; TARGET_BINNING_PT='2., 3., 4., 5.' ;

BINNINT_TAG='b-ycoarse'
SAMPLE_TAG='twodirs'

INPUTS=(
    # 'rootfiles_data-mc/0nXn-'$SAMPLE_TAG'-25_Dbdt/'$BINNINT_TAG'/draw_2025PbPb_BeamA-2025-prompt.root,rootfiles_data-mc/0nXn-'$SAMPLE_TAG'-25_Dbdt/'$BINNINT_TAG'/cook_BeamA-2025-prompt.root,rootfiles_data-mc/0nXn-'$SAMPLE_TAG'-25_Dbdt/'$BINNINT_TAG'/cook_BeamA-nonprompt.root 0nXn-'$SAMPLE_TAG'-25_Dbdt/'$BINNINT_TAG'/hists__2025PbPb_BeamA-2025-prompt_BeamA-nonprompt'
    'rootfiles_data-mc/0nXn-'$SAMPLE_TAG'-25_Dbdt/'$BINNINT_TAG'/draw_2025PbPb_BeamA-prompt.root,rootfiles_data-mc/0nXn-'$SAMPLE_TAG'-25_Dbdt/'$BINNINT_TAG'/cook_BeamA-prompt.root,rootfiles_data-mc/0nXn-'$SAMPLE_TAG'-25_Dbdt/'$BINNINT_TAG'/cook_BeamA-nonprompt.root 0nXn-'$SAMPLE_TAG'-25_Dbdt/'$BINNINT_TAG'/hists__2025PbPb_BeamA-prompt_BeamA-nonprompt'
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

    IFS=',' ; inputfilelist=($inputfiles) ; unset IFS ;
    for i in ${inputfilelist[@]} ; do ls $i ; done ;
    
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
            ./draw_fprompt.exe "rootfiles/"$itag_fit".root" $itag_draw $SAVE_PNG
        }
        
        inputs_collect=$inputs_collect',rootfiles/'$itag_draw'.root'
        output_collect=$output_collect'-'$var
    done

    inputs_collect=${inputs_collect#*,}
    echo "  "$inputs_collect
    output_collect=${tag_input/hists__/collect__}'_'${output_collect#*-}
    echo "  "$output_collect

    [[ ${4:-0} -eq 1 ]] && {
        ./collect_fprompt.exe "$inputs_collect" $output_collect $TAG_TARGET_BINNING "$TARGET_BINNING_Y" "$TARGET_BINNING_PT" $SAVE_PNG
    }
done
