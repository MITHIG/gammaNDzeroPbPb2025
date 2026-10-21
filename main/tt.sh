#!/bin/bash 

TAG_BINNING="b-default"

fitopt="3G-Peaky;;"

INPUTS_DATA=( "2025PbPb" )
CUTEVTS=(
    "0nXn-gammaN-25"
    # "0nXn-Ngamma-25"
)
INPUTS_TEMPLATE=( "2025-SoftQCD-BeamA" "2025-SoftQCD-BeamB" )
CUTDS=( "Dbdt-gammaN" "Dbdt-Ngamma" )

make toy_save.exe toy_draw.exe || exit 1

for cutevt_tag in "${CUTEVTS[@]}" ; do
    for cutd_tag in "${CUTDS[@]}" ; do
        [[ ($cutd_tag == *gammaN* && $cutevt_tag == *Ngamma*) || ($cutd_tag == *Ngamma* && $cutevt_tag == *gammaN*) ]] && continue

        cut_tag=$cutevt_tag'_'$cutd_tag ##

        for template_tag in "${INPUTS_TEMPLATE[@]}" ; do
            [[ ($cutevt_tag == *gammaN* && $template_tag == *BeamB*) || ($cutevt_tag == *Ngamma* && $template_tag == *BeamA*) || ($cutevt_tag == *0n0n* && $template_tag == *BeamB*) ]] && continue

            for data_tag in "${INPUTS_DATA[@]}" ; do
                itag_data_fit=$cut_tag'/'$TAG_BINNING'/fithist_'$data_tag'_'$template_tag$fit_tag ##

                ls 'rootfiles/'$itag_data_fit'.root'
                [[ ${1:-0} -eq 1 ]] && {
                ./toy_save.exe 'rootfiles/'$itag_data_fit'.root' $itag_data_fit'/toy_save_boostrap' 0 500 # 0: boostrap
                ./toy_save.exe 'rootfiles/'$itag_data_fit'.root' $itag_data_fit'/toy_save_closure' 1 500 # 1: closure
                }
                itag_toy=$itag_data_fit'/toy'
                [[ ${2:-0} -eq 1 ]] && {
                    ./toy_draw.exe 'rootfiles/'$itag_data_fit'/toy_save_boostrap.root,rootfiles/'$itag_data_fit'/toy_save_closure.root' $itag_toy
                }
            done
        done
    done
done 
