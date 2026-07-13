#!/bin/bash

INPUTS_DATA=( # lumi is nb-1 - directly from brilcalc
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260426-yrefmva_PbPbUPC_HIForward_Dpt-2_Dsize_24PD.root;2025 PbPb (5.36 TeV);2025PbPb;0.06036" # 2025
)
CUTEVTS=(
    "1;#gammaN (Xn0n);gammaN-0nXn-25"
    # "2;N#gamma (0nXn);Ngamma-0nXn-25"
)
INPUTS_MC=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260426-yrefmva_HiForest_260328_prompt_GNucleusToD0-PhotonBeamA_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.85]{YTHIA}8#scale[0.5]{ }#gammaN (5.36 TeV);BeamA"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260426-yrefmva_HiForest_260328_prompt_GNucleusToD0-PhotonBeamB_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.85]{YTHIA}8#scale[0.5]{ }#gammaN (5.36 TeV);BeamB"
)

make save_datasets.exe cook_datasets.exe fit_datasets.exe || exit 1

for cutevtstr in "${CUTEVTS[@]}" ; do
    IFS=';' ; cutevttags=($cutevtstr) ; unset IFS ; cutevt_tag="${cutevttags[2]}" ;
    cut_tag=$cutevt_tag
    echo -e "\033[33m"$cut_tag"\033[0m"

    ## Loop data
    for input_data in "${INPUTS_DATA[@]}" ; do

        IFS=';' ; input_data_tags=($input_data) ; unset IFS ; data_tag=${input_data_tags[2]} ; data_lumi=${input_data_tags[3]} ;
        echo -e "\033[33;2m"$cut_tag" / \033[0m\033[33m"$data_tag"\033[0m"
        
        itag_dataset_data="rootfiles/"$cut_tag"/dataset_"$data_tag
        [[ ${1:-0} -eq 1 ]] && {
            ./save_datasets.exe "$input_data" $itag_dataset_data "$cutevtstr" 0
        }

        ## Loop MC
        for input_mc in "${INPUTS_MC[@]}" ; do
            IFS=';' ; input_mc_tags=($input_mc) ; unset IFS ; mc_tag=${input_mc_tags[2]} ;
            [[ ($cutevt_tag == *gammaN* && $mc_tag == *BeamB*) || ($cutevt_tag == *Ngamma* && $mc_tag == *BeamA*) ]] && continue
            echo -e "\033[33;2m"$cut_tag" / "$data_tag" / \033[0m\033[33m"$mc_tag"\033[0m"

            itag_dataset_mc="rootfiles/"$cut_tag"/dataset_"$mc_tag
            [[ ${2:-0} -eq 1 ]] && {
                ./save_datasets.exe "$input_mc" $itag_dataset_mc "$cutevtstr" 1
            }

            itag_cook="rootfiles/"$cut_tag"/dataset-cook_"$data_tag"_"$mc_tag
            [[ ${3:-0} -eq 1 ]] && {
                ./cook_datasets.exe $itag_dataset_data".root" $itag_dataset_mc".root" $itag_cook
            }

            itag_fit=$cut_tag"/roofit_"$data_tag"_"$mc_tag
            [[ ${4:-0} -eq 1 ]] && {
                ./fit_datasets.exe $itag_cook".root" $itag_fit
            }

        done

    done

done
