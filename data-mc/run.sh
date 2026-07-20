#!/bin/bash

echo -e "usage: \033[32m./run.sh\033[0m \033[36;1m[data ds]\033[0m \033[2m(1: save, 2: cook, 3: save + cook)\033[0m \033[36;1m[mc ds]\033[0m \033[2m(1: save, 2: cook, 3: save + cook)\033[0m \033[36;1m[fit] [splot]\033[0m \033[2m(1: make, 2: draw, 3: make + draw)\033[0m"

TAG_BIN="" ; BINNING='' ; # default
# TAG_BIN="_yincl" ; BINNING='-2,2' ;

## !! do not give more than one data ##
INPUTS_DATA=( # lumi is nb-1 - directly from brilcalc
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260426-yrefmva_PbPbUPC_HIForward_Dpt-2_Dsize_24PD.root;2025 PbPb;2025PbPb" # 2025
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260426-yrefmva_HiForest_260328_prompt_GNucleusToD0-PhotonBeamA_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.85]{YTHIA}8#scale[0.5]{ }#gammaN;BeamAclose" # !! to replace by the one without gmatch filter
)
CUTEVTS=(
    "1;#gammaN (Xn0n);gammaN-0nXn-25"
    # "2;N#gamma (0nXn);Ngamma-0nXn-25"
    # "3;0nXn + Xn0n;twodirs-0nXn-25" # have not add y reflection
)
INPUTS_MC=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_prompt_GNucleusToD0-PhotonBeamA_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;Prompt P#scale[0.8]{YTHIA}8#scale[0.5]{ }#gammaN;BeamA-prompt"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_prompt_GNucleusToD0-PhotonBeamB_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;Prompt P#scale[0.8]{YTHIA}8#scale[0.5]{ }#gammaN;BeamB-prompt"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_nonprompt_GNucleusToD0-PhotonBeamA_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;Nonprompt P#scale[0.8]{YTHIA}8#scale[0.5]{ }#gammaN;BeamA-nonprompt"
)

make save_datasets.exe cook_datasets.exe fit_datasets.exe splot_datasets.exe draw_datasets.exe || exit 1

for cutevtstr in "${CUTEVTS[@]}" ; do
    IFS=';' ; cutevttags=($cutevtstr) ; unset IFS ; cutevt_tag="${cutevttags[2]}" ;
    cut_tag=$cutevt_tag
    echo -e "\033[33m"$cut_tag"\033[0m"

    ## Loop data
    for input_data in "${INPUTS_DATA[@]}" ; do

        IFS=';' ; input_data_tags=($input_data) ; unset IFS ; data_tag=${input_data_tags[2]} ;
        echo -e "\033[33;2m"$cut_tag" / \033[0m\033[33m"$data_tag"\033[0m"
        
        itag_dataset_data="rootfiles/"$cut_tag"/dataset_"$data_tag
        [[ ${1:-0} -eq 1 ||  ${1:-0} -eq 3 ]] && {
            ./save_datasets.exe "$input_data" $itag_dataset_data "$cutevtstr" 0
        }
        itag_cook_data=${itag_dataset_data/\/dataset_/\/dataset-cook_}$TAG_BIN
        [[ ${1:-0} -eq 2 || ${1:-0} -eq 3 ]] && {
            ./cook_datasets.exe $itag_dataset_data".root" $itag_cook_data "$BINNING"            
        }

        ## Loop MC
        for input_mc in "${INPUTS_MC[@]}" ; do
            IFS=';' ; input_mc_tags=($input_mc) ; unset IFS ; mc_tag=${input_mc_tags[2]} ;
            [[ ($cutevt_tag == *gammaN* && $mc_tag == *BeamB*) || ($cutevt_tag == *Ngamma* && $mc_tag == *BeamA*) ]] && continue
            echo -e "\033[33;2m"$cut_tag" / "$data_tag" / \033[0m\033[33m"$mc_tag"\033[0m"

            itag_dataset_mc="rootfiles/"$cut_tag"/dataset_"$mc_tag #
            [[ ${2:-0} -eq 1 || ${2:-0} -eq 3 ]] && {
                ./save_datasets.exe "$input_mc" $itag_dataset_mc "$cutevtstr" 1
            }

            itag_cook_mc=${itag_dataset_mc/\/dataset_/\/dataset-cook_}$TAG_BIN
            [[ ${2:-0} -eq 2 || ${2:-0} -eq 3 ]] && {
                ./cook_datasets.exe $itag_dataset_mc".root" $itag_cook_mc "$BINNING"            
            }

            itag_roofit=$cut_tag"/roofit_"$data_tag"_"$mc_tag$TAG_BIN #
            [[ ${3:-0} -eq 1 ]] && {
                ./fit_datasets.exe $itag_cook_data".root" $itag_cook_mc".root" $itag_roofit
            }

            itag_splot=${itag_roofit/\/roofit_/\/splot_}
            [[ ${4:-0} -eq 1 || ${4:-0} -eq 3 ]] && ./splot_datasets.exe "rootfiles/"$itag_roofit".root" $itag_splot
            itag_draw=${itag_splot/\/splot_/\/draw_}
            [[ ${4:-0} -eq 2 || ${4:-0} -eq 3 ]] && ./draw_datasets.exe "rootfiles/"$itag_splot".root" $itag_draw

        done

    done

done
