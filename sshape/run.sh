#!/bin/bash

SAVE_PNG=1

INPUTS_DATA=(
    '/eos/cms/store/group/phys_heavyions/wangj/Forest2025OO/crab_skim_261118_IonPhysics1-20_trkpt0p7_Dpt2.root;OO (5.36 TeV);OO;D0Skim/Tree'
)

INPUTS_TEMPLATE=(
    '/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_GNucleusToD0-BeamA_SoftQCD_KPiKKPiPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);GNucleusToD0-BeamA;Tree'
)

fitopt='3G-Peaky-Fwidth;;'

CUTS=(
    'TMath::Abs(Dtrk1Eta) < 1.5 && TMath::Abs(Dtrk2Eta) < 1.5 && Dtrk1Pt > 0.7 && Dtrk2Pt > 0.7 && (Dtrk1PixelHit+Dtrk1StripHit)>=11 && (Dtrk2PixelHit+Dtrk2StripHit)>=11 && TMath::Abs(Dtrk1PtErr/Dtrk1Pt)<0.1 && TMath::Abs(Dtrk2PtErr/Dtrk2Pt)<0.1 && Dpt > 2.5 && DsvpvDistance/DsvpvDisErr > 4 && Dalpha < 0.15 && Dchi2cl > 0.15;;Deta15;0'
    'TMath::Abs(Dtrk1Eta) < 1.5 && TMath::Abs(Dtrk2Eta) < 1.5 && Dtrk1Pt > 0.7 && Dtrk2Pt > 0.7 && Dpt > 2.5 && Ddls > 4;;Deta15;1'
)

make hist_save.exe hist_fit.exe || exit 1

for cutstr in "${CUTS[@]}" ; do
    IFS=';' ; cuttags=($cutstr) ; unset IFS ; cut_tag="${cuttags[2]}" ; cut_flag=${cuttags[3]} ;

    [[ x$cut_tag == x ]] && { echo "warning: missed cut_tag. skip." ; continue ; }

    for input_template in "${INPUTS_TEMPLATE[@]}" ; do

        IFS=';' ; input_template_tags=($input_template) ; unset IFS ; template_tag=${input_template_tags[2]} ; template_treename=${input_template_tags[3]} ;

        itag_template=$cut_tag"/template_"$template_tag #
        if [[ ${1:-0} -eq 1 ]] ; then
            [[ ($template_treename == *D0Skim* && $cut_flag -gt 0) || (($template_treename != *D0Skim* && $cut_flag -eq 0)) ]] && {
                ./hist_save.exe "$input_template" $template_treename "$cutstr" $itag_template 1 # 1: is template
            }
        fi

        for input_data in "${INPUTS_DATA[@]}" ; do

            IFS=';' ; input_data_tags=($input_data) ; unset IFS ; data_tag=${input_data_tags[2]} ; data_treename=${input_data_tags[3]} ;
            [[ ($data_treename == *D0Skim* && $cut_flag -gt 0) || (($data_treename != *D0Skim* && $cut_flag -eq 0)) ]] || continue ;
            echo -e "\033[33;2m"$cut_tag" / "$template_tag" / \033[0m\033[33m"$data_tag"\033[0m \033[2m(data)\033[0m"

            itag_data=$cut_tag'/savehist_'$data_tag #
            [[ ${2:-0} -eq 1 ]] && {
                ./hist_save.exe "$input_data" $data_treename "$cutstr" $itag_data 0 # 0: not template
            }

            IFS=';' ; fitopts=($fitopt) ; unset IFS ; fit_tag=${fitopts[2]} ;
            itag_data_fit=$cut_tag'/fithist_'$data_tag'_'$template_tag$fit_tag ##
            [[ ${3:-0} -eq 1 ]] && {
                ./hist_fit.exe 'rootfiles/'$itag_data'.root' 'rootfiles/'$itag_template'.root' $itag_data_fit "$fitopt" $SAVE_PNG
            }
        done
    done 
done
