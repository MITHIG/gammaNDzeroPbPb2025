#!/bin/bash

SAVE_PNG=0
LABELS_DIR=('gammaN' 'Ngamma')

binning=${1:-0}
isNgamma=${2:-0}
execu=${3:-00000}

evtcut_tag_DEFAULT='0nXn-'${LABELS_DIR[$isNgamma]}'-25'
evtzdc=${4:-"ZDCgammaN;;"} ## !!
IFS=';' evtzdcs=($evtzdc)
evtgap=${5:-"HFEMax_eta5 < 16;;"} ; [[ $isNgamma -eq 0 ]] && evtgap=${evtgap/HFEMax/HFEMaxPlus} || evtgap=${evtgap/HFEMax/HFEMaxMinus} ;
IFS=';' evtgaps=($evtgap)
cutdtopo=${6:-"((Dy<-1 && Dmva_BDT>0.143) || (Dy>=-1 && Dy<0 && Dmva_BDT>0.142) || (Dy>=0 && Dy<1 && Dmva_BDT>0.123) || (Dy>=1 && Dmva_BDT>0.098));;"}
dcut_tag_DEFAULT='Dbdt-'${LABELS_DIR[$isNgamma]}
IFS=';' cutdtopos=($cutdtopo)
fitopt=${7:-"3G-Peaky;;"}

runlevel=${8:-00000}
##

if [[ $binning -eq 0 ]] ; then
    TAG_BINNING="b-default" ; BINNING_Y='-2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.' ; BINNING_PT='2., 5.' ;
elif [[ $binning -eq 1 ]] ; then
    TAG_BINNING="b-ptdiff" ; BINNING_Y='-2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.' ; BINNING_PT='2., 3., 4., 5.' ;
else
    echo "bad binning option: "$binning" , abort."
    exit
fi

LUMINOSITY=0.060361

INPUT_DATA="/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260426-yrefmva_PbPbUPC_HIForward_Dpt-2_Dsize_24PD.root;2025 PbPb (5.36 TeV);2025PbPb"
input_data=$INPUT_DATA
CUTEVT_BASE="isL1ZDCOr && cscTightHalo2015Filter && selectedVtxFilter"
cutevtstr=$CUTEVT_BASE' && '${evtzdcs[0]}' && '${evtgaps[0]}';'${evtzdcs[1]}${evtgaps[1]}';'$evtcut_tag_DEFAULT${evtzdcs[2]}${evtgaps[2]}
INPUTS_TEMPLATE=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_GNucleusToD0-BeamA_SoftQCD_KPiKKPiPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamA"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_GNucleusToD0-BeamB_SoftQCD_KPiKKPiPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamB"
)
input_template=${INPUTS_TEMPLATE[$isNgamma]}
INPUTS_MC=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_prompt_GNucleusToD0-BeamA_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Deff.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamA"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_prompt_GNucleusToD0-BeamB_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Deff.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamB"
)
input_mc=${INPUTS_MC[$isNgamma]}
INPUTS_MC_EVT=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_prompt_GNucleusToD0-BeamA_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamA"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_prompt_GNucleusToD0-BeamB_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamB"
)
input_mc_evt=${INPUTS_MC_EVT[$isNgamma]}
CUT_BASE="TMath::Abs(Dtrk1PtErr/Dtrk1Pt)<0.1 && TMath::Abs(Dtrk2PtErr/Dtrk2Pt)<0.1 && TMath::Abs(Dtrk1Eta) < 2.4 && TMath::Abs(Dtrk2Eta) < 2.4 && Dtrk1Pt > 0.5 && Dtrk2Pt > 0.5 && Dchi2cl > 0.05 && (DsvpvDistance/DsvpvDisErr) > 1. && DsvpvDisErr>1.e-8 && DsvpvDisErr_2D>1.e-8"
cutdstr=${CUT_BASE}' && '${cutdtopos[0]}';'${cutdtopos[1]}';'$dcut_tag_DEFAULT${cutdtopos[2]}
CUT_SIGNALWIN="Dmass > 1.83 && Dmass < 1.9" # for event selection efficiency

# parse cuts
cut_tag_DEFAULT=$evtcut_tag_DEFAULT'_'$dcut_tag_DEFAULT
IFS=';' ; cutevttags=($cutevtstr) ; unset IFS ; cutevt="${cutevttags[0]}" ; cutevt_tex="${cutevttags[1]}" ; cutevt_tag="${cutevttags[2]}"
IFS=';' ; cutdtags=($cutdstr) ; unset IFS ; cutd="${cutdtags[0]}" ; cutd_tex="${cutdtags[1]}" ; cutd_tag="${cutdtags[2]}"
cut_tag=$cutevt_tag"_"$cutd_tag
cut_tex=${cutevt_tex} ; [[ x$cutevt_tex == x ]] || cut_tex=$cut_tex', ' ; cut_tex=$cut_tex$cutd_tex ;
cutstr=${cutevt}' && '${cutd}';'$cut_tex';'$cut_tag
echo -e "\033[33m"$cut_tag"\033[0m \033[2m(cut)\033[0m"
echo -e "\033[2m"$cutstr"\033[0m"

######################
# Event efficiency   #
######################
IFS=';' ; input_mc_evt_tags=($input_mc_evt) ; unset IFS ; mc_evt_tag=${input_mc_evt_tags[2]} ;
itag_evteff_calc=$cut_tag_DEFAULT"/"$TAG_BINNING"/evteffcalc_"$mc_evt_tag # default tag
run_level=${runlevel:0:1}
[[ ${execu:0:1} -gt 0 ]] && {
    itag_evteff=$cut_tag"/"$TAG_BINNING"/evteffsave_"$mc_evt_tag
    echo -e "\033[33;2m"$cut_tag" / \033[0m\033[33m"$mc_evt_tag"\033[0m \033[2m(event efficiency)\033[0m"
    [[ $run_level -eq 1 || $run_level -eq 3 ]] && ./evteff_save.exe "$input_mc_evt" "$cutevtstr" "${cutdstr} && $CUT_SIGNALWIN" $itag_evteff "$BINNING_Y" "$BINNING_PT" #
    itag_evteff_calc=${itag_evteff/evteffsave/evteffcalc}
    [[ $run_level -eq 2 || $run_level -eq 3 ]] && ./evteff_calc.exe "rootfiles/"$itag_evteff".root" $itag_evteff_calc $SAVE_PNG #
}

####################
# Mass template    #
####################
IFS=';' ; input_template_tags=($input_template) ; unset IFS ; template_tag=${input_template_tags[2]} ;
itag_template=$cut_tag_DEFAULT"/"$TAG_BINNING"/template_"$template_tag #
run_level=${runlevel:1:1}
[[ ${execu:1:1} -eq 1 ]] && {
    itag_template=$cut_tag"/"$TAG_BINNING"/template_"$template_tag #
    [[ $run_level -eq 1 || $run_level -eq 3 ]] && ./hist_save.exe "$input_template" "$cutstr" $itag_template "$BINNING_Y" "$BINNING_PT" 1 # 1: is_template
    echo -e "\033[33;2m"$cut_tag" / \033[0m\033[33m"$template_tag"\033[0m \033[2m(mass template)\033[0m"
}

####################
# Fill data mass   #
####################
IFS=';' ; input_data_tags=($input_data) ; unset IFS ; data_tag=${input_data_tags[2]}
itag_data=$cut_tag_DEFAULT"/"$TAG_BINNING"/savehist_"$data_tag #
run_level=${runlevel:2:1}
[[ ${execu:2:1} -eq 1 ]] && {
    itag_data=$cut_tag"/"$TAG_BINNING"/savehist_"$data_tag #
    [[ $run_level -eq 1 || $run_level -eq 3 ]] && ./hist_save.exe "$input_data" "$cutstr" $itag_data "$BINNING_Y" "$BINNING_PT" 0 # 0: not template
    echo -e "\033[33;2m"$cut_tag" / \033[0m\033[33m"$data_tag"\033[0m \033[2m(data)\033[0m"
}

####################
# Mass fitting     #
####################
IFS=';' ; fitopts=($fitopt) ; unset IFS ; fit_tag=${fitopts[2]} ;
itag_data_fit=$cut_tag_DEFAULT"/"$TAG_BINNING"/fithist_"$data_tag"_"$template_tag''
run_level=${runlevel:3:1}
[[ ${execu:3:1} -eq 1 ]] && {
    itag_data_fit=$cut_tag"/"$TAG_BINNING"/fithist_"$data_tag"_"$template_tag$fit_tag
    [[ $run_level -eq 2 || $run_level -eq 3 ]] && ./hist_fit.exe "rootfiles/"$itag_data".root" "rootfiles/"$itag_template".root" $itag_data_fit "$fitopt" $SAVE_PNG
    echo -e "\033[33;2m"$cut_tag" / "$template_tag" / \033[0m\033[33m"$data_tag"\033[0m \033[2m(fit)\033[0m"
}

####################
# D efficiency     #
####################
IFS=';' ; input_mc_tags=($input_mc) ; unset IFS ; mc_tag=${input_mc_tags[2]} ;
itag_deff_calc=$cut_tag_DEFAULT"/"$TAG_BINNING"/deffcalc_"$mc_tag
run_level=${runlevel:4:1}
[[ ${execu:4:1} -eq 1 ]] && {
    itag_deff=$cut_tag"/deffsave_"$mc_tag # no binning info
    [[ ($run_level -eq 1 || $run_level -eq 3) && $TAG_BINNING == "b-default" ]] && ./eff_save.exe "$input_mc" "$cutevtstr" "$cutdstr" $itag_deff null
    itag_deff_calc=$cut_tag"/"$TAG_BINNING"/deffcalc_"$mc_tag
    [[ $run_level -eq 2 || $run_level -eq 3 ]] && ./eff_calc.exe "rootfiles/"$itag_deff".root" $itag_deff_calc "$BINNING_Y" "$BINNING_PT" $SAVE_PNG
    echo -e "\033[33;2m"$cut_tag" / \033[0m\033[33m"$mc_tag"\033[0m \033[2m(D efficiency)\033[0m"
}

itag_fprompt='null'

####################
# Cross-section    #
####################
itag_xsec=$cut_tag'/'$TAG_BINNING'/xsec_'${itag_data_fit##*/}'_'${itag_deff_calc##*/}'_'${itag_evteff_calc##*/}'_'$itag_fprompt
echo "  itag_data_fit:     "$itag_data_fit
echo "  itag_deff_calc:    "$itag_deff_calc
echo "  itag_evteff_calc:  "$itag_evteff_calc
echo "  itag_fprompt:      "$itag_fprompt
echo "  lumi:              "$LUMINOSITY" nb-1"
echo "              ==> "$itag_xsec

run_level=$runlevel
[[ $run_level -gt 0 ]] && {
    ./xsec_calc.exe "rootfiles/"$itag_data_fit".root" "rootfiles/"$itag_deff_calc".root" "rootfiles/"$itag_evteff_calc".root" $itag_fprompt $LUMINOSITY $itag_xsec $SAVE_PNG
    [[ $isNgamma -eq 1 ]] && {
        file_Ngamma='rootfiles/'$itag_xsec'.root'
        file_gammaN=${file_Ngamma//Ngamma/gammaN} ; file_gammaN=${file_gammaN//BeamB/BeamA} ;
        outputname=${itag_xsec//-Ngamma/} ; outputname=${outputname//-BeamB/} ;
        ls $file_gammaN
        ls $file_Ngamma
        echo $outputname
        [[ -f $file_gammaN && -f $file_Ngamma ]] && {
            set -x
            ./xsec_collect.exe "${file_gammaN},${file_Ngamma}" $outputname 0
            set +x
        }
    }
}
