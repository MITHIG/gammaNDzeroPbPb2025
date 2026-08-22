#!/bin/bash

SAVE_PNG=0

TAG_BINNING="_b-default" ; BINNING_Y='-2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.' ; BINNING_PT='2., 5.' ;
# TAG_BINNING="_b-ptdiff" ; BINNING_Y='-2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.' ; BINNING_PT='2., 3., 4., 5.' ;
# TAG_BINNING="_b-ycoarse" ; BINNING_Y='-2., -1., 0., 1., 2.' ; BINNING_PT='2., 5.' ;
# TAG_BINNING="_b-incl" ; BINNING_Y='-2., 2.' ; BINNING_PT='2., 5.' ;

fitopt="3P;Triple gaus signal;"
# fitopt="P;Double gaus signal;_f-2gaus"
# fitopt="3;No KK/#pi#pi;_f-nopeaky"

INPUTS_DATA=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260426-yrefmva_PbPbUPC_HIForward_Dpt-2_Dsize_24PD.root;2025 PbPb (5.36 TeV);2025PbPb" # 2025
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260814-yinclbdt_PbPbUPC_HIForward_Trig-3_Dsize.root;2025 PbPb ZB (5.36 TeV);2025PbPb" # 2025 ZB + inclusive BDT
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_2023PbPbUPC_Jan2024ReReco_20260212Forest_HIForward_Dpt-2_Trig-2_Dsize_xbr.root;2023 PbPb (Jan2024 ReReco);2023PbPb-recoJan2024"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260426-yrefmva_2023PbPbUPC_Feb2025ReReco_20260521Forest_HIForward_Dpt-2_Trig-2_Dsize.root;2023 PbPb (Feb2025 ReReco);2023PbPb-recoFeb2025"
)
CUTEVTS=( # lumi is nb-1 - directly from brilcalc
    # # 2025 0nXn
    "isL1ZDCOr && cscTightHalo2015Filter && selectedVtxFilter && ZDCgammaN && HFEMaxPlus_eta5 < 16;Xn0n (#gammaN);gammaN-0nXn-25;0.060361"
    "isL1ZDCOr && cscTightHalo2015Filter && selectedVtxFilter && ZDCNgamma && HFEMaxMinus_eta5 < 16;0nXn (N#gamma);Ngamma-0nXn-25;0.060361"
    # "isZeroBias && cscTightHalo2015Filter && selectedVtxFilter && ZDCgammaN && HFEMaxPlus_eta5 < 16;Xn0n (#gammaN) ZB;gammaN-0nXn-25-ZB;0.0082426"
    # "isZeroBias && cscTightHalo2015Filter && selectedVtxFilter && ZDCNgamma && HFEMaxMinus_eta5 < 16;0nXn (N#gamma) ZB;Ngamma-0nXn-25-ZB;0.0082426"
    # # 2025 0nAn -> need to change BDT
    # "isZeroBias && cscTightHalo2015Filter && selectedVtxFilter && ZDCsumPlus < 1100 && HFEMaxPlus_eta5 < 16;An0n (#gammaN) ZeroBias;gammaN-0nAn-25-ZB;0.0082426"
    # "isZeroBias && cscTightHalo2015Filter && selectedVtxFilter && ZDCsumMinus < 1000 && HFEMaxMinus_eta5 < 16;0nAn (N#gamma) ZeroBias;Ngamma-0nAn-25-ZB;0.0082426"
    # # 2025 0n0n
    # "isZeroBias && cscTightHalo2015Filter && selectedVtxFilter && ZDCsumPlus < 1100 && ZDCsumMinus < 1000 && HFEMaxPlus_eta5 < 16 && HFEMaxMinus_eta5 < 16;0n0n (Both gap);gammaN-0n0n-2gap-25;0.0082426"
    # # 2023
    # "isL1ZDCOr && cscTightHalo2015Filter && selectedVtxFilter && ZDCgammaN && HFEMaxPlus_eta5 < 9.2 && ClusterCompatibilityFilter && nVtx <= 3;#gammaN (23);gammaN-0nXn-23;0.007803"
)
INPUTS_TEMPLATE=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_prompt_GNucleus-QCD-PhotonBeamA_Bin-Pthat0_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2024-SoftQCD-BeamA"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_prompt_GNucleus-QCD-PhotonBeamB_Bin-Pthat0_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 N#gamma (5.36 TeV);2024-SoftQCD-BeamB"
)
INPUTS_MC=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_prompt_GNucleusToD0-PhotonBeamA_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2024-SoftQCD-BeamA"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_prompt_GNucleusToD0-PhotonBeamB_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2.root;P#scale[0.8]{YTHIA}8 N#gamma (5.36 TeV);2024-SoftQCD-BeamB"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260814-yinclbdt_HiForest_260328_prompt_GNucleusToD0-PhotonBeamA_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2024-SoftQCD-BeamA"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260814-yinclbdt_HiForest_260328_prompt_GNucleusToD0-PhotonBeamB_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2.root;P#scale[0.8]{YTHIA}8 N#gamma (5.36 TeV);2024-SoftQCD-BeamB"
)
# CUT_BASE="Dpt>=2 && Dpt<5 && TMath::Abs(Dtrk1PtErr/Dtrk1Pt)<0.1 && TMath::Abs(Dtrk2PtErr/Dtrk2Pt)<0.1 && TMath::Abs(Dtrk1Eta) < 2.4 && TMath::Abs(Dtrk2Eta) < 2.4 && Dtrk1Pt > 0.5 && Dtrk2Pt > 0.5 && Dchi2cl > 0.05 && (DsvpvDistance/DsvpvDisErr) > 1. && DsvpvDisErr>1.e-8 && DsvpvDisErr_2D>1.e-8"
CUT_BASE="TMath::Abs(Dtrk1PtErr/Dtrk1Pt)<0.1 && TMath::Abs(Dtrk2PtErr/Dtrk2Pt)<0.1 && TMath::Abs(Dtrk1Eta) < 2.4 && TMath::Abs(Dtrk2Eta) < 2.4 && Dtrk1Pt > 0.5 && Dtrk2Pt > 0.5 && Dchi2cl > 0.05 && (DsvpvDistance/DsvpvDisErr) > 1. && DsvpvDisErr>1.e-8 && DsvpvDisErr_2D>1.e-8"
CUTDS=(
    "${CUT_BASE} && ((Dy<-1 && Dmva_BDT>0.143) || (Dy>=-1 && Dy<0 && Dmva_BDT>0.142) || (Dy>=0 && Dy<1 && Dmva_BDT>0.123) || (Dy>=1 && Dmva_BDT>0.098));BDT;Dbdt-gammaN"
    "${CUT_BASE} && ((Dy>=1 && Dmva_BDT>0.143) || (Dy<1 && Dy>=0 && Dmva_BDT>0.142) || (Dy<0 && Dy>=-1 && Dmva_BDT>0.123) || (Dy<-1 && Dmva_BDT>0.098));BDT;Dbdt-Ngamma"
    # "${CUT_BASE} && ((Dy<-1 && Dmva_BDT>0.160) || (Dy>=-1 && Dy<0 && Dmva_BDT>0.142) || (Dy>=0 && Dy<1 && Dmva_BDT>0.123) || (Dy>=1 && Dmva_BDT>0.098));BDT;Dbdt-gammaN-ZB"
    # "${CUT_BASE} && ((Dy>=1 && Dmva_BDT>0.160) || (Dy<1 && Dy>=0 && Dmva_BDT>0.142) || (Dy<0 && Dy>=-1 && Dmva_BDT>0.123) || (Dy<-1 && Dmva_BDT>0.098));BDT;Dbdt-Ngamma-ZB"
    # "${CUT_BASE} && Dmva_BDT>0.12;BDT;Dbdt-0n0n-ZB"
    # "${CUT_BASE} && Dmva_BDT>0;D pre-cuts;Dpre"
    # "${CUT_BASE} && (Dtrk1PixelHit+Dtrk1StripHit)>=11 && (Dtrk2PixelHit+Dtrk2StripHit)>=11 && DpassCut23PAS;Cut23PAS;D23pas"
    # "${CUT_BASE} && DpassCut23PAS;Cut23PAS;D23pas-nonhit"
)

##
echo "usage: ./run_main.sh [template] [save hist] [fit hist] [save eff] [calc eff] [calc xsec]"
echo "                          1          2          3           4          5           6 "

make hist_save.exe hist_fit.exe eff_save.exe eff_calc.exe xsec_calc.exe || exit 1

for cutevtstr in "${CUTEVTS[@]}" ; do
    IFS=';' ; cutevttags=($cutevtstr) ; unset IFS ; cutevt="${cutevttags[0]}" ; cutevt_tex="${cutevttags[1]}" ; cutevt_tag="${cutevttags[2]}" ; cutevt_lumi=${cutevttags[3]} ;

    [[ x$cutevt_tag == x ]] && { echo "warning: missed cutevt_tag. skip." ; continue ; }

    ####################
    # Event efficiency #
    ####################
    itag_event=$cutevt_tag"/event_eff" #

    ## Loop cuts
    for cutdstr in "${CUTDS[@]}" ; do
        IFS=';' ; cutdtags=($cutdstr) ; unset IFS ; cutd="${cutdtags[0]}" ; cutd_tex="${cutdtags[1]}" ; cutd_tag="${cutdtags[2]}"

        [[ x$cutd_tag == x ]] && { echo "warning: missed cutd_tag. skip." ; continue ; }
        [[ ($cutd_tag == *gammaN* && $cutevt_tag == *Ngamma*) || ($cutd_tag == *Ngamma* && $cutevt_tag == *gammaN*) ]] && continue 

        cut_tag=$cutevt_tag"_"$cutd_tag
        cutstr=${cutevt}" && "${cutd}" ; "${cutevt_tex}", "${cutd_tex}" ; "$cut_tag

        echo -e "\033[33m"$cut_tag"\033[0m"

        # loop template
        for input_template in "${INPUTS_TEMPLATE[@]}" ; do
            
            IFS=';' ; input_template_tags=($input_template) ; unset IFS ; template_tag=${input_template_tags[2]} ;
            [[ ($cutevt_tag == *gammaN* && $template_tag == *BeamB*) || ($cutevt_tag == *Ngamma* && $template_tag == *BeamA*) || ($cutevt_tag == *0n0n* && $template_tag == *BeamB*) ]] && continue 
            echo -e "\033[33;2m"$cut_tag" / \033[0m\033[33m"$template_tag"\033[0m"

            ####################
            # Mass template    # -> [ 3 min ]
            ####################
            itag_template=$cut_tag"/template_"$template_tag$TAG_BINNING #
            if [[ ${1:-0} -eq 1 ]] ; then
                echo "    -> generate mass templates from MC (about 3 min)"
                ./hist_save.exe "$input_template" "$cutstr" $itag_template "$BINNING_Y" "$BINNING_PT" 0 # 0: !isdata
            fi

            ## Loop data
            for input_data in "${INPUTS_DATA[@]}" ; do
                
                IFS=';' ; input_data_tags=($input_data) ; unset IFS ; data_tag=${input_data_tags[2]}
                echo -e "\033[33;2m"$cut_tag" / "$template_tag" / \033[0m\033[33m"$data_tag"\033[0m"

                ####################
                # Fill data mass   # -> [ 13 min ]
                ####################
                itag_data=$cut_tag"/savehist_"$data_tag$TAG_BINNING
                [[ ${2:-0} -eq 1 ]] && {
                    echo "    -> fill data mass (about 13 min)"
                    ./hist_save.exe "$input_data" "$cutstr" $itag_data "$BINNING_Y" "$BINNING_PT" 1 # 1: isdata
                }

                ####################
                # Mass fitting     #
                ####################
                IFS=';' ; fitopts=($fitopt) ; unset IFS ; fit_tag=${fitopts[2]} ; 
                itag_data_fit=$cut_tag"/fithist_"$data_tag"_"$template_tag$TAG_BINNING$fit_tag ## 
                [[ ${3:-0} -eq 1 ]] && {
                    echo "    -> fit invariant mass"
                    ./hist_fit.exe "rootfiles/"$itag_data".root" "rootfiles/"$itag_template".root" $itag_data_fit "$fitopt" $SAVE_PNG
                }

                ## Loop MC
                for input_mc in "${INPUTS_MC[@]}" ; do

                    IFS=';' ; input_mc_tags=($input_mc) ; unset IFS ; mc_tag=${input_mc_tags[2]} ;
                    [[ ($cutevt_tag == *gammaN* && $mc_tag == *BeamB*) || ($cutevt_tag == *Ngamma* && $mc_tag == *BeamA*) || ($cutevt_tag == *0n0n* && $mc_tag == *BeamB*) ]] && continue 
                    echo -e "\033[33;2m"$cut_tag" / "$template_tag" / "$data_tag" / \033[0m\033[33m"$mc_tag"\033[0m"

                    ####################
                    # D efficiency     # -> [ 36 min ]
                    ####################
                    itag_deff=$cut_tag"/effsave_"$mc_tag # no binning info
                    if [[ ${4:-0} -eq 1 ]] ; then
                        echo "    -> generate D efficiency table from MC (about 36 min)"
                        # ./eff_save.exe "$input_mc" "$cutevtstr" "$cutdstr" $itag_deff "$input_data"
                        ./eff_save.exe "$input_mc" "$cutevtstr" "$cutdstr" $itag_deff null
                    fi

                    itag_deff_calc=${itag_deff/effsave/effcalc}$TAG_BINNING 
                    [[ ${5:-0} -eq 1 ]] && {
                        ./eff_calc.exe "rootfiles/"$itag_deff".root" $itag_deff_calc "$BINNING_Y" "$BINNING_PT" $SAVE_PNG
                    }
                    
                    ####################
                    # Cross-section    #
                    ####################
                    # itag_xsec=$cut_tag"/xsec_"${itag_data_fit##*/}"_"${itag_deff##*/}"_null_null"$TAG_BINNING
                    echo "    itag_data_fit:  "$itag_data_fit
                    echo "    itag_deff:      "$itag_deff_calc
                    echo "    itag_evteff:    null"
                    echo "    itag_fprompt:   null"
                    echo "    lumi:           "$cutevt_lumi" nb-1"
                    # echo "                ==> "$itag_xsec
                    [[ ${6:-0} -eq 1 ]] && {
                        echo "    -> calculate cross sections"
                        ./xsec_calc.exe "rootfiles/"$itag_data_fit".root" "rootfiles/"$itag_deff_calc".root" null null $cutevt_lumi $cut_tag 
                    }
                done
            done
        done
    done 
done

# TAG_BINNING="_b-yextend" ; BINNING_Y='-2.4, -2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2., 2.4' ; BINNING_PT='2., 5.' ;
# TAG_BINNING="_b-ptlow" ; BINNING_Y='-2., 2.' ; BINNING_PT='1., 2., 3., 4., 5.' ;
