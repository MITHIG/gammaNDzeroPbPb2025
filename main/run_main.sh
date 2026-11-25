#!/bin/bash

SAVE_PNG=1
ver='v2'

TAG_BINNING="b-default" ; BINNING_Y='-2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.' ; BINNING_PT='2., 5.' ; INPUT_FPROMPT='../fprompt/rootfiles/0nXn-twodirs-25_Dbdt/b-ycoarse/collect__2025PbPb_BeamA-prompt_BeamA-nonprompt_Dip3D-Dip3Dsig/'$TAG_BINNING'.root;;fprompt-25MC' ;
# TAG_BINNING="b-ptdiff" ; BINNING_Y='-2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.' ; BINNING_PT='2., 3., 4., 5.' ; INPUT_FPROMPT='null' ;
# TAG_BINNING="b-ycoarse" ; BINNING_Y='-2., -1., 0., 1., 2.' ; BINNING_PT='2., 5.' ; INPUT_FPROMPT='null' ;

templateopt='1.07;#alpha = 1.07;'
# templateopt='1;Fixed to MC;-tmc'
fitopt='3G-Peaky-DSig;;'
# fitopt='Peaky-DSig-Exp-Fwidth;;_f-2023'
# fitopt='Peaky;;_f-2gaus'
femdpu='0.95'

INPUTS_DATA=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260928-symmva_PbPbUPC_HIForward_Dpt-2_Dsize_Golden_24PD.root;2025 PbPb (5.36 TeV);2025PbPb;0.057607"

    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260928-symmva_PbPbUPC_HIForward_Dpt-2_Dsize_24PD.root;2025 PbPb (5.36 TeV);2025PbPb-mixjson;0.060361"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260426-yrefmva_PbPbUPC_HIForward_Dpt-2_Dsize_24PD.root;2025 PbPb (5.36 TeV);2025PbPb;0.060361" # 2025
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_261002-ptdiff_PbPbUPC_HIForward_Dpt-1_xbr_Dsize_gap-17_12ePD.root;2025 PbPb (5.36 TeV);2025PbPb" # 2025 low pT
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260814-yinclbdt_PbPbUPC_HIForward_Trig-3_Dsize.root;2025 PbPb ZB (5.36 TeV);2025PbPb" # 2025 ZB + inclusive BDT
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_2023PbPbUPC_Jan2024ReReco_20260212Forest_HIForward_Dpt-2_Trig-2_Dsize_xbr.root;2023 PbPb (Jan2024 ReReco);2023PbPb-recoJan2024" # 2023 Jan24 rereco
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260426-yrefmva_2023PbPbUPC_Feb2025ReReco_20260521Forest_HIForward_Dpt-2_Trig-2_Dsize.root;2023 PbPb (Feb2025 ReReco);2023PbPb-recoFeb2025" # 2023 Feb25 rereco
)
CUTEVTS=( # lumi is nb-1 - directly from brilcalc
    # # --- 2025 0nXn ZDCOr
    "isL1ZDCOr && cscTightHalo2015Filter && selectedVtxFilter && ZDCgammaN && HFEMaxPlus_forest < 19.4;Xn0n (#gammaN);0nXn-gammaN-25;"
    "isL1ZDCOr && cscTightHalo2015Filter && selectedVtxFilter && ZDCNgamma && HFEMaxMinus_forest < 17.2;0nXn (N#gamma);0nXn-Ngamma-25;"

    # "isL1ZDCOr && cscTightHalo2015Filter && selectedVtxFilter && ZDCgammaN && HFEMaxPlus_eta5 < 17.2;Xn0n (#gammaN);0nXn-gammaN-25-gapeta5;"
    # "isL1ZDCOr && cscTightHalo2015Filter && selectedVtxFilter && ZDCNgamma && HFEMaxMinus_eta5 < 15.2;0nXn (N#gamma);0nXn-Ngamma-25-gapeta5;"
    # "isL1ZDCOr && cscTightHalo2015Filter && selectedVtxFilter && ZDCgammaN && HFEMaxPlus_eta5 < 16;Xn0n (#gammaN);0nXn-gammaN-25-preapp;"
    # "isL1ZDCOr && cscTightHalo2015Filter && selectedVtxFilter && ZDCNgamma && HFEMaxMinus_eta5 < 16;0nXn (N#gamma);0nXn-Ngamma-25-preapp;"
    # # --- 2025 0nXn ZB
    # "isZeroBias && cscTightHalo2015Filter && selectedVtxFilter && ZDCgammaN && HFEMaxPlus_forest < 19.4;Xn0n (#gammaN) ZB;0nXn-gammaN-25-ZB;0.00804756"
    # "isZeroBias && cscTightHalo2015Filter && selectedVtxFilter && ZDCNgamma && HFEMaxMinus_forest < 17.2;0nXn (N#gamma) ZB;0nXn-Ngamma-25-ZB;0.00804756"
    # # --- 2025 0nAn -> need to change inputs
    # "isZeroBias && cscTightHalo2015Filter && selectedVtxFilter && ZDCsumPlus < 1100 && HFEMaxPlus_eta5 < 16;An0n (#gammaN) ZeroBias;0nAn-gammaN-25-ZB;0.0082426"
    # "isZeroBias && cscTightHalo2015Filter && selectedVtxFilter && ZDCsumMinus < 1000 && HFEMaxMinus_eta5 < 16;0nAn (N#gamma) ZeroBias;0nAn-Ngamma-25-ZB;0.0082426"
    # # --- 2025 0n0n -> need to change inputs
    # "isZeroBias && cscTightHalo2015Filter && selectedVtxFilter && ZDCsumPlus < 1100 && ZDCsumMinus < 1000 && HFEMaxPlus_eta5 < 16 && HFEMaxMinus_eta5 < 16;0n0n (Both gap);0n0n-2gap-25;0.0082426"
    # # --- 2023
    # "isL1ZDCOr && cscTightHalo2015Filter && selectedVtxFilter && ZDCgammaN && HFEMaxPlus_eta5 < 9.2 && ClusterCompatibilityFilter && nVtx <= 3;#gammaN (23);0nXn-gammaN-23;0.007803"
)
INPUTS_TEMPLATE=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_GNucleusToD0-BeamA_SoftQCD_KPiKKPiPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamA"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_GNucleusToD0-BeamB_SoftQCD_KPiKKPiPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 N#gamma (5.36 TeV);2025-SoftQCD-BeamB"

    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_prompt_GNucleus-QCD-PhotonBeamA_Bin-Pthat0_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2024-SoftQCD-BeamA"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_prompt_GNucleus-QCD-PhotonBeamB_Bin-Pthat0_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 N#gamma (5.36 TeV);2024-SoftQCD-BeamB"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_261002-ptdiff_HiForest_260904_GNucleusToD0-BeamA_SoftQCD_KPiKKPiPi_2025_trkpt0p1_Drej-genmatched_Dpt-1_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamA"
)
INPUTS_MC=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_GNucleusToD0-BeamA_SoftQCD_KPiKKPiPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Deff.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamA;0"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_GNucleusToD0-BeamB_SoftQCD_KPiKKPiPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Deff.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamB;0"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_prompt_GNucleusToD0-BeamA_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Deff.root;Prompt P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamA-prompt;1"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_prompt_GNucleusToD0-BeamB_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Deff.root;Prompt P#scale[0.8]{YTHIA}8 N#gamma (5.36 TeV);2025-SoftQCD-BeamB-prompt;1"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260928-symmva_HiForest_260904_nonprompt_GNucleusToD0-BeamA_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Deff.root;Nonprompt P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamA-nonprompt;2"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260928-symmva_HiForest_260904_nonprompt_GNucleusToD0-BeamB_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Deff.root;Nonprompt P#scale[0.8]{YTHIA}8 N#gamma (5.36 TeV);2025-SoftQCD-BeamB-nonprompt;2"
    
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260928-symmva_HiForest_260904_prompt_GNucleusToD0-BeamA_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-1_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamA"
)
INPUTS_MC_EVT=( # only one sample work for each cut
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_GNucleusToD0-BeamA_SoftQCD_KPiKKPiPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamA"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_GNucleusToD0-BeamB_SoftQCD_KPiKKPiPi_2025_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamB"

    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_261002-ptdiff_HiForest_260904_GNucleusToD0-BeamA_SoftQCD_KPiKKPiPi_2025_trkpt0p1_Drej-genmatched_Dpt-1_Dsize.root;P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV);2025-SoftQCD-BeamA"
)
CUT_BASE="TMath::Abs(Dtrk1PtErr/Dtrk1Pt)<0.1 && TMath::Abs(Dtrk2PtErr/Dtrk2Pt)<0.1 && TMath::Abs(Dtrk1Eta) < 2.4 && TMath::Abs(Dtrk2Eta) < 2.4 && Dtrk1Pt > 0.5 && Dtrk2Pt > 0.5 && Dchi2cl > 0.05 && (DsvpvDistance/DsvpvDisErr) > 1. && DsvpvDisErr>1.e-8 && DsvpvDisErr_2D>1.e-8"
CUT_SIGNALWIN="Dmass > 1.83 && Dmass < 1.9" # for event selection efficiency
CUTDS=(
    "${CUT_BASE} && ((Dy<-1 && Dmva_BDT>0.143) || (Dy>=-1 && Dy<0 && Dmva_BDT>0.142) || (Dy>=0 && Dy<1 && Dmva_BDT>0.123) || (Dy>=1 && Dmva_BDT>0.098));BDT;Dbdt-gammaN;"
    "${CUT_BASE} && ((Dy>=1 && Dmva_BDT>0.143) || (Dy<1 && Dy>=0 && Dmva_BDT>0.142) || (Dy<0 && Dy>=-1 && Dmva_BDT>0.123) || (Dy<-1 && Dmva_BDT>0.098));BDT;Dbdt-Ngamma;"

    # "${CUT_BASE} && ((Dy<-1 && Dmva_BDT>0.143) || (Dy>=-1 && Dy<0 && Dmva_BDT>0.142) || (Dy>=0 && Dy<1 && Dmva_BDT>0.123) || (Dy>=1 && Dmva_BDT>0.098));BDT;Dbdt-gammaN-weight;(Dy<-1)*(1.27+(-0.841)*0.5*TMath::Erfc((Dmva_BDT-0.0831)/0.03))+(Dy>=-1 && Dy<0)*(0.793+1.02*0.5*TMath::Erfc((Dmva_BDT-(-0.00214))/0.0772))+(Dy>=0 && Dy<1)*(0.834+0.714*0.5*TMath::Erfc((Dmva_BDT-0.0218)/0.03))+(Dy>=1)*(1.02*(1-((Dmva_BDT-0.08)/0.556)^2))"
    # "${CUT_BASE} && ((Dy>=1 && Dmva_BDT>0.143) || (Dy<1 && Dy>=0 && Dmva_BDT>0.142) || (Dy<0 && Dy>=-1 && Dmva_BDT>0.123) || (Dy<-1 && Dmva_BDT>0.098));BDT;Dbdt-Ngamma-weight;(Dy>=1)*(1.27+(-0.841)*0.5*TMath::Erfc((Dmva_BDT-0.0831)/0.03))+(Dy<1 && Dy>=0)*(0.793+1.02*0.5*TMath::Erfc((Dmva_BDT-(-0.00214))/0.0772))+(Dy<0 && Dy>=-1)*(0.834+0.714*0.5*TMath::Erfc((Dmva_BDT-0.0218)/0.03))+(Dy<-1)*(1.02*(1-((Dmva_BDT-0.08)/0.556)^2))"
    # "${CUT_BASE} && (Dtrk1PixelHit+Dtrk1StripHit)>=11 && (Dtrk2PixelHit+Dtrk2StripHit)>=11 && DpassCut23PAS;2023 PAS cut;D23pas;"
    # "${CUT_BASE};Prefilter;prefilter;"
    # # for different mva, don't use if you don't know:
    # "${CUT_BASE} && ((Dy<-1 && Dmva_BDT>0.160) || (Dy>=-1 && Dy<0 && Dmva_BDT>0.142) || (Dy>=0 && Dy<1 && Dmva_BDT>0.123) || (Dy>=1 && Dmva_BDT>0.098));BDT;Dbdt-gammaN-ZB;"
    # "${CUT_BASE} && ((Dy>=1 && Dmva_BDT>0.160) || (Dy<1 && Dy>=0 && Dmva_BDT>0.142) || (Dy<0 && Dy>=-1 && Dmva_BDT>0.123) || (Dy<-1 && Dmva_BDT>0.098));BDT;Dbdt-Ngamma-ZB;"
    # "${CUT_BASE} && Dmva_BDT>0.12;BDT;Dbdt-0n0n-ZB;"
)

## Do not change below

##
echo "usage: ./run_main.sh [evt eff (2)] [template (1)] [hist (2)] [D eff (2)] [calc xsec]"
echo "                         1             2              3           4          5"

make evteff_save.exe evteff_calc.exe hist_save.exe hist_fit.exe eff_save.exe eff_calc.exe xsec_calc.exe xsec_collect.exe || exit 1

# loop event selections
for cutevtstr in "${CUTEVTS[@]}" ; do
    IFS=';' ; cutevttags=($cutevtstr) ; unset IFS ; cutevt="${cutevttags[0]}" ; cutevt_tex="${cutevttags[1]}" ; cutevt_tag="${cutevttags[2]}" ; cutevt_lumi=${cutevttags[3]} ;

    [[ x$cutevt_tag == x ]] && { echo "warning: missed cutevt_tag. skip." ; continue ; }

    # Loop D cuts
    for cutdstr in "${CUTDS[@]}" ; do
        IFS=';' ; cutdtags=($cutdstr) ; unset IFS ; cutd="${cutdtags[0]}" ; cutd_tex="${cutdtags[1]}" ; cutd_tag="${cutdtags[2]}" ; cutd_weight="${cutdtags[3]}"

        [[ x$cutd_tag == x ]] && { echo "warning: missed cutd_tag. skip." ; continue ; }
        [[ ($cutd_tag == *gammaN* && $cutevt_tag == *Ngamma*) || ($cutd_tag == *Ngamma* && $cutevt_tag == *gammaN*) ]] && continue 

        cut_tag=$cutevt_tag"_"$cutd_tag
        cutstr=${cutevt}" && "${cutd}" ; "${cutevt_tex}", "${cutd_tex}" ; "$cut_tag

        echo -e "\033[33m"$cut_tag"\033[0m \033[2m(cut)\033[0m"

        # loop event selection sample
        for input_mc_evt in "${INPUTS_MC_EVT[@]}" ; do
            IFS=';' ; input_mc_evt_tags=($input_mc_evt) ; unset IFS ; mc_evt_tag=${input_mc_evt_tags[2]} ;
            [[ ($cutevt_tag == *gammaN* && $mc_evt_tag == *BeamB*) || ($cutevt_tag == *Ngamma* && $mc_evt_tag == *BeamA*) || ($cutevt_tag == *0n0n* && $mc_evt_tag == *BeamB*) ]] && continue

            ######################
            # Event efficiency   #
            ######################
            itag_evteff=$ver'/'$cut_tag"/"$TAG_BINNING"/evteffsave_"$mc_evt_tag #
            echo -e "\033[33;2m"$ver'/'$cut_tag" / \033[0m\033[33m"$mc_evt_tag"\033[0m \033[2m(event efficiency)\033[0m"
            if [[ ${1:-0} -eq 1 || ${1:-0} -eq 3 ]] ; then
                ./evteff_save.exe "$input_mc_evt" "$cutevtstr" "${cutdstr} && $CUT_SIGNALWIN" $itag_evteff "$BINNING_Y" "$BINNING_PT"
            fi
            itag_evteff_calc=${itag_evteff/evteffsave/evteffcalc}
            if [[ ${1:-0} -eq 2 || ${1:-0} -eq 3 ]] ; then
                ./evteff_calc.exe "rootfiles/"$itag_evteff".root" $itag_evteff_calc $SAVE_PNG
            fi
            # echo -e "  \033[2m"$itag_evteff"\n  -> "$itag_evteff_calc"\033[0m"
        done

        # loop template
        for input_template in "${INPUTS_TEMPLATE[@]}" ; do
            
            IFS=';' ; input_template_tags=($input_template) ; unset IFS ; template_tag=${input_template_tags[2]} ;
            [[ ($cutevt_tag == *gammaN* && $template_tag == *BeamB*) || ($cutevt_tag == *Ngamma* && $template_tag == *BeamA*) || ($cutevt_tag == *0n0n* && $template_tag == *BeamB*) ]] && continue 
            IFS=';' ; templateopt_tags=($templateopt) ; unset IFS ; templateopt_tag=${templateopt_tags[2]} ;
            template_tag=$template_tag$templateopt_tag
            echo -e "\033[33;2m"$ver'/'$cut_tag" / \033[0m\033[33m"$template_tag"\033[0m \033[2m(mass template)\033[0m"

            ####################
            # Mass template    #
            ####################
            itag_template=$ver'/'$cut_tag"/"$TAG_BINNING"/template_"$template_tag #
            if [[ ${2:-0} -eq 1 ]] ; then
                ./hist_save.exe "$input_template" "$cutstr" $itag_template "$BINNING_Y" "$BINNING_PT" 1 "$templateopt" # 1: is_template
            fi
            # echo -e "  \033[2m"$itag_template"\033[0m"
            ## Loop data
            for input_data in "${INPUTS_DATA[@]}" ; do
                
                IFS=';' ; input_data_tags=($input_data) ; unset IFS ; data_tag=${input_data_tags[2]} ; data_lumi=${input_data_tags[3]} ; 
                echo -e "\033[33;2m"$ver'/'$cut_tag" / "$template_tag" / \033[0m\033[33m"$data_tag"\033[0m \033[2m(data)\033[0m"

                ####################
                # Fill data mass   # -> [ 13 min ]
                ####################
                itag_data=$ver'/'$cut_tag"/"$TAG_BINNING"/savehist_"$data_tag #
                [[ ${3:-0} -eq 1 || ${3:-0} -eq 3 ]] && {
                    ./hist_save.exe "$input_data" "$cutstr" $itag_data "$BINNING_Y" "$BINNING_PT" 0 "1" # 0: not template
                }

                ####################
                # Mass fitting     #
                ####################
                IFS=';' ; fitopts=($fitopt) ; unset IFS ; fit_tag=${fitopts[2]} ; 
                itag_data_fit=$ver'/'$cut_tag"/"$TAG_BINNING"/fithist_"$data_tag"_"$template_tag$fit_tag ## 
                [[ ${3:-0} -eq 2 || ${3:-0} -eq 3 ]] && {
                    ./hist_fit.exe "rootfiles/"$itag_data".root" "rootfiles/"$itag_template".root" $itag_data_fit "$fitopt" $SAVE_PNG
                }
                # echo -e "  \033[2m"$itag_data"\n  -> "$itag_data_fit"\033[0m"
                ## Loop MC
                for input_mc in "${INPUTS_MC[@]}" ; do

                    IFS=';' ; input_mc_tags=($input_mc) ; unset IFS ; mc_tag=${input_mc_tags[2]} ; pnp_flag=${input_mc_tags[3]}
                    [[ ($cutevt_tag == *gammaN* && $mc_tag == *BeamB*) || ($cutevt_tag == *Ngamma* && $mc_tag == *BeamA*) || ($cutevt_tag == *0n0n* && $mc_tag == *BeamB*) ]] && continue 
                    echo -e "\033[33;2m"$ver'/'$cut_tag" / "$template_tag" / "$data_tag" / \033[0m\033[33m"$mc_tag"\033[0m \033[2m(D efficiency, "$pnp_flag")\033[0m"

                    ####################
                    # D efficiency     # -> [ 36 min ]
                    ####################
                    itag_deff=$ver'/'$cut_tag"/deffsave_"$mc_tag # no binning info
                    if [[ ${4:-0} -eq 1 || ${4:-0} -eq 3 ]] ; then
                        ./eff_save.exe "$input_mc" "$cutevtstr" "$cutdstr" $itag_deff $pnp_flag "$cutd_weight"
                    fi

                    itag_deff_calc=$ver'/'$cut_tag"/"$TAG_BINNING"/deffcalc_"$mc_tag
                    [[ ${4:-0} -eq 2 || ${4:-0} -eq 3 ]] && {
                        ./eff_calc.exe "rootfiles/"$itag_deff".root" $itag_deff_calc "$BINNING_Y" "$BINNING_PT" $SAVE_PNG
                    }
                    itag_deff_calc_weight=$itag_deff_calc'_weight'

                    ####################
                    # Cross-section    #
                    ####################
                    lumi=$data_lumi
                    [[ x$cutevt_lumi != x ]] && lumi=$cutevt_lumi
                    itag_xsec=$ver'/'$cut_tag'/'$TAG_BINNING'/xsec_'${itag_data_fit##*/}'_'${itag_deff_calc##*/}'_'${itag_evteff_calc##*/}
                    # itag_xsec_w=$ver'/'$cut_tag'/'$TAG_BINNING'/xsec_'${itag_data_fit##*/}'_'${itag_deff_calc_weight##*/}'_'${itag_evteff_calc##*/}
                    echo "    itag_data_fit: "$itag_data_fit
                    echo "    itag_deff:     "$itag_deff_calc
                    echo "    itag_evteff:   "$itag_evteff_calc
                    echo "    lumi:          "$lumi" nb-1"
                    echo "                >> "$itag_xsec
                    [[ ${5:-0} -eq 1 ]] && {
                        ./xsec_calc.exe "rootfiles/"$itag_data_fit".root" "rootfiles/"$itag_deff_calc".root" "rootfiles/"$itag_evteff_calc".root" $lumi $itag_xsec $femdpu $SAVE_PNG
                        # [[ -f "rootfiles/"$itag_deff_calc_weight".root" ]] && {
                        #     ./xsec_calc.exe "rootfiles/"$itag_data_fit".root" "rootfiles/"$itag_deff_calc_weight".root" "rootfiles/"$itag_evteff_calc".root" $lumi $itag_xsec_w $SAVE_PNG                            
                        # }
                    }

                    [[ $cutevt_tag == *Ngamma* ]] || continue

                    itag_fprompt='null'
                    input_fprompt='null'
                    IFS=';' ; fprompts=($INPUT_FPROMPT) ; unset IFS ; fprompt_tag=${fprompts[2]} ; fprompt_file=${fprompts[0]}
                    [[ -f $fprompt_file ]] && { input_fprompt=$fprompt_file ; itag_fprompt=$fprompt_tag ; }
                    
                    file_Ngamma='rootfiles/'$itag_xsec'.root'
                    file_gammaN=${file_Ngamma//Ngamma/gammaN} ; file_gammaN=${file_gammaN//BeamB/BeamA} ;
                    outputname=${itag_xsec//-Ngamma/} ; outputname=${outputname//-BeamB/} ; outputname=$outputname'_'$itag_fprompt

                    # file_Ngamma_w='rootfiles/'$itag_xsec_w'.root'
                    # file_gammaN_w=${file_Ngamma_w//Ngamma/gammaN} ; file_gammaN_w=${file_gammaN_w//BeamB/BeamA} ;
                    # outputname_w=${itag_xsec_w//-Ngamma/} ; outputname_w=${outputname_w//-BeamB/} ; outputname_w=$outputname_w'_'$itag_fprompt

                    echo -e " \033[2m>> [collect]\033[0m "$outputname
                    echo "    input_fprompt:   "$input_fprompt
                    [[ -f $file_gammaN && -f $file_Ngamma ]] && {
                        [[ ${5:-0} -gt 0 ]] && {
                            ./xsec_collect.exe "${file_gammaN},${file_Ngamma}" "$input_fprompt" $outputname $SAVE_PNG
                            # [[ -f $file_Ngamma_w && -f $file_gammaN_w ]] && {
                            #     ./xsec_collect.exe "${file_gammaN_w},${file_Ngamma_w}" "$input_fprompt" $outputname_w $SAVE_PNG
                            # }
                        }
                    }
                    echo -e " \033[2m<< [collect]\033[0m"
                done
            done
        done
    done 
done
wait

# TAG_BINNING="b-yextend" ; BINNING_Y='-2.4, -2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2., 2.4' ; BINNING_PT='2., 5.' ;
# TAG_BINNING="b-ptlow" ; BINNING_Y='-2., -1., 0., 1., 2.' ; BINNING_PT='1., 2.' ;
# fitopt="Peaky;Double gaus signal;_f-2gaus"
# fitopt="3G;No KK/#pi#pi;_f-nopeaky"
