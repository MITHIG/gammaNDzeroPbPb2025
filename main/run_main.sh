#!/bin/bash

IVER=""
# IVER="-finebin"
INPUTS_DATA=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260426-yrefmva_PbPbUPC_HIForward_Dpt-2_Dsize_24PD.root;2025 PbPb (5.36 TeV);2025PbPb;0.06036"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_2023PbPbUPC_Jan2024ReReco_20260212Forest_HIForward_Dpt-2_Trig-2_Dsize_xbr.root;2023 PbPb (Jan2024 ReReco);2023PbPb-recoJan2024;0.007803"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260426-yrefmva_2023PbPbUPC_Feb2025ReReco_20260521Forest_HIForward_Dpt-2_Trig-2_Dsize.root;2023 PbPb (Feb2025 ReReco);2023PbPb-recoFeb2025;0.007803"
)
CUTEVTS=(
    "isL1ZDCOr && cscTightHalo2015Filter && ZDCgammaN && HFEMaxPlus_eta5 < 16 && selectedVtxFilter;#gammaN (25);gammaN-25"
    "isL1ZDCOr && cscTightHalo2015Filter && ZDCNgamma && HFEMaxMinus_eta5 < 16 && selectedVtxFilter;N#gamma (25);Ngamma-25"
    # "isL1ZDCOr && cscTightHalo2015Filter && ZDCgammaN && HFEMaxPlus_eta5 < 9.2 && selectedVtxFilter && ClusterCompatibilityFilter && nVtx <= 3;#gammaN (23);gammaN-23"
    # "isZeroBias && cscTightHalo2015Filter && ZDCsumPlus < 1000 && ZDCsumMinus < 1000 && (HFEMaxPlus_eta5 < 16 || HFEMaxMinus_eta5 < 16) && selectedVtxFilter;ZDC 0n0n, one-side gap;0n0n1gap-25"
)
INPUTS_TEMPLATE=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260426-yrefmva_HiForest_260328_prompt_GNucleusToD0-PhotonBeamA_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.85]{YTHIA}8#scale[0.5]{ }#gammaN (5.36 TeV);BeamA"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260426-yrefmva_HiForest_260328_prompt_GNucleusToD0-PhotonBeamB_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2_Dsize.root;P#scale[0.85]{YTHIA}8#scale[0.5]{ }#gammaN (5.36 TeV);BeamB"
)
INPUTS_MC=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260426-yrefmva_HiForest_260328_prompt_GNucleusToD0-PhotonBeamA_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2.root;P#scale[0.85]{YTHIA}8#scale[0.5]{ }#gammaN (5.36 TeV);2024-Pthat0-BeamA"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260426-yrefmva_HiForest_260328_prompt_GNucleusToD0-PhotonBeamB_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2.root;P#scale[0.85]{YTHIA}8 N#gamma (5.36 TeV);2024-Pthat0-BeamB"
)
CUT_BASE="Dpt>=2 && Dpt<5 && fabs(Dy)<2 && TMath::Abs(Dtrk1PtErr/Dtrk1Pt)<0.1 && TMath::Abs(Dtrk2PtErr/Dtrk2Pt)<0.1 && TMath::Abs(Dtrk1Eta) < 2.4 && TMath::Abs(Dtrk2Eta) < 2.4 && Dtrk1Pt > 0.5 && Dtrk2Pt > 0.5 && Dchi2cl > 0.05 && (DsvpvDistance/DsvpvDisErr) > 1. && DsvpvDisErr>1.e-8 && DsvpvDisErr_2D>1.e-8"
CUTDS=(
    # "${CUT_BASE} && (Dtrk1PixelHit+Dtrk1StripHit)>=11 && (Dtrk2PixelHit+Dtrk2StripHit)>=11 && DpassCut23PAS;Cut23PAS;D23pas"
    "${CUT_BASE} && ((Dy>=-2 && Dy<-1 && Dmva_BDT>0.143) || (Dy>=-1 && Dy<0 && Dmva_BDT>0.142) || (Dy>=0 && Dy<1 && Dmva_BDT>0.123) || (Dy>=1 && Dy<2 && Dmva_BDT>0.098));Optimized BDT;Dbdt-gammaN"
    "${CUT_BASE} && ((Dy<2 && Dy>=1 && Dmva_BDT>0.143) || (Dy<1 && Dy>=0 && Dmva_BDT>0.142) || (Dy<0 && Dy>=-1 && Dmva_BDT>0.123) || (Dy<-1 && Dy>=-2 && Dmva_BDT>0.098));Optimized BDT;Dbdt-Ngamma"
    # "${CUT_BASE} && Dmva_BDT>0;D pre-cuts;Dpre"
)

##

make save_hist.exe fit_hist.exe eff_save.exe eff_calc.exe xsec_calc.exe || exit 1

for cutevtstr in "${CUTEVTS[@]}" ; do
    IFS=';' ; cutevttags=($cutevtstr) ; unset IFS ; cutevt="${cutevttags[0]}" ; cutevt_tex="${cutevttags[1]}" ; cutevt_tag="${cutevttags[2]}"

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
        # echo -e "\033[2m"$cutstr"\033[0m"

        for input_template in "${INPUTS_TEMPLATE[@]}" ; do
            
            IFS=';' ; input_template_tags=($input_template) ; unset IFS ; template_tag=${input_template_tags[2]} ;
            [[ ($cutevt_tag == *gammaN* && $template_tag == *BeamB*) || ($cutevt_tag == *Ngamma* && $template_tag == *BeamA*) ]] && continue 
            echo -e "  \033[33m"$template_tag"\033[0m"

            ####################
            # Mass template    #
            ####################
            itag_template=$cut_tag"/mass_templates"$IVER #
            if [[ ${1:-0} -eq 1 ]] ; then
                ./save_hist.exe "$input_template" "$cutstr" $itag_template 0
                # ./save_hist.exe "$INPUT_MASS_TEMPLATE" "${CUT_BASE};D precut;Dpre" $itag_template 0
            elif [[ ${1:-0} -eq 2 ]] ; then
                template_dir="rootfiles/"$itag_template".root"
                template_dir=${template_dir%/*}
                mkdir -p $template_dir
                cp -v "rootfiles/"gammaN-25_Dbdt/mass_templates.root "rootfiles/"$itag_template".root"
            fi

            ## Loop data
            for input_data in "${INPUTS_DATA[@]}" ; do
                
                IFS=';' ; input_data_tags=($input_data) ; unset IFS ; data_tag=${input_data_tags[2]} ; data_lumi=${input_data_tags[3]} ;
                echo -e "    \033[33m"$data_tag"\033[0m "$data_lumi" nb-1"

                ####################
                # Fill data mass   #
                ####################
                itag_data=$cut_tag"/savehist_"$data_tag$IVER
                [[ ${2:-0} -eq 1 ]] && {
                    ./save_hist.exe "$input_data" "$cutstr" $itag_data 
                }

                ####################
                # Mass fitting     #
                ####################
                [[ ${3:-0} -eq 1 ]] && {
                    ./fit_hist.exe "rootfiles/"$itag_data".root" "rootfiles/"$itag_template".root"
                }
                itag_data_fit=$cut_tag"/fithist_"$data_tag$IVER ## 

                ## Loop MC
                for input_mc in "${INPUTS_MC[@]}" ; do

                    IFS=';' ; input_mc_tags=($input_mc) ; unset IFS ; mc_tag=${input_mc_tags[2]} ;
                    [[ ($cutevt_tag == *gammaN* && $mc_tag == *BeamB*) || ($cutevt_tag == *Ngamma* && $mc_tag == *BeamA*) ]] && continue 
                    echo -e "      \033[33m"$mc_tag"\033[0m"                    

                    ####################
                    # D efficiency     #
                    ####################
                    itag_deff=$cut_tag"/saveeff_"$mc_tag"_"$data_tag$IVER
                    [[ ${4:-0} -eq 1 ]] && {
                        # ./eff_save.exe "$input_mc" "$cutevtstr" "$cutdstr" $itag_deff "$input_data"
                        ./eff_save.exe "$input_mc" "$cutevtstr" "$cutdstr" $itag_deff null
                    }
                    [[ ${5:-0} -eq 1 ]] && {
                        ./eff_calc.exe "rootfiles/"$itag_deff".root"
                    }
                    itag_deff=$cut_tag"/calceff_"$mc_tag"_"$data_tag$IVER ##
                    
                    ####################
                    # Cross-section    #
                    ####################
                    echo "      itag_data_fit:     "$itag_data_fit
                    echo "      itag_data_deff:    "$itag_deff
                    echo "      itag_data_evteff:  null"
                    echo "      itag_data_fprompt: null"
                    [[ ${6:-0} -eq 1 ]] && {
                        ./xsec_calc.exe "rootfiles/"$itag_data_fit".root" "rootfiles/"$itag_deff".root" null null $data_lumi $cut_tag
                    }
                done
            done
        done
    done 
done

