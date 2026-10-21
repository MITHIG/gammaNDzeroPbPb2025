#!/bin/bash

INPUTS=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260928-symmva_PbPbUPC_HIForward_Dpt-2_Dsize_24PD.root;2025 PbPb;2025PbPb;1" # 2025
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260426-yrefmva_PbPbUPC_HIForward17_Dpt-2.root;2025 PbPb;2025PbPb;1" # 2025
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260426-yrefmva_2023PbPbUPC_Feb2025ReReco_20260521Forest_HIForward_Dpt-2_Trig-2_Dsize.root;2023 PbPb #it{Feb2025};2023PbPb-recoFeb2025;1" # 2023-Feb2025
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_2023PbPbUPC_Jan2024ReReco_20260212Forest_HIForward_Dpt-2_Trig-2_Dsize_xbr.root;2023 PbPb #it{Jan2024};2023PbPb-recoJan2024;1" # 2023-Jan2024

    '/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_prompt_GNucleusToD0-BeamA_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-2.root;P#scale[0.8]{YTHIA}8 #gammaN;2025-SoftQCD-BeamA;0'
    '/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPbMC/Dzero_260714-gen_HiForest_260904_prompt_GNucleusToD0-BeamB_SoftQCD_KPi_2025_trkpt0p1_Drej-genmatched_Dpt-2.root;P#scale[0.8]{YTHIA}8 N#gamma;2025-SoftQCD-BeamB;0'
)

EVTSELS=(
    "L1ZDCOr-gammaN-2025;L1_ZDCOr + PV filter + cscTightHalo,#gamma #rightarrow Z+ (#gammaN);"
    "L1ZDCOr-Ngamma-2025;L1_ZDCOr + PV filter + cscTightHalo,#gamma #rightarrow Z- (N#gamma);"
    # "L1ZeroBias-gammaN-2025;L1_ZDCOr + PV filter + cscTightHalo,#gamma #rightarrow Z+;"
    # "L1ZeroBias-Ngamma-2025;L1_ZDCOr + PV filter + cscTightHalo,#gamma #rightarrow Z-;"
    # "MC-gammaN-2025;PV filter,#gamma #rightarrow Z+;"
    # "MC-Ngamma-2025;PV filter,#gamma #rightarrow Z-;"
    # "L1ZDCOr-gammaN-2023;;"
    # "L1ZDCOr-Ngamma-2023;;"
)
CUT_BASE="TMath::Abs(Dtrk1PtErr/Dtrk1Pt)<0.1 && TMath::Abs(Dtrk2PtErr/Dtrk2Pt)<0.1 && TMath::Abs(Dtrk1Eta) < 2.4 && TMath::Abs(Dtrk2Eta) < 2.4 && Dtrk1Pt > 0.5 && Dtrk2Pt > 0.5 && Dchi2cl > 0.05 && (DsvpvDistance/DsvpvDisErr) > 1. && DsvpvDisErr>1.e-8 && DsvpvDisErr_2D>1.e-8"
DSELS=(
    # "(1);;Dnoreq"
    "${CUT_BASE} && ((Dy<-1 && Dmva_BDT>0.143) || (Dy>=-1 && Dy<0 && Dmva_BDT>0.142) || (Dy>=0 && Dy<1 && Dmva_BDT>0.123) || (Dy>=1 && Dmva_BDT>0.098));BDT && in signal window;Dbdt-gammaN"
    "${CUT_BASE} && ((Dy>=1 && Dmva_BDT>0.143) || (Dy<1 && Dy>=0 && Dmva_BDT>0.142) || (Dy<0 && Dy>=-1 && Dmva_BDT>0.123) || (Dy<-1 && Dmva_BDT>0.098));BDT && in signal window;Dbdt-Ngamma"
    # "Dyrefbdt_sig;Analysis BDT cut, in signal window;Dyrefbdt_sig"
)

# TAG_BINNING="_b-default" ; BINNING_Y='-2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.' ; BINNING_PT='2., 5.' ;
TAG_BINNING="_b-ycoarse" ; BINNING_Y='-2., -1., 0., 1., 2.' ; BINNING_PT='2., 5.' ;

make savehist.exe cookhist.exe || exit 1

for evtselstr in "${EVTSELS[@]}" ; do
    IFS=';' ; evtseltags=($evtselstr) ; unset IFS ; evtsel_tag="${evtseltags[0]}" ;
    [[ x$evtsel_tag == x ]] && { echo "warning: missed evtsel_tag. skip." ; continue ; }

    for dselstr in "${DSELS[@]}" ; do
        IFS=';' ; dseltags=($dselstr) ; unset IFS ; dsel_tag="${dseltags[2]}" ;
        [[ x$dsel_tag == x ]] && { echo "warning: missed dsel_tag. skip." ; continue ; }

        sel_tag=$evtsel_tag"_"$dsel_tag
        # echo -e "\033[33m"$sel_tag"\033[0m"

        for inputstr in "${INPUTS[@]}" ; do
            IFS=';' ; inputtags=($inputstr) ; unset IFS ; input_tag="${inputtags[2]}" ; isdata="${inputtags[3]}" ;
            [[ x$input_tag == x ]] && { echo "warning: missed input_tag. skip." ; continue ; }

            # common
            [[ ($evtsel_tag == *gammaN* && $dsel_tag == *Ngamma*) || ($evtsel_tag == *Ngamma* && $dsel_tag == *gammaN*) ]] && continue
            [[ $evtsel_tag == *ZeroBias* && $dsel_tag != *Dnoreq* ]] && continue
            [[ ($evtsel_tag == *gammaN* && $input_tag == *BeamB*) || ($evtsel_tag == *Ngamma* && $input_tag == *BeamA*) ]] && continue
            [[ ($evtsel_tag == *MC* && $isdata -gt 0) ]] && continue
            [[ ("$inputstr" == *_Dsize* && $dsel_tag == *Dnoreq*) || ("$inputstr" != *_Dsize* && $dsel_tag != *Dnoreq*) ]] && continue
            # year matching
            [[ ($evtsel_tag == *2023* && $input_tag == *2025PbPb*) || ($evtsel_tag == *2025* && $input_tag == *2023PbPb*) ]] && continue
            [[ $evtsel_tag == *2023* && $input_tag == *Beam* ]] && continue
            [[ $dsel_tag == *bdt* && $input_tag == *Jan2024* ]] && continue 

            echo -e "\033[33;2m"$sel_tag"\033[0m \033[33m/ "$input_tag"\033[0m"
            echo -e "\033[2m"$inputstr"\033[0m"
            itag_savehist=$sel_tag"/savehist_"$input_tag$TAG_BINNING
            [[ ${1:-0} -eq 1 ]] && {
                ./savehist.exe "$inputstr" "$evtselstr" "$dselstr" $itag_savehist "$BINNING_Y" "$BINNING_PT" &
            }
            itag_cookhist=$sel_tag"/cookhist_"$input_tag$TAG_BINNING
            [[ ${2:-0} -eq 1 ]] && {
                ./cookhist.exe "rootfiles/"$itag_savehist".root" $itag_cookhist
            }
        done
    done    
done

wait        
