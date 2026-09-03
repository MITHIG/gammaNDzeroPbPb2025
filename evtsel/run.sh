#!/bin/bash

INPUTS=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260426-yrefmva_PbPbUPC_HIForward_Dpt-2_Dsize_24PD.root;2025 PbPb;2025PbPb-Denh;0.060361" # 2025
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260426-yrefmva_PbPbUPC_HIForward17_Dpt-2.root;2025 PbPb;2025PbPb;" # 2025
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260426-yrefmva_2023PbPbUPC_Feb2025ReReco_20260521Forest_HIForward_Dpt-2_Trig-2_Dsize.root;2023 PbPb #it{Feb2025};2023PbPb-recoFeb2025;0.007803" # 2023-Feb2025
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_2023PbPbUPC_Jan2024ReReco_20260212Forest_HIForward_Dpt-2_Trig-2_Dsize_xbr.root;2023 PbPb #it{Jan2024};2023PbPb-recoJan2024;0.007803" # 2023-Jan2024
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_prompt_GNucleusToD0-PhotonBeamA_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2.root;P#scale[0.8]{YTHIA}8 #gammap;2024-SoftQCD-BeamA;1"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2024PbPb/Dzero_260714-gen_HiForest_260328_prompt_GNucleusToD0-PhotonBeamB_Bin-Pthat0_Kpi_trkpt0p1_Drej-genmatched_Dpt-2.root;P#scale[0.8]{YTHIA}8 p#gamma;2024-SoftQCD-BeamB;1"
)

EVTSELS=(
    "L1ZDCOr-gammaN-2025;;"
    "L1ZDCOr-Ngamma-2025;;"
    "L1ZeroBias-gammaN-2025;;"
    "L1ZeroBias-Ngamma-2025;;"
    # "L1ZDCOr-gammaN-2023;;"
    # "L1ZDCOr-Ngamma-2023;;"
)

DSELS=(
    "Dnoreq;;"
    "Dyrefbdt_sig;Analysis BDT cut, in signal window;"
    # "Dpascut_sig;PAS cut, in signal window;"
    # Dyinclbdt_sign
)

# TAG_BINNING="_b-default" ; BINNING_Y='-2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5, 2.' ; BINNING_PT='2., 5.' ;
TAG_BINNING="_b-ycoarse" ; BINNING_Y='-2., -1., 0., 1., 2.' ; BINNING_PT='2., 5.' ;

make savehist.exe cookhist.exe || exit 1

for evtselstr in "${EVTSELS[@]}" ; do
    IFS=';' ; evtseltags=($evtselstr) ; unset IFS ; evtsel_tag="${evtseltags[0]}" ;
    [[ x$evtsel_tag == x ]] && { echo "warning: missed evtsel_tag. skip." ; continue ; }

    for dselstr in "${DSELS[@]}" ; do
        IFS=';' ; dseltags=($dselstr) ; unset IFS ; dsel_tag="${dseltags[0]}" ;
        [[ x$dsel_tag == x ]] && { echo "warning: missed dsel_tag. skip." ; continue ; }

        sel_tag=$evtsel_tag"_"$dsel_tag
        echo -e "\033[33m"$sel_tag"\033[0m"

        for inputstr in "${INPUTS[@]}" ; do
            IFS=';' ; inputtags=($inputstr) ; unset IFS ; input_tag="${inputtags[2]}" ; lumi="${inputtags[3]}" ;
            [[ x$input_tag == x ]] && { echo "warning: missed input_tag. skip." ; continue ; }

            [[ ($evtsel_tag == *2023* && $input_tag == *2025PbPb*) || ($evtsel_tag == *2025* && $input_tag == *2023PbPb*) ]] && continue
            [[ ($evtsel_tag == *gammaN* && $input_tag == *BeamB*) || ($evtsel_tag == *Ngamma* && $input_tag == *BeamA*) ]] && continue
            [[ $evtsel_tag == *2023* && $input_tag == *Beam* ]] && continue
            [[ $dsel_tag == *bdt* && $input_tag == *Jan2024* ]] && continue
            [[ ("$inputstr" == *_Dsize* && $dsel_tag == *Dnoreq*) || ("$inputstr" != *_Dsize* && $dsel_tag != *Dnoreq*) ]] && continue
            [[ $evtsel_tag == *ZeroBias* && $dsel_tag != *Dnoreq* ]] && continue

            echo -e "\033[33;2m"$sel_tag"\033[0m \033[33m/ "$input_tag"\033[0m"
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
