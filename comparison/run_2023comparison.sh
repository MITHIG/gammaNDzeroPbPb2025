#!/bin/bash

INPUTS=(
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_2023PbPbUPC_Jan2024ReReco_20260212Forest_HIForward_Dpt-2_Trig-2.root;2023 (Jan24 Reco);23-rJan24"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260426-yrefmva_2023PbPbUPC_Feb2025ReReco_20260521Forest_HIForward_Dpt-2_Trig-2.root;2023 (Feb25 Reco);23-rFeb25"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_2023PbPbUPC_Jan2024ReReco_20260212Forest_HIForward_Dpt-2_Trig-2_Dsize_xbr.root;2023 (Jan24 Reco);23-rJan24"
    "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260426-yrefmva_2023PbPbUPC_Feb2025ReReco_20260521Forest_HIForward_Dpt-2_Trig-2_Dsize.root;2023 (Feb25 Reco);23-rFeb25"
)

VARS=(
    nTrackInAcceptanceHP
    ZDCsumPlus
    ZDCsumMinus
    HFEMaxPlusforest
    HFEMaxPlusforest-Low
    HFEMaxMinusforest
    #
    Dmass
    Dalpha
    Dalpha-Low
    Ddls
    Dtrk1Pt
    Dtrk2Pt
    Dchi2cl
    Dtrk1ptrel
    Dtrk2ptrel
    Dtrk1nhit
    Dtrk2nhit
)

CUTS=(
    # "(1);HLT_ZDCOr && PV filter;nocut"
    "(ZDCgammaN && gapgammaN && selectedBkgFilter && selectedVtxFilter);#gammaN (ZDCXOR, gap), vertex, bkgrej;gammaN"
    #
    "(ZDCgammaN && gapgammaN && selectedBkgFilter && selectedVtxFilter && Dpt>2 && Dpt<5 && fabs(Dy)<2);#gammaN (ZDCXOR, gap), vertex, bkgrej%%2 < p_{T} < 5 GeV, |y| < 2;gammaN-Dpre"
    "(ZDCgammaN && gapgammaN && selectedBkgFilter && selectedVtxFilter && Dpt>2 && Dpt<5 && fabs(Dy)<2 && Dtrk1PtErr/Dtrk1Pt<0.1 && Dtrk2PtErr/Dtrk2Pt<0.1 && DpassCut23PAS && (Dtrk1PixelHit+Dtrk1StripHit)>=11 && (Dtrk2PixelHit+Dtrk2StripHit)>=11);#gammaN (ZDCXOR, gap), vertex, bkgrej%%2<p_{T}<5 GeV, |y| < 2%%DpassCut23PAS;gammaN-D23pas"
)
##

make savehist.exe calchists.exe drawhists.exe || exit 1

for var in "${VARS[@]}" ; do
    for cutstr in "${CUTS[@]}" ; do
        IFS=';' ; cuttags=($cutstr) ; unset IFS ; cut_tag=${cuttags[2]} ; 

        [[ ($var == D* && ${cutstr} == *Dpt*) || ($var != D* && ${cutstr} != *Dpt*) ]] || { continue ; }
        [[ x$cut_tag == x ]] && { echo "warning: missed cut_tag. skip." ; continue ; }
        
        echo -e "\033[33m"$var" \033[33;2m("$cut_tag")\033[0m"

        compare_list=''
        tag_list=''
        for inputstr in "${INPUTS[@]}" ; do
            IFS=';' ; inputtags=($inputstr) ; unset IFS ; input_tag=${inputtags[2]} ; 
            [[ x$input_tag == x ]] && { echo "warning: missed input_tag. skip." ; continue ; }
            [[ ($var == D* && ${inputtags[0]} == *Dsize*) || ($var != D* && ${inputtags[0]} != *Dsize*) ]] || { continue ; }

            echo -e "    \033[33m"$var" \033[33;2m"${inputtags[0]}"\033[0m"

            itag="rootfiles/"$var"/"$cut_tag"/"$input_tag"_savehist"

            [[ ${1:-0} -eq 1 ]] && {
                ./savehist.exe "$inputstr" "$cutstr" "$var" $itag
            }

            [[ ${2:-0} -eq 1 ]] && {
                ./calchists.exe $itag".root"
            }

            itag=${itag/_savehist/_calchist}
            compare_list=$compare_list","$itag".root"
            tag_list=$tag_list"_"$input_tag
            itag=''
        done

        compare_list=${compare_list#,}
        tag_list=${tag_list#_}"__"$cut_tag

        [[ ${3:-0} -eq 1 ]] && {
            echo $compare_list
            echo $tag_list
            # ./drawhists.exe "rootfiles/${var}/calc_${tag1}.root","rootfiles/${var}/calc_${tag2}.root" "Jan2024ReReco-Feb2025ReReco__"${cut_tag}
        }
        
    done
done
