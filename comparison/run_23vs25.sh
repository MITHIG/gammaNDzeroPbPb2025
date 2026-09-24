#!/bin/bash

DRAW_PNG=1

INPUTS=(
    # ZeroBias
    '/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_HiForest_260218_HIPhysicsRawPrime0-5_HIRun2023A_ZB_374970.root;2023 ZB (374970);d23-rp'
    '/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260212-hfle_HiForest_260218_HIPhysicsRawPrime0-10_HIRun2025A_highrZB_399766.root;2025 ZB (399766);d25-rp'

    # # Empty BX
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_HiForest_260218_HIEmptyBX_HIRun2023A_PromptReco_v2.root;2023 EmptyBX;d23-rp"
    # "/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260212-hfle_HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1.root;2025 EmptyBX;d25-rp"
)

VARS=(
    'ZDCsumPlus;1'
    'ZDCsumMinus;1'

    # 'ZDCsumPlus-low;1'
    # 'ZDCsumMinus-low;1'
    # 'HFEMaxPlusforest-zoom;1'
    # 'HFEMaxMinusforest-zoom;1'

    # HFEMaxPlusforest
    # HFEMaxMinusforest
    # nTrackInAcceptanceHP
    # nVtx
    
    #
    # Dmass
    # Dalpha
    # Ddls
    # Dtrk1Pt
    # Dtrk2Pt
    # Dchi2cl
    # Dtrk1ptrel
    # Dtrk2ptrel
    # Dtrk1nhit
    # Dtrk2nhit
    # # Dalpha-zoom
)

CUTEVTS=(
    "1;Zero Bias;zerobias;"
    # "isNotBptxOR;HLT_HIL1NotBptxOR;isNotBptxOR;"
)
##

CUTDS=(
    "1;;" # for event variables 
    # "Dpt>0;;-Dnocut"
    # "TMath::Abs(Dtrk1PtErr/Dtrk1Pt)<0.1 && TMath::Abs(Dtrk2PtErr/Dtrk2Pt)<0.1 && TMath::Abs(Dtrk1Eta) < 2.4 && TMath::Abs(Dtrk2Eta) < 2.4 && Dtrk1Pt > 0.5 && Dtrk2Pt > 0.5 && Dchi2cl > 0.05 && (DsvpvDistance/DsvpvDisErr) > 1. && DsvpvDisErr>1.e-8 && DsvpvDisErr_2D>1.e-8;Precuts;-Dprecut"
    # "Dtrk1PtErr/Dtrk1Pt<0.1 && Dtrk2PtErr/Dtrk2Pt<0.1 && DpassCut23PAS && (Dtrk1PixelHit+Dtrk1StripHit)>=11 && (Dtrk2PixelHit+Dtrk2StripHit)>=11;DpassCut23PAS;-D23pas"
)

make savehist.exe calchists.exe drawhists.exe || exit 1

# D cut
for cutdstr in "${CUTDS[@]}" ; do
    IFS=';' ; cutdtags=($cutdstr) ; unset IFS ; cutd=${cutdtags[0]} ; cutd_tex=${cutdtags[1]} ; cutd_tag=${cutdtags[2]} ; 

    for varp in "${VARS[@]}" ; do
        IFS=';' ; varps=(${varp}) ; unset IFS ; var=${varps[0]} ; noratio=${varps[1]} ;
        [[ ($var == D* && ${cutd} != 1) || ($var != D* && ${cutd} == 1) ]] || { continue ; }

        # event cut
        for cutevtstr in "${CUTEVTS[@]}" ; do
            IFS=';' ; cutevttags=($cutevtstr) ; unset IFS ; cutevt=${cutevttags[0]} ; cutevt_tex=${cutevttags[1]} ; cutevt_tag=${cutevttags[2]} ; cutevt_year=${cutevttags[3]}

            cut_tag=$cutevt_tag$cutd_tag
            cutstr=$cutevt" && "$cutd";"$cutevt_tex"%%"$cutd_tex";"$cut_tag
            
            echo -e "\033[33m"$var" \033[33;2m("$cut_tag")\033[0m \033[2m"$cutstr"\033[0m"

            compare_list=''
            tag_list=''
            for inputstr in "${INPUTS[@]}" ; do
                IFS=';' ; inputtags=($inputstr) ; unset IFS ; input_tag=${inputtags[2]} ; 
                [[ x$input_tag == x ]] && { echo "warning: missed input_tag. skip." ; continue ; }
                [[ ($var == D* && ${inputtags[0]} == *Dsize*) || ($var != D* && ${inputtags[0]} != *Dsize*) ]] || { continue ; }
                [[ ($cutevt_year == 2025 && $input_tag == d23*) || ($cutevt_year == 2023 && $input_tag == d25*) ]] && { continue ; } #

                echo -e "    \033[33m"$var" \033[33;2m("$input_tag")\033[0m"

                itag="rootfiles/"$cut_tag"/"$var"_"$input_tag"_savehist" #
                echo "    "$itag

                [[ ${1:-0} -eq 1 ]] && {
                    ./savehist.exe "$inputstr" "$cutstr" "$var" $itag &
                }

                [[ ${2:-0} -eq 1 ]] && {
                    ./calchists.exe $itag".root"
                }

                itag=${itag/_savehist/_calchist}
                compare_list=$compare_list","$itag".root"
                tag_list=$tag_list"_"$cut_tag"-"$input_tag
                itag=''
            done
            compare_list=${compare_list#,}
            tag_list=${tag_list#_}

            echo $compare_list
            echo $tag_list
            [[ ${3:-0} -eq 1 ]] && {
                ./drawhists.exe "$compare_list" "$tag_list" $DRAW_PNG $noratio
            }
        done

        echo
        echo "--> manual_draw_list" 
        manual_draw_list=(
        )
        for items in "${manual_draw_list[@]}" ; do
            IFS=';' ; draw_opts=($items) ; unset IFS ; do_save_png=${draw_opts[1]}

            compare_list=
            tag_list=
            IFS=',' ; draw_tags=(${draw_opts[0]}) ; unset IFS ;
            for itag in "${draw_tags[@]}" ; do
                jtag_cut=${itag%%__*}
                jtag_input=${itag##*__}
                compare_list=$compare_list",rootfiles/"$jtag_cut"/"$var"_"$jtag_input"_calchist.root"
                tag_list=$tag_list"_"${itag/__/-}
            done
            compare_list=${compare_list#,}
            tag_list=${tag_list#_}
            echo $compare_list
            echo $tag_list
            [[ ${4:-0} -eq 1 ]] && {
                ./drawhists.exe "$compare_list" "$tag_list" $do_save_png
            }
        done
        echo
    done
    wait
done
