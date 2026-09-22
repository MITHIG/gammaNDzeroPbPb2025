#!/bin/bash

isNgamma=${1}

make evteff_save.exe evteff_calc.exe hist_save.exe hist_fit.exe eff_save.exe eff_calc.exe xsec_calc.exe || exit 1

parse_input_tag() {
    IFS=';' ; inputs=($1) ; unset IFS ;
    echo "${inputs[2]}"
}

evtgaps=(
    'HFEMax_eta5 < 16;;'
    
    'HFEMax_eta5 < 23;Rap gap 23 GeV;-gap23'
    'HFEMax_eta5 < 22;Rap gap 22 GeV;-gap22'
    'HFEMax_eta5 < 21;Rap gap 21 GeV;-gap21'
    'HFEMax_eta5 < 20;Rap gap 20 GeV;-gap20'
    'HFEMax_eta5 < 19;Rap gap 19 GeV;-gap19'
    'HFEMax_eta5 < 18;Rap gap 18 GeV;-gap18'
    'HFEMax_eta5 < 17;Rap gap 17 GeV;-gap17'
    'HFEMax_eta5 < 15;Rap gap 15 GeV;-gap15'
    'HFEMax_eta5 < 14;Rap gap 14 GeV;-gap14'
    'HFEMax_eta5 < 13;Rap gap 13 GeV;-gap13'
    'HFEMax_eta5 < 12;Rap gap 12 GeV;-gap12'
    'HFEMax_eta5 < 11;Rap gap 11 GeV;-gap11'
    'HFEMax_eta5 < 10;Rap gap 10 GeV;-gap10'
    'HFEMax_eta5 < 9;Rap gap 9 GeV;-gap9'
    'HFEMax_eta5 < 8;Rap gap 8 GeV;-gap8'
    'HFEMax_eta5 < 7;Rap gap 7 GeV;-gap7'
    'HFEMax_eta5 < 6;Rap gap 6 GeV;-gap6'
    'HFEMax_eta5 < 5;Rap gap 5 GeV;-gap5'
    'HFEMax_eta5 < 4;Rap gap 4 GeV;-gap4'
)

make_dcut_string() {
    IFS=',' ; cutvalues=($1) ; unset IFS ;
    [[ $isNgamma -eq 0 ]] && { echo '((Dy<-1 && Dmva_BDT>'${cutvalues[0]}') || (Dy>=-1 && Dy<0 && Dmva_BDT>'${cutvalues[1]}') || (Dy>=0 && Dy<1 && Dmva_BDT>'${cutvalues[2]}') || (Dy>=1 && Dmva_BDT>'${cutvalues[3]}'))' ; }
    [[ $isNgamma -eq 1 ]] && { echo '((Dy>=1 && Dmva_BDT>'${cutvalues[0]}') || (Dy<1 && Dy>=0 && Dmva_BDT>'${cutvalues[1]}') || (Dy<0 && Dy>=-1 && Dmva_BDT>'${cutvalues[2]}') || (Dy<-1 && Dmva_BDT>'${cutvalues[3]}'))' ; }
}

cutdtopos=(
    "`make_dcut_string 0.143,0.142,0.123,0.098`;;"

    "`make_dcut_string 0.123,0.122,0.103,0.078`;BDT shift -0.02;-bdtsM0p02"
    "`make_dcut_string 0.103,0.102,0.083,0.058`;BDT shift -0.04;-bdtsM0p04"
    "`make_dcut_string 0.083,0.082,0.063,0.038`;BDT shift -0.06;-bdtsM0p06"
    "`make_dcut_string 0.063,0.062,0.043,0.018`;BDT shift -0.08;-bdtsM0p08"
    "`make_dcut_string 0.043,0.042,0.023,-0.002`;BDT shift -0.1;-bdtsM0p10"
    "`make_dcut_string 0.023,0.022,0.003,-0.022`;BDT shift -0.12;-bdtsM0p12"
    "`make_dcut_string 0.003,0.002,-0.017,-0.042`;BDT shift -0.14;-bdtsM0p14"
    "`make_dcut_string -0.017,-0.018,-0.037,-0.062`;BDT shift -0.16;-bdtsM0p16"
    "`make_dcut_string -0.037,-0.038,-0.057,-0.082`;BDT shift -0.18;-bdtsM0p18"
    "`make_dcut_string -0.057,-0.058,-0.077,-0.102`;BDT shift -0.20;-bdtsM0p20"

    "`make_dcut_string 0.163,0.162,0.143,0.118`;BDT shift 0.02;-bdts0p02"
    "`make_dcut_string 0.183,0.182,0.163,0.138`;BDT shift 0.04;-bdts0p04"
    "`make_dcut_string 0.203,0.202,0.183,0.158`;BDT shift 0.06;-bdts0p06"
    "`make_dcut_string 0.223,0.222,0.203,0.178`;BDT shift 0.08;-bdts0p08"
    "`make_dcut_string 0.243,0.242,0.223,0.198`;BDT shift 0.10;-bdts0p10"
)

fitopts=(
    '3G-Peaky;;'
    # '3G-Peaky-Exp;;_f-exp'
)

make_zdccut_string() {
    IFS=',' ; cutvalues=($1) ; unset IFS ;
    [[ $isNgamma -eq 0 ]] && { echo 'ZDCsumPlus < '${cutvalues[0]}' && ZDCsumMinus > '${cutvalues[1]} ; }
    [[ $isNgamma -eq 1 ]] && { echo 'ZDCsumPlus > '${cutvalues[0]}' && ZDCsumMinus < '${cutvalues[1]} ; }
}
evtzdcs=(
    "`make_zdccut_string 1100,1000`;;"
    
    # "`make_zdccut_string 1000,1000`;ZDC 1000 GeV;-zdc1000"
    # "`make_zdccut_string 1100,1100`;ZDC 1100 GeV;-zdc1100"
    # "`make_zdccut_string 1200,1200`;ZDC 1200 GeV;-zdc1200"
    # "`make_zdccut_string 900,900`;ZDC 900 GeV;-zdc900"
)

for evtzdc in "${evtzdcs[@]}" ; do
    evtzdc_tag=`parse_input_tag "$evtzdc"`
    [[ x$evtzdc_tag == x ]] && def_evtzdc=1 || def_evtzdc=0
    # echo $evtzdc_tag' -> default : '$def_evtzdc

    for evtgap in "${evtgaps[@]}" ; do
        evtgap_tag=`parse_input_tag "$evtgap"`
        [[ x$evtgap_tag == x ]] && def_evtgap=1 || def_evtgap=0
        # echo $evtgap_tag' -> default : '$def_evtgap

        for cutdtopo in "${cutdtopos[@]}" ; do
            cutdtopo_tag=`parse_input_tag "$cutdtopo"`
            [[ x$cutdtopo_tag == x ]] && def_cutdtopo=1 || def_cutdtopo=0

            for fitopt in "${fitopts[@]}" ; do
                fitopt_tag=`parse_input_tag "$fitopt"`
                [[ x$fitopt_tag == x ]] && def_fitopt=1 || def_fitopt=0

                [[ $(($def_evtzdc+$def_evtgap+$def_cutdtopo+$def_fitopt)) -eq 3 ]] || continue

                echo -n "--"
                [[ $def_evtzdc -eq 0 ]] && echo -en " \033[1m["$evtzdc"]\033[0m" || echo -en " \033[2m["$evtzdc"]\033[0m"
                [[ $def_evtgap -eq 0 ]] && echo -en " \033[1m["$evtgap"]\033[0m" || echo -en " \033[2m["$evtgap"]\033[0m"
                [[ $def_cutdtopo -eq 0 ]] && echo -en " \033[1m["$cutdtopo"]\033[0m" || echo -en " \033[2m["$cutdtopo"]\033[0m"
                [[ $def_fitopt -eq 0 ]] && echo -en " \033[1m["$fitopt"]\033[0m" || echo -en " \033[2m["$fitopt"]\033[0m"
                echo
                
                [[ $def_evtzdc -eq 0 ]] && ./run_var.sh $isNgamma 00110 "$evtzdc" "$evtgap" "$cutdtopo" "$fitopt" ${2:-0}
                [[ $def_evtgap -eq 0 ]] && ./run_var.sh $isNgamma 10110 "$evtzdc" "$evtgap" "$cutdtopo" "$fitopt" ${2:-0} &
                [[ $def_cutdtopo -eq 0 ]] && ./run_var.sh $isNgamma 01111 "$evtzdc" "$evtgap" "$cutdtopo" "$fitopt" ${2:-0} &
                [[ $def_fitopt -eq 0 ]] && ./run_var.sh $isNgamma 00010 "$evtzdc" "$evtgap" "$cutdtopo" "$fitopt" ${2:-0}
                
            done            
        done
    done    
done

wait
