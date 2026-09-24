#!/bin/bash

BINNING_TAG='b-default'
# BINNING_TAG='b-ptdiff'

SAVE_PNG=1

FILENAME="xsec_fithist_2025PbPb_2025-SoftQCD_deffcalc_2025-SoftQCD_evteffcalc_2025-SoftQCD_null"
filename() {
    # dirname=${1//__DIR__/${tdir}}
    fname='rootfiles/'${1}'/'$BINNING_TAG'/'$FILENAME'.root'
    echo $fname
}

INPUTS=()
add_input() {
    local output_name="$1"
    local label="$2"
    local legdx="${3-}"
    local legdy="${4-}"
    local numscan="${5}"
    local varscan="${6}"
    local is_systmatics="${7}"
    shift 7

    local IFS=,
    local input_files="$*"
    INPUTS+=("${input_files}^${output_name}^${label}^${legdx}^${legdy}^${numscan}^${varscan}^${is_systmatics}")
}

add_input '0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/gap' '' 0 0 1 'Rapidity Gap Threshold' 1 \
          "`filename 0nXn-25_Dbdt`;Rap gap 16 GeV;;16" \
          "`filename 0nXn-25-gap13_Dbdt`;Rap gap 13 GeV;;13" \
          "`filename 0nXn-25-gap14_Dbdt`;Rap gap 14 GeV;;14" \
          "`filename 0nXn-25-gap15_Dbdt`;Rap gap 15 GeV;;15" \
          "`filename 0nXn-25-gap17_Dbdt`;Rap gap 17 GeV;;17" \
          "`filename 0nXn-25-gap18_Dbdt`;Rap gap 18 GeV;;18" \
          "`filename 0nXn-25-gap19_Dbdt`;Rap gap 19 GeV;;19" \
          "`filename 0nXn-25-gap20_Dbdt`;Rap gap 20 GeV;;20"

add_input '0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/bdt' '' 0 0 1 'BDT Working Point Shift' 1 \
          "`filename 0nXn-25_Dbdt`;Nominal BDT;;0" \
          "`filename 0nXn-25_Dbdt-bdtsM0p14`;BDT shift -0.14;;-0.14" \
          "`filename 0nXn-25_Dbdt-bdtsM0p12`;BDT shift -0.12;;-0.12" \
          "`filename 0nXn-25_Dbdt-bdtsM0p10`;BDT shift -0.10;;-0.10" \
          "`filename 0nXn-25_Dbdt-bdtsM0p08`;BDT shift -0.08;;-0.08" \
          "`filename 0nXn-25_Dbdt-bdtsM0p06`;BDT shift -0.06;;-0.06" \
          "`filename 0nXn-25_Dbdt-bdtsM0p04`;BDT shift -0.04;;-0.04" \
          "`filename 0nXn-25_Dbdt-bdtsM0p02`;BDT shift -0.02;;-0.02"

add_input '0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/fitting' '' 0 0 0 'Fitting function' 1 \
          "`filename 0nXn-25_Dbdt`;Nominal;;0" \
          'rootfiles/0nXn-25_Dbdt/'$BINNING_TAG'/xsec_fithist_2025PbPb_2025-SoftQCD_f-exp_deffcalc_2025-SoftQCD_evteffcalc_2025-SoftQCD_null.root;Expo background;;1' \
          'rootfiles/0nXn-25_Dbdt/'$BINNING_TAG'/xsec_fithist_2025PbPb_2025-SoftQCD_f-2gaus_deffcalc_2025-SoftQCD_evteffcalc_2025-SoftQCD_null.root;Double gaus signal;;1' \
          'rootfiles/0nXn-25_Dbdt/'$BINNING_TAG'/xsec_fithist_2025PbPb_2024-SoftQCD_deffcalc_2025-SoftQCD_evteffcalc_2025-SoftQCD_null.root;Pythia templates;;2'

##

make syst_comp.exe syst_summary.exe|| exit 1

for inputs in "${INPUTS[@]}" ; do
    IFS='^' read -r -a args <<< "$inputs"

    outputtag=${args[1]}
    tex=${args[6]}
    is_syst=${args[7]}
    # echo "${args[0]}"
    # echo "${args[1]}"
    set -x
    ./syst_comp.exe "${args[0]}" "${args[1]}" "${args[2]}" "${args[3]}" "${args[4]}" "${args[5]}" "${args[6]}" $SAVE_PNG
    set +x

    [[ $is_syst -eq 1 ]] && {
        systs=$systs',rootfiles/'$outputtag'.root;'$tex';' ;
    }
done
systs=${systs##,}
echo "$systs"
./syst_summary.exe "$systs" '0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/syst' ''


# add_input '0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/gap-extensive' '' 0 0 1 'Rapidity Gap Threshold' \
#           "`filename 0nXn-25_Dbdt`;Rap gap 16 GeV;;16" \
#           "`filename 0nXn-25-gap4_Dbdt`;Rap gap 4 GeV;;4" \
#           "`filename 0nXn-25-gap5_Dbdt`;Rap gap 5 GeV;;5" \
#           "`filename 0nXn-25-gap6_Dbdt`;Rap gap 6 GeV;;6" \
#           "`filename 0nXn-25-gap7_Dbdt`;Rap gap 7 GeV;;7" \
#           "`filename 0nXn-25-gap8_Dbdt`;Rap gap 8 GeV;;8" \
#           "`filename 0nXn-25-gap9_Dbdt`;Rap gap 9 GeV;;9" \
#           "`filename 0nXn-25-gap10_Dbdt`;Rap gap 10 GeV;;10" \
#           "`filename 0nXn-25-gap11_Dbdt`;Rap gap 11 GeV;;11" \
#           "`filename 0nXn-25-gap12_Dbdt`;Rap gap 12 GeV;;12" \
#           "`filename 0nXn-25-gap13_Dbdt`;Rap gap 13 GeV;;13" \
#           "`filename 0nXn-25-gap14_Dbdt`;Rap gap 14 GeV;;14" \
#           "`filename 0nXn-25-gap15_Dbdt`;Rap gap 15 GeV;;15" \
#           "`filename 0nXn-25-gap17_Dbdt`;Rap gap 17 GeV;;17" \
#           "`filename 0nXn-25-gap18_Dbdt`;Rap gap 18 GeV;;18" \
#           "`filename 0nXn-25-gap19_Dbdt`;Rap gap 19 GeV;;19" \
#           "`filename 0nXn-25-gap20_Dbdt`;Rap gap 20 GeV;;20" \
#           "`filename 0nXn-25-gap21_Dbdt`;Rap gap 21 GeV;;21" \
#           "`filename 0nXn-25-gap22_Dbdt`;Rap gap 22 GeV;;22" \
#           "`filename 0nXn-25-gap23_Dbdt`;Rap gap 23 GeV;;23"

# add_input '0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/bdt-extensive' '' 0 0 1 'BDT Working Point Shift' \
#           "`filename 0nXn-25_Dbdt`;Nominal BDT;;0" \
#           "`filename 0nXn-25_Dbdt-bdtsM0p20`;BDT shift -0.20;;-0.20" \
#           "`filename 0nXn-25_Dbdt-bdtsM0p18`;BDT shift -0.18;;-0.18" \
#           "`filename 0nXn-25_Dbdt-bdtsM0p16`;BDT shift -0.16;;-0.16" \
#           "`filename 0nXn-25_Dbdt-bdtsM0p14`;BDT shift -0.14;;-0.14" \
#           "`filename 0nXn-25_Dbdt-bdtsM0p12`;BDT shift -0.12;;-0.12" \
#           "`filename 0nXn-25_Dbdt-bdtsM0p10`;BDT shift -0.10;;-0.10" \
#           "`filename 0nXn-25_Dbdt-bdtsM0p08`;BDT shift -0.08;;-0.08" \
#           "`filename 0nXn-25_Dbdt-bdtsM0p06`;BDT shift -0.06;;-0.06" \
#           "`filename 0nXn-25_Dbdt-bdtsM0p04`;BDT shift -0.04;;-0.04" \
#           "`filename 0nXn-25_Dbdt-bdtsM0p02`;BDT shift -0.02;;-0.02" \
#           "`filename 0nXn-25_Dbdt-bdts0p02`;BDT shift 0.02;;0.02" \
#           "`filename 0nXn-25_Dbdt-bdts0p04`;BDT shift 0.04;;0.04" \
#           "`filename 0nXn-25_Dbdt-bdts0p06`;BDT shift 0.06;;0.06" \
#           "`filename 0nXn-25_Dbdt-bdts0p08`;BDT shift 0.08;;0.08" \
#           "`filename 0nXn-25_Dbdt-bdts0p10`;BDT shift 0.10;;0.10"
