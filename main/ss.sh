#!/bin/bash

TDIRS=("gammaN" "Ngamma") ; TEXDIRS=("Xn0n (#gammaN)" "0nXn (N#gamma)")
FILENAMES=("xsec_fithist_2025PbPb_2024-SoftQCD-BeamA_deffcalc_2025-SoftQCD-BeamA_evteffcalc_2025-SoftQCD-BeamA_null.root"
           "xsec_fithist_2025PbPb_2024-SoftQCD-BeamB_deffcalc_2025-SoftQCD-BeamB_evteffcalc_2025-SoftQCD-BeamB_null.root")
isNgamma=${1:-0} ; tdir=${TDIRS[$isNgamma]} ; texdir=${TEXDIRS[$isNgamma]} ;

filename() {
    dirname=${1//__DIR__/${tdir}}
    fname='rootfiles/'"${dirname}"'/b-default/'${FILENAMES[$isNgamma]}
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
    shift 6

    local IFS=,
    local input_files="$*"
    INPUTS+=("${input_files}^${output_name}^${label}^${legdx}^${legdy}^${numscan}^${varscan}")
}

# add_input 'MC_digi-gt_0nXn_gammaN' 'P#scale[0.8]{YTHIA}8 #gammaN, 2025 RECO, #bf{Vary Digi GT}, Fixed Digi 15_1_2 + Era 2025' 'Private MC (2024 GEN-SIM)' '.+num.+,.+den.+' -0.1 0 \
    #           'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-pri-BeamA-Reco25-Digi-cmssw1512_gt141X_era25.root;GT-141X;' \
    #           'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-pri-BeamA-Reco25-Digi-cmssw1512_gt151X_era25.root;GT-151X;' \
    #           'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-pri-BeamA-Reco25-Digi-cmssw1512-gt151X-era25-SiPixelQuality24.root;GT-151X (24 SiPixelQuality);'
add_input '0nXn-'${tdir}'-25_Dbdt-'${tdir}'/b-default/variations/gap-extensive' "$texdir" -0.1 0 1 'Rapidity Gap Threshold' \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__`;Rap gap 16 GeV;;16" \
          "`filename 0nXn-__DIR__-25-gap4_Dbdt-__DIR__`;Rap gap 4 GeV;;4" \
          "`filename 0nXn-__DIR__-25-gap5_Dbdt-__DIR__`;Rap gap 5 GeV;;5" \
          "`filename 0nXn-__DIR__-25-gap6_Dbdt-__DIR__`;Rap gap 6 GeV;;6" \
          "`filename 0nXn-__DIR__-25-gap7_Dbdt-__DIR__`;Rap gap 7 GeV;;7" \
          "`filename 0nXn-__DIR__-25-gap8_Dbdt-__DIR__`;Rap gap 8 GeV;;8" \
          "`filename 0nXn-__DIR__-25-gap9_Dbdt-__DIR__`;Rap gap 9 GeV;;9" \
          "`filename 0nXn-__DIR__-25-gap10_Dbdt-__DIR__`;Rap gap 10 GeV;;10" \
          "`filename 0nXn-__DIR__-25-gap11_Dbdt-__DIR__`;Rap gap 11 GeV;;11" \
          "`filename 0nXn-__DIR__-25-gap12_Dbdt-__DIR__`;Rap gap 12 GeV;;12" \
          "`filename 0nXn-__DIR__-25-gap13_Dbdt-__DIR__`;Rap gap 13 GeV;;13" \
          "`filename 0nXn-__DIR__-25-gap14_Dbdt-__DIR__`;Rap gap 14 GeV;;14" \
          "`filename 0nXn-__DIR__-25-gap15_Dbdt-__DIR__`;Rap gap 15 GeV;;15" \
          "`filename 0nXn-__DIR__-25-gap17_Dbdt-__DIR__`;Rap gap 17 GeV;;17" \
          "`filename 0nXn-__DIR__-25-gap18_Dbdt-__DIR__`;Rap gap 18 GeV;;18" \
          "`filename 0nXn-__DIR__-25-gap19_Dbdt-__DIR__`;Rap gap 19 GeV;;19" \
          "`filename 0nXn-__DIR__-25-gap20_Dbdt-__DIR__`;Rap gap 20 GeV;;20" \
          "`filename 0nXn-__DIR__-25-gap21_Dbdt-__DIR__`;Rap gap 21 GeV;;21" \
          "`filename 0nXn-__DIR__-25-gap22_Dbdt-__DIR__`;Rap gap 22 GeV;;22" \
          "`filename 0nXn-__DIR__-25-gap23_Dbdt-__DIR__`;Rap gap 23 GeV;;23"

add_input '0nXn-'${tdir}'-25_Dbdt-'${tdir}'/b-default/variations/bdt' "$texdir" -0.1 0 1 'BDT Working Point Shift' \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__`;Nominal BDT;;0" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdtsM0p20`;BDT shift -0.20;;-0.20" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdtsM0p18`;BDT shift -0.18;;-0.18" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdtsM0p16`;BDT shift -0.16;;-0.16" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdtsM0p14`;BDT shift -0.14;;-0.14" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdtsM0p12`;BDT shift -0.12;;-0.12" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdtsM0p10`;BDT shift -0.10;;-0.10" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdtsM0p08`;BDT shift -0.08;;-0.08" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdtsM0p06`;BDT shift -0.06;;-0.06" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdtsM0p04`;BDT shift -0.04;;-0.04" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdtsM0p02`;BDT shift -0.02;;-0.02" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdts0p02`;BDT shift 0.02;;0.02" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdts0p04`;BDT shift 0.04;;0.04" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdts0p06`;BDT shift 0.06;;0.06" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdts0p08`;BDT shift 0.08;;0.08" \
          "`filename 0nXn-__DIR__-25_Dbdt-__DIR__-bdts0p10`;BDT shift 0.10;;0.10"

make syst_comp.exe || exit 1

for inputs in "${INPUTS[@]}" ; do
    IFS='^' read -r -a args <<< "$inputs"

    # echo "${args[0]}"
    # echo "${args[1]}"
    set -x
    ./syst_comp.exe "${args[0]}" "${args[1]}" "${args[2]}" "${args[3]}" "${args[4]}" "${args[5]}" "${args[6]}"
    set +x
done
