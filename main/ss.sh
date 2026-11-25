#!/bin/bash
ver='v2'

SAVE_PNG=1

[[ $# -gt 0 ]] || {
    echo 'usage: ./ss.sh [idx_binning] ([comp syst] [collect syst])'
    exit
}

BINNINGS=(
    'b-default;fprompt-25MC'
    'b-ptdiff;null'
)

ibin=${1}
IFS=';' ; binnings=(${BINNINGS[$ibin]}) ; BINNING_TAG=${binnings[0]} ; FP_TAG=${binnings[1]} ; unset IFS;
echo 'BINNING_TAG : '$BINNING_TAG
echo '     FP_TAG : '$FP_TAG

FILENAME="xsec_fithist_2025PbPb_2025-SoftQCD_deffcalc_2025-SoftQCD_evteffcalc_2025-SoftQCD_"$FP_TAG

# FILENAME="xsec_fithist_2025PbPb_2025-SoftQCD_deffcalc_2025-SoftQCD_evteffcalc_2025-SoftQCD_null"
filename() {
    # dirname=${1//__DIR__/${tdir}}
    fname='rootfiles/'$ver'/'${1}'/'$BINNING_TAG'/'$FILENAME'.root'
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

add_input $ver'/0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/zdcscan' '' 0 0 1 'ZDC Threshold' 1 \
          "`filename 0nXn-25_Dbdt`;ZDC+/- 1100/1000 GeV;;1050" \
          "`filename 0nXn-25-zdc800_Dbdt`;ZDC 800 GeV;;800" \
          "`filename 0nXn-25-zdc900_Dbdt`;ZDC 900 GeV;;900" \
          "`filename 0nXn-25-zdc1000_Dbdt`;ZDC 1000 GeV;;1000" \
          "`filename 0nXn-25-zdc1100_Dbdt`;ZDC 1100 GeV;;1100" \
          "`filename 0nXn-25-zdc1200_Dbdt`;ZDC 1200 GeV;;1200" \
          "`filename 0nXn-25-zdc1300_Dbdt`;ZDC 1300 GeV;;1300" \
          "`filename 0nXn-25-zdc1300_Dbdt`;ZDC 1400 GeV;;1400" \
          "`filename 0nXn-25-zdc1300_Dbdt`;ZDC 1500 GeV;;1500"

add_input $ver'/0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/gapscan' '' 0 0 1 'Rapidity Gap Threshold' 1 \
          "`filename 0nXn-25_Dbdt`;Nominal;;18" \
          "`filename 0nXn-25-gap15_Dbdt`;Rap gap 15 GeV;;15" \
          "`filename 0nXn-25-gap16_Dbdt`;Rap gap 16 GeV;;16" \
          "`filename 0nXn-25-gap17_Dbdt`;Rap gap 17 GeV;;17" \
          "`filename 0nXn-25-gap19_Dbdt`;Rap gap 19 GeV;;19" \
          "`filename 0nXn-25-gap20_Dbdt`;Rap gap 20 GeV;;20" \
          "`filename 0nXn-25-gap21_Dbdt`;Rap gap 21 GeV;;21"

add_input $ver'/0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/bdtscan' '' 0 0 1 'BDT Working Point Shift' 1 \
          "`filename 0nXn-25_Dbdt`;Nominal BDT;;0" \
          "`filename 0nXn-25_Dbdt-bdtsM0p14`;BDT shift -0.14;;-0.14" \
          "`filename 0nXn-25_Dbdt-bdtsM0p12`;BDT shift -0.12;;-0.12" \
          "`filename 0nXn-25_Dbdt-bdtsM0p10`;BDT shift -0.10;;-0.10" \
          "`filename 0nXn-25_Dbdt-bdtsM0p08`;BDT shift -0.08;;-0.08" \
          "`filename 0nXn-25_Dbdt-bdtsM0p06`;BDT shift -0.06;;-0.06" \
          "`filename 0nXn-25_Dbdt-bdtsM0p04`;BDT shift -0.04;;-0.04" \
          "`filename 0nXn-25_Dbdt-bdtsM0p02`;BDT shift -0.02;;-0.02" \
          "`filename 0nXn-25_Dbdt-bdts0p02`;BDT shift +0.02;;0.02" \
          "`filename 0nXn-25_Dbdt-bdts0p04`;BDT shift +0.04;;0.04"
          # "`filename 0nXn-25_Dbdt-bdts0p06`;BDT shift +0.06;;0.06" \
          # "`filename 0nXn-25_Dbdt-bdts0p08`;BDT shift +0.08;;0.08" \
          # "`filename 0nXn-25_Dbdt-bdts0p10`;BDT shift +0.10;;0.10"


add_input $ver'/0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/fitting' '' 0 0 0 'Invariant mass fit' 1 \
          "`filename 0nXn-25_Dbdt`;Nominal;;0" \
          'rootfiles/v2/0nXn-25_Dbdt/'$BINNING_TAG'/xsec_fithist_2025PbPb_2025-SoftQCD_f-exp_deffcalc_2025-SoftQCD_evteffcalc_2025-SoftQCD_'$FP_TAG'.root;Expo background;;1' \
          'rootfiles/v2/0nXn-25_Dbdt/'$BINNING_TAG'/xsec_fithist_2025PbPb_2024-SoftQCD_deffcalc_2025-SoftQCD_evteffcalc_2025-SoftQCD_'$FP_TAG'.root;Pythia templates;;2' \
          'rootfiles/v2/0nXn-25_Dbdt/'$BINNING_TAG'/xsec_fithist_2025PbPb_2025-SoftQCD-tmc_deffcalc_2025-SoftQCD_evteffcalc_2025-SoftQCD_'$FP_TAG'.root;Fixed to MC signal widths;;3' \
          'rootfiles/v2/0nXn-25_Dbdt/'$BINNING_TAG'/xsec_fithist_2025PbPb_2025-SoftQCD-talpha2_deffcalc_2025-SoftQCD_evteffcalc_2025-SoftQCD_'$FP_TAG'.root;Alternative reso scale;;4'

# add_input $ver'/0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/gap-extensive' '' 0 0 1 'Rapidity Gap Threshold (GeV)' 0 \
#                       "`filename 0nXn-25_Dbdt`;Nominal;;18" \
#                       "`filename 0nXn-25-gap4_Dbdt`;Rap gap 4 GeV;;4" \
#                       "`filename 0nXn-25-gap5_Dbdt`;Rap gap 5 GeV;;5" \
#                       "`filename 0nXn-25-gap6_Dbdt`;Rap gap 6 GeV;;6" \
#                       "`filename 0nXn-25-gap7_Dbdt`;Rap gap 7 GeV;;7" \
#                       "`filename 0nXn-25-gap8_Dbdt`;Rap gap 8 GeV;;8" \
#                       "`filename 0nXn-25-gap9_Dbdt`;Rap gap 9 GeV;;9" \
#                       "`filename 0nXn-25-gap10_Dbdt`;Rap gap 10 GeV;;10" \
#                       "`filename 0nXn-25-gap11_Dbdt`;Rap gap 11 GeV;;11" \
#                       "`filename 0nXn-25-gap12_Dbdt`;Rap gap 12 GeV;;12" \
#                       "`filename 0nXn-25-gap13_Dbdt`;Rap gap 13 GeV;;13" \
#                       "`filename 0nXn-25-gap14_Dbdt`;Rap gap 14 GeV;;14" \
#                       "`filename 0nXn-25-gap15_Dbdt`;Rap gap 15 GeV;;15" \
#                       "`filename 0nXn-25-gap16_Dbdt`;Rap gap 16 GeV;;16" \
#                       "`filename 0nXn-25-gap17_Dbdt`;Rap gap 17 GeV;;17" \
#                       "`filename 0nXn-25-gap19_Dbdt`;Rap gap 19 GeV;;19" \
#                       "`filename 0nXn-25-gap20_Dbdt`;Rap gap 20 GeV;;20" \
#                       "`filename 0nXn-25-gap21_Dbdt`;Rap gap 21 GeV;;21" \
#                       "`filename 0nXn-25-gap22_Dbdt`;Rap gap 22 GeV;;22" \
#                       "`filename 0nXn-25-gap23_Dbdt`;Rap gap 23 GeV;;23" \
#                       "`filename 0nXn-25-gap24_Dbdt`;Rap gap 24 GeV;;24"

# add_input $ver'/0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/zdc-extensive' '' 0 0 1 'ZDC Threshold (GeV)' 0 \
    #           "`filename 0nXn-25_Dbdt`;ZDC+/- 1100/1000 GeV;;1050" \
    #           "`filename 0nXn-25-zdc600_Dbdt`;ZDC 600 GeV;;600" \
    #           "`filename 0nXn-25-zdc700_Dbdt`;ZDC 700 GeV;;700" \
    #           "`filename 0nXn-25-zdc800_Dbdt`;ZDC 800 GeV;;800" \
    #           "`filename 0nXn-25-zdc900_Dbdt`;ZDC 900 GeV;;900" \
    #           "`filename 0nXn-25-zdc1000_Dbdt`;ZDC 1000 GeV;;1000" \
    #           "`filename 0nXn-25-zdc1100_Dbdt`;ZDC 1100 GeV;;1100" \
    #           "`filename 0nXn-25-zdc1200_Dbdt`;ZDC 1200 GeV;;1200" \
    #           "`filename 0nXn-25-zdc1300_Dbdt`;ZDC 1300 GeV;;1300" \
    #           "`filename 0nXn-25-zdc1400_Dbdt`;ZDC 1400 GeV;;1400" \
    #           "`filename 0nXn-25-zdc1500_Dbdt`;ZDC 1500 GeV;;1500"

# add_input $ver'/0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/gap-range' '' 0 0 0 'Rapidity Gap Threshold' 0 \
#               "`filename 0nXn-25_Dbdt`;HF E_{max} in 3 < |#eta| < 5.2;;" \
#               "`filename 0nXn-25-gapeta5_Dbdt`;HF E_{max} in 3 < |#eta| < 5;;" \

##

make syst_comp.exe syst_summary.exe syst_bdt.exe || exit 1

for inputs in "${INPUTS[@]}" ; do
    IFS='^' read -r -a args <<< "$inputs"

    outputtag=${args[1]}
    tex=${args[6]}
    is_syst=${args[7]}
    # echo "${args[0]}"
    # echo "${args[1]}"

    [[ ${2:-0} -eq 1 || ${2:-0} -eq 3 ]] && {
        ./syst_comp.exe "${args[0]}" "${args[1]}" "${args[2]}" "${args[3]}" "${args[4]}" "${args[5]}" "${args[6]}" $SAVE_PNG
    }
    # echo "${args[1]}"
    [[ "${args[1]}" == *bdtscan* ]] && {
        [[ ${2:-0} -eq 2 || ${2:-0} -eq 3 ]] && {
            inputname_bdt="${args[1]}"
            outputname_bdt=${inputname_bdt/bdtscan/bdtexp}
            inputname_bdt='rootfiles/'$inputname_bdt'.root;BDT Working Point Shift;'
            ./syst_bdt.exe "$inputname_bdt" $outputname_bdt $SAVE_PNG
        }
    }

    [[ $is_syst -eq 1 ]] && {
        outputtag_syst=$outputtag
        [[ $outputtag == *bdtscan ]] && outputtag_syst=${outputtag/bdtscan/bdtexp} ;
        systs=$systs',rootfiles/'$outputtag_syst'.root;'$tex';' ;
    }
done
systs=${systs##,}
echo -e " \033[2m>> [ syst summary ]\033[0m"
echo "$systs"
[[ ${3:-0} -eq 1 && x$systs != x ]] && ./syst_summary.exe "$systs" "`filename 0nXn-25_Dbdt`" $ver'/0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/syst' $SAVE_PNG
echo -e " \033[2m<<\033[0m"

# add_input $ver'/0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/bdt-extensive' '' 0 0 1 'BDT Working Point Shift' \
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

# add_input $ver'/0nXn-25_Dbdt/'$BINNING_TAG'/'$FILENAME'/bdt_weight' '' 0 0 1 'BDT Weight' 1 \
    #           "`filename 0nXn-25_Dbdt`;No weight;;0" \
    #           'rootfiles/0nXn-25_Dbdt/'$BINNING_TAG'/xsec_fithist_2025PbPb_2025-SoftQCD_deffcalc_2025-SoftQCD_weight_evteffcalc_2025-SoftQCD_'$FP_TAG'.root;BDT weight;;1'

