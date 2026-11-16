#!/bin/bash

SAVE_PNG=0

INPUTS=()

add_input() {
    local output_name="$1"
    local label="$2"
    local title="${3-}"
    local include_patterns="${4-}"
    local legdx="${5-}"
    local legdy="${6-}"
    shift 6

    local IFS=,
    local input_files="$*"
    INPUTS+=("${input_files}^${output_name}^${label}^${title}^${include_patterns}^${legdx}^${legdy}")
}

# add_input '0nXn-vs-0nAn_gammaN' 'ZeroBias (Same D cut)' '' '.+mass_mc-.+,.+_width.+' 0 0 \
    #           '../rootfiles/gammaN-0nXn-25-ZB_Dbdt-gammaN-ZB/fithist_2025PbPb_2024-SoftQCD-BeamA_ycoarse.root;Xn0n;Xn0n' \
    #           '../rootfiles/gammaN-0nAn-25-ZB_Dbdt-gammaN-ZB/fithist_2025PbPb_2024-SoftQCD-BeamA_ycoarse.root;An0n;An0n'
# add_input '0nXn-vs-0nAn_Ngamma' 'ZeroBias (Same D cut)' '' '.+mass_mc-.+,.+_width.+' 0 0 \
    #           '../rootfiles/Ngamma-0nXn-25-ZB_Dbdt-Ngamma-ZB/fithist_2025PbPb_2024-SoftQCD-BeamB_ycoarse.root;0nXn;0nXn' \
    #           '../rootfiles/Ngamma-0nAn-25-ZB_Dbdt-Ngamma-ZB/fithist_2025PbPb_2024-SoftQCD-BeamB_ycoarse.root;0nAn;0nAn'
# add_input '0nAn-vs-0nXn_gammaN' 'Analysis D cut' '2025 PbPb (5.36 TeV)' '.+mass_mc-.+,.+_width.+' 0 0 \
    #           '../rootfiles/0nXn-gammaN-25_Dbdt-gammaN/b-ycoarse/fithist_2025PbPb_2025-SoftQCD-BeamA.root;Xn0n (ZDCOr);Xn0n' \
    #           '../rootfiles/0nAn-gammaN-25-ZB_Dbdt-gammaN/b-ycoarse/fithist_2025PbPb_2025-SoftQCD-BeamA.root;An0n (ZB);An0n'
# add_input '0nAn-vs-0nXn_0nXn_gammaN' '#gammaN,Same cuts' '2025 PbPb (5.36 TeV)' '.+yield+.,.+corr.+' 0 0 \
    #               '../rootfiles/0nXn-gammaN-25_Dbdt-gammaN/b-ycoarse/xsec_fithist_2025PbPb_2025-SoftQCD-BeamA_deffcalc_2025-SoftQCD-BeamA_evteffcalc_2025-SoftQCD-BeamA.root;Xn0n (ZDCOr);0nXn' \
    #               '../rootfiles/0nAn-gammaN-25-ZB_Dbdt-gammaN/b-ycoarse/xsec_fithist_2025PbPb_2025-SoftQCD-BeamA_deffcalc_2025-SoftQCD-BeamA_evteffcalc_2025-SoftQCD-BeamA.root;An0n (ZeroBias);0nAn'
# add_input 'ZDCOr-vs-ZB_0nXn_gammaN' 'Xn0n (#gammaN),Same cuts' '2025 PbPb (5.36 TeV)' '.+yield+.,.+corr.+' 0 0 \
    #               '../rootfiles/0nXn-gammaN-25_Dbdt-gammaN/b-ycoarse/xsec_fithist_2025PbPb_2025-SoftQCD-BeamA_deffcalc_2025-SoftQCD-BeamA_evteffcalc_2025-SoftQCD-BeamA.root;L1 ZDCOr;ZDCOr' \
    #               '../rootfiles/0nXn-gammaN-25-ZB_Dbdt-gammaN/b-ycoarse/xsec_fithist_2025PbPb_2025-SoftQCD-BeamA_deffcalc_2025-SoftQCD-BeamA_evteffcalc_2025-SoftQCD-BeamA.root;L1 ZeroBias;ZB'
# add_input 'ZDCOr-vs-ZB_0nXn_Ngamma' '0nXn (N#gamma),Same cuts' '2025 PbPb (5.36 TeV)' '.+yield.+,.+corr.+' 0 0 \
    #               '../rootfiles/0nXn-Ngamma-25_Dbdt-Ngamma/b-ycoarse/xsec_fithist_2025PbPb_2025-SoftQCD-BeamB_deffcalc_2025-SoftQCD-BeamB_evteffcalc_2025-SoftQCD-BeamB.root;L1 ZDCOr;ZDCOr' \
    #               '../rootfiles/0nXn-Ngamma-25-ZB_Dbdt-Ngamma/b-ycoarse/xsec_fithist_2025PbPb_2025-SoftQCD-BeamB_deffcalc_2025-SoftQCD-BeamB_evteffcalc_2025-SoftQCD-BeamB.root;L1 ZeroBias;ZB'
add_input 'eff_prompt-vs-inclusive_gammaN' 'Xn0n (#gammaN),Same analysis cut' 'P#scale[0.8]{YTHIA}8 #gammaN (5.36 TeV)' '.+_num.+,.+_den.+' 0 0 \
          'rootfiles/0nXn-gammaN-25_Dbdt-gammaN/b-default/deffcalc_2025-SoftQCD-BeamA.root;Prompt D;prompt' \
          'rootfiles/0nXn-gammaN-25_Dbdt-gammaN/b-default/deffcalc_2025-SoftQCD-BeamA-inclusive.root;Inclusive D;inclusive'
add_input 'eff_prompt-vs-inclusive_Ngamma' 'Xn0n (N#gamma),Same analysis cut' 'P#scale[0.8]{YTHIA}8 N#gamma (5.36 TeV)' '.+_num.+,.+_den.+' 0 0 \
          'rootfiles/0nXn-Ngamma-25_Dbdt-Ngamma/b-default/deffcalc_2025-SoftQCD-BeamB.root;Prompt D;prompt' \
          'rootfiles/0nXn-Ngamma-25_Dbdt-Ngamma/b-default/deffcalc_2025-SoftQCD-BeamB-inclusive.root;Inclusive D;inclusive'

make compare.exe || exit 1

for inputs in "${INPUTS[@]}" ; do
    IFS='^' read -r -a args <<< "$inputs"

    set -x
    ./compare.exe "${args[0]}" "${args[1]}" "${args[2]}" "${args[3]}" "${args[4]}" "${args[5]}" "${args[6]}" $SAVE_PNG
    set +x
done
