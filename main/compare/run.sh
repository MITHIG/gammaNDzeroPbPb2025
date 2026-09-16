#!/bin/bash

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
# add_input 'ZDCOr-vs-ZB_0nXn_gammaN' 'Xn0n (#gammaN),Analysis D cut' '' '' 0 0 \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/xsec_fithist_2025PbPb_2024-SoftQCD-BeamA_ycoarse_calceff_2024-SoftQCD-BeamA_ycoarse_null_null.root;L1 ZDCOr;ZDCOr' \
    #           '../rootfiles/gammaN-0nXn-25-ZB_Dbdt-gammaN/xsec_fithist_2025PbPb_2024-SoftQCD-BeamA_ycoarse_calceff_2024-SoftQCD-BeamA_ycoarse_null_null.root;L1 ZeroBias;ZB'
# add_input 'ZDCOr-vs-ZB_0nXn_Ngamma' '0nXn (N#gamma),Analysis D cut' '' '' 0 0 \
    #           '../rootfiles/Ngamma-0nXn-25_Dbdt-Ngamma/xsec_fithist_2025PbPb_2024-SoftQCD-BeamB_ycoarse_calceff_2024-SoftQCD-BeamB_ycoarse_null_null.root;L1 ZDCOr;ZDCOr' \
    #           '../rootfiles/Ngamma-0nXn-25-ZB_Dbdt-Ngamma/xsec_fithist_2025PbPb_2024-SoftQCD-BeamB_ycoarse_calceff_2024-SoftQCD-BeamB_ycoarse_null_null.root;L1 ZeroBias;ZB'

make compare.exe || exit 1

for inputs in "${INPUTS[@]}" ; do
    IFS='^' read -r -a args <<< "$inputs"

    set -x
    ./compare.exe "${args[0]}" "${args[1]}" "${args[2]}" "${args[3]}" "${args[4]}" "${args[5]}" "${args[6]}"
    set +x
done
