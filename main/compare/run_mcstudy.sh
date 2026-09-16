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

# add_input 'MC_25-vs-24_0nXn_gammaN' 'Xn0n (#gammaN), Same analysis D cut' '' '.+num.+,.+den.+' 0 0 \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-SoftQCD-BeamA.root;2024 MC;2024-SoftQCD' \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_SoftQCD-BeamA.root;2025 MC;2025-SoftQCD'
# add_input 'MC_checkforest_0nXn_gammaN' 'Xn0n (#gammaN),Same analysis D cut' '2025 Official MC' '.+num.+,.+den.+' 0 0 \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_SoftQCD-BeamA.root;0904 Forest;2025-SoftQCD' \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2025-SoftQCD-BeamA-0328forest.root;0328 Forest;2025-SoftQCD-0328forest'
# add_input 'MC_checkskim_0nXn_gammaN' 'Xn0n (#gammaN),Same analysis D cut' '2024 Official MC' '.+num.+,.+den.+' 0 0 \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-SoftQCD-BeamA.root;0714 skim;2024-SoftQCD' \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-SoftQCD-BeamA-rpskim.root;0906 skim;2024-SoftQCD-rpskim'

add_input 'MC_reco_0nXn_gammaN' 'Xn0n (#gammaN), 2024 Digi, #bf{Vary RECO}, Fixed RECO GT 141X' 'Private MC (2024 GEN-SIM)' '.+num.+,.+den.+' 0 0 \
          'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-SoftQCD-BeamA.root;Official 24 MC;' \
          'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-pri-BeamA-Digi24-Reco-cmssw1417-gt141X-era24.root;14_1_7 Era-24' \
          'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-pri-BeamA-Digi24-Reco-cmssw1512_gt141X_era24.root;15_1_2 Era-24' \
          'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-pri-BeamA-Digi24-Reco-cmssw1512_gt141X_era25.root;15_1_2 Era-25'
add_input 'MC_reco-gt_0nXn_gammaN' 'Xn0n (#gammaN), 2024 Digi, #bf{Vary RECO GT}, Fixed RECO 15_1_2 + Era 2025' 'Private MC (2024 GEN-SIM)' '.+num.+,.+den.+' 0 0 \
          'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-pri-BeamA-Digi24-Reco-cmssw1512_gt141X_era25.root;GT-141X;' \
          'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-pri-BeamA-Digi24-Reco-cmssw1512_gt151X_era25.root;GT-151X;'
add_input 'MC_digi_0nXn_gammaN' 'P#scale[0.8]{YTHIA}8 #gammaN, 2025 RECO, #bf{Vary Digi}, Fixed Digi GT 141X' 'Private MC (2024 GEN-SIM)' '.+num.+,.+den.+' 0 0 \
          'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-pri-BeamA-Reco25-Digi-cmssw1417-gt141X-era24.root;14_1_7 Era-24;' \
          'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-pri-BeamA-Reco25-Digi-cmssw1512_gt141X_era25.root;15_1_2 Era-25;'
add_input 'MC_digi-gt_0nXn_gammaN' 'P#scale[0.8]{YTHIA}8 #gammaN, 2025 RECO, #bf{Vary Digi GT}, Fixed Digi 15_1_2 + Era 2025' 'Private MC (2024 GEN-SIM)' '.+num.+,.+den.+' 0 0 \
          'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-pri-BeamA-Reco25-Digi-cmssw1512_gt141X_era25.root;GT-141X;' \
          'rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-pri-BeamA-Reco25-Digi-cmssw1512_gt151X_era25.root;GT-151X;'
# add_input 'MC_sum_0nXn_gammaN' 'P#scale[0.8]{YTHIA}8 #gammaN' 'Private MC (2024 GEN-SIM)' '.+num.+,.+den.+' -0.15 0 \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-private-BeamA-Digi24-Reco-cmssw1417-gt141X-era24.root;24 Digi + 24 RECO;Private-Digi24-Reco-cmssw1417-gt141X-era24' \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-private-BeamA-Digi24-Reco-cmssw1512-gt141X-era24.root;24 Digi + 25 RECO/GT-141X;Private-Digi24-Reco-cmssw1512-gt141X-era24' \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-private-BeamA-Digi24-Reco-cmssw1512-gt151X-era25.root;24 Digi + 25 RECO;Private-Reco25-Digi-cmssw1417-gt141X-era24' \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-private-BeamA-Digi25-gt141X-Reco-cmssw1512-gt151X-era25.root;25 Digi/GT-141X + 25 RECO;Private-Reco25-Digi-cmssw1512-gt141X-era25' \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-private-BeamA-Digi25-Reco-cmssw1512-gt151X-era25.root;2025 Digi + 25 RECO;Private-Reco25-Digi-cmssw1512-gt151X-era25'
# add_input 'MC_sumbrief_0nXn_gammaN' 'P#scale[0.8]{YTHIA}8 #gammaN, #bf{Vary GT}' 'Private MC (2024 GEN-SIM)' '.+num.+,.+den.+' -0.15 0 \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-SoftQCD-BeamA.root;Official 2024;2024-SoftQCD' \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-private-BeamA-Digi24-Reco-cmssw1512-gt141X-era24.root;Digi 141X + RECO 141X;Private-Digi24-Reco-cmssw1512-gt141X-era24' \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-private-BeamA-Digi25-gt141X-Reco-cmssw1512-gt151X-era25.root;Digi 141X + RECO 151X;Private-Reco25-Digi-cmssw1512-gt141X-era25' \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-private-BeamA-Digi25-Reco-cmssw1512-gt151X-era25.root;Digi 151X + RECO 151X;Private-Reco25-Digi-cmssw1512-gt151X-era25'
# add_input 'MC_gen_0nXn_gammaN' 'P#scale[0.8]{YTHIA}8 #gammaN, 2025 Digi, 2025 RECO' '' '.+num.+,.+den.+' -0.1 0 \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-private-BeamA-Digi25-Reco-cmssw1512-gt151X-era25.root;Private (2024 GEN-SIM);Private-Reco25-Digi-cmssw1512-gt151X-era25' \
    #           '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_SoftQCD-BeamA.root;Official (2025 GEN);2025-SoftQCD'

make compare.exe || exit 1

for inputs in "${INPUTS[@]}" ; do
    IFS='^' read -r -a args <<< "$inputs"

    set -x
    ./compare.exe "${args[0]}" "${args[1]}" "${args[2]}" "${args[3]}" "${args[4]}" "${args[5]}" "${args[6]}"
    set +x
done
