#!/bin/bash

INPUTS=(
    # '../rootfiles/gammaN-0nXn-25-ZB_Dbdt-gammaN-ZB/fithist_2025PbPb_2024-SoftQCD-BeamA_ycoarse.root;Xn0n;Xn0n,../rootfiles/gammaN-0nAn-25-ZB_Dbdt-gammaN-ZB/fithist_2025PbPb_2024-SoftQCD-BeamA_ycoarse.root;An0n;An0n'^'0nXn-vs-0nAn_gammaN'^'ZeroBias (Same D cut)'^'.+mass_mc-.+,.+_width.+'
    # '../rootfiles/Ngamma-0nXn-25-ZB_Dbdt-Ngamma-ZB/fithist_2025PbPb_2024-SoftQCD-BeamB_ycoarse.root;0nXn;0nXn,../rootfiles/Ngamma-0nAn-25-ZB_Dbdt-Ngamma-ZB/fithist_2025PbPb_2024-SoftQCD-BeamB_ycoarse.root;0nAn;0nAn'^'0nXn-vs-0nAn_Ngamma'^'ZeroBias (Same D cut)'^'.+mass_mc-.+,.+_width.+'
    # '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/xsec_fithist_2025PbPb_2024-SoftQCD-BeamA_ycoarse_calceff_2024-SoftQCD-BeamA_ycoarse_null_null.root;L1 ZDCOr;ZDCOr,../rootfiles/gammaN-0nXn-25-ZB_Dbdt-gammaN/xsec_fithist_2025PbPb_2024-SoftQCD-BeamA_ycoarse_calceff_2024-SoftQCD-BeamA_ycoarse_null_null.root;L1 ZeroBias;ZB'^'ZDCOr-vs-ZB_0nXn_gammaN'^'Xn0n (#gammaN),Analysis D cut'^
    # '../rootfiles/Ngamma-0nXn-25_Dbdt-Ngamma/xsec_fithist_2025PbPb_2024-SoftQCD-BeamB_ycoarse_calceff_2024-SoftQCD-BeamB_ycoarse_null_null.root;L1 ZDCOr;ZDCOr,../rootfiles/Ngamma-0nXn-25-ZB_Dbdt-Ngamma/xsec_fithist_2025PbPb_2024-SoftQCD-BeamB_ycoarse_calceff_2024-SoftQCD-BeamB_ycoarse_null_null.root;L1 ZeroBias;ZB'^'ZDCOr-vs-ZB_0nXn_Ngamma'^'0nXn (N#gamma),Analysis D cut'^
    # '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-SoftQCD-BeamA.root;2024 MC;2024-SoftQCD,../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_SoftQCD-BeamA.root;2025 MC;2025-SoftQCD'^'MC25-vs-24_0nXn_gammaN'^'Xn0n (#gammaN),Same analysis D cut'^'.+num.+,.+den.+'
    # '../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_SoftQCD-BeamA.root;25 MC (0904Forest);2025-SoftQCD,../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-SoftQCD-BeamA.root;24 MC (0328Forest);2024-SoftQCD,../rootfiles/gammaN-0nXn-25_Dbdt-gammaN/b-default/deffcalc_2024-SoftQCD-BeamA-0328forest.root;25 MC (0328Forest);2025-SoftQCD-0328forest'^'MC_checkforest_0nXn_gammaN'^'Xn0n (#gammaN),Same analysis D cut'^'.+num.+,.+den.+'
)

make compare.exe || exit 1

for inputs in "${INPUTS[@]}" ; do
    IFS='^' ; args=($inputs) ; unset IFS ;
    # echo $inputs
    # for ii in "${arguments[@]}" ; do echo "    "$ii ; done ;
    
    set -x
    ./compare.exe "${args[0]}" ${args[1]} "${args[2]}" "${args[3]}"
    set +x
done
