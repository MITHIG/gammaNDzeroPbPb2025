#!/bin/bash

Vars=(
    'HFEMax_eta5;(3 < |#eta| < 5);HFEMax_eta5'
    'HFEMax_forest;(3 < |#eta| < 5.2);HFEMax_forest'
)

Inputs=(
    '/eos/cms/store/group/phys_heavyions/wangj/Forest2025PbPb/Dzero_260212-hfle_HiForest_260218_HIEmptyBX_HIRun2025A_PromptReco_v1.root;2025 EmptyBX;emptybx-2025'
    # '/eos/cms/store/group/phys_heavyions/wangj/Forest2023PbPb/Dzero_260212-hfle_HiForest_260218_HIEmptyBX_HIRun2023A_PromptReco_v2.root;2023 EmptyBX;emptybx-2023'
)

make percentile.exe draw.exe || exit 1

for inputstr in "${Inputs[@]}" ; do
    IFS=';' ; pinput=($inputstr) ; unset IFS ; itag_input=${pinput[2]}
    for varstr in "${Vars[@]}" ; do
        IFS=';' ; pvar=($varstr) ; unset IFS ; itag_var=${pvar[2]}

        itag=$itag_input'_'$itag_var

        [[ ${1:-0} -eq 1 ]] && ./percentile.exe "$inputstr" "$varstr" $itag &

        [[ ${2:-0} -eq 1 ]] && ./draw.exe 'rootfiles/'$itag'.root' $itag
    done
done
wait
