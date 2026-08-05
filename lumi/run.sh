#!/bin/bash

JSONS=(
    jsons/Collisions25HI_399465_400007_Golden_400059_400426_DCSOnlyTkPc.json:HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v15,HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v16,HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v15,HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v16
    # /eos/user/c/cmsdqm/www/CAF/certification/Collisions25HI/Cert_Collisions2025_HI_399465_400426_Golden.json:HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v15,HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v16,HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v15,HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v16
)

for jj in "${JSONS[@]}" ; do
    IFS=':' ; pars=($jj) ; unset IFS ;
    fjson=${pars[0]}
    IFS=',' ; hlts=(${pars[1]}) ; unset IFS ;
    for ll in ${hlts[@]} ; do
        echo $fjson" --> "$ll
        jtag=${fjson##*/} ; jtag=${jtag%%.*} ;
        outdir=results/$jtag
        mkdir -p $outdir
        brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json -u /nb -i $fjson --hltpath $ll > $outdir"/"$ll".txt"
    done
    # brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json -u /nb -i $jj
done

