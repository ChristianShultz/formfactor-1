#!/bin/tcsh

# run_dist_src_harom-jitptx

set srcmod="srca.mod  srcf.mod srcb.mod"
set srcxmlout="distillution.ini.xml.out"

/bin/rm -f $srcmod $srcxmlout


# run_prop_dist_chroma-jitptx

set propmod="propa.mod propb.mod propf.mod"
set propdistxmlout="prop_dist.ini.xml.out"

/bin/rm -f $propmod $propdistxmlout


# run_peram_dist_harom-jitptx

set perammod="perama.mod peramb.mod peramf.mod"
set peramdistxmlout="peram_dist.ini.xml.out"

/bin/rm -f $perammod $peramdistxmlout
