#!/bin/tcsh -x

# make the sources for all quark lines
~/qcd/git/harom/scalar-Nd3/mainprogs/main/harom -i ./make_source_dist.ini.xml

# invert and save solutions
/dist/scidac/chroma/master/scalar-9q/bin/chroma -i ./distillution.ini.xml

# generate perambulators from some of the solutions
~/qcd/git/harom/scalar-Nd3/mainprogs/main/harom -i ./peram_dist.ini.xml

# Run redstar using harom to produce the hadron_node sdbs
./run_harom.csh pi_pi_3pt.xml


