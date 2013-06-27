#!/usr/bin/perl
#
# Run meson_hadron_node

$, = ' ';               # set output field separator
##$\ = "\n";              # set output record separator
$\ = "";              # set output record separator

die "Usage: $0 <proj list file> <output xml file>\n" unless @ARGV == 2;

$list_file = shift(@ARGV);
$xml_file = shift(@ARGV);

die "Input parameter file= $list_file  does not exist\n" unless -f $list_file;

open(FILE, "< $list_file") || die "Error opening $list_file\n";
open(XML, "> $xml_file");

printf XML <<EOF;
<?xml version="1.0"?>
<ProjectedOps>
EOF

$key = "";
$prev_key = "";

while($_ = <FILE>)
{
  chomp $_;

  ($key, $op, $weight) = split(' ', $_);
  if ($key ne $prev_key)
  {
    if ($prev_key ne "")
    {
print XML <<EOF;
    </Weights>
  </elem>
EOF
    }

print XML <<EOF;
  <elem>
    <Key>
      <id>$key</id> <smear></smear>
    </Key>
    <Weights>
EOF
    
   $prev_key = $key;
  }

  print XML "      <elem><id>$op</id> <weight>$weight</weight></elem>\n";
}

printf XML "    </Weights>\n  </elem>\n</ProjectedOps>\n";

exit(0);

