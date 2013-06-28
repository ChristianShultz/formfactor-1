#!/usr/bin/perl



use strict;
use warnings; 


die("usage: $0 <time/directory> <space/directory> <FF #> ") unless ($#ARGV == 2);


my $time_base = $ARGV[0]; 
my $space_base = $ARGV[1];
my $ff = $ARGV[2];

my $time = "${time_base}/t_ins_fits";
my $space = "${space_base}/t_ins_fits";
my $ff_parse = "F_${ff}";

die("$time does not exist") unless (-d $time);
die("$space does not exist") unless (-d $space); 



my @timeQ2 =  `ls $time | grep $ff_parse | grep jack | grep fit`;
my @spaceQ2 = `ls $space | grep $ff_parse | grep jack | grep fit`;

chomp @timeQ2;
chomp @spaceQ2; 

my $rn = int(rand(1000)); 
my $tt = make_gnu_file("time",$time,$rn,@timeQ2);
my $ss = make_gnu_file("space",$space,$rn,@spaceQ2);
my $rat = make_rat_file("rat",$time,$space,$rn,\@timeQ2,\@spaceQ2); 

# print $tt . " " . $ss . "\n";


my $plot = "tmp/plot_${rn}.gp";
my $rat_plot = "tmp/rat_plot_${rn}.gp";

open PLOT , ">" . $plot; 
open RATPLOT , ">" . $rat_plot; 

print PLOT "plot \'$tt\' using 1:2:3 w yerr title \'temporal\' , \\\n";
print PLOT "\'$ss\' using 1:2:3 w yerr title \'spatal\' \n";
print RATPLOT "plot \'$rat\' using 1:2:3 w yerr title \'ratio\' \n";

system("gnuplot -persist $plot");
system("gnuplot -persist $rat_plot"); 


close PLOT;
close RATPLOT;

sub make_gnu_file
{
  my ($outname,$stem,$rn,@files) = @_; 

  mkdir "tmp"; 

  my $fout = "tmp/${outname}_${rn}.dat";

  open OUT , ">" , $fout; 

  foreach my $file (@files)
  {

    my ($a,$b,$q,@trash) = split (/_/,$file);

    $q =~ s/p/./g;

    my $f2 = "${stem}/${file}";


    my $calc = `calc $f2 | awk \'{ print \$2 \" \" \$3 }\' `;

#    print "$file $q $f2 " . "\n" . "$calc" . "\n";

    print OUT "$q $calc"; 

  }

  close OUT; 

  return $fout;  
}


sub make_rat_file
{
  my ($outname,$tstem,$sstem,$rn,$tref,$sref) = @_;

  mkdir "tmp";

  my $fout = "/tmp/${outname}_${rn}.dat";

  my @one = @{$tref};
  my @two = @{$sref};
  my %one = map{$_=>1} @one;
  my %two = map{$_=>2} @two; 
  my @intersection = grep($one{$_},@two); 

  open OUT , ">" , $fout; 


  foreach my $elem (@intersection)
  {
    my $tf = $tstem . "/${elem}";
    my $sf = $sstem . "/${elem}";
    my ($a,$b,$q,@trash) = split (/_/,$elem);
    $q =~ s/p/./g;

    my $calc = `calcbc \' $tf / $sf \' |  awk \'{ print \$2 \" \" \$3 }\' `;

    print OUT "$q $calc"; 

  }

  close OUT; 

  return $fout;
}


