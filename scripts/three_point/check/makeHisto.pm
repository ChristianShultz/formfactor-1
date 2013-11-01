#!/urs/bin/perl
#

package makeHisto;

use POSIX;

use strict; 

sub new
{
  my $invocant = shift; 
  my $class = ref($invocant) || $invocant; 

  my $self= {
      MAX => undef, 
      MIN => undef, 
      BIN => undef, 
      DATA => [],
      @_, 
  };

  return bless $self , $class; 
}


sub max{
  my $self = shift; 
  if (@_) {$self->{MAX} = shift;}
  return $self->{MAX};
}

sub min{
  my $self = shift; 
  if (@_) {$self->{MIN} = shift;}
  return $self->{MIN};
}

sub bin{
  my $self = shift; 
  if (@_) {$self->{BIN} = shift;}
  return $self->{BIN};
}

sub data{
  my $self = shift; 
  if (@_) {$self->{DATA} = shift;}
  return $self->{DATA};
}

sub update_min{
  my $self = shift; 
  my $min = shift; 
  $self->min($min) unless $self->min(); 
  $self->min($min) unless $min > $self->min(); 
}

sub update_max{
  my $self = shift; 
  my $max = shift; 

  $self->max($max) unless $self->max();
  $self->max($max) unless $max < $self->max(); 
}

sub push_data{
  my $self = shift; 
  my $d = shift; 
  my @dd = @{ $self->data() }; 
  
  push @dd , $d; 
  $self->data(\@dd); 
}

sub add_data{
 my $self = shift; 
 my $d = shift; 
  
 $self->update_max($d);
 $self->update_min($d);  
 $self->push_data($d);
}

sub histo
{
  my $self = shift; 
  my %histo = (); 

  foreach my $datum ( @{ $self->data() } ) 
  { 
    $histo{ceil( $datum  / $self->bin() ) }++; 
  }
  
  my $min;
  my $max;

  while ( my ($k,$v) = each (%histo) )
  {
    $min = $k if ! defined($min) || $k < $min; 
    $max = $k if ! defined($max) || $k > $max; 
  }

  print "range [ " . $self->min .  " ,  "  . $self->max() . " ]\n";

  print $min . " " . $max . "\n";

  open OUT , ">" , "histo.out" ;
  for (my $i = $min; $i <= $max; $i++)
  {
    my $h = 0; 
    $h = $histo{$i} if exists $histo{$i}; 
    print OUT $i * $self->bin() . " " . $h . "\n";
  }
  close OUT;

  open GNU , ">" , "histo.gp";
  print GNU "plot \'histo.out\' u 1:2 smooth freq w boxes \n";
  close GNU;

  system ( "gnuplot -p histo.gp" ) ;  
}



1;
