#!/usr/bin/perl


package findPhase;
use strict; 

sub new
{
  my $invocant = shift; 
  my $class = ref($invocant) || $invocant; 

  my $self = {

    MF => 0.21616,
    MI => 0.1483,
    XI => 3.45,
    LS => 16.,
    ELEM => -1,
    ID => "uninit",
    STEM => "lattice_mat_elem_",
    MAT => "RhoPi_",
    JMU => 400,
    RADRL => 1000,
    RADIM => 1000,
    LATRL => 1000, 
    LATIM => 1000,
    PSRC => [],
    PSNK => [],
    HSRC => 1000,
    HSNK => 1000,
    PHIRAD => "uninit",
    PHILAT => "uninit",    
    @_, 
  };

  return bless $self, $class; 
}



sub mf {
  my $self = shift; 
  if(@_) {$self->{MF} = shift;}
  return $self->{MF};
}

sub mi {
  my $self = shift; 
  if(@_) {$self->{MI} = shift;}
  return $self->{MI};
}

sub xi {
  my $self = shift; 
  if(@_) {$self->{XI} = shift;}
  return $self->{XI};
}

sub ls {
  my $self = shift; 
  if(@_) {$self->{LS} = shift;}
  return $self->{LS};
}

sub elem {
  my $self = shift; 
  if(@_) {$self->{ELEM} = shift;}
  return $self->{ELEM};
}

sub id {
  my $self = shift; 
  if(@_) {$self->{ID} = shift;}
  return $self->{ID};
}

sub stem {
  my $self = shift; 
  if(@_) {$self->{STEM} = shift;}
  return $self->{STEM};
}

sub matrix_element_stem {
  my $self = shift;
  if(@_) {$self->{MAT} = shift;}
  return $self->{MAT};
}

sub jmu{
  my $self = shift;
  if(@_) {$self->{JMU} = shift;}
  return $self->{JMU};

}

sub rad_real {
  my $self = shift; 
  if(@_) {$self->{RADRL} = shift;}
  return $self->{RADRL};
}


sub rad_imag {
  my $self = shift; 
  if(@_) {$self->{RADIM} = shift;}
  return $self->{RADIM};
}


sub lat_real {
  my $self = shift; 
  if(@_) {$self->{LATRL} = shift;}
  return $self->{LATRL};
}


sub lat_imag {
  my $self = shift; 
  if(@_) {$self->{LATIM} = shift;}
  return $self->{LATIM};
}


sub lat_elem{
  my $self = shift; 
  my $foo = $self->stem();
  my $bar =$self->elem();
  return $foo . $bar . ".jack"; 
}


sub psrc{
  my $self = shift;
  if(@_) {$self->{PSRC} = shift;}
  return $self->{PSRC};
}

sub psnk{
  my $self = shift;
  if(@_) {$self->{PSNK} = shift;}
  return $self->{PSNK};
}

sub hsrc{
  my $self = shift;
  if(@_) {$self->{HSRC} = shift;}
  return $self->{HSRC};
}


sub hsnk{
  my $self = shift;
  if(@_) {$self->{HSNK} = shift;}
  return $self->{HSNK};
}

sub phirad{
   my $self = shift;
  if(@_) {$self->{PHIRAD} = shift;}
  return $self->{PHIRAD};
}

sub philat{
   my $self = shift;
  if(@_) {$self->{PHILAT} = shift;}
  return $self->{PHILAT};
}



sub rad_str
{
  my $self = shift;

  my $id = $self->matrix_element_stem . $self->hsnk() . "_" . $self->hsrc();

  my $snk = " " . $self->mf() . " " . $self->psnk_str();
  my $src = " " . $self->mi() . " " . $self->psrc_str(); 
  my $fac = " " . $self->momentum_factor(); 

  return $id . $snk . $src . $fac . " " . $self->jmu; 
}


sub psnk_str
{
  my $self = shift; 
  return $self->p_str( $self->psnk() );
}


sub psrc_str
{
  my $self = shift; 
  return $self->p_str( $self->psrc() );
}

sub momentum_factor
{
  my $self = shift; 
  return 2.*3.1415 / $self->xi() / $self->ls();
}


sub p_str
{
  my $self = shift;
  my $ref = shift; 
  my @mom = @{ $ref };
  my $ret = "";

  foreach my $p (@mom)
  {
    $ret = $ret . $p . " "; 
  }

  return $ret; 
}

sub q_str
{ 
  my $self = shift; 
  my @psnk = @{ $self->psnk() };
  my @psrc = @{ $self->psrc() };
  
  my $str = "";

  for my $i (0..2)
  {
    my $q = $psnk[$i] - $psrc[$i]; 
    $str = $str . " " . $q;   
  }
  return $str; 
}

sub parse_id
{
  my $self = shift; 
  my $str = $self->id(); 

  #  print $self->id(); 

  my ($snk,$ins,$src) = split(/\./,$str); 

  # print $src . " " . $snk . "\n";

  # source params
  my @src = @{ $self->parse_state($src) };  
  $self->psrc( $src[0] );
  $self->hsrc( $src[1] );

  # sink params
  my @snk = @{ $self->parse_state($snk) };
  $self->psnk( $snk[0] );
  $self->hsnk( $snk[1] ); 

  # ins params
  my ($trash,$mu) = split(/_/,$ins); 
  $self->jmu($mu); 
}



sub parse_state
{
  my $self = shift;
  my $str = shift;

  # split the name 
  my($op,$sp,$h,$p,$I,$t) = split(/,/,$str); 

  # return array
  my @params; 

  # get the momentum
  my $tmp;
  my $ct = 0; 
  my $minus = "-";
  $p = substr($p,1,length($p)); 
  my @m; 

  for my $i (0..2)
  {
    $tmp = substr($p,$ct,1);

    if($tmp =~ m/$minus/)
    {
      $tmp = substr($p,$ct,2);

      $ct++; 
    }

    $ct++;
    push @m, $tmp; 
  }

  push @params , \@m;

# get the helicity
$h = substr($h,1,length($h));
push @params, $h; 


return \@params; 
}


sub out_str
{
  my $self = shift;

  my $str = $self->id() . "   " ;
  $str = $str . $self->philat() . " - " . $self->phirad() . " = ";
  
  my $diff =  $self->philat() - $self->phirad();
  
  my $diffbyPi = sprintf("%.3f", $diff / 3.14159 );

  $str = $str . $diff . "  divide by pi " . $diffbyPi;
  return $str; 
}

sub elem_out_str
{
  my $self = shift;

  my $str = $self->stem() . $self->elem() .  "   " . $self->id() . "   ";
  $str = $str . $self->philat() . " - " . $self->phirad() . " = ";
  
  my $diff =  $self->philat() - $self->phirad();
  
  my $diffbyPi = sprintf("%.3f", $diff / 3.14159 );

  $str = $str . $diff . "  divide by pi " . $diffbyPi;
  return $str; 

}


sub phase_diff 
{
  my $self = shift; 
  return ( $self->philat() - $self->phirad() );
}

1;  
