#!/usr/bin/perl


# a parameter package for optimized operators with a 
# tie in for radmat operator state overlaps 
#
# 

package OPparams;
use strict;

sub new
{
  my $invocant = shift; 
  my $class = ref($invocant) || $invocant; 
  my $self = {
    T0 => undef,
    STATE => undef,
    TZ => undef, 
    REP => undef,
    REPSTEM => undef, 
    MOM => undef, 
    PID => undef, 
    NCFG => undef,
    ENSEMBLE => undef,  
    RECONDIR => undef,
    @_,
  }; ## self 

  return bless $self , $class; 
}

#
# WITH an arg sets the variable
#
# WITHOUT and arg fetches the variable
#

sub t0 {
  my $self = shift; 
  if (@_) {$self->{T0} = shift}
  return $self->{T0}
}

sub state {
  my $self = shift; 
  if (@_) {$self->{STATE} = shift}
  return $self->{STATE}
}

sub tz {
  my $self = shift; 
  if (@_) {$self->{TZ} = shift}
  return $self->{TZ}
}

sub irrep {
  my $self = shift; 
  if (@_) {$self->{REP} = shift}
  return $self->{REP}
}

sub irrep_stem {
  my $self = shift; 
  if (@_) {$self->{REPSTEM} = shift}
  return $self->{REPSTEM}
}

sub mom {
  my $self = shift; 
  if (@_) {$self->{MOM} = shift}
  return $self->{MOM}
}

sub pid {
  my $self = shift; 
  if (@_) {$self->{PID} = shift}
  return $self->{PID}
}


sub ncfg {
  my $self = shift; 
  if (@_) {$self->{NCFG} = shift}
  return $self->{NCFG}
}
sub ensemble {
  my $self = shift; 
  if (@_) {$self->{ENSEMBLE} = shift}
  return $self->{ENSEMBLE}
}

sub recon_dir {
  my $self = shift; 
  if (@_) {$self->{RECONDIR} = shift}
  return $self->{RECONDIR};
}


#
#
# parameterless methods
#
#

sub op_name
{
  my $self = shift; 

  my $s = $self->pid() . "_" . $self->mom() . "_" . $self->irrep_stem(); 

  return $s; 
}





1;
