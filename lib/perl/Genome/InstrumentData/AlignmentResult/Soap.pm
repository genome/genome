package Genome::InstrumentData::AlignmentResult::Soap;

use strict;
use warnings;
use File::Basename;

use Genome;

class Genome::InstrumentData::AlignmentResult::Soap {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'soap', is_param=>1 },
    ],
};

sub required_arch_os { 'x86_64' }

sub required_rusage { 
    "-R 'select[tmp>50000 && mem>10000] span[hosts=1] rusage[tmp=50000, mem=10000]' -M 10000000 -n 4";
}

sub _run_aligner {
    my $self = shift;

    $self->fatal_message('Generating new soap results is no longer supported.');
}

sub aligner_params_for_sam_header {
    my $self = shift;
    return 'soap ' . $self->aligner_params;
}

sub fillmd_for_sam {
    # soap appears to already fill in the MD string... so...
    return 0;
}

#                                      _______
###################################   /______/|
### Taken from the Soap package ###   |     | |
### Author: lh3   Version 0.1.1 ###   |trash| |
###################################   |_____|/

sub __mating {
  my ($s1, $s2) = @_;
  my $isize = 0;
  if ($s1->[2] ne '*' && $s1->[2] eq $s2->[2]) { # then calculate $isize
    my $x1 = ($s1->[1] & 0x10)? $s1->[3] + length($s1->[9]) : $s1->[3];
    my $x2 = ($s2->[1] & 0x10)? $s2->[3] + length($s2->[9]) : $s2->[3];
    $isize = $x2 - $x1;
  }
  # update mate coordinate
  if ($s2->[2] ne '*') {
    @$s1[6..8] = (($s2->[2] eq $s1->[2])? "=" : $s2->[2], $s2->[3], $isize);
    $s1->[1] |= 0x20 if ($s2->[1] & 0x10);
  } else {
    $s1->[1] |= 0x8;
  }
  if ($s1->[2] ne '*') {
    @$s2[6..8] = (($s1->[2] eq $s2->[2])? "=" : $s1->[2], $s1->[3], -$isize);
    $s2->[1] |= 0x20 if ($s1->[1] & 0x10);
  } else {
    $s2->[1] |= 0x8;
  }
}

sub _soap2sam {
  my $input_file = shift;
  my $output_file = shift;
  my $is_paired = shift;
  my $in_fh = IO::File->new($input_file);
  my $out_fh = IO::File->new(">>$output_file");
  # core loop
  my @s1 = ();
  my @s2 = ();
  my ($s_last, $s_curr) = (\@s1, \@s2);
  while (<$in_fh>) {
    s/[\177-\377]|[\000-\010]|[\012-\040]//g;
    next if (&__soap2sam_aux($_, $s_curr, $is_paired) < 0);
    if (@$s_last != 0 && $s_last->[0] eq $s_curr->[0]) {
      &__mating($s_last, $s_curr);
      $out_fh->print(join("\t", @$s_last), "\n");
      $out_fh->print(join("\t", @$s_curr), "\n");
      @$s_last = (); @$s_curr = ();
    } else {
      $out_fh->print(join("\t", @$s_last), "\n") if (@$s_last != 0);
      my $s = $s_last; $s_last = $s_curr; $s_curr = $s;
    }
  }
  $out_fh->print(join("\t", @$s_last), "\n") if (@$s_last != 0);
  $in_fh->close();
  $out_fh->close();
}

sub __soap2sam_aux {
  my ($line, $s, $is_paired) = @_;
  chomp($line);
  my @t = split(/\s+/, $line);
  return -1 if (@t < 9 || $line =~ /^\s/ || !$t[0]);
  @$s = ();
  # fix SOAP-2.1.x bugs
  @t = @t[0..2,4..$#t] unless ($t[3] =~ /^\d+$/);
  # read name
  $s->[0] = $t[0];
  $s->[0] =~ s/\/[12]$//g;
  # initial flag (will be updated later)
  $s->[1] = 0;
  $s->[1] |= 1 | 1<<($t[4] eq 'a'? 6 : 7);
  $s->[1] |= 2 if ($is_paired);
  # read & quality
  $s->[9] = $t[1];
  $s->[10] = (length($t[2]) > length($t[1]))? substr($t[2], 0, length($t[1])) : $t[2];
  # cigar
  $s->[5] = length($s->[9]) . "M";
  # coor
  $s->[2] = $t[7]; $s->[3] = $t[8];
  $s->[1] |= 0x10 if ($t[6] eq '-');
  # mapQ
  $s->[4] = $t[3] == 1? 30 : 0;
  # mate coordinate
  $s->[6] = '*'; $s->[7] = $s->[8] = 0;
  # aux
  push(@$s, "NM:i:$t[9]");
  my $md = '';
  if ($t[9]) {
    my @x;
    for (10 .. $#t) {
      push(@x, sprintf("%.3d,$1", $2)) if ($t[$_] =~ /^([ACGT])->(\d+)/i);
    }
    @x = sort(@x);
    my $a = 0;
    for (@x) {
      my ($y, $z) = split(",");
      $md .= (int($y)-$a) . $z;
      $a += $y - $a + 1;
    }
    $md .= length($t[1]) - $a;
  } else {
    $md = length($t[1]);
  }
  push(@$s, "MD:Z:$md");
  return 0;
}
