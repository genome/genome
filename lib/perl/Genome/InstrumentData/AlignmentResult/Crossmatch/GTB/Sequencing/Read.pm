package Genome::InstrumentData::AlignmentResult::Crossmatch::GTB::Sequencing::Read;

### Adapted from code by Nancy Hansen <nhansen@mail.nih.gov>


############################################################

=head1 NAME

Read.pm - A Perl module to contain a single Sequencing read.

=head1 DESCRIPTION

  A Perl module for handling Sequencing reads: quality metrics,
  trimming, etc.

=head1 DATE

 July 19, 2007

=head1 AUTHORS

 Nancy Hansen <nhansen@mail.nih.gov>

=head1 PUBLIC METHODS

=cut

############################################################
use strict;
use warnings;

use Carp;

our $VERSION  = '0.01';

use vars qw( );

###########################################################

=over 4

=item new()

  This method creates a Read object.

  Input:  -name - the name of the read (optional) 
          -seq - the base sequence of the read
          -quals - a reference to an array of phred-style
             quality scores, one for each base OR
          -qual_string - fastq-style string of base 
             quality scores for the read, one character
             per base (no requirement for 0 qual first)
             Either quals or qual_string can be specified
             in the new method, but not both.
          -qual_offset - option to pass a character to 
             act as the offset when calculating quality
             scores from fastq quality string characters
             (default is '!')

  Output: New Read object

=cut

###########################################################
sub new {
    my $class = shift;
    my %params = @_;
    my $name = $params{-name};
    my $seq = $params{-seq};
    my $quals = $params{-quals};
    my $qual_string = $params{-qual_string};
    my $qual_offset = $params{-qual_offset} || '!';

    if ($quals && $qual_string)
    {
        croak "Please specify -quals or -qual_string in $class constructor, but not both!\n";
    }
 
    my $self = {name => $name,
                seq => $seq,
                quals => $quals,
                qual_string => $qual_string,
                qual_offset => $qual_offset };

    bless $self, $class;
    return $self;
}

###########################################################

=item name()

  This method gets or sets the value of name for this
  Read object.  (optional)

  Input: Optional argument sets value.
  Output: Value of "name" (scalar string).

=cut

###########################################################
sub name {
    my $self  = shift;

    if (defined (my $new_name = shift))
    {
        $self->{name} = $new_name;
    }

    return $self->{name};

} ## end name

###########################################################

=item seq()

  This method gets or sets the value of seq for this
  Read object.

  Input: Optional argument sets value.
  Output: Value of "seq" (scalar string).

=cut

###########################################################
sub seq {
    my $self  = shift;

    if (defined (my $new_seq = shift))
    {
        $self->{seq} = $new_seq;
    }

    return $self->{seq};

} ## end seq

###########################################################

=item seq_length()

  This method gets or sets the value of seq_length for this
  Read object.  If it is undefined, it calculates the length
  as the length of the string value of the "seq" field.

  Input: Optional argument sets value.
  Output: Value of "seq_length" (scalar number).

=cut

###########################################################
sub seq_length {
    my $self  = shift;

    if (defined (my $new_seq_length = shift))
    {
        $self->{seq_length} = $new_seq_length;
    }

    if (!defined ($self->{seq_length}))
    {
        $self->{seq_length} = length $self->seq();
    }

    return $self->{seq_length};

} ## end seq_length

###########################################################

=item quals()

  This method gets or sets the value of quals, which
  is a reference to an array of quality scores for the 
  read, one for each base.

  Input: None.  If no value, object will attempt to 
     populate it by converting the "qual_string" field.
  Output: Value of "quals" (reference to an array of numbers)

=cut

###########################################################
sub quals {
    my $self  = shift;

    if (!defined ($self->{quals}) && defined ($self->{qual_string})) 
    {
        $self->{quals} = [];
        foreach my $char (split //, $self->{qual_string})
        {
            my $qual_offset_char = $self->{qual_offset};
            my $qual = ord($char) - ord($qual_offset_char);
            push @{$self->{quals}}, $qual;
        }
    }

    return $self->{quals};

} ## end quals

###########################################################

=item qual_string()

  This method gets or sets the value of qual_string, which
  is a fastq-style quality string, with one character per
  base.

  Input: None.  If qual_string has no value, the object will
     attempt to populate it by constructing it from the
     quals array.
  Output: Value of "qual_string" (scalar string).

=cut

###########################################################
sub qual_string {
    my $self  = shift;

    if (!defined ($self->{qual_string}) && defined ($self->{quals})) 
    {
        $self->{qual_string} = '';
        foreach my $qual (@{$self->{quals}})
        {
            my $qual_offset_char = $self->{qual_offset};
            my $qual_char = chr(ord($qual_offset_char) + $qual);
            $self->{qual_string} .= $qual_char;
        }
    }
    return $self->{qual_string};

} ## end qual_string

###########################################################

=item numerical_qual_string()

  This method constructs a numerical quality string by 
  concatenating the numerical quality scores.

  Input: None.
  Output: Value of "numerical_qual_string" (scalar string).

=cut

###########################################################
sub numerical_qual_string {
    my $self  = shift;

    if (!defined ($self->{numerical_qual_string}))
    {
        my $ra_quals = $self->quals();
        $self->{numerical_qual_string} = ($ra_quals) ? 
                  join ' ', @{$ra_quals} : undef;

    }
    return $self->{numerical_qual_string};

} ## end numerical_qual_string

###########################################################

=item pig_trim_endpoints()

  This method trims the object according to the "PIG"
  algorithm ("Peak Interval Greedy").

  Input: -kmer_size => size of kmers to use as smoothing windows 
                (default 14)
         -max_expected_errors => maximum acceptable expected
                errors per base (default 0.01)
         -min_window => minimum length of read window to return.
                (default 0)
  Output: Array containing the start and endpoints, or (0,-1) if
         no suitable window could be found.

=cut

###########################################################
sub pig_trim_endpoints {
    my $self  = shift;
    my %params = @_;

    my $k = $params{-kmer_size} || 12;
    my $max_errors = $params{-max_expected_errors} || .1;
    my $min_window = $params{-min_window} || 20;

    my $MILLION = 1e6 + 1e-6;
    my $maxEperM = $MILLION * $max_errors;

    # here I've inserted Paul's implementation of fpktrim, with mods to fit module:
    my @ePerMs = map { int(10**(-$_/10) * $MILLION) } @{ $self->quals() };

    my @FAIL = (0, -1); # zero-length window
    my @ePerM_k = ();
    my $minEperM_k = 9e9;
    my $minEperM_kid = 0;
    my $totEperM = 0;
    my $currEperM_k = 0;

    my $len = @ePerMs;
    #print "len $len\n";
    my $i = 0;
    foreach my $ePerM (@ePerMs) {
        $totEperM += $ePerM;
        #print "$ePerM $i $totEperM\n";
        $currEperM_k += $ePerM;
        if ($i >= $k-1) {
            $currEperM_k -= $ePerMs[$i - $k] 
                unless $i < $k;
            $ePerM_k[$i] = $currEperM_k;
            if ($currEperM_k < $minEperM_k) {
                $minEperM_k = $currEperM_k;
                $minEperM_kid = $i;
            }
        }
        $i++;
    }
    $ePerMs[$len] = 9e9; # convenient to have boundary value after last base position
    if ($minEperM_k <= $maxEperM) { # i.e, if we can start with a kmer of few enough errors
        # then grow until it would have too many
        my $currEperM = $minEperM_k;
        my $winstop  = $minEperM_kid;
        my $winstart = $winstop + 1 - $k;
        for (;;) {
            if ($winstart > 0 and $ePerMs[$winstart - 1] <= $ePerMs[$winstop + 1]) {
                my $newstart = $winstart - 1;
                my $newEperM = $currEperM + $ePerMs[$newstart];
                if ($newEperM < $maxEperM) {
                    $winstart = $newstart;
                    $currEperM = $newEperM;
                    next;
                }
                else {
                    last;
                }
            }
            else {
                my $newstop = $winstop + 1;
                my $newEperM = $currEperM + $ePerMs[$newstop];
                if ($newEperM < $maxEperM) {
                    $winstop = $newstop;
                    $currEperM = $newEperM;
                    next;
                }
                else {
                    last;
                }
            }
        }
        # We've expanded as far as we can, are we long enough?
        if ($winstop + 1 - $winstart >= $min_window) {
            return ($winstart + 1, $winstop + 1);
        }
    }
    return @FAIL;

} ## end pig_trim_endpoints

###########################################################

=item fpk_trim_endpoints()

  This method trims the object according to the "FPK"
  algorithm ("First Passing KMERs").

  Input: -kmer_size => size of kmers to use as smoothing windows 
                (default 14)
         -max_expected_errors => maximum acceptable expected
                errors per base (default 0.01)
         -min_window => minimum length of read window to return.
                (default 0)
  Output: Array containing the start and endpoints, or (0,-1) if
         no suitable window could be found.

=cut

###########################################################
sub fpk_trim_endpoints {
    my $self  = shift;
    my %params = @_;

    my $k = $params{-kmer_size} || 12;
    my $max_errors = $params{-max_expected_errors} || .1; # per base
    my $min_window = $params{-min_window} || 20;

    my $MILLION = 1e6 + 1e-6;
    my $maxEperM = $MILLION * $max_errors;

    # here I've inserted Paul's implementation of fpktrim, with mods to fit module:
    my @ePerMs = map { int(10**(-$_/10) * $MILLION) } @{ $self->quals() };

    my @FAIL = (0, -1); # zero-length window
    my @ePerM_k = ();
    my $totEperM = 0;
    my $currEperM_k = 0;

    my $len = @ePerMs;
    my $i = 0;
    foreach my $ePerM (@ePerMs) {
        $totEperM += $ePerM;
        $currEperM_k += $ePerM;
        if ($i >= $k-1) {
            $currEperM_k -= $ePerMs[$i - $k] 
                unless $i < $k;
            $ePerM_k[$i] = $currEperM_k;
        }
        $i++;
    }
    $ePerM_k[$len] = 9e9; # convenient to have boundary value after last base position
    my ($winstart, $winstop) = @FAIL;

    for ($i = $k - 1; $i < $len; $i++) { # find first passing kmer
        if ($ePerM_k[$i] <= $maxEperM) {
            $winstart = $i + 1 - $k;
            $winstop = $i;
            last;
        }
    }
    until ($ePerM_k[$i] > $maxEperM) {
        $winstop = $i++;
    }
    # We've expanded as far as we can, are we long enough?
    if ($winstop + 1 - $winstart >= $min_window) {
        return ($winstart + 1, $winstop + 1);
    }
    else {
        return @FAIL;
    }

} ## end fpk_trim_endpoints

###########################################################

=item qtrim_endpoints()

  This method trims the object according to Paul's "qtrim"
  algorithm (which finds the longest size stretch of quals
  >= to some minimum quality "min_q"--in case of tie, first
  stretch in read wins.

  Input: -min_q => minimum quality score (phred style, 
               default 10).
  Output: Array containing the start and endpoints, or (0,-1) if
         no suitable window could be found.

=cut

###########################################################
sub qtrim_endpoints {
    my $self  = shift;
    my %params = @_;
    my $min_q = defined($params{-min_q}) ? $params{-min_q} : 10;

    my ($best, $beststart) = (0, 0);

    my @quals = @{$self->quals()};
    my $length = @quals;

    my $i = 0;
    my $count = 0;
    foreach my $q (@quals) {
        if ($q >= $min_q) {
            $count++;
        }
        else {
            ($best, $beststart) = ($count, $i - $count) if ($count > $best);
            $count = 0;
        }
        $i++;
    }
    ($best, $beststart) = ($count, $i - $count) if ($count > $best);

    return ($beststart, $beststart + $best - 1);

} ## end qtrim_endpoints

###########################################################

=item subseq()

  This method returns a new GTB::Sequencing::Read object
  consisting only of the bases and qualities from the
  specified starting point (-start) to the specified end
  (-end), both 1-based.  The new object's name field will
  be the string from the original read with ".$start-$end"
  appended.

  Input: -start => $start (default 1)
         -end => $end (default last base of read)
  Output: GTB::Sequencing::Read object

=cut

###########################################################
sub subseq {
    my $self  = shift;
    my %params = @_;

    my $start = $params{-start} || 1;
    my $end = $params{-end} || $self->seq_length();

    if (($start !~ /^\d+$/) || ($start < 1) || ($start > $self->seq_length()))
    {
        croak "Starting point $start must be numerical and within the bounds of the read!\n";
    }
    
    if (($end !~ /^\d+$/) || ($end < 1) || ($end > $self->seq_length()))
    {
        croak "Ending point $end must be numerical and within the bounds of the read!\n";
    }

    my $offset = $start - 1;
    my $new_length = $end - $start + 1;

    my $name = $self->name();
    $name .= ".$start-$end";

    my $seq = $self->seq();
    $seq = substr($seq, $offset, $new_length);
    my $qual_string = $self->qual_string();
    $qual_string = substr($qual_string, $offset, $new_length);

    my $class = ref $self;
    return $class->new(-name => $name,
                       -seq => $seq,
                       -qual_string => $qual_string);

}

###########################################################

=item comp()

  This method returns a new GTB::Sequencing::Read object
  that is the reverse complement of the object passed.

  Input: None.
  Output: GTB::Sequencing::Read object

=cut

###########################################################
sub comp {
    my $self  = shift;

    my $name = $self->name();
    $name .= ".comp";

    my $seq = $self->seq();
    $seq = reverse $seq;
    $seq =~ tr/ATGCatgc/TACGtacg/;

    my $qual_string = $self->qual_string();
    $qual_string = reverse $qual_string;

    my $class = ref $self;
    return $class->new(-name => $name,
                       -seq => $seq,
                       -qual_string => $qual_string);

} # end comp

###########################################################

=item qual_mask()

  This method returns a new GTB::Sequencing::Read object
  for which bases in the sequence whose qualities are
  less than a specified minimum quality (-min_qual) are
  replaced with 'N' (by default, or -replacement_char if
  specified in the parameters). The new object's name field 
  will be the string from the original read with 
  ".qual_mask_<min_qual>" appended.  Quality scores will
  be converted to 0 for these bases unless the 
  "-no_qual_conversion" parameter is true

  Input: -min_qual => $min (default 20)
         -replacement_char => $char (default 'N')
         -no_qual_conversion => true will suppress conversion
             or qualities less than min_qual to 0 (default false).
  Output: GTB::Sequencing::Read object

=cut

###########################################################
sub qual_mask {
    my $self  = shift;
    my %params = @_;

    my $min_qual = $params{-min_qual} || 20;
    my $replacement_char = $params{-replacement_char} || 'N'; 
    my $no_qual_conversion = $params{-no_qual_conversion} || 0; 

    my $name = $self->name();
    $name .= ".qual_mask_$min_qual";

    my $seq = $self->seq();
    my @quals = @{$self->quals()};

    # move through read, creating new read and quality scores:

    my $new_seq = '';
    my @new_quals = ();
    while (my $last_base = chop $seq)
    {
        my $last_qual = pop @quals;
        my $new_base = ($last_qual < $min_qual) ? $replacement_char : $last_base;
        $new_seq = $new_base.$new_seq;
        my $new_qual = ($no_qual_conversion) ? $last_qual : 
                                ($last_qual < $min_qual) ? 0 : $last_qual;
        unshift @new_quals, $new_qual;
    }

    my $class = ref $self;
    return $class->new(-name => $name,
                       -seq => $new_seq,
                       -quals => \@new_quals);


} ## end qual_mask
    
###########################################################

=item read_reads_from_sqr_file()

  This method reads in the contents of the specified sqr
  file, and outputs a reference to an array of Read 
  objects.

  Input: Class name, -sqr_file => path to .sqr file
  Output: Reference to an array of GTB::Sequencing::Read
     objects.

=cut

###########################################################
sub read_reads_from_sqr_file {
    my $class  = shift;
    my %params = @_;

    my $sqr_file = $params{-sqr_file}
        or croak "Must pass a valid path to a .sqr file as -sqr_file to read_reads_from_sqr_file!\n";

    my $ra_return_seqs = [];
    open SQR, $sqr_file
        or croak "Couldn\'t open file $sqr_file: $!\n";
    while (<SQR>)
    {
        chomp;
        my ($seq, $qual, $name, $rev_qual, $rev_seq) = split /\t/, $_;
        push @{$ra_return_seqs}, $class->new(
                                        -name => $name,
                                        -seq => $seq,
                                        -qual_string => $qual);
    }

    return $ra_return_seqs;
}

###########################################################

1;
__END__

=back
