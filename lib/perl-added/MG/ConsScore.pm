package MG::ConsScore;

use strict;
use warnings;
use Carp;

sub new
{
    my $class = shift;
    my %args = @_;
    my $self = { _location => undef,
                 _file_handles => undef,
                };
    if(exists($args{'-location'}))
    {
        # Temporary until disks come back up and we are sure that nothing is defaulting to temp202
        if ($args{'-location'} eq '/gscmnt/temp202/info/medseq/josborne/conservation/b36/fixed-records/') {
            $args{'-location'} = '/gscmnt/sata835/info/medseq/model_data/2741951221/v36-build93636924/ucsc_conservation/';
        }
        $self->{_location} = $args{'-location'};
    }
    bless($self, $class);
    return $self;
}

=head1 NAME

MG::ConsScore

=head1 SYNOPSIS

  use MG::ConsScore;
  my $c = new MG::ConsScore;
  my $array_ref = $c->get_scores(2,[10000,10001,10002]);

=cut


sub get_scores
{
    use IO::File;
#    use Set::IntSpan;
    my $self = shift;
    my $chromosome = shift;
    # should check if $coordinates is an array ref or a string
    # with integer ranges (to use with Set::IntSpan?)
    my $coordinates = shift;
    my $coordref = [ ];
    if(ref($coordinates)  eq 'ARRAY')
    {
        # do stuff for array ref
        $coordref = $coordinates;        
    }
    elsif($coordinates =~ /\d+([,\-]\d+){0,}/x )
    {
        #turn this into an array ref of the values.
        # somehow.
    }
    else
    {
        carp "non-useful form for coordinates (use either array ref or string)...";
        return;
    }
    my @results;
    # build file path/name
    my $file = $self->{_location} . "/chr".$chromosome."-rec";
    # open file, 
    my $fh = new IO::File;
    if(!exists($self->{_file_handles}->{$chromosome} ) )
    {
        $self->{_file_handles}->{$chromosome} = $fh;
        $fh->open($file) or croak "can't open $file : $!";
    }
    else
    {
        $fh = $self->{_file_handles}->{$chromosome};
    }
    @$coordref = sort { $a <=> $b } @$coordref;
    foreach my $pos (@{$coordref})
    {
        # start seeking until the first position
        $fh->seek(($pos -1)*2,0);
        my $tmpval;
        $fh->read($tmpval,2);
        # unpack value, divide by 1000 and push onto array.
        my $score = unpack("n",$tmpval)/1000;
        push(@results, [ $pos, $score ]);
    }
    # return array ref of array refs of positions, scores?
    return \@results;
}

=head1 DESCRIPTION

B<MG::ConsScore> is a module for accessing the conservation scores out of fixed record
files.  This is a quick 'n' dirty module that probably could be improved upon.

=head1 OPTIONS

=over 4

=back

=head1 FILES

Requires the fixed record conservation files.

=over 4


=back

=head1 DIAGNOSTICS

=over 4

=back

=head1 REQUIRES


=head1 SEE ALSO


=head1 BUGS


=head1 AUTHOR

John Osborne B<josborne@watson.wustl.edu>

=head1 COPYRIGHT AND LICENSE

Use modify distribute under the same terms as perl
Copyright 2008 John Osborne

=cut


1;


