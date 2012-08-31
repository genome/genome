package Genome::Model::Tools::ChimeraSlayer::Reader;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::ChimeraSlayer::Reader {
    is => 'Genome::Utility::IO::Reader',
};

# 0      ChimeraSlayer
# 1      chimera_AJ007403            # the accession of the chimera query
# 2      S000387216                  # reference parent A
# 3      S000001688                  # reference parent B
# 4      0.9422                      # divergence ratio of query to chimera (left_A, right_B)
# 5      90.00                       # percent identity between query and chimera(left_A, right_B)
# 6      0                           # confidence in query as a chimera related to (left_A, right_B)
# 7      1.0419                      # divergence ratio of query to chimera (right_A, left_B)
# 8      99.52                       # percent identity between query and chimera(right_A, left_B)
# 9      100                         # confidence in query as a chimera related to (right_A, left_B)
# 10     YES                         # ** verdict as a chimera or not **
# 11     NAST:4032-4033              # estimated approximate chimera breakpoint in NAST coordinates
# 12     ECO:767-768                 # estimated approximate chimera breakpoint according to the E. coli unaligned reference seq coordinates
 
my @fields = (qw/ 
    id parent_A parent_B
    divergence_ratio_left_A_right_B percent_identity_left_A_right_B confidence_left_A_right_B
    divergence_ratio_right_A_left_B percent_identity_right_A_left_B confidence_right_A_left_B
    verdict est_breakpoint_NAST est_breakpoint_ECOLI
    /);
sub fields {
    return @fields;
}

sub read {
    my $self = shift;

    while ( my $line = $self->input->getline ) {
        next if not $line =~ s/^ChimeraSlayer\t//;
        chomp $line;
        my @tokens = split(/\t/, $line);
        Carp::confess( # need 10 or 12 fields
            'Malformed chimera slayer line! Got '.@tokens.' fields, but expected 10 or '.@fields.' from chimera slayer line: '.$line
        ) if @tokens != 10 and @tokens != @fields; # 10 or 12
        my %chimera;
        @chimera{@fields} = @tokens;
        Carp::confess( # check verdict
            'Malformed chimera slayer line! Verdict is expected to be YES, NO or UNKNOWN but is '.$chimera{verdict}.'. Line: '.$line
        ) if not $chimera{verdict} or not grep { $chimera{verdict} eq $_ } (qw/ YES NO UNKNOWN /);
        return \%chimera;
    }

    return;
}

1;

