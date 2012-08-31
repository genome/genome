package Genome::Model::Tools::PhredPhrap;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::PhredPhrap {
    is => 'Command',
    is_abstract => 1,
    has_optional => [ __PACKAGE__->properties_hash ],
};

sub properties_hash {
    return (
        version  => {
            type => 'String',
            default => __PACKAGE__->default_version,#=> $versions[0],
            doc => "Version to use: " . join(', ', __PACKAGE__->versions)
        },
        forcelevel  => {
            type => 'Integer',
            default => 1,
            doc => 'Relaxes stringency to varying degree during final contig merge pass.  Allowed values are integers from 0 (most stringent, to 10 (least stringent). (1)',
        },
        minmatch  => {
            type => 'Integer',
            default => 17,
            doc => 'Minimum length of matching word to nucleate SWAT comparison. if minmatch = 0, a full (non-banded, comparison is done [N.B. NOT PERMITTED IN CURRENT VERSION]. Increasing -minmatch can dramatically decrease the time required for the pairwise sequence comparisons; in phrap, it also tends to have the effect of increasing assembly stringency. However it may cause some significant matches to be missed, and it may increase the risk of incorrect joins in phrap in certain situations (by causing implied overlaps between reads with high-quality discrepancies to be missed). (17)',
        },
        minscore  => {
            type => 'Integer',
            default => 30,
            doc => 'Minimum alignment score. (30)',
        },
        revise_greedy  => {
            type => 'Boolean',
            default => 0,
            doc => 'Splits initial greedy assembly into pieces at "weak joins", and then tries to reattach them to give higher overall score.  Use of this option should correct some types of missassembly.',
        },
        shatter_greedy  => {
            type => 'Boolean',
            default => 0,
            doc => 'Breaks assembly at weak joins (as with revise-greedy, but does not try to reattach pieces.',
        },
        view  => {
            type => 'Boolean',
            default => 0,
            doc => 'Create ".view" file suitable for input to phrapview. (0)',
        },
        new_ace  => {
            type => 'Boolean',
            default => 1,
            doc => 'Create ".ace" file for viewing in consed. Default is to create an acefile. (1)',
        },
        trim_qual => {
            type => 'Integer',
            default => 13,
            doc => 'Quality value used in to define the "high-quality" part of a read, (the part which should overlap; this is used to adjust qualities at ends of reads. (13)'
        },
        vector_bound => {
            type => 'Integer',
            default => 80,
            doc => 'Number of potential vector bases at beginning of each read.  Matches that lie entirely within this region are assumed to represent vector matches and are ignored. (80)',
        },
    );
}

#- COMMAND -#
sub phrap_command_name {
    my $self = shift;

    # TODO get full path of phrap executable?
    return sprintf('phrap%s', (( $self->version and $self->version ne 'phrap') ? ('.' . $self->version) : ''));
}

#- VERSIONS -#
my @versions = (qw/ phrap manyreads longreads /);
sub versions {
    return @versions;
}

sub default_version {
    return $versions[0];
}

1;

=pod

=head1 Name

=head1 Synopsis

=head1 Methods

=head1 Disclaimer

 Copyright (C) 2006 Washington University Genome Sequencing Center

 This module is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

=head1 Author(s)

 Eddie Belter <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
