#!/usr/bin/env genome-perl

use strict;
use warnings;

use Data::Dumper;
use Finfo::Logging 'info_msg';
use Finishing::Assembly::Factory;

my $factory = Finishing::Assembly::Factory->connect
(
    'cmap_admin',
    #'dbi:mysql:cmap_sorting:mysql1', 'sorting_admin', 's_admin', { AutoCommit => 1},
);
my $org = $factory->get_organism('pan_troglodytes');
my $assem = $org->get_assembly(name => '2.1_051011');
my %chr_params;
$chr_params{name} = \@ARGV if @ARGV;
my $chr_iterator = $org->chromosome_iterator(%chr_params);
while ( my $chr = $chr_iterator->next )
{
    $factory->txn_do
    (
        sub
        {
            my $chr_name = $chr->name;
            main->info_msg(sprintf('Working on chr %s', $chr_name));

            my $is_ordered;
            if ( $chr_name =~ s/_random//i or $chr_name eq 'Un' )
            {
                $is_ordered = 0;
            }
            else
            {
                $is_ordered = 1;
            }
            
            my $ctgs_chr;
            if ( $chr_name eq $chr->name )
            {
                $ctgs_chr = $chr;
            }
            else
            {
                $ctgs_chr = $org->get_chromosome(name => $chr_name);
            }
            
            my $first_scaff = $chr->first_scaffold_for_assembly($assem);
            my $ctg = $first_scaff->scaffold->first_contig;
            $factory->txn_do( sub{ _update_contig($ctg, $ctgs_chr, $is_ordered); } );

            my $chr_last_ctg = $ctg;
            
            while ( $ctg = _next_contig($ctg) )
            {
                $factory->txn_do
                (
                    sub
                    {
                        _update_contig($ctg, $chr, $is_ordered);
                    } 
                );

                $chr_last_ctg = $ctg;
            }

            unless ( $is_ordered )
            {
                $factory->txn_do
                (
                    sub
                    {
                        foreach my $gap ( $chr_last_ctg->right_gaps )
                        {
                            $gap->delete;
                        }
                        $first_scaff->delete;
                        $chr->delete;
                        return 1;
                    }
                );
            }
        }
    );
}

exit 0;

sub _next_contig
{
    my ($ctg) = @_;

    return unless $ctg;

    return $ctg->right_contig;
}

sub _update_contig
{
    my ($ctg, $chr, $is_ordered) = @_;

    $ctg->chromosome($chr, $is_ordered);

    unless ( $is_ordered )
    {
        foreach my $gap ( $ctg->left_gaps )
        {
            $gap->delete;
        }
    }

    my $ctg_chr = $ctg->chromosome;
    main->fatal_msg("Did not set contigs chromosome") unless $ctg_chr
        and $ctg_chr->name eq $chr->name;
    
    return 1;
}

=pod

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

