package Finishing::Assembly::Scaffold;

use strict;
use warnings;

use base 'Finishing::Assembly::Item';

sub contig_count
{
    return shift->contigs->count;
}

sub reverse_orientation{
    my $self = shift;
    my @contigs_to_complement_ace;
    if ($self->orientation eq '+'){
        $self->fatal_msg("No point in reversing a scaffold already in positive orientation!"); #no soup for you
    }
    eval{
        $self->schema->txn_do(
            sub{
                my @contigs = $self->contigs->all;

                while (@contigs){
                    my $contig_a = shift @contigs;
                    $contig_a->complement;

                    my $contig_b = pop @contigs;
                    next unless $contig_b;
                    $contig_b->complement;

                    my $temp_a = $contig_a->contig_num;
                    my $temp_b = $contig_b->contig_num;
                    $contig_a->set_contig_num(0);
                    $contig_b->set_contig_num($temp_a);
                    $contig_a->set_contig_num($temp_b);

                    push @contigs_to_complement_ace, $contig_a;
                    push @contigs_to_complement_ace, $contig_b;
                }

                $self->orientation($self->orientation eq '+'? '-' : '+');
            }
        );
    };

    unless ($@){
        foreach my $contig (@contigs_to_complement_ace){
            my $acefile = $contig->acefile;
            copy($acefile, $acefile.".bak");
            my $ace_factory = Finishing::Assembly::Factory->connect('ace', $acefile);
            $ace_factory->commit;
        }
    }
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Scaffold.pm $
#$Id: Scaffold.pm 31362 2007-12-28 18:17:08Z adukes $
