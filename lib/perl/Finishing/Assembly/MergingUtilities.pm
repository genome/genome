package Finishing::Assembly::MergingUtilities;

use strict;
use warnings;
no warnings 'reserved';
use File::Copy 'cp';

#use Finishing::Assembly::Factory;
use Finfo::Std;
use Finishing::Assembly::DBIx::AssemblyImporter;

my %schema :name(schema:p) :isa(object);


sub remove_and_replace{  
    my ($self, %params) = @_;
    my @missing_params;
    
    my $remove_start = delete $params{remove_start} or push @missing_params, 'remove_start';
    my $remove_stop = delete $params{remove_stop} or push @missing_params, 'remove_stop';
    my $rp_contig = delete $params{contig} or push @missing_params, 'contig';
    my $db_asm = delete $params{db_asm} or push @missing_params, 'db_asm';
    my $scaf_num = delete $params{scaffold_num}; 
    push @missing_params,'scaffold_num' unless defined $scaf_num;
    

    $self->fatal_msg("Missing params: ".join(" ", @missing_params)) if @missing_params;
    $self->fatal_msg("Replacement contig Finishing::Assembly::Contig!") unless $rp_contig->isa("Finishing::Assembly::Contig");
    $self->fatal_msg("Start contig num must be less than end contig num!") unless $remove_start <= $remove_stop;
    
    $self->schema($db_asm->result_source->schema);
    my $scaffold = $db_asm->find_related('scaffolds', {scaffold_num => $scaf_num});


    if ($scaffold->orientation and $scaffold->orientation eq '-') {
        $scaffold->reverse_orientation;
        #This allows us to preserve the left and right gap info easier by aligning left and right to high low num
    }

    my $contig_subtract = $remove_stop - $remove_start;
    #my @contigs = sort {$a->contig_num <=> $b->contig_num} $scaffold->contigs->all; fully load to memory
    my $ctg_itr = $scaffold->contigs;

    my $removed_read_count=0;
    my $removed_contig_count=0;
    my @removed_contig_ids;
    my @left_gap_infos;
    my @right_gap_infos;

    #$scaffold->schema->txn_do(
    $self->schema->txn_do(
        sub{

            my $importer = Finishing::Assembly::DBIx::AssemblyImporter->new(
                store_reads => 0, 
                store_ace   => 1,
                store_tags  => 0,
            );
            #Importing the new contig here as contig 0, so we can attach the gaps to the new contig, then delete the old one
            my $new_contig = $importer->import_single_contig(
                contig     => $rp_contig, 
                scaffold   => $scaffold, 
                contig_num => 0,
                dbix_assembly => $db_asm,
                #schema     => $schema,
            );

            #foreach my $contig (@contigs){
            while (my $contig = $ctg_itr->next) {
                my $num = $contig->contig_num;

                next if $num < $remove_start;
                
                if ($num <= $remove_stop){
                    if ($num == $remove_start){
                        my @gaps = $contig->left_gaps->all;
                        for my $gap (@gaps){
                            $gap->set_right_contig($new_contig);
                            $gap->update;
                        }
                    }
                    elsif ($num == $remove_stop){
                        my @gaps = $contig->right_gaps->all;
                        for my $gap (@gaps){
                            $gap->set_left_contig($new_contig);
                            $gap->update;
                        }
                    }
                    $removed_contig_count++;
                    $removed_read_count += $contig->read_count;
                    push @removed_contig_ids, $contig->id;

                    $contig->delete;

                    $self->info_msg("Removing ".$contig->name.", id ".$contig->id);

                    next;
                }

                #renumber higher contigs
                my $new_num = $num-$contig_subtract;
               
                $contig->contig_num($new_num) 
                    and $self->debug_msg("Changing Contig$num to Contig$new_num");
                $contig->update;
            }
            $self->info_msg("$removed_contig_count Contigs removed with $removed_read_count reads");

            #now set the new contig to the proper number
            $new_contig->contig_num($remove_start);
            $new_contig->update;


            my $replacement_contig = $scaffold->get_contig($remove_start);
            $self->create_neighbor_gaps($replacement_contig, \@left_gap_infos, \@right_gap_infos);
            $self->create_replaced_contig_events(replacement_contig => $replacement_contig, removed_ids => \@removed_contig_ids);
            $self->info_msg($self->assembly->read_iterator->count." reads in assembly after replace");
        }
    ); #end txn_do
}

sub create_replaced_contig_events{
    my ($self, @ids) = @_;
    my %params;
    my $replacement_id;

    my $removed_ids = delete $params{'removed_ids'};
    my $schema = $self->schema;
    $schema->txn_do(
        sub{
            foreach (@$removed_ids){
                my $event = $schema->resultset('Event')->create(
                    {
                        type => "replaced_contig_event",
                        created_by => "contiguous_bac_merge",
                    }
                );
                my $replaced_contig = $schema->resultset('ReplacedContig')->create(
                    {
                        id => $_,
                        replacement_id => $replacement_id,
                    }
                );
                my $replaced_contig_event = $schema->resultset('ReplacedContigEvent')->create(
                    {
                        old_contig_id => $_,
                        new_contig_id => $replacement_id,
                        event_id => $event->id,
                    }
                );
            }
        }
    ); #end txn_do
}
1;
