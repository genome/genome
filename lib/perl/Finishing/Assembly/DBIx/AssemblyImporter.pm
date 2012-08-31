package Finishing::Assembly::DBIx::AssemblyImporter;

=pod

=head1 Finishing::Assembly::DBIx::AssemblyImporter

=head1 SYNOPSIS

This module is used to import an ace file into the DBIX db/filesystem schema for the Assembly Group.

=head1 DESCRIPTION

my $importer = Finishing::Assembly::AceImporter->new
    (
    root_path   => <path_for_single_contig_acefiles>,   ( optional.  default => '/gscmnt/843/finishing/assembly/WholeGenomeImprovement/' )
    store_tags  => 1, ( optional, default 1 )
    store_reads => 1, ( optional, default 1 )
    store_ace   => 1, ( optional, default 1 ) 
    );

$importer->import_acefile
    (
    ace_assembly => <Finishing::Assembly::Assembly object>, assembly with acefile source
    dbix_assembly => <Finishing::Assembly::Assembly object>, assembly with dbix database source
    );


=head1 METHODS

=over 4

=item import_acefile 

=item add_item_callback

=back

=cut

use Data::Dumper;
use strict;
use warnings;
no warnings 'reserved';

use DateTime;
use Finfo::Std;
use Finishing::Assembly::Ace::Exporter;

my %root_path :name(root_path:o) :isa(string) :default('/gscmnt/853/finishing/assembly/WholeGenomeImprovement/');
my %store_tags :name(store_tags:o) :default(1);
my %store_reads :name(store_reads:o) :default(1);
my %store_ace :name(store_ace:o) :default(1);

my %dbix_assembly :name(dbix_assembly:p); #DBIx::Class object
my %ace_assembly :name(ace_assembly:p); #Behavior assembly
my %schema :name(schema:p);

sub import_ace_assembly {
    my ($self, %params) = @_;
    #get and check for required params

    my @missing_params;
    my $dbix_assembly = delete $params{dbix_assembly} or push @missing_params, 'dbix_assembly';
    my $ace_assembly = delete $params{ace_assembly} or push @missing_params, 'ace_assembly';
    $self->fatal_msg("Needed parameters missing:".join(' ', @missing_params)) if @missing_params;

    $self->schema($dbix_assembly->result_source->schema);
    $self->dbix_assembly($dbix_assembly); 
    $self->ace_assembly($ace_assembly);

    my $contigs = $ace_assembly->contigs;

    my @new_contigs;

    my $finished = 0;
    my $counter = 0;
    my $total = $contigs->count;
    while (!$finished){
        $self->schema->txn_do ( 
            sub{
                for (1..500){
                    my $contig = $contigs->next;
                    unless ($contig){
                        $finished = 1;
                        last;
                    }
                    my $new_contig=$self->_import_contig($contig);
                    $self->fatal_msg("Couldn't import contig: ".$contig->name."!") unless $new_contig; 
                    $counter++;
                    print $counter.'...' if $counter % 20 == 0;
                }
            }
        );
    }
    print "\n$counter of $total contigs imported\n";
    $self->schema->txn_do(
        sub{

            if ($self->store_tags) {
                foreach my $tag (@{$ace_assembly->tags}){ 
                    $self->fatal_msg( "couldn't add tag to ". $self->dbix_assembly->name) 
                        unless $dbix_assembly->add_to_tags(convert_tag($tag, 'a'));
                }
            }
            return \@new_contigs;
        }
    );
}

sub _import_contig {
    my ($self, $contig) = @_;
    my $scaffold;
    my $new_contig;
    my $scaffold_num;
    my $contig_num;
    my $already = 0;
    eval{
        if ($contig->name =~ /Contig\d+\.\d+/) {
            ($scaffold_num, $contig_num) = $contig->name =~ /Contig(\d+)\.(\d+)/;
        }
        else {
            #TODO what to do with orphans?
            $self->fatal_msg("Contig name must be in the form Contig<scaf_num>.<contig_num>!");
        }
        $self->schema->txn_do(
            sub{
                $scaffold = $self->dbix_assembly->find_or_create_related( 'scaffolds', {
                        scaffold_num => $scaffold_num 
                    }
                );
            }
        );
        $self->fatal_message("Can't find or create new scaffold for assembly!") unless $scaffold;

        $self->schema->txn_do(
            sub{
                if ( $new_contig = $scaffold->find_related('contigs', { contig_num => $contig->contig_num } ) ){
                    #if ( $new_contig = $scaffold->get_contig($contig->contig_num) ){
                    #$self->fatal_msg("Contig $scaffold_num.$contig_num already exists in assembly: ".$contig->assembly->name."\n")reads.placed.out;
                    $self->info_msg("Contig $scaffold_num.$contig_num already exists in assembly: ".$new_contig->assembly->name."\n");
                    $already = 1;
                }
                else {
                    $self->schema->txn_do(
                        sub{
                            $new_contig = $scaffold->add_to_contigs(
                                {
                                    contig_num => $contig_num,
                                    scaffold_id => $scaffold->id,
                                    length => $contig->base_count,
                                    assembly_id => $scaffold->assembly_id,
                                }
                            );
                        }
                    );
                }
                unless ($already){
                    my $quals = $contig->qualities;
                    my $qual_string = "@$quals";

                    $self->schema->txn_do(
                        sub{
                            $new_contig->create_related('consensus',
                                {
                                    bases => $contig->base_string,
                                    qualities => $qual_string,
                                }
                            );
                        }
                    );
                    if ($self->store_tags){
                        foreach my $tag(@{$contig->tags}){
                            my $arg = convert_tag($tag, 'c');
                            $self->fatal_msg("couldn't add tag to ".$new_contig->name) 
                                unless $new_contig->add_to_tags($arg);
                        }
                    }
                    if ($self->store_reads){
                        foreach my $read ($contig->reads->all){
                            $self->_import_read($read,$new_contig);
                        }
                    }
                }
            }
        );
    };

    if ($@){
        $self->fatal_msg("Error in importing Contig $scaffold_num.$contig_num!\n $@");
    }else{
        return 1 if $already;
        if ($self->store_ace){
            print 'writing ace contig: '.$new_contig->name.' '.$new_contig->id."\n";
            $self->_write_acefile_for_contig($contig, $new_contig->id);
        }
    }

    return $new_contig;
}

sub _import_read {
    my ($self, $read, $new_contig) = @_;

    my $name = $read->name;

    my $new_read;
    $self->schema->txn_do(
        sub {

            if ($new_read = $new_contig->find_related('assembled_reads',{name => $name})){
                $self->info_msg("Read $name is already in assembly!" );
            }
            else{
                $self->schema->txn_do(
                    sub{
                        $new_read = $new_contig->add_to_assembled_reads(
                            {
                                name            => $name,
                                template_name   => $read->template,
                                direction       => $read->direction, 
                                complemented    => $read->complemented,
                                length          => $read->length,
                                start_position  => $read->position, 
                                stop_position   => $read->position + $read->length - 1,
                                assembly_id     => $new_contig->assembly_id,
                                contig_id       => $new_contig->id,
                            }
                        );
                    }
                );

                $self->fatal_msg("Couldn't create read $name!") unless $new_read;
            }
            if ($self->store_tags){
                foreach my $tag (@{$read->tags}){
                    $self->fatal_msg("Couldn't add tag to ".$new_read->name) 
                        unless $new_read->add_to_tags(convert_tag($tag, 'r'));
                    }
            }
        }
    );
    return $read;
}

sub _write_acefile_for_contig{
    my ($self, $contig, $contig_id) = @_;
    my $path = $self->dbix_assembly->file_path;
    my $ace_path = "$path/$contig_id.ace";
    $self->fatal_msg("Can't generate assembly file path for writing acefile!") unless defined $path;
    if(!-e $path) { 
        system "mkdir -p $path/"; }
    if (-e $ace_path){
        $self->fatal_msg("error contig acefile found! $ace_path\n");
    }else{
        my $exp = Finishing::Assembly::Ace::Exporter->new(file => $ace_path);
        $exp->export_contig(contig=>$contig, new_name=>$contig_id);
        $exp->close;
    }
}


#TODO - not up to date, still using api, which is probably ok
sub import_single_contig {
    my ($self, %params) = @_;
    #takes in a contig object, adds it to the sorting database, and writes the contig out in ace format
    #get and check for required params
    my @missing_params;
    my $dbx_asm = delete $params{dbix_assembly} or push @missing_params, 'dbix_assembly';
    my $contig = delete $params{contig} or push @missing_params, 'contig';
    my $scaffold = delete $params{scaffold} or push @missing_params, 'scaffold';
    my $contig_num = delete $params{contig_num};
    push @missing_params, 'contig_num' unless defined $contig_num;
    
    $self->fatal_msg("Needed parameters missing:".join(' ', @missing_params)) if @missing_params;

    #$self->schema($scaffold->schema) or $self->fatal_msg("Can't get schema from scaffold:".ref($scaffold));
    $self->dbix_assembly($dbx_asm);
    $self->schema($dbx_asm->result_source->schema) or $self->fatal_msg("Can't get schema");
    my $new_contig;
    eval{
        $self->schema->txn_do(
            sub{
                #$new_contig = $scaffold->get_contig($contig->contig_num);
                $new_contig = $scaffold->get_contig($contig_num);
                $self->fatal_msg("Found ".$new_contig->name." already in assembly") if $new_contig; 
                #$new_contig = $scaffold->add_contig($contig, $contig_num);
                $new_contig = $scaffold->add_to_contigs(
                    {
                        contig_num => $contig_num,
                        scaffold_id => $scaffold->id,
                        length => $contig->base_count,
                        assembly_id => $scaffold->assembly_id,
                    }
                );

                my $quals = $contig->qualities;
                my $qual_string = "@$quals";

                #$new_contig->create_sequence(
                #    bases => $contig->base_string,
                #    qualities => $qual_string,
                #);

                $new_contig->create_related('consensus',
                    {
                        bases     => $contig->base_string,
                        qualities => $qual_string,
                    }
                );

                if ($self->store_tags){
                    foreach my $tag (@{$contig->tags}){
                        #$self->fatal_msg("Couldn't add read to ".$new_contig->name) unless $new_contig->add_tag($tag);
                        my $arg = convert_tag($tag, 'c');
                        $self->fatal_msg("couldn't add tag to ".$new_contig->name) unless $new_contig->add_to_tags($arg);
                    }
                }
                if ($self->store_reads){
                    foreach my $read ($contig->reads->all){
                        $self->_import_read($read,$new_contig);
                    }
                }
            }
        );
    };
    if ($@){
        $self->fatal_msg("Couldn't import contig! \n $@");
    }else{
        if ($self->store_ace){
            $self->_write_acefile_for_contig($contig, $new_contig->id);
        }
    }

    return $new_contig;
}


sub convert_tag {
    my ($tag, $type) = @_;

    my $date = sprintf("TO_DATE('%s', 'YYYY-MM-DD HH24:MI:SS')", $tag->date);

    my $arg = {
        tag_type      => $tag->type,
        source        => $tag->source,
        creation_date => \$date, 
    };
    
    unless ($type eq 'a') {
        $arg->{start_position} = $tag->start;
        $arg->{stop_position}  = $tag->stop;
        $arg->{tag_comment}    = $tag->comment  if $tag->comment;
    }
    
    $arg->{no_trans} = $tag->no_trans ? 1 : 0 if $type eq 'c';
    $arg->{text} = $tag->text if $tag->text;

    return $arg;
}

1;

#$Header$
