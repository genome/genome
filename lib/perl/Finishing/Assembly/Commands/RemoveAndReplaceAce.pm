package Finishing::Assembly::Commands::RemoveAndReplaceAce;

use strict;
use warnings;
no warnings 'reserved';
use Data::Dumper;
use IO::File;
use Finishing::Assembly::Ace::Exporter;
use Finishing::Assembly::Factory;

use Finfo::Std;

#attributes

my %replacement_ace_file :name(replacement_ace_file:r) :isa(string) :desc('Acefile containing the replacement contig');
my %source_ace_file :name(source_ace_file:r) :isa(string) :desc('Acefile containing contigs to be replaced');
my %scaffold_num :name(scaffold_num:r) :isa('int non_neg') :desc('scaffold num where replacement contig goes');
my %start :name(start:r) :desc('Contig number to start replacing at');
my %stop :name(stop:r) :desc('Conting number to stop replacing at');
my %output_file :name(output_file:r) :isa(string) :desc('output file name for source_ace_file with replaced contigs');
my %log_file :name(log_file:r) :isa(string) :desc('log file for outputting contigs replaced, renamed');

my %source_assembly :name(source_assembly:p);
my %source_contigs :name(source_contigs:p);
my %replacement_contig :name(replacement_contig:p);
my %exporter :name(exporter:p);

sub START{
    my $self = shift;
    
    Finfo::Logging::add_appender(
        class => 'Finishing::Assembly::Commands::RemoveAndReplaceAce',
        type => 'file',
        params => {filename=>$self->log_file, mode =>'clobber'},
    );

    $self->info_msg("Replacing contigs Contig".$self->scaffold_num.".".$self->start." to Contig".$self->scaffold_num.".".$self->stop." in file \n".$self->source_ace_file."\n from file\n".$self->replacement_ace_file);
    $self->source_assembly(Finishing::Assembly::Factory->connect('ace', $self->source_ace_file)->get_assembly);
    $self->source_contigs($self->source_assembly->contigs);
    
    $self->replacement_contig( Finishing::Assembly::Factory->connect('ace', $self->replacement_ace_file)->get_assembly->contigs->first );
    
    $self->exporter( Finishing::Assembly::Ace::Exporter->new(file=>$self->output_file) );
    

    $self->fatal_msg("stop must be larger than start!") unless $self->stop >= $self->start;
}

sub execute{
    my $self = shift;

    my $start = $self->start;
    my $stop = $self->stop;
    my $subtract = $stop - $start;
    my $assembly_tags = $self->source_assembly->tags;
    $self->exporter->export_assembly_tags($assembly_tags);
    my $acefile_for_deletes = $self->output_file.".deletes";
    #for (1..100){  #dont overwrite for multi inserts on a single scaffold
    #    next if -e $acefile_for_deletes."$_";
    #    $acefile_for_deletes .= "$_";
    #    last;
    #}
    my $ex_deletes = Finishing::Assembly::Ace::Exporter->new(file =>$acefile_for_deletes);
    my $total_ctgs=0;
    my $deletes = 0;
    my $deletes_written = 0;
    my $ctgs_written = 0;

    my $replacement_written = 0;
    print ">>writing new ace file : ".$self->replacement_ace_file;
    while (my $contig = $self->source_contigs->next){
        $total_ctgs++;
        my ($scaffold_num, $contig_num) = $self->get_nums($contig);
        unless($scaffold_num == $self->scaffold_num){
            $self->exporter->export_contig(contig => $contig) and $ctgs_written++;
            print "-";
            next;
        }
        if ($contig_num < $self->start){
            $self->exporter->export_contig(contig => $contig) and $ctgs_written++;
            print ".";
            next;
        }elsif ($contig_num <= $self->stop){
            print "+";
            $deletes++;
            $ex_deletes->export_contig(contig => $contig) and $deletes_written++; 
            $self->exporter->export_contig(contig => $self->replacement_contig, new_name => 'Contig'.$self->scaffold_num.'.'.$self->start) and $ctgs_written++ and $replacement_written++ unless $replacement_written;
            next;

        }elsif ($contig_num > $self->stop){
            print ".";
            $self->exporter->export_contig(contig => $contig) and $ctgs_written++;

            #my $new_contig_num = $contig_num - $subtract;
            #TODO consider this, is it necessary?, easier to maintain gaps on fs without it
            #$self->exporter->export_contig(contig => $contig, new_name => "Contig$scaffold_num.$new_contig_num") and $ctgs_written++;
        }
    }
    print "\n";
    my $scaf = $self->scaffold_num;
    $self->info_msg(">>Replaced contigs in assembly ace file \n".$self->source_ace_file."\n:replaced $deletes of $total_ctgs original contigs between Contig$scaf.$start to Contig$scaf.$stop inclusive.\n Wrote a total of $ctgs_written in new assembly file \n".$self->replacement_ace_file."\nWrote $deletes_written replaced contigs to \n$acefile_for_deletes");
    $self->exporter->close;
    $ex_deletes->close;
}

sub get_nums{
    my ($self, $contig) = @_;
    my ($s, $c) = $contig->name =~ /Contig(\d+)\.(\d+)/;
    return ($s, $c);
}

1;
