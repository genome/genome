package Genome::Model::RnaSeq::DetectFusionsResult::BreakFusionResult::Index;

use strict;
use warnings;

class Genome::Model::RnaSeq::DetectFusionsResult::BreakFusionResult::Index {
    is => 'Genome::Model::RnaSeq::DetectFusionsResult::Index::Base',
    has => [
        twobit_file => {
            is => "Text",
            is_calculated => 1,
            calculate_from => ["output_dir"],
            calculate => q|
                $output_dir . "/all_sequences.2bit"
            |
        },
    ],
    has_input => [
        #TODO This should be moved up to the parent class
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            doc => 'The annotation build from which to derive the gene file',
        },
    ],
};

sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/breakfusion-index/' . $self->id;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_) or return;
    $self->_prepare_staging_directory;

    $self->prepare_gene_file($self->temp_staging_directory . '/refGene.txt');

    $self->_download_ucsc_table('chainSelf.txt.gz', $self->temp_staging_directory);

    $self->_split_chain_annotation_file();
    unlink($self->temp_staging_directory . "/chainSelf.txt");

    $self->_get_2bit_for_reference();

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub prepare_gene_file {
    my $self = shift;
    my $gene_file = shift;

    my $annotation_build = $self->annotation_build;

    my $itr = $annotation_build->transcript_iterator;

    my $output_file = Genome::Sys->open_file_for_writing($gene_file);

    while (my $tx = $itr->next()){
        next if $tx->coding_region_start eq "NULL";
        next unless $tx->gene_name;

        my @line = ();
        push @line, $tx->transcript_name, $tx->chrom_name;

        #strand needs to be +/- instead of +1/-1
        (my $strand_val = $tx->strand) =~ s/\+1/\+/;
        $strand_val =~ s/-1/-/;

        push @line, $strand_val, $tx->transcript_start, $tx->transcript_stop, $tx->coding_region_start, $tx->coding_region_stop;

        my @start, my @stop = ();
        for (grep {$_->structure_type =~/cds_exon|utr_exon/ } $tx->ordered_sub_structures){
            push @start, $_->structure_start;
            push @stop, $_->structure_stop;
        }

        #yes, these trailing commas are needed. no, I don't like it either.
        push @line, scalar(@start), join(",",@start) . ",", join(",",@stop) . ",", $tx->transcript_id, $tx->gene_name;

        $output_file->say(join("\t", @line));

        Genome::Transcript->unload();
        Genome::TranscriptStructure->unload();
    }

    $output_file->close();

    return 1;
}
sub _get_2bit_for_reference {
    my $self = shift;
    if(my $path = $self->reference_build->full_consensus_path("2bit")){
        Genome::Sys->create_symlink($path, $self->temp_staging_directory . "/all_sequences.2bit");
    }else{
        my $path = $self->reference_build->full_consensus_path("fa");
        my $output_path = $self->temp_staging_directory . "/all_sequences.2bit";
        Genome::Sys->shellcmd(
            cmd => "faToTwoBit $path $output_path",
            input_files => [$path],
            output_files => [$output_path]
        );
    }
}

sub _split_chain_annotation_file {
    my $self = shift;

    my $fh = Genome::Sys->open_file_for_reading($self->temp_staging_directory . "/chainSelf.txt");
    my %out_fhs = ();

    while (<$fh>){
        my $chr = (split(/\s+/, $_))[2];
        unless($out_fhs{$chr}){
            $out_fhs{$chr} = Genome::Sys->open_file_for_writing($self->temp_staging_directory . "/chainSelf.chr" . $chr . ".tab" );
        }
        $out_fhs{$chr}->print($_);
    }

    for(values  %out_fhs){
        $_->close();
    }

    $fh->close();
    return 1;
}

1;

