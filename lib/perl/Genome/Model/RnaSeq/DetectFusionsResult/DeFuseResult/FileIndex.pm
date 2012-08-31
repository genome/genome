package Genome::Model::RnaSeq::DetectFusionsResult::DeFuseResult::FileIndex;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::DetectFusionsResult::DeFuseResult::FileIndex {
    is => 'Genome::Model::RnaSeq::DetectFusionsResult::Index::Base',
    has => [
        genome_fasta => {
            is => 'Text',
            is_calculated => 1,
            calculate =>  q{ "$output_dir/all_sequences.fa" },
            calculate_from => ['output_dir'],
        },
        gene_models => {
            is => 'Text',
            is_calculated => 1,
            calculate => q{ "$output_dir/ensemble.gtf" },
            calculate_from => ['output_dir'],
        },
        repeats_filename => {
            is => 'Text',
            is_calculated => 1,
            calculate => q{ "$output_dir/rmsk.txt" },
            calculate_from => ['output_dir'],
        },
        est_fasta => {
            is => 'Text',
            is_calculated => 1,
            calculate => q{ "$output_dir/est.fa" },
            calculate_from => ["output_dir"],
        },
        est_alignments => {
            is => 'Text',
            is_calculated => 1,
            calculate => q{ "$output_dir/intronEst.txt"},
            calculate_from => ["output_dir"],
        },
        unigene_fasta => {
            is => 'Text',
            is_calculated => 1,
            calculate => q{ "$output_dir/Hs.seq.uniq"},
            calculate_from => ["output_dir"],
        },
    ],
    has_param => [
        ensembl_release => {
            is => 'Text',
            doc => 'the ensembl release number, used to determine the gtf file to download for use in deFuse'
        },
    ],
    doc => 'create all the files needed for a successful run of deFuse'
};


sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/defuse-index/' . $self->id;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_) or return;

    unless ($self->reference_build->species_name eq "human"){
        die($self->error_message("deFuse only supports human fusion detection"));
    }

    $self->_prepare_staging_directory;

    $self->_create_modified_fasta();
    $self->_download_ensembl_file($self->ensembl_release);
    $self->_download_est_file($self->_ucsc_build_name);

    for(qw(rmsk.txt.gz intronEst.txt.gz)){
        $self->_download_ucsc_table($_, $self->temp_staging_directory);
    }

    $self->_prepare_output_directory;

    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub _create_modified_fasta {
    my $self = shift;
    my $fasta_in = Genome::Sys->open_file_for_reading($self->reference_build->full_consensus_path("fa"));
    my $fasta_out = Genome::Sys->open_file_for_writing($self->temp_staging_directory . "/all_sequences.fa");

    while(my $line = <$fasta_in>){
        if($line =~ /^>/){
            #remove everything past the chromosome name
            $line =~ s/(^>\S*)(\s.+)/$1\n/;
        }
        print $fasta_out $line;
    }

    $fasta_in->close();
    $fasta_out->close();
}

#this file doesn't appear to be versioned in any way...
sub _download_unigene_clusters{
    my $self = shift;

    my $cmd = "wget -O ".$self->temp_staging_directory .  "/Hs.seq.uniq.gz ftp://ftp.ncbi.nih.gov/repository/UniGene/Homo_sapiens/Hs.seq.uniq.gz";
    $cmd .= " && gunzip " . $self->temp_staging_directory . "/Hs.seq.uniq.gz";

    Genome::Sys->shellcmd(
        cmd => $cmd,
        output_files => [$self->temp_staging_directory . "/Hs.seq.uniq"],
    );
}

sub _download_est_file {
    my ($self, $release) = @_;
    my $path = "ftp://hgdownload.cse.ucsc.edu/goldenPath/$release/bigZips/est.fa.gz";


    my $cmd = "wget -O ". $self->temp_staging_directory ."/est.fa.gz $path";
    $cmd .= " && gunzip " . $self->temp_staging_directory . "/est.fa.gz";

    Genome::Sys->shellcmd(
        cmd => $cmd,
        output_files => [$self->temp_staging_directory . "/est.fa"],
    );
}

sub _download_ensembl_file {
    my ($self, $release) = @_;
    #my $header = '#' . join("\t",('bin','swScore','milliDiv','milliDel','milliIns','genoName','genoStart','genoEnd','genoLeft','strand','repName','repClass','repFamily','repStart','repEnd','repLeft','id'));
    my $path = "ftp://ftp.ensembl.org/pub/release-$release/gtf/homo_sapiens/*";

    my $cmd = "wget -O " . $self->temp_staging_directory . "/ensembl.gtf.gz $path";
    $cmd .= " && gunzip " . $self->temp_staging_directory . "/ensembl.gtf.gz";

    Genome::Sys->shellcmd(
        cmd => $cmd,
        output_files => [$self->temp_staging_directory . "/ensembl.gtf"],
    );
}

1;
