package Genome::Model::Tools::Somatic::ReadCounts;

use warnings;
use strict;

use Genome;
use Workflow;
use Carp;
use IO::File;
use Data::Dumper;

class Genome::Model::Tools::Somatic::ReadCounts {

    is  => ['Command'],
    has => [
        bam_readcount_version => {
            is => 'Version',
            doc => 'Version of bam readcount to utilize',
        },
        tumor_bam => {
            is => 'String',
            doc =>'Path to the tumor bam file',
        },
        normal_bam => {
            is => 'String',
            doc =>'Path to the normal bam file',
        },
        sites_file => {
            is => 'String',
            doc =>'the sites of interest in annotation format.'
        },
        minimum_mapping_quality => {
            is => 'String',
            default => 1,
            doc =>'the minimum mapping quality of a read to include it in the counts'
        },
        reference_sequence => {
            is => 'String',
            doc =>'the reference sequence to use, defaults to NCBI-human-build36',
            example_values => [Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa'],
        },
        output_file => {
            is => 'String',
            doc =>'path to output file',
        },
    ],
};

sub help_brief {
    return "generate read count statistics";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    gmt somatic read-counts --tumor /path/to/tumor.bam --normal /path/to/normal.bam --sites /path/to/sites.file --output /path/to/output.out
EOS
}

sub help_detail {
    return <<EOS 
    Produces a tab-delimited file of statistics from the 'bam-readcount' command.
EOS
}

sub execute {
    my ($self) = @_;
    my ($readcount_regions_fh, $readcount_regions_file)  =Genome::Sys->create_temp_file();
    my $anno_fh = IO::File->new($self->sites_file);
    my $output_fh = IO::File->new($self->output_file, ">");
    unless($output_fh) {
        $self->error_message("Couldn't open output file " . $self->output_file);
        return 0;
    }
    unless($anno_fh) {
        $self->error_message("Couldn't open sites file " . $self->sites_file);
        return 0;
    }

    unless (`uname -a` =~ /x86_64/) {
        $self->error_message("Must run on a 64 bit machine");
        return 0;
    }
    
    while (my $line = $anno_fh->getline) {
        chomp $line;
        my ($chr, $pos,) = split /\t/, $line;
        $readcount_regions_fh->print("$chr\t$pos\t$pos\n");
    }
    $readcount_regions_fh->close;
    my $min_mapping_quality = $self->minimum_mapping_quality;
    my $normal_temp = Genome::Sys->create_temp_file_path;
    my $tumor_temp = Genome::Sys->create_temp_file_path;
    my $normal_rc = Genome::Model::Tools::Sam::Readcount->execute(
        bam_file => $self->normal_bam,
        minimum_mapping_quality => $min_mapping_quality,
        output_file => $normal_temp,
        reference_fasta => $self->reference_sequence,
        region_list => $readcount_regions_file,
        use_version => $self->bam_readcount_version,
    );
    my $tumor_rc = Genome::Model::Tools::Sam::Readcount->execute(
        bam_file => $self->tumor_bam,
        minimum_mapping_quality => $min_mapping_quality,
        output_file => $tumor_temp,
        reference_fasta => $self->reference_sequence,
        region_list => $readcount_regions_file,
        use_version => $self->bam_readcount_version,
    );
    unless ($normal_rc->result and $tumor_rc->result) {
        $self->error_message("readcount failed\n");
        return;
    }

    my @normal_lines = `cat $normal_temp`;
    my @tumor_lines  = `cat $tumor_temp`;
    my %hash_of_arrays;
    $hash_of_arrays{'Normal'}=\@normal_lines;
    $hash_of_arrays{'Tumor'}=\@tumor_lines;
    $self->make_excel_friendly_output_sheet($output_fh, \%hash_of_arrays);
    return 1;
}

sub make_excel_friendly_output_sheet {
    my ($self, $output_fh, $tumor_normal_hash_ref) = @_;
    my %hash_of_arrays = %{$tumor_normal_hash_ref};
    $output_fh->print("CHROM\tPOS\tREF\tVAR\tGENE\tMUTATION\t" . join("\t\t\t", sort keys %hash_of_arrays) . "\n");
    my $anno_fh = IO::File->new($self->sites_file);
    while (my $line = $anno_fh->getline) {
        chomp $line;
        my ($chrom, $pos, $foo, $ref, $var, $bstype, $gene, $transcript,
        $organism, $source, $version, $some_shit, $transcipt_status, $type, @foobar)= split /\t/, $line; #add one field before  second type to get this crap to work correctyly with new annotator
        $output_fh->print("$chrom\t$pos\t$ref\t$var\t$gene\t$type\t");

        for my $flowcell_id(sort keys %hash_of_arrays) {
            if(my ($line1, $line2) = grep( /$chrom\t$pos/, @{$hash_of_arrays{$flowcell_id}})) {
                my ($stats_for_var1) = ($line1 =~ m/($var\S+)/);
                my ($stats_for_ref1) = ($line1 =~ m/($ref\S+)/);
                my ($base_var1, $count_var1, ) = split /:/, $stats_for_var1;
                my ($base_ref1, $count_ref1, ) = split /:/, $stats_for_ref1;
                my $percent = 0;
                if($count_ref1 > 0 || $count_var1 > 0) {
                    $percent = $count_var1 / ($count_ref1 + $count_var1);
                }
                $output_fh->print("\t$count_ref1\t$count_var1\t$percent");
            }
            else {
                $output_fh->print("\t0\t0\t0");
            }
        }
        $output_fh->print("\n");
    }
}

1;
