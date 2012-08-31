package Genome::Model::Tools::Sam::UniqueReadNames;

use strict;
use warnings;

use Genome;
use IO::File;
use File::Basename;

class Genome::Model::Tools::Sam::UniqueReadNames {
    is  => 'Genome::Model::Tools::Sam',
    has => [
	input_bam_file => {
	    is  => 'String',
	    doc => 'Input File',
	},
        output_bam_file => {
            is  => 'String',
            doc => 'Output file, defaults to file_name.unique.bam',
	    is_optional => 1,
        },
        reference_faidx => {
            is => 'String',
            doc => 'There reference fasta file',
            is_optional => 1,
            default_value => Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa.fai',
        },
    ],
};

sub help_brief {
    'Tool to sort and unique-ify read names in BAM files';
}

sub help_detail {
    return <<EOS
    Tool to sort and unique-ify read names in BAM files.
EOS
}



sub execute {
    my $self = shift;

    #Inputs/Outputs
    my ($input_basename,$input_dirname,$input_suffix) = File::Basename::fileparse($self->input_bam_file,qw/bam/);
    $input_basename =~ s/\.$//;

    my $real_output_file;
    unless ($self->output_bam_file) {
	$self->output_bam_file($input_dirname .'/'. $input_basename .'_unique.bam');
    }

    my ($output_basename,$output_dirname,$output_suffix) = File::Basename::fileparse($self->output_bam_file,qw/bam/);
    my $tmp_dir = File::Temp::tempdir( DIR => $output_dirname, CLEANUP => 1 );

    #Sort by name
    my $name_sorted_bam_file = $tmp_dir .'/'. $input_basename .'_name_sorted.bam';
    unless (Genome::Model::Tools::Sam::SortBam->execute(
        file_name => $self->input_bam_file,
        name_sort => 1,
        output_file => $name_sorted_bam_file,
        use_version => $self->use_version,
    )) {
        die('Failed to name sort bam file '. $self->input_bam_file .' to '. $name_sorted_bam_file);
    }

    #Convert to SAM
    my $name_sorted_sam_file = $tmp_dir .'/'. $input_basename .'_name_sorted.sam';
    my $cmd = $self->samtools_path .' view -o '. $name_sorted_sam_file .' '. $name_sorted_bam_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$name_sorted_bam_file],
        output_files => [$name_sorted_sam_file],
    );

    #Unique-ify read names
    my $unique_name_sorted_sam_file = $tmp_dir .'/'. $input_basename .'_name_sorted_unique.sam';
    my $in_sam_fh = IO::File->new($name_sorted_sam_file,'r');
    my $out_sam_fh = IO::File->new($unique_name_sorted_sam_file,'w');
    my $prior_read_name;
    my $i = 1;
    while (my $line = $in_sam_fh->getline) {
        chomp($line);
        my @entry = split("\t",$line);
        my $read_name = $entry[0];
        if ($prior_read_name) {
            if ($prior_read_name eq $read_name) {
                my $new_read_name = $read_name .'-'. $i++;
                $line =~ s/$read_name/$new_read_name/;
            } else {
                $i = 1;
            }
        } else {
            $prior_read_name = $read_name;
        }
        print $out_sam_fh $line ."\n";
    }
    $in_sam_fh->close;
    $out_sam_fh->close;

    #Convert back to BAM
    unless (Genome::Model::Tools::Sam::SamToBam->execute(
        sam_file => $unique_name_sorted_sam_file,
        bam_file => $self->output_bam_file,
        ref_list => $self->reference_faidx,
        is_sorted => 0,
        index_bam => 1,
        use_version => $self->use_version,
    )) {
        die('Failed to convert sam file '. $unique_name_sorted_sam_file .' to bam file '. $self->output_bam_file);
    }
    return 1;
}


1;
