package Genome::Model::Tools::DetectVariants2::Filter::LibrarySupport;

use warnings;
use strict;

use Genome;

class Genome::Model::Tools::DetectVariants2::Filter::LibrarySupport{
    is => 'Genome::Model::Tools::DetectVariants2::Filter',
    doc => "Outputs list of input indels from Sniper that have high library support.",
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants2 filter library-support --input-directory /path/to/sniper/out/ --output-directory /some/path
EOS
}

sub help_detail {                           
    return <<EOS
Requires an input list of indels from somatic sniper.
Outputs list of indels that have high library support. The output columns
are the same as the input indel_file with the addition of two new columns:
the number of libraries containing the indel and a "score" for the indel.
The score is calculated as the product of "strength", size, and the number
of libraries where strength is a function of the number of reads supporting
the indel.
EOS
}

sub _variant_type { 'indels' };

sub _filter_variants {

    # for each row in the indel output file of somatic sniper,
    # count how many distinct libraries have the same cigar string,
    # regardless of position, etc.

    my ($self) = @_;

    # This should be a file from sniper only
    my $input_fh = Genome::Sys->open_file_for_reading($self->input_directory . "/indels.hq");

    # TODO note: write these both to the scratch dir, then at the end of the method -- if multi lib has size, it is the HQ and single lib is the LQ, otherwise single lib is the HQ and we have no LQ. TODO ... this may not be necessary or something we want or need to do
    my $multi_lib_output_file = $self->_temp_scratch_directory . "/multi_lib_output";
    my $multi_lib_output_fh = Genome::Sys->open_file_for_writing($multi_lib_output_file);
    my $single_lib_output_file = $self->_temp_scratch_directory . "/single_lib_output";
    my $single_lib_output_fh = Genome::Sys->open_file_for_writing($single_lib_output_file);
    my $lq_output_file = $self->_temp_staging_directory . "/indels.lq";
    my $lq_output_fh = Genome::Sys->open_file_for_writing($lq_output_file);
    my $hq_output_file = $self->_temp_staging_directory . "/indels.hq";

    while ( my $line = $input_fh->getline() ) {
        chomp $line;
        my @fields = split /\t/, $line;

        my $chr    = $fields[0];
        my $pos    = $fields[1];
        my $indel1 = $fields[4];
        my $indel2 = $fields[5];
        my $indel1_size = $fields[6];
        my $indel2_size = $fields[7];
        my $tumor_reads_indel1 = $fields[13];
        my $tumor_reads_indel2 = $fields[14];
        my $normal_reads_indel1 = $fields[26]; # number of reads that support indel 1 in normal
        my $normal_reads_indel2 = $fields[27]; # number of reads that support indel 2 in normal
        my $num_libs_match=0;
        my $indel_size=0;
        my $indel_strength=0;
        unless($indel1 eq '*' || $indel2 eq '*') {
            #skip double indel case. probably ref error.
            $lq_output_fh->print($line . "\n");
            next;
        }

        if($indel1 ne '*') {
            if($normal_reads_indel1 > 0) {
                $lq_output_fh->print($line . "\n");
                next;
            }
            if($indel1_size <= 2 && $indel1_size >= -2) {
                if($tumor_reads_indel2> 0) {
                    if($tumor_reads_indel1/$tumor_reads_indel2 < (.2/(abs $indel1_size))) {
                        $lq_output_fh->print($line . "\n");
                        next;
                    }
                }
            }
            $indel_strength = $tumor_reads_indel1/($tumor_reads_indel1 + $tumor_reads_indel2) * $tumor_reads_indel1 * $normal_reads_indel2;
            $indel_size=$indel1_size;
            if(defined $fields[34]) {
                $num_libs_match = $fields[34];
            }

        }
        elsif($indel2 ne '*') {
            if($normal_reads_indel2 > 0) {
                $lq_output_fh->print($line . "\n");
                next;
            }
            if($indel2_size <= 2 && $indel2_size >= -2) {
                if($tumor_reads_indel1 > 0) {
                    if($tumor_reads_indel2/$tumor_reads_indel1 < (.2 / (abs $indel2_size))) {
                        $lq_output_fh->print($line . "\n");
                        next;
                    }
                }
            }
            $indel_strength = $tumor_reads_indel2/($tumor_reads_indel2 + $tumor_reads_indel1)  * $tumor_reads_indel2 * $normal_reads_indel1;
            $indel_size = $indel2_size;

            if(defined $fields[35]) {
                $num_libs_match = $fields[35];
            }
        }
        else { 
            die "indel1 and indel2 ... *\n";
        }

        # number of distinct libraries that contain this indel
        if($num_libs_match< 2) {
            my $indel_score = sprintf( "%.3f", abs($indel_strength * $indel_size ));
            if($indel_score > 0) {
                print $single_lib_output_fh "$line\t$num_libs_match\t$indel_score\n";
            }    
            else {
                print $lq_output_fh $line . "\n";
            }
        }
        else {
            my $indel_score = sprintf( "%.3f", abs($indel_strength * $indel_size * $num_libs_match));
            print $multi_lib_output_fh "$line\t$num_libs_match\t$indel_score\n";
        }
    }

    $multi_lib_output_fh->close;
    $single_lib_output_fh->close;
    $lq_output_fh->close;

    # If multiple libraries existed, use only those with multiple library support for HQ, if only one library existed, the single library indels are HQ
    if (-s $multi_lib_output_file) {
        Genome::Sys->copy_file($multi_lib_output_file, $hq_output_file);
        # FIXME Since we already have lq indels, we must add the single library indels to LQ should do this in a smarter way
        unless (system("cat $single_lib_output_file >> $lq_output_file") == 0) {
            $self->error_message("Failed to add single library output to lq output");
            die $self->error_message;
        }
    } else {
        Genome::Sys->copy_file($single_lib_output_file, $hq_output_file);
    }

    $self->_generate_standard_files;

    return 1; 
}

sub _generate_standard_files {
    my $self = shift;

    my $convert = Genome::Model::Tools::Bed::Convert::Indel::SniperToBed->create( 
                source => $self->_temp_staging_directory . "/indels.hq",
                output => $self->_temp_staging_directory . "/indels.hq.bed");

    unless($convert->execute){
        $self->error_message("Failed to convert hq output to bed.");
        die $self->error_message;
    }

    my $convert_lq = Genome::Model::Tools::Bed::Convert::Indel::SniperToBed->create( 
                source => $self->_temp_staging_directory . "/indels.lq",
                output => $self->_temp_staging_directory . "/indels.lq.bed");

    unless($convert_lq->execute){
        $self->error_message("Failed to convert lq output to bed.");
        die $self->error_message;
    }

    return 1;
}

1;
