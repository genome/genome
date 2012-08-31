package Genome::Model::Tools::Gatk::SomaticIndel;

use strict;
use warnings;
use FileHandle;
use Genome;
use File::Basename qw/dirname/;
use File::Spec::Functions;

class Genome::Model::Tools::Gatk::SomaticIndel {
    is => 'Genome::Model::Tools::Gatk',
    has => [
        normal_bam => {
            is => 'Text',
            doc => "BAM File for Normal Sample",
            is_optional => 0,
            is_input => 1,
        },
        tumor_bam => {
            is => 'Text',
            doc => "BAM File for Tumor Sample",
            is_optional => 0,
            is_input => 1
        },
        output_file => {
            is => 'Text',
            doc => "Output file to receive formatted lines",
            is_optional => 0,
            is_input => 1,
            is_output => 1,
        },
        bed_output_file => {
            is => 'Text',
            doc => "Optional abbreviated output in BED format",
            is_optional => 1,
            is_input => 1,
            is_output => 1,
        },
        formatted_file => {
            is => 'Text',
            doc => "Optional output file of indels in annotation format",
            is_optional => 1,
            is_input => 1,
            is_output => 1,
        },
        somatic_file => {
            is => 'Text',
            doc => "Optional output file for Somatic indels parsed from formatted file",
            is_optional => 1,
            is_input => 1,
            is_output => 1,
        },
        gatk_params => {
            is => 'Text',
            doc => "Parameters for GATK",
            is_optional => 1,
            is_input => 1,
            is_output => 1,
            default => "-T IndelGenotyperV2 --somatic --window_size 300 -et NO_ET",
        },
        reference => {
            is => 'Text',
            doc => "Parameters for GATK",
            is_optional => 1,
            is_input => 1,
            default => "/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa",
        },
        mb_of_ram => {
            is => 'Text',
            doc => 'The amount of RAM to use, in megabytes',
            default => 5000,
        },
        skip_if_output_present => {
            is => 'Text',
            doc => "Skip if output is present",
            is_optional => 1,
            is_input => 1,
        },
    ],
    has_param => [
        lsf_resource => {
            default => "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>1000] span[hosts=1] rusage[mem=6000, tmp=1000]' -M 6000000",
        }
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Runs the GATK somatic indel detection pipeline"
}

sub help_synopsis {
return <<EOS
        This command runs the GATK somatic indel detection pipeline
        EXAMPLE:	gmt gatk somatic-indel --normal-bam Normal.bam --tumor-bam Tumor.bam --output-file GATK.indel --bed-output GATK.indel.bed
EOS
}

sub help_detail {
return <<EOS

EOS
}

sub execute {
    my $self = shift;

    ## Run GATK ##
    my $path_to_gatk = $self->gatk_path;
    my $gatk_params = $self->gatk_params;
    my $reference = $self->reference;

    ## Add reference to GATK params ##
    $gatk_params = "-R $reference " . $gatk_params;
    #-I /gscmnt/sata905/info/model_data/2858219475/build103084961/alignments/103084961_merged_rmdup.bam
    #-I /gscmnt/sata871/info/model_data/2858334303/build103084933/alignments/103084933_merged_rmdup.bam
    #-O gatk_testing/indels.GATK.H_GP-13-0890-01A-01-1.tsv -o gatk_testing/indels.GATK.H_GP-13-0890-01A-01-1.out 

    my $output_file = $self->output_file;
    my $vcf_output_file = $output_file . ".vcf";
    my $ram = $self->mb_of_ram;
    my $cmd = 'java -Xms'.$ram.'m -Xmx'.$ram.'m -jar ';
    my @args = ($path_to_gatk, $gatk_params, "-I:normal", $self->normal_bam, "-I:tumor", $self->tumor_bam, "--verboseOutput", $output_file, "--out", $vcf_output_file);
    my $whitelist_args = _infer_whitelist_args($reference);
    push(@args, $whitelist_args) if $whitelist_args;
    $cmd .= join(" ", @args);

    ## Optionally append BED output file ##

    my $bed_output_file = $self->output_file.".bed";
    if ($self->bed_output_file) {
        $bed_output_file = $self->bed_output_file;
    }

    $cmd .= " --bedOutput $bed_output_file";

    ## Run GATK Command ##
    my $return;
    if ($self->skip_if_output_present && -e $output_file) {

    }
    else {
        system("touch $output_file"); # This will create an empty output file to help prevent GATK from crashing 
        system("touch $bed_output_file"); # This will create an empty output file to help prevent GATK from crashing 
        if ($self->somatic_file) {
            system("touch " . $self->somatic_file);
        }

        $return = Genome::Sys->shellcmd(
                cmd => "$cmd",
                output_files => [$output_file],
                skip_if_output_is_present => 0,
                allow_zero_size_output_files => 1,
        );
        unless($return) {
            $self->error_message("Failed to execute GATK: GATK Returned $return");
            die $self->error_message;
        }

        unless(-s $output_file) {
            $self->warning_message('No indels returned by GATK.');
        }
    }

    if ($self->formatted_file) {
        my $formatted_output_file = $self->formatted_file;

        ## Format GATK Indels ##

        if ($self->skip_if_output_present && -s $formatted_output_file) {
        }
        else {
            print "Formatting indels for annotation...\n";

            my $cmd_obj = Genome::Model::Tools::Gatk::FormatIndels->create(
                variants_file => $output_file,
                output_file => $formatted_output_file,
            );
            $cmd_obj->execute;
        }
        if ($self->somatic_file) {
            if ($self->skip_if_output_present && -s $self->somatic_file) {

            }
            else {
                print "Parsing out Somatic indels...\n";
                ## Parse the results to the somatic output file ##
                parse_somatic($formatted_output_file, $self->somatic_file);
            }
        }
    }
    return 1;
}


sub _infer_whitelist_args {
    my $reference = shift;

    my $seqdict = catfile(dirname($reference), "seqdict", "seqdict.sam");
    return unless -f $seqdict;
    my @sequences;
    open(FH, "<$seqdict") or return;
    while (<FH>) {
        next unless /^\@SQ/;
        chomp;
        my @fields = split(/[\t:]/);
        push(@sequences, $fields[2] . ':1-' . $fields[4] . "\n");
    }
    close(FH);

    my $whitelist_file = Genome::Sys->create_temp_file_path() . '.interval_list';
    Genome::Sys->write_file($whitelist_file, @sequences);
    return "-L $whitelist_file";
}


################################################################################################
# Parse_Somatic - isolate somatic indels 
#
################################################################################################

sub parse_somatic
{
    my $FileName = shift(@_);
    my $OutFileName = shift(@_);

    open(OUTFILE, ">$OutFileName") or die "Can't open outfile: $!\n";

    ## Parse the variants file ##

    my $input = new FileHandle ($FileName);
    my $lineCounter = 0;

    while (<$input>) {
        chomp;
        my $line = $_;
        $lineCounter++;

        my @lineContents = split(/\t/, $line);
        #my $somatic_status = $lineContents[17];

        if (($lineContents[16] && $lineContents[16] =~ 'SOMATIC') || ($lineContents[17] && $lineContents[17] =~ 'SOMATIC')) {
            print OUTFILE "$line\n";
        }
        else {

        }
    }

    close($input);
    close(OUTFILE);
}

1;

