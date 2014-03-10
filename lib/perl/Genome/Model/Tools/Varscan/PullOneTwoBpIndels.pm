package Genome::Model::Tools::Varscan::PullOneTwoBpIndels;

use strict;
use warnings;

use FileHandle;
use File::Basename;
use Genome;
use Genome::Info::IUB;

class Genome::Model::Tools::Varscan::PullOneTwoBpIndels {
    is => 'Command',
    doc => "Generate lists of 1-2 bp and 3+ indels, run GATK recalibration, and then sort and index bams",
    has => [
        varscan_version => {
            is => 'Text',
            doc => "Varscan version to use when running varscan validation",
            is_input => 1,
            is_optional => 1,
        },
        project_name => { 
            is => 'Text',
            doc => "Name of the project i.e. ASMS",
            is_optional => 1
        },
        list_of_indel_files_to_validate => { 
            is => 'Text',
            doc => "A file listing indel files to include in the validation process, 1 filename per line, no headers. Files should either be named ending in '.bed' or be in annotation format: chr start stop ref var... If you specify a somatic validation build and the variants/indels.hq.bed file is not in your file list, those indels will be included in the validation process in addition to your file list. If any file on your list ends with '.bed', then it will be converted to annotation format for processing via Varscan, etc, which expect 1-based start positions (not bed).",
            is_optional => 1,
        },
        small_indel_output_bed => {
            is => 'Text',
            doc => "File of small indels to be realigned. BED format - must be named *.bed! Note that the '.bed' output version of this file has padded start and stop to allow more robust local realignment. The '.annotation_format' output version has true coordinates.",
            is_optional => 1,
        },
        large_indel_output_bed => {
            is => 'Text',
            doc => "File of large indels to be processed using other tools. BED format - must be named *.bed!",
            is_optional => 1,
        },
        tumor_bam => {
            is => 'Text',
            doc => "Tumor Bam File (Validation Bam)",
            is_optional => 1,
        },
        normal_bam => {
            is => 'Text',
            doc => "Normal Bam File (Validation Bam)",
            is_optional => 1,
        },
        relapse_bam => {
            is => 'Text',
            doc => "(Optional) Relapse Bam File (Validation Bam)",
            is_optional => 1,
        },
        reference_fasta => {
            is => 'Text',
            doc => "Reference Fasta",
            is_optional => 1,
        },
        varscan_indel_output => {
            is => 'Text',
            doc => "gmt varscan validation output run on realigned bams",
            is_optional => 1,
        },
        varscan_snp_output => {
            is => 'Text',
            doc => "gmt varscan validation output run on realigned bams",
            is_optional => 1,
        },
        final_output_file => {
            is => 'Text',
            doc => "gmt varscan process-validation-indels final output file labeling indels as Somatic or otherwise",
            is_output => 1,
            is_optional => 1,
        },
        skip_if_output_present => {
            is => 'Boolean',
            doc => "Skip Creating new Bam Files if they exist",
            is_optional => 1,
            default => "",
        },
        realigned_bam_file_directory => {
            is => 'Text',
            doc => "Where to dump the realigned bam file",
            is_optional => 1,
        },
        normal_purity => {
            is => 'Float',
            doc => "Normal purity param to pass to varscan",
            is_optional => 0,
            default => 1,
        },
        min_var_frequency => {
            is => 'Float',
            doc => "Minimum variant frequency to pass to varscan",
            is_optional => 0,
            default => 0.08,
        },
        somatic_validation_build => {
            is => 'Text',
            doc => "Somatic validation build to use. If this is set, defaults will be set inside that build directory for final_output_file,realigned_bam_file_directory,small_indel_output_bed,large_indel_output_bed.",
            is_optional => 1,
            is_input => 1,
            is_output => 1,
        },
    ],
};

sub help_synopsis {
    return <<EOS
        Generate lists of 1-2 bp and 3+ indels, run GATK recalibration, and then sort and index bams. Run varscan validation on the realigned bams, and then run varscan process-validation-indels on the output.
EOS
}

sub help_detail {
    return <<EOS
                Generate lists of 1-2 bp and 3+ indels, run GATK recalibration, and then sort and index bams. Run varscan validation on the realigned bams, and then run varscan process-validation-indels on the output.
EOS
}

##############################################################################################
# convert_indel_bed_file - for converting to anno format and returning a temp file path
#
##############################################################################################
sub convert_anno_indel_file_to_bed {
    my $self = shift;
    my $anno_file = shift;
    my $bed_file = Genome::Sys->create_temp_file_path;
    my %convert_params = (source => $anno_file, output => $bed_file);
    my $convert_class = "Genome::Model::Tools::Bed::Convert::AnnotationToBed";
    $convert_class->execute(%convert_params) || ($self->error_message("Could not convert annotation file to BED format: $anno_file") and return);
    return $bed_file;
}

################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {
    my $self = shift;

    # Make sure that if output paths arent set that the somatic variation build is, and set good defaults
    if ($self->somatic_validation_build) {
        my $build = Genome::Model::Build->get($self->somatic_validation_build);
        die $self->error_message("Could not get a build for id " . $self->somatic_validation_build) unless ($build);

        unless($build->normal_sample) {
            $self->status_message('No normal sample found.  Skipping.');
            return 1;
        }

        my $base_dir = $build->data_directory . "/validation/small_indel";
        Genome::Sys->create_directory($base_dir);
        unless ($self->final_output_file) {
            $self->final_output_file($base_dir."/final_output");
        }
        unless ($self->realigned_bam_file_directory) {
            my $realigned_dir = "$base_dir/realigned_bams";
            $self->realigned_bam_file_directory($realigned_dir);
        }
        unless ($self->small_indel_output_bed) {
            $self->small_indel_output_bed("$base_dir/small_indels.bed");
        }
        unless ($self->large_indel_output_bed) {
            $self->large_indel_output_bed("$base_dir/large_indels.bed");
        }
        unless ($self->list_of_indel_files_to_validate) {
            $self->list_of_indel_files_to_validate("$base_dir/indel_files_to_validate");
        }
        unless ($self->varscan_indel_output) {
            $self->varscan_indel_output("$base_dir/varscan_indels");
        }
        unless ($self->varscan_snp_output) {
            $self->varscan_snp_output("$base_dir/varscan_snps");
        }

        # Set the varscan version
        if ($self->varscan_version) {
            die $self->error_message("Please set varscan_version or somatic_variation_build but not both.");
        } else {
            my $pp = $build->processing_profile;
            my $varscan_version = $pp->varscan_validation_version;
            unless ($varscan_version) {
                die $self->error_message("Could not get varscan version from processing profile " . $pp->id);
            }
            $self->varscan_version($varscan_version);
        }
    } else {
        my @required_properties = qw(final_output_file realigned_bam_file_directory small_indel_output_bed list_of_indel_files_to_validate varscan_indel_output varscan_snp_output tumor_bam normal_bam reference_fasta varscan_version);
        my $fail = 0;
        for my $property (@required_properties) {
            unless (defined $self->$property) {
                $fail = 1;
                $self->error_message("$property is not set and must be if somatic_validation_build is not set");
            }
        }
        die $self->error_message("All of the above properties must be set unless somatic_validation_build is set.") if $fail;
    }

    my $project_name = $self->project_name;
    my $small_indel_list = $self->small_indel_output_bed;
    my $large_indel_list = $self->large_indel_output_bed;
    my $file_list_file = $self->list_of_indel_files_to_validate;
    my $normal_bam = $self->normal_bam;
    my $tumor_bam = $self->tumor_bam;
    my $reference = $self->reference_fasta;
    my $output_indel = $self->varscan_indel_output;
    my $output_snp = $self->varscan_snp_output;
    my $final_output_file = $self->final_output_file;
    my $skip_if_output_present = $self->skip_if_output_present;
    my $somatic_validation_build = $self->somatic_validation_build;
    Genome::Sys->create_directory($self->realigned_bam_file_directory);

    # check small/large indel list filename for bed nomenclature (to eliminate future confusion) #
    unless ($small_indel_list =~ m/\.bed$/i && $large_indel_list =~ m/\.bed$/i) {
        $self->error_message("Both small and large indel files must end in .bed");
        return;
    }

    # process somatic validation build #
    if (defined($somatic_validation_build)) {

        #get some parameters from the build #
        my $build = Genome::Model::Build->get($somatic_validation_build);
        my $model = Genome::Model->get($build->model_id);
        my $ref_seq_build_id = $model->reference_sequence_build->build_id;
        my $ref_seq_build = Genome::Model::Build->get($ref_seq_build_id);
        $reference = $ref_seq_build->full_consensus_path('fa');
        $normal_bam = $build->normal_bam;
        $tumor_bam = $build->tumor_bam;
        my $build_dir = $build->data_directory;
        my $val_build_indels_file = $build_dir . "/variants/indels.hq.bed";

        # check to see if the build's indel file is in the file list #
        if (-s $file_list_file) {
            my $contains_som_val_build_file = 0;
            my $file_list_fh = new IO::File $file_list_file,"r";
            while (my $line = $file_list_fh->getline) {
                chomp $line;
                if ($line eq $val_build_indels_file) { $contains_som_val_build_file++; }
            }
            $file_list_fh->close;

            # add the build indel file to file list if needed
            unless ($contains_som_val_build_file) {
                my $file_list_append_fh = new IO::File $file_list_file,">>";
                print $file_list_append_fh $val_build_indels_file . "\n";
                $file_list_append_fh->close;
            }
        } else {
            my $file_list_append_fh = new IO::File $file_list_file,">>";
            if (-s $val_build_indels_file) {
                print $file_list_append_fh $val_build_indels_file . "\n";
            }
            if ($build->indel_variant_list) {
                my $indel_variant_list = $build->indel_variant_list->output_dir . "/indels.hq.bed";
                if (-s $indel_variant_list) {
                    print $file_list_append_fh $indel_variant_list . "\n";
                }
            }
            $file_list_append_fh->close;
        }
    }

    # else if there is no som-val build defined - require other params #
    else {
        unless (defined($file_list_file) && defined($tumor_bam) && defined($normal_bam) && defined($reference)) {
            $self->error_message("If a somatic-validation build is not defined, you must provide the following parameters:\n  file_list\n  tumor-bam\n  normal-bam\n  reference-fasta\n  output-indel\n  output-snp");
            return;
        }                    
    }

    # define output filenames #
    (my $small_indel_list_nobed = $small_indel_list) =~ s/\.bed$/\.annotation_format/;
    $small_indel_list =~ s/\.bed$/.padded1bp.bed/;
    my $realigned_bam_file_directory = $self->realigned_bam_file_directory;

    my $realigned_normal_bam_file = basename($normal_bam,qr{\.bam});
    $realigned_normal_bam_file = "$realigned_bam_file_directory/$realigned_normal_bam_file.realigned.bam";

    my $realigned_tumor_bam_file = basename($tumor_bam,qr{\.bam});
    $realigned_tumor_bam_file = "$realigned_bam_file_directory/$realigned_tumor_bam_file.realigned.bam";

    my $relapse_bam;
    my $realigned_relapse_bam_file;
    if ($self->relapse_bam){
        $relapse_bam = $self->relapse_bam;
        $realigned_relapse_bam_file = basename($relapse_bam,qr{\.bam});
        $realigned_relapse_bam_file = "$realigned_bam_file_directory/$realigned_relapse_bam_file.realigned.bam";

    }

    # open output filehandles #
    open(INDELS_OUT_BED, ">$small_indel_list") or die "Can't open small indel output bed file: $!\n";
    open(INDELS_OUT_ANNO, ">$small_indel_list_nobed") or die "Can't open small indel output annotation_format file: $!\n";
    if($large_indel_list) { open(LARGE_INDELS_OUT_BED, ">$large_indel_list") or die "Can't open large indel output file: $!\n"; }

    my $file_input = new FileHandle ($file_list_file);
    unless($file_input) {
        $self->error_message("Unable to open file list: $file_list_file");
        return;
    }

    # read through the file list and print each files' indels into the appropriate size-specific output file #
    while (my $file = <$file_input>) {
        chomp($file);

        # unless file is already in .bed format, convert it to bed format before further processing #
        unless ($file =~ m/\.bed$/i) { $file = $self->convert_anno_indel_file_to_bed($file); }

        # open filehandle for processing indels in $file #
        my $indel_input = new IO::File $file,"r";
        unless($indel_input) {
            $self->error_message("Unable to open indel-containing file: $file");
            return;
        }

        # process indels in indel file #
        while (my $line = <$indel_input>) {
            chomp($line);
            my ($chr, $start, $stop, $ref_var, @everything_else) = split(/\t/, $line);
            my $size;
            my $annostart;
            my $annostop;
            my $type;
            my ($ref,$var) = split(/\//,$ref_var);
            if ($ref eq '-' || $ref eq '0' || $ref eq '*') { #ins
                #count number of bases inserted
                $ref = "-";
                $size = length($var);
                $annostart = ($start);
                $annostop = ($stop + 1);
                $type = 'INS';
            }
            elsif ($var eq '-' || $var eq '0' || $var eq '*') { #del
                $var = "-";
                $size = length($ref);
                $annostart = ($start + 1);
                $annostop = ($stop);
                $type = 'DEL';
            }
            else {
                print "ERROR: Line $line in file $file has wrong insertion or deletion nomenclature. Either ref or var should be 0 or -\n";
                $size = 0;  #this will include this indel despite its wrongness
            }
            if ( $size > 0 && $size <= 2) {
                #Add 1 bp padding to bed because we just want to look at regions
                $start-=1;
                $stop+=1;
                print INDELS_OUT_BED "$chr\t$start\t$stop\t$ref/$var\n";
                print INDELS_OUT_ANNO "$chr\t$annostart\t$annostop\t$ref\t$var\n";
            }
            elsif ( $size > 2 && $large_indel_list) {
                print LARGE_INDELS_OUT_BED "$chr\t$start\t$stop\t$ref/$var\t$type\n";
            }
        }
        close($indel_input);
    }
    close($file_input);

    # set up realignment jobs and varscan indel calling jobs using input files defined above #
    my $min_freq = $self->min_var_frequency;
    my $normal_purity = $self->normal_purity;
    my $varscan_params = "--validation 1 --somatic-p-value 1.0e-02 --p-value 0.10 --min-coverage 8 --min-var-freq $min_freq --normal-purity $normal_purity";
    my $default_varscan_params = "--validation 1 --somatic-p-value 1.0e-02 --p-value 0.10 --min-coverage 8 --min-var-freq 0.08 --normal-purity 1";
    my $bsub = qq(bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -R "select[type==LINUX64 && mem>16000 && tmp>10000] rusage[mem=16000, tmp=10000]" -M 16000000 );

    # put several jobs into array @cmds
    my @cmds;
    my $user = $ENV{USER};
    if ($skip_if_output_present && -e $realigned_normal_bam_file && -e $realigned_tumor_bam_file) {
        if ($self->relapse_bam){
            push(@cmds,"$bsub -J varscan_validation_tumnor \'gmt varscan validation --normal-bam $realigned_normal_bam_file --tumor-bam $realigned_tumor_bam_file --output-indel $output_indel.tumnor --output-snp $output_snp.tumnor --reference $reference --varscan-params \"$varscan_params\"\'");
            push(@cmds, "$bsub -J varscan_validation_relnor \'gmt varscan validation --normal-bam $realigned_normal_bam_file --tumor-bam $realigned_relapse_bam_file --output-indel $output_indel.relnor --output-snp $output_snp.relnor --reference $reference --varscan-params \"$varscan_params\"\'");
            push(@cmds,"$bsub -J varscan_validation_reltum \'gmt varscan validation --normal-bam $realigned_tumor_bam_file --tumor-bam $realigned_relapse_bam_file --output-indel $output_indel.reltum --output-snp $output_snp.reltum --reference $reference --varscan-params \"$varscan_params\"\'");
            push(@cmds,"$bsub -N -u $user\@$ENV{GENOME_EMAIL_DOMAIN} -J varscan_process_validation_tumnor -w \'ended(JOB0))\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel.tumnor --validation-snp-file $output_snp.tumnor --variants-file $small_indel_list_nobed --output-file $final_output_file.tumnor\'");
            push(@cmds,"$bsub -N -u $user\@$ENV{GENOME_EMAIL_DOMAIN} -J varscan_process_validation_relnor -w \'ended(JOB1)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel.relnor --validation-snp-file $output_snp.relnor --variants-file $small_indel_list_nobed --output-file $final_output_file.relnor\'");
            push(@cmds,"$bsub -N -u $user\@$ENV{GENOME_EMAIL_DOMAIN} -J varscan_process_validation_reltum -w \'ended(JOB2)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel.reltum --validation-snp-file $output_snp.reltum --variants-file $small_indel_list_nobed --output-file $final_output_file.reltum\'");
        }
        else {
            push(@cmds,"$bsub -J varscan_validation \'gmt varscan validation --normal-bam $realigned_normal_bam_file --tumor-bam $realigned_tumor_bam_file --output-indel $output_indel --output-snp $output_snp --reference $reference --varscan-params \"$varscan_params\"\'");

            push(@cmds,"$bsub -N -u $user\@$ENV{GENOME_EMAIL_DOMAIN} -J varscan_process_validation -w \'ended(JOB0)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel --validation-snp-file $output_snp --variants-file $small_indel_list_nobed --output-file $final_output_file\'");
        }
    }
    elsif ($self->relapse_bam) {
        my $bsub_normal_output = "$realigned_bam_file_directory/realignment_normal.out";
        my $bsub_normal_error = "$realigned_bam_file_directory/realignment_normal.err";
        push(@cmds,"$bsub -J $realigned_normal_bam_file -o $bsub_normal_output -e $bsub_normal_error \'java -Xmx16g -Djava.io.tmpdir=/tmp -jar $ENV{GENOME_SW}/gatk/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -et NO_ET -T IndelRealigner -targetIntervals $small_indel_list -o $realigned_normal_bam_file -I $normal_bam -R $reference  --targetIntervalsAreNotSorted\'");

        my $bsub_tumor_output = "$realigned_bam_file_directory/realignment_tumor.out";
        my $bsub_tumor_error = "$realigned_bam_file_directory/realignment_tumor.err";
        push(@cmds,"$bsub -J $realigned_tumor_bam_file -o $bsub_tumor_output -e $bsub_tumor_error \'java -Xmx16g -Djava.io.tmpdir=/tmp -jar $ENV{GENOME_SW}/gatk/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -et NO_ET -T IndelRealigner -targetIntervals $small_indel_list -o $realigned_tumor_bam_file -I $tumor_bam -R $reference --targetIntervalsAreNotSorted\'");

        my $bsub_relapse_output = "$realigned_bam_file_directory/realignment_relapse.out";
        my $bsub_relapse_error = "$realigned_bam_file_directory/realignment_relapse.err";
        push(@cmds,"$bsub -J $realigned_relapse_bam_file -o $bsub_relapse_output -e $bsub_relapse_error \'java -Xmx16g -Djava.io.tmpdir=/tmp -jar $ENV{GENOME_SW}/gatk/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -et NO_ET -T IndelRealigner -targetIntervals $small_indel_list -o $realigned_relapse_bam_file -I $relapse_bam -R $reference --targetIntervalsAreNotSorted\'");
        push(@cmds,"$bsub -J bamindex_normal -w \'ended(JOB0)\' \'samtools index $realigned_normal_bam_file\'");
        push(@cmds,"$bsub -J bamindex_tumor -w \'ended(JOB1)\' \'samtools index $realigned_tumor_bam_file\'");
        push(@cmds,"$bsub -J bamindex_relapse -w \'ended(JOB2)\' \'samtools index $realigned_relapse_bam_file\'");

        push(@cmds,"$bsub -J varscan_validation_tumnor -w \'ended(JOB3) && ended(JOB4)\' \'gmt varscan validation --normal-bam $realigned_normal_bam_file --tumor-bam $realigned_tumor_bam_file --output-indel $output_indel.tumnor --output-snp $output_snp.tumnor --reference $reference --varscan-params \"$varscan_params\"\'");
        push(@cmds,"$bsub -J varscan_validation_relnor -w \'ended(JOB3) && ended(JOB5)\' \'gmt varscan validation --normal-bam $realigned_normal_bam_file --tumor-bam $realigned_relapse_bam_file --output-indel $output_indel.relnor --output-snp $output_snp.relnor --reference $reference --varscan-params \"$default_varscan_params\"\'");
        push(@cmds,"$bsub -J varscan_validation_reltum -w \'ended(JOB4) && ended(JOB5)\' \'gmt varscan validation --normal-bam $realigned_tumor_bam_file --tumor-bam $realigned_relapse_bam_file --output-indel $output_indel.reltum --output-snp $output_snp.reltum --reference $reference --varscan-params \"$default_varscan_params\"\'");

        push(@cmds,"$bsub -N -u $user\@$ENV{GENOME_EMAIL_DOMAIN} -J varscan_process_validation_tumnor -w \'ended(JOB6)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel.tumnor --validation-snp-file $output_snp.tumnor --variants-file $small_indel_list_nobed --output-file $final_output_file.tumnor\'");
        push(@cmds,"$bsub -N -u $user\@$ENV{GENOME_EMAIL_DOMAIN} -J varscan_process_validation_relnor -w \'ended(JOB7)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel.relnor --validation-snp-file $output_snp.relnor --variants-file $small_indel_list_nobed --output-file $final_output_file.relnor\'");
        push(@cmds,"$bsub -N -u $user\@$ENV{GENOME_EMAIL_DOMAIN} -J varscan_process_validation_reltum -w \'ended(JOB8)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel.reltum --validation-snp-file $output_snp.reltum --variants-file $small_indel_list_nobed --output-file $final_output_file.reltum\'");
    }
    else{
        # This was added so things are testable for now
        my $cmd = Genome::Model::SomaticValidation::Command::ValidateSmallIndels->create(
            varscan_version => $self->varscan_version,
            final_output_file => $self->final_output_file,
            realigned_bam_file_directory => $self->realigned_bam_file_directory,
            small_indel_output_bed => $small_indel_list,
            varscan_indel_output => $self->varscan_indel_output,
            varscan_snp_output => $self->varscan_snp_output,
            tumor_bam => $tumor_bam,
            normal_bam => $normal_bam,
            reference_fasta => $reference,
        );
        unless ($cmd->execute) {
            die $self->error_message("Failed to execute Genome::Model::SomaticValidation::Command::ValidateSmallIndels");
        }
    }

    if (@cmds) {
        # run jobs in @cmds #
        my @jobids;
        foreach my $cmd (@cmds){        
            # waiting on two jobs #
            if($cmd =~ /JOB(\d).+JOB(\d)/){
                my $j1 = $1;
                my $j2 = $2;
                my $jobid1 = $jobids[$j1];
                my $jobid2 = $jobids[$j2];
                $cmd =~ s/JOB$j1/$jobid1/g;
                $cmd =~ s/JOB$j2/$jobid2/g;            

                # waiting on one job #
            } elsif($cmd =~ /JOB(\d)/){
                my $job = $1;
                my $jobid = $jobids[$job];
                $cmd =~ s/JOB$job/$jobid/g;
            }           

            print STDERR "\nRunning command:  $cmd\n";
            my $id = `$cmd`;
            if($id =~ /<(\d+)>/){
                $id = $1;
            } else {
                die "job not submitted correctly\n";
            }
            push(@jobids,$id);
        }    

        print STDERR "\njob ids: " . join(" ",@jobids) . "\n";
    }

    return 1;
}
