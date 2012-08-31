package Genome::Model::Tools::Varscan::PullOneTwoBpIndels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# PullOneTwoBpIndels - Generate lists of 1-2 bp and 3+ indels, run GATK recalibration, and then sort and index bams
#
#       AUTHOR:         Will Schierding (wschierd@genome.wustl.edu)
#
#       CREATED:        11/29/2010 by W.S.
#       MODIFIED:       11/29/2010 by W.S.
#
#       NOTES:
#
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;
use File::Basename; #for file name parsing
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it
use Genome::Info::IUB;

class Genome::Model::Tools::Varscan::PullOneTwoBpIndels {
    is => 'Command',

    has => [                               
# specify the command's single-value properties (parameters) <---
        project_name => { 
            is => 'Text',
            doc => "Name of the project i.e. ASMS" ,
            is_optional => 1},
        
        file_list => { 
            is => 'Text',
            doc => "File of indel files to include,  1 name per line, tab delim, no headers. Should be chr start stop ref var blahhhhhhhhhhhh." ,
            is_optional => 1},

        small_indel_outfile => {
            is => 'Text',
            doc => "File of small indels to be realigned" ,
            is_optional => 0},

        large_indel_outfile => {
            is => 'Text',
            doc => "File of large indels to be realigned" ,
            is_optional => 1},

        tumor_bam   => {
            is => 'Text',
            doc => "Tumor Bam File (Validation Bam)" ,
            is_optional => 1},

        normal_bam  => {
            is => 'Text',
            doc => "Normal Bam File (Validation Bam)" ,
            is_optional => 1},

        relapse_bam => {
            is => 'Text',
            doc => "(Optional) Relapse Bam File (Validation Bam)" ,
            is_optional => 1},

        reference_fasta => {
            is => 'Text',
            doc => "Reference Fasta" ,
            is_optional => 1,
            #default => "/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/all_sequences.fa"},
        },

        output_indel => {
            is => 'Text',
            doc => "gmt varscan validate input" ,
            is_optional => 1},

        output_snp  => {
            is => 'Text',
            doc => "gmt varscan validate input" ,
            is_optional => 1},

        final_output_file   => {
            is => 'Text',
            doc => "process-validation-indels output file" ,
            is_optional => 0},

        skip_if_output_present      => {
            is => 'Boolean',
            doc => "Skip Creating new Bam Files if they exist" ,
            is_optional => 1,
            default => ""},

        realigned_bam_file_directory => {
            is => 'Text',
            doc => "Where to dump the realigned bam file",
            is_optional => 0},

        normal_purity => {
            is => 'Float',
            doc => "Normal purity param to pass to varscan",
            is_optional => 0,
            default => 1},

        min_var_frequency => {
            is => 'Float',
            doc => "Minimum variant frequency to pass to varscan",
            is_optional => 0,
            default => 0.08},
        
        somatic_validation_build => {
            is => 'Integer',
            doc => "directory for ",
            is_optional => 1,
        }

        ],    
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Generate lists of 1-2 bp and 3+ indels, run GATK recalibration, and then sort and index bams"
}

sub help_synopsis {
    return <<EOS
        Generate lists of 1-2 bp and 3+ indels, run GATK recalibration, and then sort and index bams
EXAMPLE:        gmt varscan pull-one-two-bp-indels
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS
                Generate lists of 1-2 bp and 3+ indels, run GATK recalibration, and then sort and index bams
EXAMPLE:        gmt varscan pull-one-two-bp-indels
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
    my $self = shift;
    my $project_name = $self->project_name;
    my $small_indel_list = $self->small_indel_outfile;
    my $large_indel_list = $self->large_indel_outfile;
    my $file_list_file = $self->file_list;
    my $normal_bam = $self->normal_bam;
    my $tumor_bam = $self->tumor_bam;
    my $reference = $self->reference_fasta;
    my $output_indel = $self->output_indel;
    my $output_snp = $self->output_snp;
    my $final_output_file = $self->final_output_file;
    my $skip_if_output_present = $self->skip_if_output_present;
    my $somatic_validation_build = $self->somatic_validation_build;

    if(defined($somatic_validation_build)){

        #get some parameters from the build
        my $build = Genome::Model::Build->get($somatic_validation_build);
        my $build_dir = $build->data_directory;
        my $model = Genome::Model->get($build->model_id);
        my $refseq_build_id = $model->reference_sequence_build->build_id;
        my $refseq_dir = Genome::Model::Build::ReferenceSequence->get($refseq_build_id)->data_directory;        
        $reference = $refseq_dir . "/all_sequences.fa";

        $normal_bam = $build->normal_bam;
        $tumor_bam = $build->tumor_bam;


        #convert the bed file to unsplit the ref/var alleles and convert to annotation format

        #create a tmp file for the converted file
        my ($tfh,$newfile) = Genome::Sys->create_temp_file;
        unless($tfh) {
            $self->error_message("Unable to create temporary file $!");
            die;
        }
        open(OUTFILE,">$newfile") || die "can't open temp file for writing ($newfile)\n";

        #open the bed file, do the conversion
        my $inFh = IO::File->new( $build_dir . "/variants/indels.hq.bed" ) || die "can't open indels file\n";
        while( my $line = $inFh->getline )
        {
            chomp($line);
            my @F = split("\t",$line);
            $F[3] =~ s/\*/-/g;
            $F[3] =~ s/0/-/g;
            my @a = split(/\//,$F[3]);

            if (($F[3] =~ /^0/) || ($F[3] =~ /^\-/)){ #indel INS
                $F[2] = $F[2]+1;
                print OUTFILE join("\t",($F[0],$F[1],$F[2],$a[0],$a[1]));;
            } elsif (($F[3] =~ /0$/) || ($F[3] =~ /\-$/)){ #indel DEL
                $F[1] = $F[1]+1;
                print OUTFILE join("\t",($F[0],$F[1],$F[2],$a[0],$a[1]));;
            } else { #SNV
                $F[1] = $F[1]+1;
                print OUTFILE join("\t",($F[0],$F[1],$F[2],$a[0],Genome::Info::IUB::variant_alleles_for_iub($a[0],$a[1])));;
            }
            
            if(@F > 3){
                print OUTFILE "\t" . join("\t",@F[4..$#F])
            }
            print OUTFILE "\n";
        }
        close(OUTFILE);
        

        #create a tmp file for the file list
        my ($tfh2,$newfile2) = Genome::Sys->create_temp_file;
        unless($tfh2) {
            $self->error_message("Unable to create temporary file $!");
            die;
        }
        open(OUTFILE,">$newfile2") || die "can't open temp file for writing ($newfile)\n";
        print OUTFILE $newfile . "\n";
        close(OUTFILE);    
        $file_list_file = $newfile2;
        
        
    } else {#no som-val build - require all these other params
        if(!(defined($file_list_file)) || 
           !(defined($tumor_bam)) || 
           !(defined($normal_bam)) || 
           !(defined($reference))){
            
            die("if a somatic-validation build is not defined, you must provide the following parameters:\n  file_list\n  tumor-bam\n  normal-bam\n  reference-fasta\n  output-indel\n  output-snp");
        }                    
    }



    unless ($small_indel_list =~ m/\.bed/i) {
        die "Indel File Must end in .bed (because will says so)";
    }

    my $small_indel_list_nobed = $small_indel_list;
    $small_indel_list_nobed =~ s/\.bed//;
    $small_indel_list_nobed = "$small_indel_list_nobed.txt";

    my $large_indel_list_nobed = $large_indel_list;
    if($large_indel_list) {
        $large_indel_list_nobed =~ s/\.bed//;
        $large_indel_list_nobed = "$large_indel_list_nobed.txt";
    }

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

    ## Open the outfiles ##
    my $bed_indel_outfile = $small_indel_list;
    my $nobed_indel_outfile = $small_indel_list_nobed;
    open(INDELS_OUT, ">$bed_indel_outfile") or die "Can't open output file1: $!\n";
    open(NOBED_INDELS_OUT, ">$nobed_indel_outfile") or die "Can't open output file2: $!\n";
    if($large_indel_list) {
        my $large_bed_indel_outfile = $large_indel_list;
        my $large_nobed_indel_outfile = $large_indel_list_nobed;
        open(LARGE_INDELS_OUT, ">$large_bed_indel_outfile") or die "Can't open output file3: $!\n";
        open(LARGE_NOBED_INDELS_OUT, ">$large_nobed_indel_outfile") or die "Can't open output file4: $!\n";
    }
    my $file_input = new FileHandle ($file_list_file);
    unless($file_input) {
        $self->error_message("Unable to open $file_list_file");
        return;
    }

    while (my $file = <$file_input>) {
        chomp($file);
        my $indel_input = new FileHandle ($file);
        unless($indel_input) {
            $self->error_message("Unable to open $file");
            return;
        }

        while (my $line = <$indel_input>) {
            chomp($line);
            my ($chr, $start, $stop, $ref, $var, @everything_else) = split(/\t/, $line);
            my $size;
            my $bedstart;
            my $bedstop;
            my $type;
            if ($ref =~ m/\//) {
                my $split = $ref;
                ($ref, $var) = split(/\//, $split);
            }
            if ($ref eq '-' || $ref eq '0') { #ins
                #count number of bases inserted
                $size = length($var);
                $bedstart = ($start);
                $bedstop = ($stop - 1);
                $type = 'INS';
            }
            elsif ($var eq '-' || $var eq '0') { #del
                $size = length($ref);
                $bedstart = ($start - 1);
                $bedstop = ($stop);
                $type = 'DEL';
            }
            else {
                print "ERROR: Line $line in file $file has wrong insertion or deletion nomenclature. Either ref or var should be 0 or -\n";
                $size = 0;  #this will include this indel despite its wrongness
            }
            if ( $size > 0 && $size <= 2) {
                #Add 1 bp padding to bed because we just want to look at regions
                $bedstart--;
                $bedstop++;
                print INDELS_OUT "$chr\t$bedstart\t$bedstop\t$ref\t$var\n";
                print NOBED_INDELS_OUT "$chr\t$start\t$stop\t$ref\t$var\n";
            }
            elsif ( $size > 2 && $large_indel_list) {
                print LARGE_INDELS_OUT "$chr\t$bedstart\t$bedstop\t$ref\t$var\t$type\n";
                print LARGE_NOBED_INDELS_OUT "$chr\t$start\t$stop\t$ref\t$var\t$type\n";
            }
        }
        close($indel_input);
    }
    close($file_input);

    my $min_freq = $self->min_var_frequency;
    my $normal_purity = $self->normal_purity;
    my $varscan_params = "--validation 1 --somatic-p-value 1.0e-02 --p-value 0.10 --min-coverage 8 --min-var-freq $min_freq --normal-purity $normal_purity";
    my $default_varscan_params = "--validation 1 --somatic-p-value 1.0e-02 --p-value 0.10 --min-coverage 8 --min-var-freq 0.08 --normal-purity 1";

    my $bsub = 'bsub -q long -R "select[type==LINUX64 && mem>16000 && tmp>10000] rusage[mem=16000, tmp=10000]" -M 16000000 ';


    my @cmds;
    my $user = $ENV{USER};
    if ($skip_if_output_present && -e $realigned_normal_bam_file && -e $realigned_tumor_bam_file) {
        if ($self->relapse_bam){
            push(@cmds,"$bsub -J varscan_validation_tumnor \'gmt varscan validation --normal-bam $realigned_normal_bam_file --tumor-bam $realigned_tumor_bam_file --output-indel $output_indel.tumnor --output-snp $output_snp.tumnor --varscan-params \"$varscan_params\"\'");
            push(@cmds, "$bsub -J varscan_validation_relnor \'gmt varscan validation --normal-bam $realigned_normal_bam_file --tumor-bam $realigned_relapse_bam_file --output-indel $output_indel.relnor --output-snp $output_snp.relnor --varscan-params \"$varscan_params\"\'");
            push(@cmds,"$bsub -J varscan_validation_reltum \'gmt varscan validation --normal-bam $realigned_tumor_bam_file --tumor-bam $realigned_relapse_bam_file --output-indel $output_indel.reltum --output-snp $output_snp.reltum --varscan-params \"$varscan_params\"\'");
            push(@cmds,"$bsub -N -u $user\@genome.wustl.edu -J varscan_process_validation_tumnor -w \'ended(JOB0))\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel.tumnor --validation-snp-file $output_snp.tumnor --variants-file $small_indel_list_nobed --output-file $final_output_file.tumnor\'");
            push(@cmds,"$bsub -N -u $user\@genome.wustl.edu -J varscan_process_validation_relnor -w \'ended(JOB1)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel.relnor --validation-snp-file $output_snp.relnor --variants-file $small_indel_list_nobed --output-file $final_output_file.relnor\'");
            push(@cmds,"$bsub -N -u $user\@genome.wustl.edu -J varscan_process_validation_reltum -w \'ended(JOB2)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel.reltum --validation-snp-file $output_snp.reltum --variants-file $small_indel_list_nobed --output-file $final_output_file.reltum\'");
        }
        else {
            push(@cmds,"$bsub -J varscan_validation \'gmt varscan validation --normal-bam $realigned_normal_bam_file --tumor-bam $realigned_tumor_bam_file --output-indel $output_indel --output-snp $output_snp --varscan-params \"$varscan_params\"\'");

            push(@cmds,"$bsub -N -u $user\@genome.wustl.edu -J varscan_process_validation -w \'ended(JOB0)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel --validation-snp-file $output_snp --variants-file $small_indel_list_nobed --output-file $final_output_file\'");
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

        push(@cmds,"$bsub -J varscan_validation_tumnor -w \'ended(JOB3) && ended(JOB4)\' \'gmt varscan validation --normal-bam $realigned_normal_bam_file --tumor-bam $realigned_tumor_bam_file --output-indel $output_indel.tumnor --output-snp $output_snp.tumnor --varscan-params \"$varscan_params\"\'");
        push(@cmds,"$bsub -J varscan_validation_relnor -w \'ended(JOB3) && ended(JOB5)\' \'gmt varscan validation --normal-bam $realigned_normal_bam_file --tumor-bam $realigned_relapse_bam_file --output-indel $output_indel.relnor --output-snp $output_snp.relnor --varscan-params \"$default_varscan_params\"\'");
        push(@cmds,"$bsub -J varscan_validation_reltum -w \'ended(JOB4) && ended(JOB5)\' \'gmt varscan validation --normal-bam $realigned_tumor_bam_file --tumor-bam $realigned_relapse_bam_file --output-indel $output_indel.reltum --output-snp $output_snp.reltum --varscan-params \"$default_varscan_params\"\'");

        push(@cmds,"$bsub -N -u $user\@genome.wustl.edu -J varscan_process_validation_tumnor -w \'ended(JOB6)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel.tumnor --validation-snp-file $output_snp.tumnor --variants-file $small_indel_list_nobed --output-file $final_output_file.tumnor\'");
        push(@cmds,"$bsub -N -u $user\@genome.wustl.edu -J varscan_process_validation_relnor -w \'ended(JOB7)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel.relnor --validation-snp-file $output_snp.relnor --variants-file $small_indel_list_nobed --output-file $final_output_file.relnor\'");
        push(@cmds,"$bsub -N -u $user\@genome.wustl.edu -J varscan_process_validation_reltum -w \'ended(JOB8)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel.reltum --validation-snp-file $output_snp.reltum --variants-file $small_indel_list_nobed --output-file $final_output_file.reltum\'");
    }
    else{
#/gscuser/dkoboldt/Software/GATK/GenomeAnalysisTK-1.0.4418/GenomeAnalysisTK.jar /gsc/scripts/pkg/bio/gatk/GenomeAnalysisTK-1.0.5336/GenomeAnalysisTK.jar
        my $bsub_normal_output = "$realigned_bam_file_directory/realignment_normal.out";
        my $bsub_normal_error = "$realigned_bam_file_directory/realignment_normal.err";
        push(@cmds,"$bsub -J $realigned_normal_bam_file -o $bsub_normal_output -e $bsub_normal_error \'java -Xmx16g -Djava.io.tmpdir=/tmp -jar $ENV{GENOME_SW}/gatk/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -et NO_ET -T IndelRealigner -targetIntervals $small_indel_list -o $realigned_normal_bam_file -I $normal_bam -R $reference  --targetIntervalsAreNotSorted\'");

        my $bsub_tumor_output = "$realigned_bam_file_directory/realignment_tumor.out";
        my $bsub_tumor_error = "$realigned_bam_file_directory/realignment_tumor.err";
        push(@cmds,"$bsub -J $realigned_tumor_bam_file -o $bsub_tumor_output -e $bsub_tumor_error \'java -Xmx16g -Djava.io.tmpdir=/tmp -jar $ENV{GENOME_SW}/gatk/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -et NO_ET -T IndelRealigner -targetIntervals $small_indel_list -o $realigned_tumor_bam_file -I $tumor_bam -R $reference --targetIntervalsAreNotSorted\'");

        push(@cmds, "$bsub -J bamindex_normal -w \'ended(JOB0)\' \'samtools index $realigned_normal_bam_file\'");
        push(@cmds, "$bsub -J bamindex_tumor -w \'ended(JOB1)\' \'samtools index $realigned_tumor_bam_file\'");
        push(@cmds, "$bsub -J varscan_validation -w \'ended(JOB2) && ended(JOB3)\' \'gmt varscan validation --normal-bam $realigned_normal_bam_file --tumor-bam $realigned_tumor_bam_file --output-indel $output_indel --output-snp $output_snp --varscan-params \"$varscan_params\"\'");

        push(@cmds,"$bsub -N -u $user\@genome.wustl.edu -J varscan_process_validation -w \'ended(JOB4)\' \'gmt varscan process-validation-indels --validation-indel-file $output_indel --validation-snp-file $output_snp --variants-file $small_indel_list_nobed --output-file $final_output_file\'");
    }

    my @jobids;
    foreach my $cmd (@cmds){        
        #waiting on two jobs
        if($cmd =~ /JOB(\d).+JOB(\d)/){
            my $j1 = $1;
            my $j2 = $2;
            my $jobid1 = $jobids[$j1];
            my $jobid2 = $jobids[$j2];
            $cmd =~ s/JOB$j1/$jobid1/g;
            $cmd =~ s/JOB$j2/$jobid2/g;            

        #waiting on one job
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
    return 1;
}

















