#ReAlign BreakDancer SV supporting reads using novoalign and produce a bam file
package Genome::Model::Tools::DetectVariants2::Filter::NovoRealign;

use strict;
use warnings;
use Genome;
use File::Copy;
use File::Basename;

class Genome::Model::Tools::DetectVariants2::Filter::NovoRealign {
    is  => 'Genome::Model::Tools::DetectVariants2::Filter',
    has_optional => [
        config_file => {
            calculate_from => 'detector_directory',
            calculate => q{ return $detector_directory.'/breakdancer_config';},
            doc  => 'breakdancer config file',
        },
        pass_staging_output => {
            is => 'FilePath',
            calculate_from => '_temp_staging_directory',
            calculate => q{ return $_temp_staging_directory . '/svs.hq'; },
        },
        fail_staging_output => {
            is => 'FilePath',
            calculate_from => '_temp_staging_directory',
            calculate => q{ return $_temp_staging_directory . '/svs.lq'; },
        },
        novoalign_version => {
            type => 'String',
            doc  => 'novoalign version to use in this process',
            default_value =>  '2.05.13',
            valid_values  => [Genome::Model::Tools::Novocraft->available_novocraft_versions],
        },
        novoalign_path => {
            type => 'String',
            doc  => 'novoalign executeable path to use',
            calculate_from => 'novoalign_version',
            calculate => q{ return Genome::Model::Tools::Novocraft->path_for_novocraft_version($novoalign_version); },
        },
        novo2sam_path => {
            type => 'String',
            doc  => 'Path to novosam.pl',
            calculate_from => 'novoalign_version',
            calculate => q{ return Genome::Model::Tools::Novocraft->path_for_novosam_version($novoalign_version); },
        },
        platform => {
            type => 'String',
            doc  => 'Path to novoalign reference sequence index',
            default_value => 'SLX',
        },
        samtools_version => {
            type => 'String',
            doc  => 'samtools version to use in this process',
            default_value =>  Genome::Model::Tools::Sam->default_samtools_version,
            valid_values  => [Genome::Model::Tools::Sam->available_samtools_versions],
        },
        samtools_path => {
            type => 'String',
            calculate_from => 'samtools_version',
            calculate => q{ return Genome::Model::Tools::Sam->path_for_samtools_version($samtools_version); },
            doc => 'path to samtools executable',
        },
        breakdancer_path => {
            type => 'String',
            calculate_from => 'detector_version',
            calculate => q{ return Genome::Model::Tools::Breakdancer->breakdancer_max_command_for_version($detector_version); },
            doc => 'path to breakdancer executable',
        },

    ],
    has_param => [
        lsf_resource => {
            default_value => "-R 'select[mem>10000] rusage[mem=10000]' -M 10000000",
        },
    ],
};

sub _variant_type { 'svs' };

sub _filter_variants {
    my $self     = shift;
    my $cfg_file = $self->config_file;

    #Allow 0 size of output
    if (-z $cfg_file) {
        $self->warning_message('0 size of breakdancer config file. Probably it is for testing of samll bam files');
        my $output_file = $self->pass_staging_output;
        `touch $output_file`;
        return 1;
    }

    my (%mean_insertsize, %std_insertsize, %readlens);

    my $fh = Genome::Sys->open_file_for_reading($cfg_file) or die "unable to open config file: $cfg_file";
    while (my $line = $fh->getline) {
        next unless $line =~ /\S+/;
        chomp $line;
        my ($mean)   = $line =~ /mean\w*\:(\S+)\b/i;
        my ($std)    = $line =~ /std\w*\:(\S+)\b/i;
        my ($lib)    = $line =~ /lib\w*\:(\S+)\b/i;
        my ($rd_len) = $line =~ /readlen\w*\:(\S+)\b/i;

        ($lib) = $line =~ /samp\w*\:(\S+)\b/i unless defined $lib;
        $mean_insertsize{$lib} = int($mean + 0.5);
        $std_insertsize{$lib}  = int($std  + 0.5);
        $readlens{$lib}        = $rd_len;
    }
    $fh->close;

    my %fastqs;
    my $prefix;
    my %conv_libs;

    my $dir = $self->detector_directory;
    opendir (DIR, $dir) or die "Failed to open directory $dir\n";

    for my $fastq (grep{/\.fastq/} readdir(DIR)){
        for my $lib (keys %mean_insertsize) {
            my $conv_lib = $lib;

            if ($lib =~ /[\(\)]/) {  #sometimes the library name is like H_KU-15901-D108132(2)-lib2
                $self->warning_message("$conv_lib contains parentesis");
                $conv_lib =~ s{\(}{\\(}g;
                $conv_lib =~ s{\)}{\\)}g;
            }

            $conv_libs{$lib} = $conv_lib;

            if ($fastq =~/^(\S+)\.${conv_lib}\.\S*([12])\.fastq/) {
                $prefix = $1;
                my $id  = $2;
                if ($lib =~ /[\(\)]/) {
                    $fastq =~ s{\(}{\\(}g;
                    $fastq =~ s{\)}{\\)}g;
                }
                push @{$fastqs{$lib}{$id}}, $dir.'/'.$fastq if defined $id;
                last;
            }
        }
    }
    closedir(DIR);

    unless (scalar keys %mean_insertsize == scalar keys %fastqs) {
        $self->warning_message("library list does not match in mean_insertsize and fastqs");
    }
    #Move breakdancer_config to output_directory so TigraValidation
    #can use it to parse out skip_libraries
    if (-s $cfg_file) {
        copy $cfg_file, $self->_temp_staging_directory;
    }
    else {
        $self->warning_message("Failed to find breakdancer_config from detector_directory: $dir");
    }

    $prefix = $self->_temp_staging_directory . "/$prefix";

    my @bams2remove; 
    my @librmdupbams;
    my @novoaligns;
    my %headerline;

    my $novo_path     = $self->novoalign_path;
    my $novosam_path  = $self->novo2sam_path;
    my $samtools_path = $self->samtools_path;

    my $ref_seq     = $self->reference_sequence_input;
    my $ref_seq_idx = $ref_seq . '.fai';
    unless (-s $ref_seq_idx) {
        $self->error_message("Failed to find ref seq fasta index file: $ref_seq_idx");
        die;
    }

    my $build_id = $self->reference_build_id;
    $self->debug_message("The reference sequence build id is : $build_id");

    my $novo_idx_obj = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_or_create(
        reference_build_id => $build_id,
        aligner_version    => $self->novoalign_version,
        aligner_name       => 'novocraft',
        aligner_params     => '-k 14 -s 3',
    );
    unless (defined $novo_idx_obj) {
        die "Could not retrieve novocraft index for reference build $build_id and aligner version " . $self->novoalign_version;
    }

    my $novo_idx = $novo_idx_obj->full_consensus_path('fa.novocraft');
    unless (-e $novo_idx) {
        die "Found no novocraft index file at $novo_idx for reference build $build_id and aligner version " . $self->novoalign_version;
    }
    $self->debug_message("Found novocraft index file: $novo_idx");

    my $skip_count = 0;
    for my $lib (keys %fastqs) {
        my @read1s = @{$fastqs{$lib}{1}};
        my @read2s = @{$fastqs{$lib}{2}};
        my $line   = sprintf "\@RG\tID:%s\tPU:%s\tLB:%s", $lib, $self->platform, $lib;

        unless(grep { -s $_ } (@read1s, @read2s)) {
            $self->warning_message('No FASTQS with data for library ' . $lib . '. Skipping library.');
            $skip_count++;
            next;
        }

        my $conv_lib = $conv_libs{$lib};
        $headerline{$line} = 1;

        my @bams;
        my $cmd;
        for (my $i=0; $i<=$#read1s; $i++) {
            unless(-s $read1s[$i] or -s $read2s[$i]) {
                $self->warning_message('FASTQs are empty: ' . $read1s[$i] . ',' . $read2s[$i] . '. Skipping this pair.');
                next;
            }

            my $fout_novo = "$prefix.$conv_lib.$i.novo";
            $cmd = $novo_path . ' -d '. $novo_idx . " -f $read1s[$i] $read2s[$i] -i $mean_insertsize{$lib} $std_insertsize{$lib} > $fout_novo";

            $self->_run_cmd($cmd);
            push @novoaligns,$fout_novo;
            
            my $sort_prefix = "$prefix.$conv_lib.$i";
            $cmd = $novosam_path . " -g $conv_lib -f ".$self->platform." -l $conv_lib $fout_novo | ". $samtools_path. " view -b -S - -t ". $ref_seq_idx .' | ' . $samtools_path." sort - $sort_prefix";
            $self->_run_cmd($cmd);
            push @bams, $sort_prefix.'.bam';
            push @bams2remove, $sort_prefix.'.bam';
        }
    
        if ($#bams>0) {
            #TODO using gmt command modules
            $cmd = $samtools_path ." merge $prefix.$conv_lib.bam ". join(' ', @bams);
            $self->_run_cmd($cmd);
            push @bams2remove, "$prefix.$conv_lib.bam";
        }
        else {
            `mv $bams[0] $prefix.$conv_lib.bam`;  #rename behaves strange to interpolate/escape, use mv for now.
        }

        $cmd = $samtools_path." rmdup $prefix.$conv_lib.bam $prefix.$conv_lib.rmdup.bam"; #using $conv_lib here will properly parse ()
        $self->_run_cmd($cmd);
        push @librmdupbams, "$prefix.$conv_lib.rmdup.bam";
    }

    my $merge_bam = "$prefix.novo.rmdup.bam";

    if (@librmdupbams) {
        my $cmd;
        if (@librmdupbams == 1) {
            if ($self->control_aligned_reads_input and $skip_count ne 1) {
                $self->error_message('It is impossible for somatic case to have only 1 per lib rmdup bam');
                die;
            }
            $self->warning_message('There is only 1 per library rmdup bam. Probably for germline purpose');
            `mv $librmdupbams[0] $merge_bam`;  #rg header already made during novo2sam step
        }
        else {
            $self->debug_message('Now merge per library rmdup bam files');
            my $header_file = $prefix . '.header';
            my $header = Genome::Sys->open_file_for_writing($header_file) or die "fail to open $header_file for writing\n";

            for my $line (keys %headerline) {
                $header->print("$line\n");
            }
            $header->close;

            $cmd = $samtools_path . " merge -h $header_file $merge_bam ". join(' ', @librmdupbams);
            $self->_run_cmd($cmd);
            unlink $header_file;
        }
    }
    else {
        if($skip_count eq scalar(keys %fastqs)) {
            my $svs_hq_file = join("/", $self->detector_directory, "svs.hq");
            my $svs_hq_fh = Genome::Sys->open_file_for_reading($svs_hq_file);

            while(<$svs_hq_fh>) {
                next if $_ =~ /^#/; #skip header lines

                #found a variant
                $self->error_message('All libraries were skipped but there appear to be variants in the svs.hq file.');
                $svs_hq_fh->close;
                die;
            }

            #if we get here we didn't find any variants
            $svs_hq_fh->close;
            $self->debug_message('No SVs in svs.hq, skipping run');
            my $output_file = $self->pass_staging_output;
            `touch $output_file`;
            return 1;
        }


        $self->error_message('There is no per library rmdup bam');
        die;
    }

    $self->debug_message('Now prepare breakdancer config file');
    my $novo_cfg    = "$prefix.novo.cfg";
    my $novo_cfg_fh = Genome::Sys->open_file_for_writing($novo_cfg) or die "failed to open $novo_cfg for writing\n";

    for my $lib (keys %fastqs) {
        #Somehow platform:illumina is a required field in config file for breakdancer1.4.2, not for 1.3
        #Hardcoded for now and remove it once it is not needed in new version of bd
        $novo_cfg_fh->printf("platform:illumina\tmap:$merge_bam\tmean:%s\tstd:%s\treadlen:%s\tsample:%s\texe:samtools view\n",$mean_insertsize{$lib},$std_insertsize{$lib},$readlens{$lib},$lib);
    }
    $novo_cfg_fh->close;

    unless (-s $novo_cfg) {
        $self->error_message("novo.cfg file $novo_cfg is not valid");
        die;
    }

    unlink (@bams2remove, @librmdupbams, @novoaligns);
    $self->debug_message('Now run breakdancer');

    my $bd_out_hq_filtered = $self->pass_staging_output;
    my $bd_out_lq_filtered = $self->fail_staging_output;
    my $bd_in_hq           = $self->detector_directory .'/svs.hq';  #DV2::Filter does not have _sv_base_name preset
    my $bd_path            = $self->breakdancer_path;

    my $cmd = $bd_path . ' -t -a '. $novo_cfg .' > '. $bd_out_hq_filtered; #-a to print out copy number per lib so tigra can skip normal lib
    $self->_run_cmd($cmd);

    my $bd_in_hq_fh  = Genome::Sys->open_file_for_reading($bd_in_hq) or die "Failed to open $bd_in_hq for reading\n";
    my $bd_out_hq_fh = Genome::Sys->open_file_for_reading($bd_out_hq_filtered) or die "Failed to open $bd_out_hq_filtered for reading\n";
    my $bd_out_lq_fh = Genome::Sys->open_file_for_writing($bd_out_lq_filtered) or die "Failed to open $bd_out_lq_filtered for writing\n";

    my %filter_match;

    while (my $line = $bd_out_hq_fh->getline) {
        next if $line =~ /^#/;
        my $match = _get_match_key($line);
        $filter_match{$match} = 1;
    }

    while (my $l = $bd_in_hq_fh->getline) {
        next if $l =~ /^#/;
        my $match = _get_match_key($l);
        $bd_out_lq_fh->print($l) unless exists $filter_match{$match};
    }

    $bd_in_hq_fh->close;
    $bd_out_hq_fh->close;
    $bd_out_lq_fh->close;

    return 1;
}


sub _validate_output {
    my $self = shift;

    unless(-d $self->output_directory){
        die $self->error_message("Could not validate the existence of output_directory");
    }
    
    my @files = glob($self->output_directory."/svs.hq");
    unless (@files) {
        die $self->error_message("Failed to get svs.hq");
    }
    return 1;
}


sub _get_match_key {
    my $line = shift;
    my @columns = split /\s+/, $line;
    #compare chr1 pos1 chr2 pos2 sv_type 5 columns
    my $match = join '-', $columns[0], $columns[1], $columns[3], $columns[4], $columns[6];
    return $match;
}


sub _run_cmd {
    my ($self, $cmd) = @_;
    
    unless (Genome::Sys->shellcmd(cmd => $cmd)) {
        $self->error_message("Failed to run $cmd");
        die $self->error_message;
    }
    return 1;
}

sub _create_bed_file {
    return 1;
}


1;
