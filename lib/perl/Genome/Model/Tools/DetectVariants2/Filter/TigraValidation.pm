#extract reads from a set of bams in breakdancer predicted regions
package Genome::Model::Tools::DetectVariants2::Filter::TigraValidation;

use strict;
use warnings;
use Genome;
use Bio::Seq;
use Bio::SeqIO;
use File::Temp;
use File::Basename;
use Carp 'confess';

#my %opts = (l=>500,p=>1000,s=>0,q=>0,n=>0,m=>10,x=>3,P=>-10,G=>-10,S=>0.02,A=>500,Q=>0);
#getopts('l:d:c:p:r:s:Q:q:t:n:a:b:f:km:MRzv:hD:x:i:P:G:I:A:S:L:',\%opts);
#die("
#Usage:   AssemblyValidation.pl <SV file, default in BreakDancer format> <bam files ... >
#Options:
#         -d DIR     Directory where all intermediate files are saved
#         -I DIR     Read intermediate files from DIR instead of creating them
#         -f FILE    Save Breakpoint sequences in file
#         -r FILE    Save relevant cross_match alignment results to file
#         -z         Customized SV format, interpret column from header that must start with # and contain chr1, start, chr2, end, type, size
#         -v FILE    Save unconfirmed predictions in FILE
#         -h         Make homo variants het by adding same amount of randomly selected wiletype reads
#         -l INT     Flanking size [$opts{l}]
#         -A INT     Esimated maximal insert size [$opts{A}]
#         -q INT     Only assemble reads with mapping quality > [$opts{q}]
#         -m INT     Minimal size (bp) for an assembled SV to be called confirmed [$opts{m}]
#         -i INT     invalidate indels are -i bp bigger or smaller than the predicted size, usually 1 std insert size
#         -a INT     Get reads with start position bp into the left breakpoint, default [50,100,150]
#         -b INT     Get reads with start position bp into the right breakpoint, default [50,100,150]
#         -S FLOAT   Maximally allowed polymorphism rate in the flanking region [$opts{S}];
#         -P INT     Substitution penalty in cross_match alignment [$opts{P}]
#         -G INT     Gap Initialization penalty in cross_match alignment [$opts{G}]
#         -M         Prefix reference name with \'chr\' when fetching reads using samtools view
#         -R         Assemble Mouse calls, NCBI reference build 37
#
#Filtering:
#         -p INT     Ignore cases that have average read depth greater than [$opts{p}]
#         -c STRING  Specify a single chromosome
#         -s INT     Minimal size of the region to analysis [$opts{s}]
#         -Q INT     minimal BreakDancer score required for analysis [$opts{Q}]
#         -t STRING  type of SV
#         -n INT     minimal number of supporting reads [$opts{n}]
#         -L STRING  Ingore calls supported by libraries that contains (comma separated) STRING
#         -k         Attach assembly results as additional columns in the input file
#         -D DIR     A directory that contains a set of supplementary reads (when they are missing from the main data stream)
#         -x INT     Duplicate supplementary reads [$opts{x}] times
#\n") unless ($#ARGV>=1);

my @FULL_CHR_LIST = (1..22, 'X', 'Y', 'MT');

class Genome::Model::Tools::DetectVariants2::Filter::TigraValidation {
    is  => 'Genome::Model::Tools::DetectVariants2::Filter',
    has_optional => [
        pass_output => {
            is => 'FilePath',
            calculate_from => '_temp_staging_directory',
            calculate => q{ return $_temp_staging_directory . '/svs.hq'; },
        },
        fail_output => {
            is => 'FilePath',
            calculate_from => '_temp_staging_directory',
            calculate => q{ return $_temp_staging_directory . '/svs.lq'; },
        },
        sv_output_name => { 
            is => 'Text',
            default_value => 'svs.out',
        },
        tigra_version => {
            type => 'String',
            doc  => 'tigra_sv version to use in this process',
            default_value => '0.1', 
            valid_values  => [Genome::Model::Tools::TigraSv->available_tigrasv_versions],
        },
        tigra_path => {
            is  => 'FilePath',
            doc => 'tigra_sv executable path to use',
            calculate_from => 'tigra_version',
            calculate      => q{return Genome::Model::Tools::TigraSv->path_for_tigrasv_version($tigra_version);},
        },
        sv_merge_path => {
            is => 'FilePath',
            default => '/gsc/scripts/opt/genome-stable/lib/perl/Genome/Model/Tools/Sv/MergeAssembledCallsets.pl',
        },
        sv_annot_path => {
            is => 'FilePath',
            default => '/gsc/scripts/opt/genome-stable/lib/perl/Genome/Model/Tools/Sv/BreakAnnot.pl',
        },
        # TODO Either point to a specific version of phrap or (even better) use the crossmatch tool
        crossmatch_path => {
            is => 'FilePath',
            default => Genome::Sys->sw_path('phrap','1.080721','cross_match'),
        },
        workflow_log_directory => {
            calculate_from => 'output_directory',
            calculate => q{ return $output_directory . '/tigrasv_by_chromosome_log'; },
        },
        breakpoint_seq_file => {
            calculate_from => '_temp_staging_directory',
            calculate => q{ return $_temp_staging_directory . '/breakpoint_seq'; },
        },
        specify_chr => {
            is => 'String',
            is_input => 1,
            valid_values => [@FULL_CHR_LIST, 'all'],
            default => 'all',
        },
        cm_aln_file => {
            calculate_from => '_temp_staging_directory',
            calculate => q{ return $_temp_staging_directory . '/cm_aln.out'; },
            doc  => 'Save relevant cross_match alignment results to file',
        },
        # tigra options, can be overridden via params
        maximum_node => {
            is => 'Number',
            default => 300,
        },
        dump_reads => {
            is => 'Boolean',
            default => 0,
        },
        write_local_ref => {
            is => 'Boolean',
            default => 1,
        },
        custom_sv_format => {
            type => 'Boolean',
            doc  => 'Customized SV format, interpret column from header that must start with # and contain chr1, start, chr2, end, type, size',
            default => 0,
        },
        asm_high_coverage => {
            type => 'Boolean',
            doc  => 'Assemble high coverage data (such as capture)',
            default => 0,
        },
        pad_local_ref => {
            type => 'Integer',
            doc  => 'Pad local reference by additional bp on both ends',
            default_value => 200,
        },
        mismatch_limit => {
            type => 'Integer',
            doc  => 'Number of mismatches required to be tagged as poorly mapped',
            default_value => 5,
        },
        flank_size => {
            type => 'Integer',
            doc  => 'Flanking size',
            default_value => 500,
        },
        est_max_ins_size => {
            type => 'Integer',
            doc  => 'Esimated maximal insert size',
            default_value => 500,
        },
        map_qual_to_asm => {
            type => 'Integer',
            doc  => 'Only assemble reads with mapping quality >',
        },
        min_size_of_confirm_asm_sv => {
            type => 'Integer',
            doc  => 'Minimal size (bp) for an assembled SV to be called confirmed',
            default_value => 10, #original value is 3, changed based on Xian's script
        },
        invalid_indel_range => {
            type => 'Integer',
            doc  => 'invalidate indels are -i bp bigger or smaller than the predicted size, usually 1 std insert size',
        },
        avg_read_depth_limit => {
            type => 'Integer',
            doc  => 'Ignore cases that have average read depth greater than',
            default_value => 1000,
        },
        min_breakdancer_score => {
            type => 'Number',
            doc  => 'minimal BreakDancer score required for analysis',
            default_value => 40, #original value is 0, changed based on Xian's script
        },
        skip_libraries => {
            type => 'String',
            doc  => 'Ingore calls supported by libraries that contains (comma separated)',
            is_input => 1,  #for now maybe the skip_libs stuff should be moved into params 
        },
        skip_call => {
            type => 'Integer',
            doc => 'Ignore tigra_sv for the calls before this number',
        },
        # Crossmatch parameters
        cm_bandwidth => {
            is => 'Number',
            default => 20,
        },
        cm_minimum_match => {
            is => 'Number',
            default => 20,
        },
        cm_minimum_score => {
            is => 'Number',
            default => 25,
        },
        cm_penalty => {
            is => 'Number',
            default => -10,
        },
        cm_gap_init => {
            is => 'Number',
            default => -10,
        },
        cm_gap_extension => {
            is => 'Number',
            default => -1,
        },
        cm_discrepency_lists => {
            is => 'Boolean',
            default => 1,
        },
        cm_tags => {
            is => 'Boolean',
            default => 1,
        },
        cm_max_polymorphism_rate => {
            type => 'Number',
            doc  => 'Maximally allowed polymorphism rate in the flanking region',
            default_value => 0.02,
        },
        _maxSV => {
            type => 'HASH',
        },
        _N50size => {
            type => 'Number',
        },
        _WeightAvgSize => {
            type => 'Number',
        },
        _tigra_data_dir => {
            type => 'String',
        },
        _run_by_workflow => {
            is => 'Boolean',
            is_input => 1,
            default => 0,
        },
        _base_output_directory => {
            is => 'Text',
            is_optional => 1,
            doc => 'Store the base output directory when using per-chromosome output dirs',
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => "-R 'select[mem>8000] rusage[mem=8000]' -M 8000000", 
        },
    ],
};

sub _variant_type { 'svs' };

my %TIGRA_PARAMS_LIST = (
    l => 'flank_size',
    c => 'specify_chr',
    R => 'reference_sequence_input',
    q => 'map_qual_to_asm',
    N => 'mismatch_limit',
    p => 'avg_read_depth_limit',
    #I => 'output_directory',   we don't want all intermediate files dumpping to output dir
    w => 'pad_local_ref',
    Q => 'min_breakdancer_score',
    M => 'min_size_of_confirm_asm_sv',
    h => 'maximum_node',
    A => 'est_max_ins_size',
    L => 'skip_libraries',
);
my %TIGRA_FLAGS_LIST = (
    r => 'write_local_ref',
    d => 'dump_reads',
    b => 'custom_sv_format',
);
my %TIGRA_PROCESSING_PARAMS = (
    i => 'invalid_indel_range',
);
my %TIGRA_PROCESSING_FLAGS = (
    z => 'skip_call',
    a => 'asm_high_coverage',
);

my %CROSSMATCH_PARAMS_LIST = (
    bandwidth => 'cm_bandwidth',
    minmatch  => 'cm_minimum_match',
    minscore  => 'cm_minimum_score',
    penalty   => 'cm_penalty',
    gap_init  => 'cm_gap_init',
    gap_ext   => 'cm_gap_extension',
);
my %CROSSMATCH_FLAGS_LIST = (
    discrep_lists => 'cm_discrepency_lists',
    tags => 'cm_tags',
);

sub _create_temp_directories {
    my $self = shift;
    local %ENV = %ENV;
    $ENV{TMPDIR} = $self->output_directory;
    return $self->SUPER::_create_temp_directories(@_);
}

sub _promote_staged_data {
    my $self = shift;
    my $output_dir = $self->SUPER::_promote_staged_data(@_);

    #since the temp dir is under the output dir, we need to remove it
    #before the final reallocation or everything is double-counted
    my $staging_dir = $self->_temp_staging_directory;
    Genome::Sys->remove_directory_tree($staging_dir);

    return $output_dir;
}

sub _resolve_output_directory {
    my $self = shift;

    if ($self->_base_output_directory and $self->_base_output_directory ne $self->output_directory) {
        return 1;
    }

    $self->_base_output_directory($self->output_directory);

    if ($self->_run_by_workflow) {
        my $output_dir = $self->output_directory;
        unless (-d $output_dir) {
            #This should only happen if a single chromosome was executed directly
            Genome::Sys->create_directory($output_dir);
        }
        #Put per-chromosome outputs in subdirectories to avoid collisions in SoftwareResults
        $self->output_directory($output_dir . '/' . $self->specify_chr);
    }

    return 1;
}


sub _filter_variants {
    my $self = shift;
    my $variant_file = $self->_breakdancer_input;

    #Allow 0 size of output
    if (-z $variant_file) {
        $self->warning_message('0 size of breakdancer input : '.$variant_file.'. Probably it is for testing of small bams');
        my $pass_out = $self->pass_output;
        `touch $pass_out`;
        my @output_files = map{$self->_temp_staging_directory .'/'.$self->_variant_type.'.merge.'.$_}qw(file out fasta);
        `touch @output_files`;
        return 1;
    }

    if ($self->specify_chr eq 'all') {
        $self->status_message("Splitting breakdancer input file by chromosome");

        # Split up breakdancer file by chromosome so tigra can be run in parallel
        my $split_obj = $self->_get_split_object; 
        
        my $rv = $split_obj->execute;
        Carp::confess 'Could not execute breakdancer split file command!' unless defined $rv and $rv == 1;

        my @use_chr_list = $self->_use_chr_list;
        unless(@use_chr_list) {
            #squaredancer includes a header even when no results, so the -z above doesn't catch this case
            $self->warning_message('0 size of breakdancer input (excluding header): '.$variant_file.'.');
            my $pass_out = $self->pass_output;
            `touch $pass_out`;
            my @output_files = map{$self->_temp_staging_directory .'/'.$self->_variant_type.'.merge.'.$_}qw(file out fasta);
            `touch @output_files`;
            return 1;
        }
        
        
        my $skip_libs    = $self->skip_libraries || $self->_get_skip_libs;

        $self->status_message("Creating workflow to parallelize by chromosome");

        # Create and execute workflow
        require Workflow::Simple;
        my $op = Workflow::Operation->create(
            name => 'Tigra by chromosome',
            operation_type => Workflow::OperationType::Command->get(ref($self)),
        );
        $op->parallel_by('specify_chr');

        if(Workflow::Model->parent_workflow_log_dir) {
            $op->log_dir(Workflow::Model->parent_workflow_log_dir);
        } elsif ($self->workflow_log_directory) {
            unless (-d $self->workflow_log_directory) {
                unless (Genome::Sys->create_directory($self->workflow_log_directory)) {
                    $self->error_message('Failed to create workflow_log_directory: '.$self->workflow_log_directory);
                    die;
                }
            }
            $op->log_dir($self->workflow_log_directory);
        }

        $self->status_message("Running workflow");

        my %options = (
            previous_result_id => $self->previous_result_id,
            output_directory   => $self->_temp_staging_directory,
            specify_chr        => \@use_chr_list,
            _run_by_workflow   => 1,
        );

        $options{skip_libraries} = $skip_libs if $skip_libs;

        my $output = Workflow::Simple::run_workflow_lsf($op, %options);

        unless (defined $output) {
            my @error;
            for (@Workflow::Simple::ERROR) {
                push @error, $_->error;
            }
            $self->error_message(join("\n", @error));
            die $self->error_message;
        }

        # Now merge together all the pass/fail files produced for each chromosome
        $self->status_message("Merging output files together");

        # use outputs in software-result output_dir to generate final outputs
        my $sr_dirs = $self->_get_sr_dirs(@use_chr_list);
        my @sr_dirs;
        map{push @sr_dirs, $sr_dirs->{$_}}@use_chr_list; # make the same order as before
                
        for my $file ($self->pass_output, $self->fail_output) {
            my $merge_obj = Genome::Model::Tools::Breakdancer::MergeFiles->create(
                input_files => join(',', map { $_ . '/' . basename($file) } @sr_dirs),
                output_file => $file,
            );
            my $merge_rv = $merge_obj->execute;
            Carp::confess 'Could not execute breakdancer merge command!' unless defined $merge_rv and $merge_rv == 1;
        }

        $self->status_message("Running MergeCallSet");

        my ($merge_index, $merge_file, $merge_annot, $merge_out, $merge_fa) = 
            map{$self->_temp_staging_directory .'/'.$self->_variant_type.'.merge.'.$_}qw(index file file.annot out fasta);

        my $idx_fh = IO::File->new(">$merge_index") or die "Failed to open $merge_index for writing\n";
        my %valid_bams = $self->_check_bam; #normal bam, tumor bam

        for my $chr (@use_chr_list) {
            my $chr_dir = $sr_dirs->{$chr};
            for my $type (sort keys %valid_bams) {
                my $name   = $type.$chr;
                my $sv_fa  = $chr_dir .'/'. basename($self->breakpoint_seq_file) .".$type.fa";
                my $sv_out = $chr_dir .'/'. $self->sv_output_name .'.'. $type;

                if (-e $sv_fa and -e $sv_out) {
                    $idx_fh->print($name.' '.$sv_out.' '.$sv_fa."\n");
                }
                else {
                    $self->warning_message("$name: $sv_fa and $sv_out are not both valid");
                }
            }
        }
        $idx_fh->close;

        my $merge_cmd = $self->sv_merge_path . ' -c -f '.$merge_fa.' -d 200 -h '.$merge_index.' 1> '.$merge_file.' 2> '.$merge_out;
        $rv = system $merge_cmd;
        unless ($rv == 0) {
            $self->warning_message("Sv merge command: $merge_cmd probably did not finish ok");
            return 1;
        }

        $self->status_message('Running SvAnnot');
        my $ref_build_id = $self->_get_ref_build_id;

        if ($ref_build_id) {
            my %annot_params = (
                sv_file     => $merge_file,
                output_file => $merge_annot,
                sv_format   => 'merged',
                annot_build => $ref_build_id,
            );
            $annot_params{repeat_mask} = 1 if $ref_build_id eq '36';

            my $annot = Genome::Model::Tools::Sv::SvAnnot->create(%annot_params);
            my $rv = $annot->execute;    
            $self->warning_message("SvAnnot probably did not finished ok") unless $rv == 1;     
        }
        else {
            $self->warning_message('No ref_build_id available. Skip SvAnnot');
            `touch $merge_annot`;
        }
        return 1;
    }

    # If running as part of a workflow, need to update certain properties to contain the chromosome name
    if ($self->_run_by_workflow) {
        $variant_file = dirname($self->output_directory) . '/svs.hq.tigra.' . $self->specify_chr; #since now svs.hq.tigra.chr is on upper level dir
        my $out_dir = $self->output_directory;
        unless (-d $out_dir) {
            unless (Genome::Sys->create_directory($out_dir)) {
                $self->error_message("Failed to create directroy: $out_dir under _run_by_workflow block");
                die;
            }
        }
    }

    die "No valid variant file existing\n" unless $variant_file and -e $variant_file;
    $self->status_message("Parsing params");
    $self->parse_params; # goes through the param string, sets object properties

    my %bam_files = $self->_check_bam;
    my @tmp_tigra_files = ();

    #MG need tigra_sv run on each single bam separately not all bams
    #for my $i (0..$#bam_files) {
    for my $type (keys %bam_files) {
        my $bam_file = $bam_files{$type};

        my $tmp_tigra_dir = File::Temp::tempdir('tigra_sv_out_'.$self->specify_chr.'_'.$type.'_XXXXXX', DIR => '/tmp', CLEANUP => 1);
        $self->_tigra_data_dir($tmp_tigra_dir); 

        # Construct tigra command and execute
        $self->status_message("Making tigra command and executing");

        my $sv_output     = $self->_temp_staging_directory . '/' .$self->sv_output_name .'.'. $type;
        my $tigra_options = $self->_get_tigra_options;
        $tigra_options = '-I ' . $tmp_tigra_dir . ' ' . $tigra_options; #hack for now, need separate tigra dump dirs

        my $tigra_sv_cmd = join(' ', $self->tigra_path, $tigra_options, $variant_file, $bam_file, '>', $sv_output);

        my $rv = Genome::Sys->shellcmd(
            cmd         => $tigra_sv_cmd,
            input_files => [$variant_file],
        );
        unless ($rv) {
            $self->error_message("Running tigra_sv failed.\nCommand: $tigra_sv_cmd");
            die;
        }
        $self->status_message("Done running tigra, now parsing");

        my $bp_file = $self->breakpoint_seq_file . '.' . $type .'.fa';
        my $bp_io;
        if ($bp_file) {
            if (-s $bp_file) {
                $self->warning_message('breakpoint seq file: '.$bp_file. ' existing, Now remove it');
                unlink $bp_file;
            }
            $bp_io = Bio::SeqIO->new(
                -file   => ">>$bp_file",
                -format => 'Fasta',
            );
        }

        my $cm_aln_file = $self->cm_aln_file . '.' . $type;
        my $cm_aln_fh   = IO::File->new(">$cm_aln_file") or die "Failed to open $cm_aln_file for writing";

        my @tigra_sv_fas = glob($tmp_tigra_dir . "/*.fa.contigs.fa"); #get tigra homo ctg list
        @tigra_sv_fas = sort{(basename ($a)=~/^\S+?\.(\d+)\./)[0]<=> (basename ($b)=~/^\S+?\.(\d+)\./)[0]}@tigra_sv_fas; #sort ctg file by chr pos
    
        # open output file for the following resume
        my $out_fh = IO::File->new(">>". $sv_output) or die "Failed to open $sv_output for writing";
    
        for my $tigra_sv_fa (@tigra_sv_fas) {
            my ($tigra_sv_name) = basename $tigra_sv_fa =~ /^(\S+)\.fa\.contigs\.fa/;
            my ($chr1,$start,$chr2,$end,$type,$size,$ori,undef) = split /\./, $tigra_sv_name; # you get the $size from $prefix        
            next if($chr2 ne $self->specify_chr && defined $self->specify_chr); 
            my $prefix = join('.',$chr1,$start,$chr2,$end,$type,$size,$ori);

            $self->_N50size(_ComputeTigraN50($tigra_sv_fa));
            $self->_WeightAvgSize(_ComputeTigraWeightedAvgSize($tigra_sv_fa));
        
            #test homo, het contigs
            for my $ctg_type ('homo', 'het') {
                $self->_cross_match_validation($ctg_type, $tigra_sv_name);
            }
        
            my $maxSV = $self->_maxSV;

            # The if statement from hell
            if (defined $maxSV && ($type eq 'CTX' && $maxSV->{type} eq $type ||
	            $type eq 'INV' && $maxSV->{type} eq $type ||
		        (($type eq $maxSV->{type} && $type eq 'DEL') ||
		        ($type eq 'ITX' && ($maxSV->{type} eq 'ITX' || $maxSV->{type} eq 'INS')) ||
		        ($type eq 'INS' && ($maxSV->{type} eq 'ITX' || $maxSV->{type} eq 'INS'))) &&
                $size >= $self->min_size_of_confirm_asm_sv && (!defined $self->invalid_indel_range || abs($maxSV->{size}-$size)<=$self->invalid_indel_range))) {

                my $scarstr = $maxSV->{scarsize}>0 ? substr($maxSV->{contig},$maxSV->{bkstart}-1,$maxSV->{bkend}-$maxSV->{bkstart}+1) : '-';

                $out_fh->printf("%s\t%d(%d)\t%s\t%d(%d)\t%s\t%d(%d)\t%s(%s)\t%s\t%d\t%d\t%d\%\t%d\t%d\t%d\t%d\t%d\t%s\t%s\ta%d.b%d\t%s\t%s\t%s\n",$maxSV->{chr1},$maxSV->{start1},$start,$maxSV->{chr2},$maxSV->{start2},$end,$maxSV->{ori},$maxSV->{size},$size,$maxSV->{type},$type,$maxSV->{het},$maxSV->{weightedsize},$maxSV->{read_len},$maxSV->{fraction_aligned}*100,$maxSV->{n_seg},$maxSV->{n_sub},$maxSV->{n_indel},$maxSV->{nbp_indel},$maxSV->{microhomology},$scarstr,$prefix,50,100, 'NA', 'NA', 'NA');
        
                if ($bp_io) {  #save breakpoint sequence
                    my $coord = join(".",$maxSV->{chr1},$maxSV->{start1},$maxSV->{chr2},$maxSV->{start2},$maxSV->{type},$maxSV->{size},$maxSV->{ori});
                    my $contigsize = $maxSV->{contiglens};
                    my $seqobj = Bio::Seq->new( 
                        -display_id => "ID:$prefix,Var:$coord,Ins:$maxSV->{bkstart}\-$maxSV->{bkend},Length:$contigsize,KmerCoverage:$maxSV->{contigcovs},Strand:$maxSV->{strand},Assembly_Score:$maxSV->{weightedsize},PercNonRefKmerUtil:$maxSV->{kmerutil},Ref_start:$maxSV->{refpos1},Ref_end:$maxSV->{refpos2},Contig_start:$maxSV->{rpos1},Contig_end:$maxSV->{rpos2},TIGRA",
                        -seq => $maxSV->{contig}, 
                    );
                    $bp_io->write_seq($seqobj);
                }

                if ($cm_aln_fh) {
                    $cm_aln_fh->printf("%s\t%d(%d)\t%s\t%d(%d)\t%s\t%d(%d)\t%s(%s)\t%s\t%d\t%d\t%d\%\t%d\t%d\t%d\t%d\t%d\t%s\t%s\ta%d.b%d\n",$maxSV->{chr1},$maxSV->{start1},$start,$maxSV->{chr2},$maxSV->{start2},$end,$maxSV->{ori},$maxSV->{size},$size,$maxSV->{type},$type,$maxSV->{het},$maxSV->{weightedsize},$maxSV->{read_len},$maxSV->{fraction_aligned}*100,$maxSV->{n_seg},$maxSV->{n_sub},$maxSV->{n_indel},$maxSV->{nbp_indel},$maxSV->{microhomology},$scarstr,$prefix,'50','100');
                    for my $aln (split /\,/, $maxSV->{alnstrs}) {
	                    $cm_aln_fh->printf("%s\n", join("\t", split /\|/, $aln));
                    }
                    $cm_aln_fh->print("\n");
                }
            }
            $self->_maxSV(undef)         if $self->_maxSV; #reset for each SV
            $self->_N50size(undef)       if $self->_N50size;
            $self->_WeightAvgSize(undef) if $self->_WeightAvgSize;
        }
        $out_fh->close;
        $cm_aln_fh->close if $cm_aln_fh;
        push @tmp_tigra_files, $sv_output;
    }
    $self->status_message("Done validating, now merging multiple tigra outputs into one");

    my $sv_out_file = $self->_temp_staging_directory . '/' . $self->sv_output_name;
    my $input_files = join ',', @tmp_tigra_files;

    my $tigra_merge = Genome::Model::Tools::Breakdancer::MergeFiles->create(
        input_files => $input_files,
        output_file => $sv_out_file,
    );
    my $merge_rv = $tigra_merge->execute;
    confess 'Could not merge multiple tigra validation lists' unless defined $merge_rv and $merge_rv == 1;

    $self->status_message("Now mapping tigra output to breakdancer output");

    my $tigra_adaptor_obj = Genome::Model::Tools::Breakdancer::TigraToBreakdancer->create(
        original_breakdancer_file => $variant_file,
        tigra_output_file         => $sv_out_file,
        pass_filter_file          => $self->pass_output,
        fail_filter_file          => $self->fail_output,
    );
    my $adaptor_rv = $tigra_adaptor_obj->execute;
    confess 'Could not produce filtered breakdancer files from tigra output!' unless defined $adaptor_rv and $adaptor_rv == 1;
    
    $self->status_message('TigraValidation finished ok.');
    return 1;
}


#overwrite this base method because there is no svs.bed out for now
sub _validate_output {
    my $self = shift;
    #my $name = $self->sv_output_name;
    my $name = "sv*out";

    unless(-d $self->output_directory){
        die $self->error_message("Could not validate the existence of output_directory");
    }
    my @files = glob($self->output_directory."/*$name*");
    unless (@files) {
        die $self->error_message("Failed to get $name from ". $self->output_directory);
    }
    return 1;
}

#Need figure out the source of breakdancer input, directly from
#breakdancer run or from novorealign filter
sub _breakdancer_input {
    my $self     = shift;
    my $bd_input = $self->input_directory . '/svs.hq';

    unless (-e $bd_input) {
        $self->error_message('Failed to find breakdancer input file from input directory: '. $self->input_directory);
        die;
    }
    
    $self->status_message("Find breakdancer input: $bd_input");
    return $bd_input;
}


sub _check_bam {
    my $self = shift;
    my %valid_bams;

    my ($tumor_bam, $normal_bam) = ($self->aligned_reads_input, $self->control_aligned_reads_input);

    if ($tumor_bam) {
        if ($self->_validate_bam('tumor', $tumor_bam)) {
            $valid_bams{tumor} = $tumor_bam;
        }
    }

    if ($normal_bam) {
        if ($self->_validate_bam('normal', $normal_bam)) {
            $valid_bams{normal} = $normal_bam;
        }
    }

    unless (%valid_bams and $valid_bams{tumor}) {
        die $self->error_message('There is no valid bam files for tigra validation');
    }

    return %valid_bams;
}


sub _validate_bam {
    my ($self, $type, $bam) = @_;

    my $bai = $bam . '.bai';
    my $noexist = '%s %s file: %s does not exist';
    my  $nosize = '%s %s file: %s has no size';
    unless (-e $bam) {
        die $self->error_message(sprintf($noexist, $type, 'bam', $bam));
    }
    unless (-s $bam) {
        die $self->error_message(sprintf( $nosize, $type, 'bam', $bam));
    }
    unless (-e $bai) {
        die $self->error_message(sprintf($noexist, $type, 'bam index', $bai));
    }
    unless (-s $bai) {
        die $self->error_message(sprintf( $nosize, $type, 'bam index', $bai));
    }
    return 1;
}


sub parse_params {
    my $self = shift;
    my $param_string = $self->params;
    return 1 unless defined $param_string;
    for my $param (sort keys %TIGRA_PARAMS_LIST) {
        if ($param_string =~ /\-$param(=|\s)(\S+)/) { #must be key/value pair like "-a=12 -t=1" or "-a 12 -t 1" not "-a 12 -t -d"
            my $value = $2;
            my $attribute = $TIGRA_PARAMS_LIST{$param};
            $self->$attribute($value);
        }
    }
    return 1;
}

#for inter-chr breakdancer chr2 column does not have 1 (chr).
sub _use_chr_list {
    my $self = shift;
    my @chr_list = ();

    for my $chr ($self->_full_chromosome_list) {
        if (-s $self->_temp_staging_directory .'/svs.hq.tigra.'.$chr) {  #what if it is not 0 size but only contains breakdancer header ?
            push @chr_list, $chr;
        }
        else {
            $self->warning_message("No valid svs.hq.tigra.$chr. Skip $chr for tigra validation");
        }
    }
    return @chr_list;
}

sub _full_chromosome_list {
     my ($self) = @_;
     return @FULL_CHR_LIST;
}

sub _get_split_object {
    my ($self) = @_;
    return Genome::Model::Tools::Breakdancer::SplitFiles->create(
            input_file           => $self->_breakdancer_input,
            output_directory     => $self->_temp_staging_directory,
            output_file_template => 'svs.hq.tigra.CHR',
        );
}


sub _get_sr_dirs {
    my ($self, @use_chr_list) = @_;
    my %sr_dirs;

    for my $chr_name (@use_chr_list) {
        my $sr_params = $self->params_for_filter_result;
        $sr_params->{chromosome_list} = $chr_name;
        delete $sr_params->{filter_version};

        my $sr = Genome::Model::Tools::DetectVariants2::Result::Filter->get_with_lock(%$sr_params);
        unless ($sr) {
            $self->error_message('Failed to find software result for chromosome '.$chr_name);
            die;
        }
        my $sr_dir = $sr->output_dir;
        unless ($sr_dir and -d $sr_dir) {
            $self->error_message("Software result output_dir for chromosome $chr_name does not exist");
            die;
        }
        $self->status_message("Chromosome $chr_name gets software_result output_dir: $sr_dir");

        #Just check whether software_result output dir is linked to
        #local or not, shortcut or exceute
        my $link_dir = $self->_temp_staging_directory . '/' . $chr_name;
        my $real_dir = readlink($link_dir);
        unless (-d $link_dir) {
            $self->error_message("SoftwareResult output dir for chromosome $chr_name is not linked to local temp_staging_dir as $link_dir");
            die;
        }
        unless ($real_dir eq $sr_dir) {
            $self->error_message("Target link dir for Chr $chr_name : $real_dir is not software_result output dir $sr_dir");
            die;
        }
        $self->status_message("Software_result output dir is correctly linked for Chr $chr_name");
        $sr_dirs{$chr_name} = $sr_dir;
    }
    return \%sr_dirs;
}


sub _get_skip_libs {
    my $self   = shift;
    my $bd_cfg = $self->input_directory . '/breakdancer_config';

    unless (-s $bd_cfg) {
        if ($self->detector_directory) {
            $bd_cfg = $self->detector_directory . '/breakdancer_config';
        }
        unless (-s $bd_cfg) {
            $self->warning_message("Failed to find valid breakdancer_config. Not using skip_libraries for TigraValidation");
            return;
        }
    }
    $self->status_message("Find breakdancer config: $bd_cfg to get skip libraries");

    my $normal_bam = $self->control_aligned_reads_input;
    
    if ($normal_bam) {
        my %libs = ();
        my $fh = Genome::Sys->open_file_for_reading($bd_cfg) or die "Failed to open $bd_cfg\n";
        while (my $line = $fh->getline) {
            if ($line =~ /$normal_bam/) {
                my @columns = split /\s+/, $line;
                my ($lib) = $columns[4] =~ /lib\:(\S+)/;
                if ($lib =~ /[\(\)]/) {  #sometimes the library name is like H_KU-15901-D108132(2)-lib2
                    $self->warning_message("$lib contains parentesis");
                    $lib =~ s{\(}{\\(}g;
                    $lib =~ s{\)}{\\)}g;
                }
                $libs{$lib} = 1;
            }
        }
        $fh->close;

        my $libs = join ",", keys %libs;
        $self->status_message("Skip libraries : $libs for tigra validation");
        return $libs;
    }
    else {
        $self->warning_message('No normal bam given. Run tigra-sv without skip_libraries');
        return;
    }
}


sub _get_tigra_options {
    my $self = shift;
    my $tigra_opts;

    for my $option (sort keys %TIGRA_PARAMS_LIST) {
        my $option_name = $TIGRA_PARAMS_LIST{$option};
        if ($option_name eq 'avg_read_depth_limit' or $option_name eq 'maximum_node') {
            next unless $self->asm_high_coverage;
        }
        if ($option_name eq 'specify_chr') {
            next if $self->$option_name eq 'all';
        }
        $tigra_opts .= '-' . $option . ' ' . $self->$option_name . ' ' if defined $self->$option_name;
    }
    for my $option (sort keys %TIGRA_FLAGS_LIST) {
        my $option_name = $TIGRA_FLAGS_LIST{$option};
        if ($option_name eq 'custom_sv_format') {
            $tigra_opts .= '-' . $option . ' ' unless $self->$option_name;
        }
        else {
            $tigra_opts .= '-' . $option . ' ' if $self->$option_name;
        }
    }
    return $tigra_opts;
}


sub _get_ref_build_id {
    my $self   = shift;
    my $ref_id = $self->reference_build_id;
    return unless $ref_id;

    my %refs = (
        106942997 => 37,  #GRCh37-lite-build37
        102671028 => 37,  #g1k-human-build37
        101947881 => 36,  #NCBI-human-build36
        109104543 => 36,  #fdu_human36_chr16_17_for_novo_test-build for Jenkins apipe test build
        107494762 => 'mouse_37',  #UCSC-mouse build37, mouse ref seq used in production 
    );
    
    if ($refs{$ref_id}) {
        return $refs{$ref_id};
    }
    else {
        return;
    }
}


sub crossmatch_options {
    my $self = shift;
    my $cm_opts;
    for my $option (sort keys %CROSSMATCH_PARAMS_LIST) {
        my $option_name = $CROSSMATCH_PARAMS_LIST{$option};
        $cm_opts .= '-' . $option . ' ' . $self->$option_name . ' ' if defined $self->$option_name;
    }
    for my $option (sort keys %CROSSMATCH_FLAGS_LIST) {
        my $option_name = $CROSSMATCH_FLAGS_LIST{$option};
        $cm_opts .= '-' . $option . ' ' if $self->$option_name;
    }
    return $cm_opts;
}

sub _cross_match_validation {
    my ($self, $ctg_type, $tigra_sv_name) = @_;
    my ($chr1,$start,$chr2,$end,$type,$size,$ori,$regionsize,$pos1,$pos2) = split /\./, $tigra_sv_name;

    my $posstr = join '_', $chr1, $pos1, $chr2, $pos2, $type, $size, $ori;
    my $head   = join '.', $chr1, $start, $chr2, $end, $type, $size, $ori;

    my $datadir = $self->_tigra_data_dir;
    my $ref_fa  = $datadir . "/$head.ref.fa"; 
    my ($tigra_sv_fa, $cm_out);

    if ($ctg_type eq 'homo') {
        $tigra_sv_fa = $datadir . "/$tigra_sv_name.fa.contigs.fa";
        $cm_out      = $datadir . "/$tigra_sv_name.stat";
    }
    elsif ($ctg_type eq 'het') {
        $tigra_sv_fa = $datadir . "/$tigra_sv_name.fa.contigs.het.fa";
        $cm_out      = $datadir . "/$tigra_sv_name.het.stat";
    }
    else {
        $self->error_message("Wrong type: $ctg_type");
        die $self->error_message;
    }
    
    unless (-s $tigra_sv_fa) {
        $self->warning_message("tigra sv fasta: $tigra_sv_fa is not valid. Skip this $ctg_type cross_match run");
        return;
    }

    my $cm_cmd_opt = $self->crossmatch_options;
    my $cm_cmd = join(' ', $self->crossmatch_path, $tigra_sv_fa, $ref_fa, $cm_cmd_opt, '>', $cm_out, '2>/dev/null');
    my $rv = Genome::Sys->shellcmd (
        cmd           => $cm_cmd,
        input_files   => [$tigra_sv_fa, $ref_fa],
    );
    
    unless ($rv) {
        $self->error_message("Running cross_match for $tigra_sv_name homo failed.\nCommand: $cm_cmd");
        die $self->error_message;
    }
    $self->status_message("Cross_match for $type contigs: $tigra_sv_name Done");
        
    my $makeup_size      = 0; # by default they are zero
    my $concatenated_pos = 0; # by default they are zero

    if ($size > 99999) {
        $makeup_size      = ($end - $self->flank_size) - ($start + $self->flank_size) - 1;
        $concatenated_pos = 2 * $self->flank_size + $self->pad_local_ref;
    }

    $self->status_message("GetCrossMatchIndel for $tigra_sv_name $type");
    my $cm_indel = Genome::Model::Tools::Sv::CrossMatchForIndel->create(
        cross_match_file     => $cm_out,
        local_ref_seq_file   => $ref_fa,
        assembly_contig_file => $tigra_sv_fa,
        per_sub_rate         => $self->cm_max_polymorphism_rate,
        ref_start_pos        => $posstr,
        ref_cat_pos          => $concatenated_pos,
        variant_size         => $makeup_size,
    );
    my $result = $cm_indel->execute;

    # unless explicitly deleted, this rather large object will stay in the UR
    # cache forever
    $cm_indel->delete;
    $cm_indel = undef;

    if ($result && $result =~ /\S+/) {
	    $self->_UpdateSVs($result,$makeup_size,$regionsize,$tigra_sv_fa,$ctg_type, $cm_out);
    }
    
    return 1;
}


sub _UpdateSVs{
    my ($self,$result,$makeup_size,$regionsize,$tigra_sv_fa,$type, $cm_out) = @_;
    #my $datadir = $self->output_directory;
    my $maxSV   = $self->_maxSV;
    my $N50size = $self->_N50size;
    my $depthWeightedAvgSize = $self->_WeightAvgSize;

    if (defined $result) {
        print STDERR "Result: $result\n";
        my ($pre_chr1,$pre_start1,$pre_chr2,$pre_start2,$ori,$pre_bkstart,$pre_bkend,$pre_size,$pre_type,$pre_contigid,$alnscore,$scar_size,$read_len,$fraction_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$strand,$microhomology,$alnstrs) = split /\s+/, $result;
        #$pre_size += $makeup_size if $n_seg >= 2;
        if (defined $pre_size && defined $pre_start1 && defined $pre_start2) {
            my ($contigseq,$contiglens,$contigcovs,$kmerutil) = _GetContig($tigra_sv_fa, $pre_contigid);
            my ($refpos1, $refpos2, $rpos1, $rpos2) = _GetRefPos($cm_out, $pre_contigid,$pre_size,$pre_type);
            $alnscore = int($alnscore*100/$regionsize); 
            $alnscore = $alnscore>100 ? 100 : $alnscore;
            if (!defined $maxSV || $maxSV->{size}<$pre_size || $maxSV->{alnscore} < $alnscore) {
	            my $N50score = int($N50size*100/$regionsize); 
                $N50score = $N50score>100 ? 100 : $N50score;
	            ($maxSV->{chr1},$maxSV->{start1},$maxSV->{chr2},$maxSV->{start2},$maxSV->{bkstart},$maxSV->{bkend},$maxSV->{size},$maxSV->{type},$maxSV->{contigid},$maxSV->{contig},$maxSV->{contiglens},$maxSV->{contigcovs},$maxSV->{kmerutil},$maxSV->{N50},$maxSV->{weightedsize},$maxSV->{alnscore},$maxSV->{scarsize},$maxSV->{a},$maxSV->{b},$maxSV->{read_len},$maxSV->{fraction_aligned},$maxSV->{n_seg},$maxSV->{n_sub},$maxSV->{n_indel},$maxSV->{nbp_indel},$maxSV->{strand},$maxSV->{microhomology},$maxSV->{refpos1},$maxSV->{refpos2},$maxSV->{rpos1},$maxSV->{rpos2}) = ($pre_chr1,$pre_start1,$pre_chr2,$pre_start2,$pre_bkstart,$pre_bkend,$pre_size,$pre_type,$pre_contigid,$contigseq,$contiglens,$contigcovs,$kmerutil,$N50score,$depthWeightedAvgSize,$alnscore,$scar_size,'50','100',$read_len,$fraction_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$strand,$microhomology,$refpos1,$refpos2,$rpos1,$rpos2);
                print STDERR "\n\nWeightedsize: $maxSV->{weightedsize}\n";
	            $maxSV->{het}     = $type;
	            $maxSV->{ori}     = $ori;
	            $maxSV->{alnstrs} = $alnstrs;

                    # add according to Ken's requirement, now cut and add it in CrossMatchIndel.pm for only INDEL
#                    if($maxSV->{strand} =~ /-/ && $maxSV->{type} =~ /DEL/i){
#                        $maxSV->{start2} -= 1;
#                        $maxSV->{bkend} -= 1;
#                    }
            }
        }
    }
    $self->_maxSV($maxSV);
    return 1;
}

# only good for small INDELs
sub _GetRefPos{
    my ($fin, $contigid, $pre_size, $pre_type) = @_;
    my ($ref_start, $ref_end, $r_start, $r_end);
    my ($ref_start1, $ref_end1, $r_start1, $r_end1);
    my ($refpos1, $refpos2, $rpos1, $rpos2);
    my ($chr, $ref_1, $ref_2);
    my $num = 0;
    my ($r1, $r2, $r3, $type, $str, $g1, $g2, $g3);
    my $check = 0;
    my $check_results = 0;
    open(CM,"<$fin") || die "unable to open $fin\n";
    while(<CM>){
        # go ahead to check DISCREPANCY contains I or D directly after find a matched ALIGNMENT
        if($check == 1 && $_ =~ /^DISCREPANCY/){
            if($pre_type =~ /DEL/){
                if($_ =~ /D-$pre_size/){
                    $check_results = 1;
                    last;
                }
            }
            elsif($pre_type =~ /INS/){
                if($_ =~ /I-$pre_size/){
                    $check_results = 1;
                    last;
                }
            }
        }
        else{
            $check = 0;
        }

        if($_ =~ /^ALIGNMENT/ && $_ =~ /$contigid/){
            my @a = split(/\s+/, $_);
            next if($a[5] ne $contigid);

            if($#a == 13){
                ($r1, $r2, $r3, $type, $str, $g1, $g2, $g3) = ($a[$#a - 7], $a[$#a - 6], $a[$#a - 5], $a[$#a - 4], $a[$#a - 3], $a[$#a - 2], $a[$#a - 1], $a[$#a]);
            }
            else{
                ($r1, $r2, $r3, $str, $g1, $g2, $g3) = ($a[$#a - 6], $a[$#a - 5], $a[$#a - 4], $a[$#a - 3], $a[$#a - 2], $a[$#a - 1], $a[$#a]);
            }
            
            ($chr, $ref_1, $ref_2) = ($str =~ /(\S+):(\d+)-(\d+)/);
            $check = 1;
        }
    }

    if($check_results == 1){
        ($refpos1, $refpos2) = ($g2, $g3) if($g1 =~ /\(/);
        ($refpos1, $refpos2) = ($g1, $g2) if($g3 =~ /\(/);
        ($rpos1, $rpos2) = ($r1, $r2) if($r3 =~ /\(/);
        ($rpos1, $rpos2) = ($r2, $r3) if($r1 =~ /\(/);
    }
    else{
        print STDERR "Error, didn't find $contigid of $pre_type and $pre_size in $fin";
    }
    $refpos1 += $ref_1;
    $refpos2 += $ref_1;

    return ($refpos1, $refpos2, $rpos1, $rpos2);
}


sub _GetContig{
    my ($fin, $contigid) = @_;
    my $in = Bio::SeqIO->newFh(-file => $fin, '-format' => 'Fasta');
    my $sequence;
    my @info;
    while (my $seq = <$in>) {
    # do something with $seq
        next unless $seq->id eq $contigid;
        @info = split /\s+/, $seq->desc;
        $sequence = $seq->seq();
        last;
    }
    return ($sequence, $info[0], $info[1], $info[$#info]);
}


sub _ComputeTigraN50{
    my ($contigfile) = @_;
    my @sizes;
    my $totalsize = 0;

    my $fh = Genome::Sys->open_file_for_reading($contigfile) or return;
    while (my $l = $fh->getline){
        chomp $l;
        next unless $l =~ /^\>/;
        next if $l =~ /\*/;
        my ($id, $size, $depth, $ioc, @extra) = split /\s+/, $l;
        next if $size <= 50 && $depth < 3 && $ioc =~ /0/;
        push @sizes, $size;
        $totalsize += $size;
    }
    $fh->close;
    my $cum_size = 0;
    my $halfnuc_size = $totalsize/2;
    my $N50size;
    @sizes = sort {$a<=>$b} @sizes;
    while ($cum_size < $halfnuc_size) {
        $N50size = pop @sizes;
        $cum_size += $N50size;
    }
    return $N50size;
}


sub _ComputeTigraWeightedAvgSize{
    my $contigfile = shift;
    my $totalsize  = 0;
    my $totaldepth = 0;

    my $fh = Genome::Sys->open_file_for_reading($contigfile) or return;
    while (my $l = $fh->getline) {
        chomp $l;
        next unless $l =~ /^\>/;
        next if $l =~ /\*/;
        my ($id,$size,$depth,$ioc,@extra) = split /\s+/, $l;
        next if $size <= 50 && (($depth<3 && $ioc=~/0/) || $depth>500);  #skip error tips or extremely short and repetitive contigs
        $l = $fh->getline;
        chomp $l;
        #$_=<CF>; chomp;
        #next if($size<=50 && (/A{10}/ || /T{10}/ || /C{10}/ || /G{10}/));  #ignore homopolymer contig
        next if $size <= 50 && $l =~ /A{10}|T{10}|C{10}|G{10}/;
        $totalsize += $size*$depth;
        $totaldepth+= $depth;
    }
    $fh->close;
    return $totaldepth > 0 ? int($totalsize*100/$totaldepth)/100 : 0;
}

sub _create_bed_file {
    return 1;
}

sub params_for_filter_result {
    my $self = shift;
    my ($params) = $self->SUPER::params_for_filter_result;

    $params->{chromosome_list} = $self->specify_chr;
    return $params;
}

1;
