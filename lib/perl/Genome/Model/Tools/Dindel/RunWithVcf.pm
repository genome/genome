package Genome::Model::Tools::Dindel::RunWithVcf;

use strict;
use warnings;

use Genome;
use Workflow;
use Workflow::Simple;
use File::Basename;
class Genome::Model::Tools::Dindel::RunWithVcf {
    is => 'Command',
    has_input => [
        input_vcf => {
            is => 'Path',
            doc => 'The Vcf file which has indels called in it',
        },
        ref_fasta => {
            is => 'Path',
            doc => 'The reference sequence that the bam_file is aligned to',
        },
        bam_file => {
            is => 'Path',
            doc => 'The BAM file which contains reads',
        },
        output_directory => {
            is => 'Path',
            doc => 'Where the output should be stored',
        },
    ],
    has_optional_input => [
        roi_bed => {
            is => 'Path',
            doc => "If you'd like to restrict dindel to variants in certain regions, supply a bed file",
        },
        include_dindel_sites => {
            is => 'Boolean',
            default => 0,
            doc => "If you'd like dindel to use all gapped reads as potential indel sites, set this to 1",
        },
        num_windows_per_file => {
            is => 'Number',
            default => '5000',
            doc => 'attempt to control the parallelization of dindel by setting max chunk size.'
        },
        make_realigned_bam => {
            is => 'Boolean',
            default => 0,
            doc => 'Set this to 1 if you wish to output a realigned bam for further counting or manual review',
        },
    ],
};

sub help_brief {
    'Run Dindel on a BAM file and have it realign to the indel sites called in the input-vcf.'
}

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {
    return <<EOS
EOS
}


sub execute {
    my $self = shift;

    unless(-d $self->output_directory) {
        Genome::Sys->create_directory($self->output_directory);
    }

    $DB::single=1;
    my $vcf_in_dindel_format = $self->convert_vcf_to_dindel_and_left_shift($self->output_directory, $self->ref_fasta, $self->input_vcf);
    my ($dindel_var_file, $library_file) = $self->get_cigar_indels($self->output_directory, $self->ref_fasta, $self->bam_file);
#    my $vcf_in_dindel_format = "/gscmnt/gc2146/info/medseq/dindel/test_2368_ALS17_kid/dindel_formatted_vcf_input.left_shifted.variants.txt";
#    my ($dindel_var_file, $library_file) = ("/gscmnt/gc2146/info/medseq/dindel/test_2368_ALS17_kid/cigar_generated_indels.variants.txt", "/gscmnt/gc2146/info/medseq/dindel/test_2368_ALS17_kid/cigar_generated_indels.libraries.txt");
    if($self->include_dindel_sites) {
         $vcf_in_dindel_format = $self->merge_callsets($self->output_directory, $vcf_in_dindel_format, $dindel_var_file);
    }
    if($self->roi_bed) {
        $vcf_in_dindel_format = $self->limit_callset($self->output_directory, $vcf_in_dindel_format, $self->roi_bed);
    }
    my @windows_files = $self->make_windows($self->output_directory, $vcf_in_dindel_format, $self->num_windows_per_file);
    my $results_dir = $self->run_parallel_analysis($self->output_directory, $self->ref_fasta, $self->bam_file, $library_file, $self->make_realigned_bam,\@windows_files);
    my $file_of_results = $self->make_fof($self->output_directory, $results_dir);
    $self->generate_final_vcf($self->output_directory, $file_of_results, $self->ref_fasta);
    if($self->make_realigned_bam) {
        $self->generate_final_bam($self->output_directory, $results_dir, $self->bam_file);
    }
    return 1;
}

sub limit_callset {
    my ($self, $output_dir, $dindel_variant_file, $roi_bed) = @_;
    $DB::single=1;
    my ($file, $path, $suffix) = fileparse($dindel_variant_file, ".variants.txt");
    my $bed_file = $self->dindel_to_bed($dindel_variant_file); #returns a temp file
    my $temp_out = Genome::Sys->create_temp_file_path;
    my $cmd ="intersectBed -a $bed_file -b $roi_bed -wa > $temp_out";
    Genome::Sys->shellcmd(cmd=>$cmd);
    my $dindel_file = $self->bed_to_dindel($temp_out); #returns a temp file
    my $output_file = "$path/$file.region_limited.variants.txt";
    Genome::Sys->shellcmd(cmd=>"cp $dindel_file $output_file");
    return $output_file;
}

sub dindel_to_bed {
    my ($self, $dindel_variant_file) = @_;
    my ($temp_fh, $temp_path) = Genome::Sys->create_temp_file();
    my $fh =  IO::File->new($dindel_variant_file);
    while(my $line = $fh->getline) {
        my ($chr, $pos, $allele, $pound, $number) = split /\s+/, $line;
        my $start = $pos -1; #FIXME should be actual bed, not this hack.  It could conceivably include some additional deletions on the edges of the target region.  need to fix bed to dindel at same time
        $temp_fh->print("$chr\t$start\t$pos\t$allele $pound $number\n");
    }
    $temp_fh->close;
    return $temp_path;
}

sub bed_to_dindel {
    my ($self, $bed_file) = @_;
     my ($temp_fh, $temp_path) = Genome::Sys->create_temp_file();
    my $fh =  IO::File->new($bed_file);
    while(my $line = $fh->getline) {
        my ($chr, $start, $stop, $string) = split /\t/, $line;
        $temp_fh->print("$chr $stop $string"); ##FIXME: if you fixed dindel to bed then this should also be fixed
    }
    $temp_fh->close;
    return $temp_path;
}


sub merge_callsets {
    my ($self, $output_dir, $dindel_format_vcf, $dindel_variants) = @_;
    my $output_file = $output_dir . "/merged_cigar_vcf.variants.txt";
    my $cmd = "cat $dindel_format_vcf $dindel_variants > $output_file";
    Genome::Sys->shellcmd(cmd=>"$cmd"); #larson says theres a Genome::sys->cat, but f that guy and that method
    return $output_file;
}



sub generate_final_bam {
    my ($self, $output_dir, $results_dir, $bam_file) = @_;
    my @sams = glob ("$results_dir/*.merged.sam");
    my $output_sam = $output_dir . "/all_events_realigned.merged.sam";
    my $output_bam =  $output_dir . "/all_events_realigned.merged.bam";
    my $header_cmd = "samtools view -H $bam_file > $output_sam";
    print "$header_cmd\n";
    Genome::Sys->shellcmd(cmd=> $header_cmd);
    for my $sam (@sams) {
        my $dump_cmd = "cat $sam >> $output_sam";
        print "$dump_cmd\n";
        Genome::Sys->shellcmd(cmd=>$dump_cmd);
    }
    my $compress_cmd = "samtools view -S $output_sam -1 > $output_bam";
    print "$compress_cmd\n";
    Genome::Sys->shellcmd(cmd=>$compress_cmd);
    my $final_bam = "$output_dir/all_events_realigned.merged.sorted";
    my $sort_cmd = "samtools sort $output_bam $final_bam";
    Genome::Sys->shellcmd(cmd=>$sort_cmd);
    unlink($output_sam);
    unlink($output_bam);
    unlink(@sams);
    return $final_bam;
}


sub generate_final_vcf {
    my ($self, $output_dir, $file_of_results, $ref_fasta) = @_;
    my $output_vcf = $output_dir . "/final_result.vcf";
    my $merger = Genome::Model::Tools::Dindel::MergeDindelOutput->create(
        dindel_file_output_list=>$file_of_results,
        output_file=>$output_vcf,
        ref_fasta=>$ref_fasta,
    );
    $merger->execute();
}

sub make_fof {
    my ($self, $output_dir, $results_dir) = @_;
    my $fof = $output_dir . "/file_of_result_files";
    my $fof_fh = IO::File->new($fof, ">");
    my @files = glob("$results_dir/result*txt");
    for my $file (@files) {
        $fof_fh->print($file ."\n");
    }
    $fof_fh->close;
    return $fof;
}

sub run_parallel_analysis {
    my ($self, $output_dir, $ref_fasta, $bam_file, $library_file, $make_realigned_bams, $window_files) = @_;
    my $results_dir = $output_dir . "/results/";
    $DB::single=1;
    unless(-d $results_dir) {
        Genome::Sys->create_directory($results_dir);
    }
    my %inputs;
    $inputs{bam_file}=$bam_file;
    $inputs{ref_fasta}=$ref_fasta;
    $inputs{library_metrics_file}=$library_file;
    $inputs{output_bam}=0;
    if($make_realigned_bams) {
        $inputs{output_bam}=1;
    }


    my @inputs;
    my @prefixes;
    for my $window (@$window_files) {
        $DB::single=1;
        my ($number) = ($window =~m/(\d+)\.txt/);
        push @inputs, "window_file_$number";
        push @inputs, "result_prefix_$number";
        $inputs{"window_file_$number"}= $window;
        $inputs{"result_prefix_$number"}= $results_dir . "/result_$number";
    }
    my $workflow = Workflow::Model->create(
        name=>"highly parallel dindel :-(",
        input_properties=> [
        'bam_file',
        'ref_fasta',
        'library_metrics_file',
        'output_bam',
        @inputs,
        ],
        output_properties=> [
        'output',
        ],
    );
    $workflow->log_dir($output_dir);
    for my $window (@$window_files) {
        my ($number) = ($window =~m/(\d+)\.txt/);

        my $analyze_op = $workflow->add_operation(
            name=>"dindel analyze window $number",
            operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Dindel::AnalyzeWindowFile"),
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"bam_file",
            right_operation=>$analyze_op,
            right_property=>"bam_file",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"output_bam",
            right_operation=>$analyze_op,
            right_property=>"output_bam",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"ref_fasta",
            right_operation=>$analyze_op,
            right_property=>"ref_fasta",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"library_metrics_file",
            right_operation=>$analyze_op,
            right_property=>"library_metrics_file",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"window_file_$number",
            right_operation=>$analyze_op,
            right_property=>"window_file",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"result_prefix_$number",
            right_operation=>$analyze_op,
            right_property=>"output_prefix",
        );
        $workflow->add_link(
            left_operation=>$analyze_op,
            left_property=>"output_prefix",
            right_operation=>$workflow->get_output_connector,
            right_property=>"output",
        );
    }
    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }   
    $self->status_message("Now launching a butt-ton of dindel jobs");
    $DB::single=1;
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %inputs);
    unless($result) {
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("parallel mpileup workflow did not return correctly.");
    } 
    return $results_dir; 
} 

sub make_windows {
    my ($self, $output_dir, $candidate_indel_file, $num_windows) =@_;
    my $window_dir = $output_dir . "/windows/";
    unless(-d $window_dir) {
        Genome::Sys->create_directory($window_dir);
    }
    my $window_prefix = $window_dir ."/dindel_window_file";
    my $make_dindel_windows = Genome::Model::Tools::Dindel::MakeDindelWindows->create(
        input_dindel_file=>$candidate_indel_file,
        output_prefix=>$window_prefix,
        num_windows_per_file=>$num_windows,
    );
    if($make_dindel_windows->execute()) {
        return glob("$window_prefix*");
    }
    else {
        $self->error_message("Dindel Window maker failed for some reason. Exiting.");
        die;
    }
}

sub get_cigar_indels {
    my ($self, $output_dir, $ref_fasta, $bam_file) = @_;
    my $get_cigar_indels = Genome::Model::Tools::Dindel::GetCigarIndels->create(
        input_bam=>$bam_file,
        output_directory=>$output_dir,
        ref_fasta=>$ref_fasta
    );
    if($get_cigar_indels->execute()) {
        my $output_variants = $get_cigar_indels->output_variants;
        my $output_libs = $get_cigar_indels->output_libraries;
        return ($output_variants, $output_libs);
    }
    else {
        $self->error_message("fail from GetCigarIndels...exiting");
        die;
    }
}

sub convert_vcf_to_dindel_and_left_shift {
    my ($self, $output_dir, $ref_fasta, $input_vcf) = @_;
    my $output_dindel_file = $output_dir ."/dindel_formatted_vcf_input";
    my $vcf_to_dindel = Genome::Model::Tools::Dindel::VcfToDindel->create(
        input_vcf=> $input_vcf,
        output_dindel_file=>$output_dindel_file,
        ref_fasta=>$ref_fasta
    );
    $vcf_to_dindel->execute();
    my $left_shifted_output = $output_dindel_file . ".left_shifted";
    my $left_shifter = Genome::Model::Tools::Dindel::RealignCandidates->create( 
        ref_fasta=>$ref_fasta,
        variant_file=>$output_dindel_file,
        output_prefix=>$left_shifted_output,
    );
    if($left_shifter->execute()) {
        return $left_shifter->output_file;
    }
    else {
        $self->error_message("the left shifter has caused the datacenter to burn to the ground. great job. was it worth not getting those indels?");
        die;
    }
}
1;
