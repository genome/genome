package Genome::Model::Tools::DetectVariants2::Polymutt;

use strict;
use warnings;

use FileHandle;

use Genome;
use Workflow;
use Workflow::Simple;
use File::Basename;
class Genome::Model::Tools::DetectVariants2::Polymutt {
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    has_param => [
        lsf_queue=> {
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKFLOW},
        },
    ],
};

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt detect-variants2 polymutt \\
    --version 0.01 \\
    --reference-build NCBI-human-build36 \\
    --pedigree-file-path /tmp/my.ped \\
    --alignment-results id:116553088/116553238/116553281 \\
    --output-directory /tmp/myresults 

gmt detect-variants2 polymutt \\
    --version 0.01 \\
    --reference-build GRCh37-lite-build37 \\
    --pedigree-file-path /gscmnt/gc2146/info/medseq/test_polymutt_pipeline/13211.ped \\
    --alignment-results id:116159911/116159901/116160033/116159884 \\
    --output-directory /gscmnt/gc4095/info/polymutt-test/t1
    
EOS
}

sub help_detail {
    return <<EOS 
This tool runs Polymutt for detection of SNPs and/or indels.
EOS
}

sub _supports_cross_sample_detection {
    my ($class, $version, $vtype, $params) = @_;
    return 1;
};

#sub output_dir {
#    # Q: "where am i getting output directory from";
#    # A: the base class 2 up ::Base has output_directory
#}

#sub ped_file {
#    Q: "implement passthrough pedfile";
#    A: pedigree_file_path is in the ::Base now
#}

sub _detect_variants {
    my $self = shift;

    $DB::single = 1;

    my $version = $self->version;
    unless ($version) {
        die $self->error_message("A version of Polymutt must be specified");
    }
    my @alignments = $self->alignment_results;
    my @glfs = $self->generate_glfs(@alignments);
    my $dat_file = $self->generate_dat();
    my $glf_index = $self->generate_glfindex(@glfs);
    $self->run_polymutt($dat_file, $glf_index);
    return 1;
}

# Parse the chr2process parameter out of the params string
sub chr2process {
    my $self = shift;
    my $params = $self->params;

    my ($chrom_string) = $params =~ m/chr2process (\S+)/;

    # TODO do some checking?

    return $chrom_string;
}

# Parse the skip_denovo parameter out of the params string
sub skip_denovo {
    my $self = shift;
    my $params = $self->params;
    my ($skip_denovo) = $params =~ m/skip_denovo/;

    return defined $skip_denovo; # Treat skip_denovo as a boolean. If it was mentioned, it is set.
}

sub run_polymutt {
    my($self, $dat_file, $glf_index) = @_;
    my $ped_file = $self->pedigree_file_path;
    my %inputs;
    $inputs{version}=$self->version;
    $inputs{dat_file}=$dat_file;
    $inputs{ped_file}=$ped_file;
    $inputs{denovo}=1;
    $inputs{glf_index}=$glf_index;
    #FIXME: if we intend to make one ped per project and just run subsets with it, this could would need to be more complex
    #i.e. if 10 families reside in one pedfile, and we supply a glfindex for just one of those families, this code is bad
    chomp(my $family_id = `head -n 1 $ped_file | cut -f 1`); 
    #FIXME:
    $inputs{output_denovo} = $self->output_directory . "/snvs.denovo.vcf.gz";
    $inputs{output_standard} = $self->output_directory . "/snvs.standard.vcf.gz";
    $inputs{output_merged}= $self->output_directory . "/snvs.vcf.gz";
    $inputs{chr2process} = $self->chr2process;
    my $workflow = Workflow::Model->create(
        name=> "Run polymutt standard and denov",
        input_properties => [
        'dat_file',
        'glf_index',
        'ped_file',
        'output_denovo',
        'output_standard',
        'output_merged',
        'denovo',
        'version',
        'chr2process',
        ],
        output_properties => [
        'output',
        ],
    );

    my $denovo_op;
    my @polymutt_operations;
    if (!$self->skip_denovo) {
        $denovo_op = $workflow->add_operation(
            name=>"denovo polymutt",
            operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Relationship::RunPolymutt"),
        );
        push @polymutt_operations, $denovo_op;
    }

    my $standard_op = $workflow->add_operation(
        name=>"standard polymutt",
        operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Relationship::RunPolymutt"),
    );
    push @polymutt_operations, $standard_op;

    for my $op (@polymutt_operations) {
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"dat_file",
            right_operation=>$op,
            right_property=>"dat_file",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"glf_index",
            right_operation=>$op,
            right_property=>"glf_index",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"ped_file",
            right_operation=>$op,
            right_property=>"ped_file",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"version",
            right_operation=>$op,
            right_property=>"version",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"chr2process",
            right_operation=>$op,
            right_property=>"chr2process",
        );

    }

    $workflow->add_link(
        left_operation=>$workflow->get_input_connector,
        left_property=>"output_standard",
        right_operation=>$standard_op,
        right_property=>"output_vcf",
    );
    my $merge_op = $workflow->add_operation(
        name=>"merge standard and denovo vcfs",
        operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Relationship::MergeAndFixVcfs"),
    );

    $workflow->add_link(
        left_operation=>$standard_op,
        left_property=>"output_vcf",
        right_operation=>$merge_op,
        right_property=>"standard_vcf",
    );
    $workflow->add_link(
        left_operation=>$workflow->get_input_connector,
        left_property=>'output_merged',
        right_operation=>$merge_op,
        right_property=>'output_vcf',
    );
     $workflow->add_link(
        left_operation=>$merge_op,
        left_property=>'output_vcf',
        right_operation=>$workflow->get_output_connector,
        right_property=>'output',
    );

    if (!$self->skip_denovo) {
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"output_denovo",
            right_operation=>$denovo_op,
            right_property=>"output_vcf",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"denovo",
            right_operation=>$denovo_op,
            right_property=>"denovo",
        );
        $workflow->add_link(
            left_operation=>$denovo_op,
            left_property=>"output_vcf",
            right_operation=>$merge_op,
            right_property=>"denovo_vcf",
        );
    }
   

    my $log_dir = $self->output_directory;
    if(Workflow::Model->parent_workflow_log_dir) {
        $log_dir = Workflow::Model->parent_workflow_log_dir;
    }
    $workflow->log_dir($log_dir);

    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }
    $self->status_message("Now launching 2 polymutt jobs");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %inputs);
    unless($result) {
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        $self->error_message("parallel polymutt did not return correctly.");
        die;
    }

}

sub generate_dat {
    my $self = shift;
    my $out_file_name = $self->output_directory . "/" . "polymutt.dat";
    my $dat_fh = IO::File->new($out_file_name, ">");
    $dat_fh->print("T\tGLF_Index\n"); #this is required because the pedfile may have arbitrarily many attributes after the first 5 columns, which you can label and polymutt will ignore. but we never do that, so it always looks like this
    $dat_fh->close;
    return $out_file_name;
}

sub generate_glfindex {
    my $self=shift;
    my @glfs = @_;
    my $glf_name = $self->output_directory ."/" . "polymutt.glfindex";
    my %sample_index = $self->parse_ped();
    my $glf_fh = IO::File->new($glf_name, ">");
    for my $glf (sort @glfs) {  #again we just assume the incoming ped file is alpha sorted, so we alpha sort
        my ($file, $path, $suffix) = fileparse($glf, ".glf");
        my $index = $sample_index{$file}; 
        unless($index) {
            $self->error_message("I am unable to match up the generated glfs with the supplied ped file. Please assist.  I was searching for a sample name match for this file: $glf\n");
            die;
        }
        $glf_fh->print("$index\t$glf\n");
    }
    $glf_fh->close;
    return $glf_name;
}

sub parse_ped {
    my $self = shift;
    my $ped_fh = Genome::Sys->open_file_for_reading($self->pedigree_file_path);
    my %sample_index;
    while(my $ped_line = $ped_fh->getline) {
        chomp($ped_line);
        my ($family_id, $individual, $mother, $father, $sex, $glf_index) = split "\t", $ped_line;
        $sample_index{"$individual"}=$glf_index; 
    }
    $DB::single=1;
    return %sample_index;
}



sub generate_glfs {
    my $self = shift;
    my @alignments = @_;
    my %inputs;
    my (@outputs, @inputs);
    $inputs{ref_fasta} = $alignments[0]->reference_build->full_consensus_path("fa");

#    my $bam_path = $a->merged_alignment_bam_path;
    for (my $i =0; $i < scalar(@alignments); $i++) {
        $DB::single=1;
        my @instrument_data = $alignments[$i]->instrument_data;
        my $output_name = $self->output_directory . "/" . $instrument_data[0]->sample_name . ".glf";
        push @outputs, $output_name;
        $inputs{"bam_$i"}=$alignments[$i]->merged_alignment_bam_path;
        $inputs{"output_glf_$i"}=$output_name;
        push @inputs, ("bam_$i", "output_glf_$i");
    }
    my $workflow = Workflow::Model->create(
        name=> "polymutt parallel glf file creation",
        input_properties => [
        'ref_fasta',
        @inputs,
        ],
        output_properties => [
        'output',
        ],
    );
    for(my $i=0; $i< scalar(@alignments); $i++) {
        my $hybridview_op = $workflow->add_operation(
            name=>"glf creation $i",
            operation_type=>Workflow::OperationType::Command->get("Genome::Model::Tools::Samtools::HybridView"),
        );

        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"ref_fasta",
            right_operation=>$hybridview_op,
            right_property=>"ref_fasta",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"bam_$i",
            right_operation=>$hybridview_op,
            right_property=>"bam",
        );
        $workflow->add_link(
            left_operation=>$workflow->get_input_connector,
            left_property=>"output_glf_$i",
            right_operation=>$hybridview_op,
            right_property=>"output_glf",
        );
        $workflow->add_link(
            left_operation=>$hybridview_op,
            left_property=>"output_glf",
            right_operation=>$workflow->get_output_connector,
            right_property=>"output",
        );
    }
    my @errors = $workflow->validate;
    $workflow->log_dir($self->output_directory);
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }
    $self->status_message("Now launching glf generation jobs");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %inputs);
    unless($result) {
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("parallel glf generation workflow did not return correctly.");
    }

    return @outputs;
}


sub generate_metrics {
    my $self = shift;

    my $metrics = {};
    
    if($self->detect_snvs) {
        my $snp_count      = 0;
        
        my $snv_output = $self->_snv_staging_output;
        my $snv_fh = Genome::Sys->open_file_for_reading($snv_output);
        while (my $row = $snv_fh->getline) {
            $snp_count++;
        }
        $metrics->{'total_snp_count'} = $snp_count;
    }

    if($self->detect_indels) {
        my $indel_count    = 0;
        
        my $indel_output = $self->_indel_staging_output;
        my $indel_fh = Genome::Sys->open_file_for_reading($indel_output);
        while (my $row = $indel_fh->getline) {
            $indel_count++;
        }
        $metrics->{'total indel count'} = $indel_count;
    }

    return $metrics;
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }

    my @versions = Genome::Model::Tools::Relationship::RunPolymutt->available_versions;
    for my $v (@versions){
        if($v eq $version){
            return 1;
        }
    }
    return 0; 
}

sub parse_line_for_bed_intersection {
    my $class = shift;
    my $line = shift;

    unless ($line) {
        die $class->error_message("No line provided to parse_line_for_bed_intersection");
    }

    my ($chromosome, $position, $_reference, $consensus) = split "\t",  $line;

    if ($consensus =~ /\-|\+/) {
        return $class->_parse_indel_for_bed_intersection($line);
    } else {
        return $class->_parse_snv_for_bed_intersection($line);
    }
}

sub _parse_indel_for_bed_intersection {
    my $class = shift;
    my $line = shift;

    my ($chromosome, $position, $_reference, $consensus, @extra) = split "\t",  $line;
    
    my @variants;
    my @indels = Genome::Model::Tools::Bed::Convert::Indel::PolymuttToBed->convert_indel($line);

    for my $indel (@indels) {
        my ($reference, $variant, $start, $stop) = @$indel;
        if (defined $chromosome && defined $position && defined $reference && defined $variant) {
            push @variants, [$chromosome, $stop, $reference, $variant];
        }
    }

    unless(@variants){
        die $class->error_message("Could not get chromosome, position, reference, or variant for line: $line");
    }

    return @variants;
}

sub _parse_snv_for_bed_intersection {
    my $class = shift;
    my $line = shift;

    my ($chromosome, $position, $reference, $consensus, @extra) = split("\t", $line);

    return [$chromosome, $position, $reference, $consensus];
}

# this is a no-op, because we're not using local tmp. 
# this runs as a workflow on LSF and we need network scratch space
sub _promote_staged_data {
    return 1;
}

1;
