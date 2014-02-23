package Genome::Model::Tools::Validation::LongIndelsGenerateMergedAssemblies;

use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::Model::Tools::Validation::LongIndelsGenerateMergedAssemblies {
    is => 'Command::V2',
    has_input => [
        long_indel_bed_file => {
            is => 'String',
            doc => 'unsorted, unannotated 3bp indel file in BED format!!! BED format!!!',
        },
        output_dir => {
            is => 'String',
            doc => 'directory for output files',
        },
        reference_transcripts => {
            is => 'String',
            doc => 'reference transcripts plus version to be used to annotate input indel file',
        },
        transcript_variant_annotator_version => {
            is => 'Number',
            doc => 'Version of the "annotate transcript-variants" tool to run during the annotation step',
        },
        normal_bam => {
            is => 'String',
        },
        tumor_bam => {
            is => 'String',
        },
        reference_fasta => {
            is => 'String',
        },
    ],
    has_transient_optional_output => [
        contigs_fasta => {
            is => 'String',
        },
    ],
};

sub execute {
    my $self = shift;

    #parse input params, declare vars that need to be in wide scope
    my $indels_full_path = $self->long_indel_bed_file;
    my $output_dir = $self->output_dir;
    my $normal_bam = $self->normal_bam;
    my $tumor_bam = $self->tumor_bam;
    my $ref_seq_fasta = $self->reference_fasta;

    #dedup indels - assumes indel file isn't massive (a million indels would probably be an issue)
    my %indelhash;
    open(OUTFILE,">" . $indels_full_path . ".dedup");
    my $inFh = IO::File->new( $indels_full_path ) || die "can't open file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        my @F = split("\t",$line);
        $F[3] =~ s/0/-/g;
        $F[4] =~ s/0/-/g;
        $line = join("\t",@F);
        unless(exists($indelhash{$line})){
            print OUTFILE $line . "\n";
        }
        $indelhash{$line} = 0;        
    }
    close($inFh);
    
    #sort indels
    my ($indels_filename_only) = fileparse($indels_full_path) .  ".dedup";
    my $sort_output = $output_dir . "/" . $indels_filename_only . ".sorted";
    my $sort_cmd = Genome::Model::Tools::Snp::Sort->create(
        output_file => $sort_output,
        snp_file => $indels_full_path,
        force => 1,
    );
    unless ($sort_cmd->execute) {
        die "Sort of indels failed.\n";
    }
    $sort_cmd->delete;

    #annotate indels
    my $anno_output = $sort_output . ".anno";
    my $anno_cmd = Genome::Model::Tools::Annotate::TranscriptVariants->create(
        output_file => $anno_output,
        annotation_filter => "top",
        variant_bed_file => $sort_output,
        use_version => $self->transcript_variant_annotator_version,
        #variant_file => $sort_output,
        reference_transcripts => $self->reference_transcripts,
    );
    unless ($anno_cmd->execute) {
        die "Annotation of sorted indels failed.\n";
    }
    $anno_cmd->delete;

    #prepare assembly inputs
    my $assembly_input = $anno_output . ".assembly_input";
    my $prepare_ass_input_cmd = Genome::Model::Tools::Validation::AnnotationToAssemblyInput->create(
        annotation_file => $anno_output,
        output_file => $assembly_input,
        minimum_size => 1,
    );
    unless ($prepare_ass_input_cmd->execute) {
        die "annotation-to-assembly-input failed.\n";
    }

    my $assembly_input_walleles = $anno_output . ".assembly_input_walleles";
    $prepare_ass_input_cmd = Genome::Model::Tools::Validation::AnnotationToAssemblyInput->create(
        annotation_file => $anno_output,
        output_file => $assembly_input_walleles,
        add_indel_alleles => 1,
        minimum_size => 1,
    );
    unless ($prepare_ass_input_cmd->execute) {
        die "annotation-to-assembly-input w/alleles failed.\n";
    }
    $prepare_ass_input_cmd->delete;

    #--------------------------------------------------------------------
    #run tigra on the list of predicted indels in the normal BAM
    my $normal_output_file = $output_dir . "/normal.csv";
    my $normal_breakpoint_file = $output_dir . "/normal.bkpt.fa";
    print "Command for gmt sv assembly-validation\nbam-files $normal_bam\noutput-file $normal_output_file\nsv-file $assembly_input\nmin-size-of-confirm-asm-sv 3\nflank_size 200\nbreakpoint_seq_file $normal_breakpoint_file\n--asm-high-coverage\n--reference-file $ref_seq_fasta\n"; ##TEST
    $self->debug_message("Doing assembly validation for normal");
    my $normal_assembly_cmd = Genome::Model::Tools::Sv::AssemblyValidation->create(
        bam_files => $normal_bam,
        output_file =>  $normal_output_file,
        sv_file => $assembly_input,
        min_size_of_confirm_asm_sv => '3',
        flank_size => '200',
        breakpoint_seq_file => $normal_breakpoint_file,
        asm_high_coverage => '1',
        reference_file => $ref_seq_fasta,
    );
    #print STDERR "Normal: gmt sv assembly-validation --bam-files $normal_bam --output_file $normal_output_file --sv_file $assembly_input --min-size-of-confirm-asm-sv 3 --flank-size 200 --breakpoint-seq-file $normal_breakpoint_file --asm-high-coverage --reference-file $ref_seq_fasta\n";
    unless ($normal_assembly_cmd->execute) {
        die "Normal SV assembly-validation failed (normal.bkpt.fa compromised).\n";
    }
    $normal_assembly_cmd->delete;

    $self->debug_message("Doing assembly validation for tumor");
    #run tigra on the list of predicted indels in the tumor BAM
    my $tumor_output_file = $output_dir . "/tumor.csv";
    my $tumor_breakpoint_file = $output_dir . "/tumor.bkpt.fa";
    my $tumor_assembly_cmd = Genome::Model::Tools::Sv::AssemblyValidation->create(
        bam_files => $tumor_bam,
        output_file =>  $tumor_output_file,
        sv_file => $assembly_input,
        min_size_of_confirm_asm_sv => '3',
        flank_size => '200',
        breakpoint_seq_file => $tumor_breakpoint_file,
        asm_high_coverage => '1',
        reference_file => $ref_seq_fasta,
    );
    #print STDERR "Tumor: gmt sv assembly-validation --bam-files $tumor_bam --output_file $tumor_output_file --sv_file $assembly_input --min-size-of-confirm-asm-sv 3 --flank-size 200 --breakpoint-seq-file $tumor_breakpoint_file --asm-high-coverage --reference-file $ref_seq_fasta\n";
    unless ($tumor_assembly_cmd->execute) {
        die "Tumor SV assembly-validation failed (tumor.bkpt.fa compromised).\n";
    }
    $tumor_assembly_cmd->delete;

    #build contigs for remapping based on the assembly results
    my $contigs_file = $output_dir . "/contigs.fa";
    $self->debug_message("Writing contigs to $contigs_file");
    my $contig_cmd = Genome::Model::Tools::Validation::BuildRemappingContigs->create(
        input_file => $assembly_input_walleles,
        normal_assembly_file => $normal_output_file,
        tumor_assembly_file => $tumor_output_file,
        normal_assembly_breakpoints_file => $normal_breakpoint_file,
        tumor_assembly_breakpoints_file => $tumor_breakpoint_file,
        output_file => $contigs_file,
        contig_size => '500',
        append_indel_alleles => 1,
        reference_sequence => $ref_seq_fasta,
    );
    unless ($contig_cmd->execute) {
        die "Failed to build contigs for remapping.\n";
    }

    $self->contigs_fasta($contigs_file);
    return 1;
}

1;

