package Genome::Model::Tools::Predictor::Ber;

use strict;
use warnings;
use Genome;
use Workflow::Simple;
use Command;
use Getopt::Long;

use Carp;
use File::Path qw(remove_tree);
use File::Spec;
use Cwd;
use English;

BEGIN {
    $ENV{WF_USE_FLOW} = 1;
}

class Genome::Model::Tools::Predictor::Ber {
    is => 'Genome::Model::Tools::Predictor::Base',
    doc => 'execute BER gene predictor',
    has => [
        ber_version => {
            is => 'Text',
            value => 'v2.5',
            doc => 'Version of Ber to use',
        },
        version => {
            is => 'Text',
            value => 'v2.5',
            doc => 'Version of Ber to use',
        },
        cleanup => {
            is => 'Boolean',
            default_value => 1,
            doc => 'cleanup intermediate files',
        },
        ber_source_path => {
            is => 'DirectoryPath',
            value => '/gscmnt/gc9002/info/annotation/BER/autoannotate_v2.5',
            doc => 'Source of BER',
        },
        lsf_queue => {
            is_param => 1,
            default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
        lsf_resource => {
            is => 'Text',
            default => "-M 4000000 -R 'select[type==LINUX64 && mem>4000] rusage[mem=4000]'",
        },
    ],
    has_transient_optional => [
        _locus_id => {
            is => 'Text',
        },
    ],

};

sub help_brief
{
    "Runs the entire BER product naming pipeline";
}

sub help_synopsis
{
    return <<"EOS"
Runs the entire BER product naming pipeline. This tool pulls data from biosql, runs blastp, hmmpfam, btab, hmmtab and anno-sqlite.bash.
EOS
}

sub help_detail
{
    return <<"EOS"
Runs the entire BER product naming pipeline. This tool pulls data from biosql, runs blastp, hmmpfam, btab, hmmtab and anno-sqlite.bash.  In addition, this script will also create the BER_product_naming.ace file, writes the tace script, parses all ace files into acedb, gathers QC stats from the phase5 file and queries acedb to make sure gene counts match, then writes the rt ticket blurb and mails the user when finished.
EOS
}

sub requires_chunking {
    return 0;
}

sub run_predictor {
    my $self = shift;

    if(-e $self->final_report_file) {
        $self->status_message("Final report already exists\n");
        return 0;
    }

    $self->status_message("Starting BER");
   
    my ($blastp_fasta_files, $hmmpfam_fasta_files) = $self->setup or croak "failed to setup BER config";

    my %inputs = (
        'output_directory' => $self->output_directory,
        'blastp_fasta_files'  => $blastp_fasta_files,
        'hmmpfam_fasta_files' => $hmmpfam_fasta_files,
        'ber_source_path'  => $self->ber_source_path,
        'gram_stain'       => $self->gram_stain ? $self->gram_stain : 'negative',
        'locus_id'         => $self->locus_id,
        'lsf_queue'        => $self->lsf_queue,
        'lsf_resource'     => $self->lsf_resource,

    );

    $self->status_message("found ".@$blastp_fasta_files." fasta files for blastp\n");
    $self->status_message("found ".@$hmmpfam_fasta_files." fasta files for hmmpfam\n");
    my $workflow = $self->generate_work_flow( (@$blastp_fasta_files != 0 ?  1 : 0), @$hmmpfam_fasta_files != 0 ? 1 : 0);
    $self->status_message("starting workflow\n");
    my $result = Workflow::Simple::run_workflow_lsf($workflow, %inputs);
    unless ($result) {
        croak "Error running BER workflow\n" . join("\n", map { $_->name . ": " . $_->error } @Workflow::Simple::ERROR);
    }
    unlink $self->raw_output_path if(-e $self->raw_output_path);
    Genome::Sys->create_symlink($result->{output_file},$self->raw_output_path);
    $self->status_message("finished running BER workflow\n");

    $self->write_final_report;
    $self->cleanup_files if($self->cleanup);

    return 1;
}

sub generate_work_flow {
    my $self = shift;
    my ($blastp, $hmmpfam) = @_;

    my $workflow = Workflow::Model->create(
        name => 'Predictor::Ber',
        input_properties => [ 'output_directory', 
                              'blastp_fasta_files', 
                              'hmmpfam_fasta_files', 
                              'gram_stain', 
                              'ber_source_path', 
                              'locus_id',
                              'lsf_queue',
                              'lsf_resource',
                          ],
        output_properties => ['output_file'],
    );

    unless (-d $self->log_directory) {
        Genome::Sys->create_directory($self->log_directory);
    }
    $workflow->log_dir($self->log_directory);
    my $blastp_operation;
    if ($blastp) {
        $self->status_message("create workflow operation for blastp\n");
    
        $blastp_operation = $workflow->add_operation(
            name => 'BER blastp',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::Tools::Predictor::Ber::Blastp',
            ),
            parallel_by => 'input_fasta_file'
        );
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => 'blastp_fasta_files',
            right_operation => $blastp_operation,
            right_property => 'input_fasta_file',
        );
        
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => 'output_directory',
            right_operation => $blastp_operation,
            right_property => 'output_directory',
        );
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => 'ber_source_path',
            right_operation => $blastp_operation,
            right_property => 'ber_source_path',
        );
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => 'lsf_queue',
            right_operation => $blastp_operation,
            right_property => 'lsf_queue',
        );
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => 'lsf_resource',
            right_operation => $blastp_operation,
            right_property => 'lsf_resource',
        );
    }

    my $hmmpfam_operation;
    if ($hmmpfam) {
        $self->status_message("create workflow operation for hmmpfam\n");
    
        $hmmpfam_operation = $workflow->add_operation(
            name => 'BER hmmpfam',
            operation_type => Workflow::OperationType::Command->create(
                command_class_name => 'Genome::Model::Tools::Predictor::Ber::Hmmpfam',
            ),
            parallel_by => 'input_fasta_file'
        );
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => 'hmmpfam_fasta_files',
            right_operation => $hmmpfam_operation,
            right_property => 'input_fasta_file',
        );
        
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => 'output_directory',
            right_operation => $hmmpfam_operation,
            right_property => 'output_directory',
        );
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => 'ber_source_path',
            right_operation => $hmmpfam_operation,
            right_property => 'ber_source_path',
        );
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => 'lsf_queue',
            right_operation => $hmmpfam_operation,
            right_property => 'lsf_queue',
        );
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => 'lsf_resource',
            right_operation => $hmmpfam_operation,
            right_property => 'lsf_resource',
        );
    }
    my $converge_operation;
    if ($blastp_operation and $hmmpfam_operation) {
        $converge_operation = $workflow->add_operation(
            name => 'converge',
            operation_type => Workflow::OperationType::Converge->create(
                input_properties => ['btab_file', 'htab_file'],
                output_properties => ['result'],
            ),
        );
        $workflow->add_link(
            left_operation => $blastp_operation,
            left_property => 'btab_file',
            right_operation => $converge_operation,
            right_property => 'btab_file',
        );
        $workflow->add_link(
            left_operation => $hmmpfam_operation,
            left_property => 'htab_file',
            right_operation => $converge_operation,
            right_property => 'htab_file',
        );
    }
    elsif ($blastp_operation) {
        $converge_operation = $blastp_operation
    }
    elsif ($hmmpfam_operation) {
        $converge_operation = $hmmpfam_operation
    }

    $self->status_message("create workflow operation for annotation\n");
    my $annotate_operation = $workflow->add_operation(
        name => 'BER annotation',
        operation_type => Workflow::OperationType::Command->create(
            command_class_name => 'Genome::Model::Tools::Predictor::Ber::Annotate',
        ),
    );
    $workflow->add_link(
        left_operation => $converge_operation,
        left_property => 'result',
        right_operation => $annotate_operation,
        right_property => 'converge_result',
    ) if($converge_operation);
    $workflow->add_link(
        left_operation => $workflow->get_input_connector,
        left_property => 'output_directory',
        right_operation => $annotate_operation,
        right_property => 'output_directory',
    );
    $workflow->add_link(
        left_operation => $workflow->get_input_connector,
        left_property => 'ber_source_path',
        right_operation => $annotate_operation,
        right_property => 'ber_source_path',
    );
    $workflow->add_link(
        left_operation => $workflow->get_input_connector,
        left_property => 'ber_source_path',
        right_operation => $annotate_operation,
        right_property => 'ber_source_path',
    );
    $workflow->add_link(
        left_operation => $workflow->get_input_connector,
        left_property => 'locus_id',
        right_operation => $annotate_operation,
        right_property => 'locus_id',
    );
    $workflow->add_link(
        left_operation => $workflow->get_input_connector,
        left_property => 'gram_stain',
        right_operation => $annotate_operation,
        right_property => 'gram_stain',
    );
    $workflow->add_link(
        left_operation => $annotate_operation,
        left_property => 'output_file',
        right_operation => $workflow->get_output_connector,
        right_property => 'output_file',
    );

    my @errors = $workflow->validate;
    if (@errors) {
        die "Could not validate workflow:\n" . join("\n", @errors);
    }
    return $workflow;
}

sub locus_id {
    my $self=shift;
    return defined $self->_locus_id ? $self->_locus_id : $self->_get_locus_from_fasta_header;
}

sub _get_locus_from_fasta_header {
    my $self = shift;
    
    #get headers
    my $output_fh = Genome::Sys->open_file_for_reading($self->input_fasta_file) 
        or die "Could not get file handle for " . $self->input_fasta_file;
    
    my $fasta_header = $output_fh->getline;
    chomp $fasta_header;
    
    #check header format 
    my ($locus_id) = $fasta_header =~ /^\>\w*?\-?(\w+)_Contig/;

    $self->_locus_id($locus_id);
    #get locus id and ensure only one exists
    return $locus_id;
}


sub setup {
    my $self = shift;
    my $locusid = $self->locus_id;
    my $config_dir =$self->output_directory.'/';

    $self->status_message("setting up inputs, creating fasta files and checking hmm and ber files\n");
    #to satisfy the Ber software it expects a directory in its data path
    my $ber_genomes_dir = $self->ber_source_path.'/data/genomes/'.$locusid;
    Genome::Sys->create_symlink($config_dir, $ber_genomes_dir) unless(-l $ber_genomes_dir);
    #and cleanup any previous run of this locus
    unlink $self->ber_source_path.'/data/db/SQLite/'.$locusid if(-e $self->ber_source_path.'/data/db/SQLite/'.$locusid);

    Genome::Sys->create_directory($self->fasta_dir);
    Genome::Sys->create_directory($self->hmm_dir);
    Genome::Sys->create_directory($self->ber_dir);

    my $asm_file = $locusid.'_asm_feature';
    unlink $config_dir.$asm_file if(-e $config_dir.$asm_file);
    my $asm_fh = Genome::Sys->open_file_for_writing($config_dir.$asm_file) 
        or croak "cann't open the file:$!.";
    my $asmbl_file = $locusid.'_asmbl_data';
    unlink $config_dir.$asmbl_file if(-e $config_dir.$asmbl_file);
    my $asmbl_fh = Genome::Sys->open_file_for_writing($config_dir.$asmbl_file) 
        or croak "cann't open the file:$!.";
    my $ident2_file = $locusid.'_ident2';
    unlink $config_dir.$ident2_file if(-e $config_dir.$ident2_file);
    my $ident2_fh = Genome::Sys->open_file_for_writing($config_dir.$ident2_file) 
        or croak "cann't open the file:$!.";
    my $stan_file = $locusid.'_stan';    
    unlink $config_dir.$stan_file if(-e $config_dir.$stan_file);
    my $stan_fh = Genome::Sys->open_file_for_writing($config_dir.$stan_file) 
        or croak "cann't open the file:$!.";

    $asm_fh->print(join("\t", qw(asmbl_id end3 end5 feat_name feat_type))."\r\n");
    $asmbl_fh->print(join("\t", qw(id name type))."\r\n");
    $ident2_fh->print(join("\t", qw(complete feat_name locus))."\r\n");
    $stan_fh->print(join("\t", qw(asmbl_data_id asmbl_id iscurrent))."\r\n");

    my $count=0;
    my $asmblid=0;
    my @blastp_files;
    my @hmmpfam_files;

    my $seq_in = Bio::SeqIO->new(-file => $self->input_fasta_file, -format => 'Fasta')
        or croak "failed to open: ".$self->input_fasta_file;
    while (my $seq = $seq_in->next_seq()) {
        my $contig = $seq->primary_id;
        my ($start, $stop) = split(/\s+/, $seq->desc);
        my $locus_count = sprintf("%05d", $count);
        $asm_fh->print(join("\t", ($count, $start, $stop, $contig, 'ORF'))."\r\n");
        $asmbl_fh->print(join("\t", ($count, 'Contig','contig'))."\r\n");
        $ident2_fh->print(join("\t", (' ', $contig, $locusid.$locus_count))."\r\n");
        $stan_fh->print(join("\t", ($count, $count, 1))."\r\n");
        $count++;

        my $gene_file = $self->fasta_dir.$contig.'.fasta';
        push(@blastp_files, $gene_file) unless(Genome::Model::Tools::Predictor::Ber::Blastp->validate_output(
            Genome::Model::Tools::Predictor::Ber::Blastp->blastp_file_name($config_dir, $contig),
            Genome::Model::Tools::Predictor::Ber::Blastp->btab_file_name($config_dir, $contig)));
        push(@hmmpfam_files, $gene_file) unless(Genome::Model::Tools::Predictor::Ber::Hmmpfam->validate_output(
            Genome::Model::Tools::Predictor::Ber::Hmmpfam->hmmpfam_file_name($config_dir, $contig),
            Genome::Model::Tools::Predictor::Ber::Hmmpfam->htab_file_name($config_dir, $contig)));

        next if(-e $gene_file and -s $gene_file);

        my $seq_out = Bio::SeqIO->new(-file => ">$gene_file", -format => 'Fasta');
        $seq_out->write_seq($seq);
        $seq_out->close;
    }

    $asm_fh->close;
    $asmbl_fh->close;
    $ident2_fh->close;
    $stan_fh->close;

    #link these files into the BER data directories also, BLAH!
    my $ber_db_csv_dir = $self->ber_source_path.'/data/db/CSV/';
    Genome::Sys->create_symlink($config_dir.$asm_file, $ber_db_csv_dir.$asm_file) unless(-l $ber_db_csv_dir.$asm_file);
    Genome::Sys->create_symlink($config_dir.$asmbl_file, $ber_db_csv_dir.$asmbl_file) unless(-l $ber_db_csv_dir.$asmbl_file);
    Genome::Sys->create_symlink($config_dir.$ident2_file, $ber_db_csv_dir.$ident2_file) unless(-l $ber_db_csv_dir.$ident2_file);
    Genome::Sys->create_symlink($config_dir.$stan_file, $ber_db_csv_dir.$stan_file) unless(-l $ber_db_csv_dir.$stan_file); 
    
    return (\@blastp_files, \@hmmpfam_files);
}

sub final_report_file {
    my $self = shift;

    return $self->output_directory.'/'.$self->locus_id.'-final-report.txt';
}

sub write_final_report {
    my $self = shift;

    my $config_dir =$self->output_directory.'/';
    my $dat_line_count = Genome::Sys->line_count($self->raw_output_path);
    my $no_blast_hits=0;
    my $no_domain_hits=0;
    my $total_genes=0;
    my $seq_in = Bio::SeqIO->new(-file => $self->input_fasta_file, -format => 'Fasta')
        or croak "failed to open: ".$self->input_fasta_file;
    while (my $seq = $seq_in->next_seq()) {
        my $contig = $seq->primary_id;

        $no_blast_hits++ if(-z Genome::Model::Tools::Predictor::Ber::Blastp->btab_file_name($config_dir, $contig) and
                                Genome::Model::Tools::Predictor::Ber::Blastp->no_blast_hits(
                                    Genome::Model::Tools::Predictor::Ber::Blastp->blastp_file_name($config_dir, $contig),
                                    Genome::Model::Tools::Predictor::Ber::Blastp->btab_file_name($config_dir, $contig)
                                    )
                                );
        $no_domain_hits++ if(-z Genome::Model::Tools::Predictor::Ber::Hmmpfam->htab_file_name($config_dir, $contig) and
                                 Genome::Model::Tools::Predictor::Ber::Hmmpfam->no_domain_hits(
                                     Genome::Model::Tools::Predictor::Ber::Hmmpfam->hmmpfam_file_name($config_dir, $contig),
                                     Genome::Model::Tools::Predictor::Ber::Hmmpfam->htab_file_name($config_dir, $contig)
                                     )
                                 );
        $total_genes++;
    }
    
    my $final_fh = Genome::Sys->open_file_for_writing($self->final_report_file);
    $final_fh->print("BER Final Report for ".$self->locus_id."\n");
    $final_fh->print("number of input sequences: $total_genes\n");
    $final_fh->print("number of sequences with no blastp hits: $no_blast_hits or ". sprintf("%.2f", $no_blast_hits/$total_genes *100)."\n");
    $final_fh->print("number of sequences with no hmmpfam hits: $no_domain_hits or ". sprintf("%.2f", $no_domain_hits/$total_genes *100)."\n");
    $final_fh->print("number of sequences in BER annotation output: $dat_line_count or ". sprintf("%.2f", $dat_line_count/$total_genes *100)."\n");
    $final_fh->close;
    
    return;
}



# Should contain all code necessary to parse the raw output of the predictor.
sub parse_output {
    return 1;
}

# Any filtering logic should go here.
sub filter_results {
    return 1;
}

# Should create an ace file from the raw output of the predictor
sub create_ace_file {
    my $self = shift;
    unlink $self->ace_file_path if(-e $self->ace_file_path);
    my $fh = Genome::Sys->open_file_for_reading($self->raw_output_path);
    my $ace_fh = Genome::Sys->open_file_for_writing($self->ace_file_path);

    while (my $line = $fh->getline) {
        chomp $line;
        my ($gene, $desc) = split(/\t/, $line);
        $ace_fh->print("Sequence \"$gene\"\nBER_product   \"$desc\"\n\n");
    }

    $fh->close;
    $ace_fh->close;
    return 1;
}

sub cleanup_files {
    my $self = shift;
    my $error = remove_tree($self->fasta_dir,
                             $self->ber_dir,
                             $self->hmm_dir,
                             $self->log_directory);
    croak "failed to remove directories" unless($error);
    
    return 1;
}

sub fasta_dir {
    shift->output_directory.'/fasta/';
}

sub ber_dir {
    shift->output_directory.'/ber/';
}

sub hmm_dir {
    shift->output_directory.'/hmm/';
}


sub log_directory {
    my $self = shift;
    return join('/', $self->output_directory, 'logs');
}

sub debug_output_path {
    my $self = shift;
    return join('/', $self->output_directory, 'ber.debug');
}

sub dump_output_path {
    my $self = shift;
    return join('/', $self->output_directory, 'ber.output');
}

sub ace_file_path {
    my $self = shift;
    return join('/', $self->output_directory, 'ber.ace');
}

sub raw_output_path {
    my $self = shift;
    return join('/', $self->output_directory, 'ber.dat');
}

1;

