package Genome::Model::Tools::Allpaths::Metrics;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::Allpaths::Metrics {
    #is => 'Genome::Model::Tools::Allpaths::Base',
    is => 'Command::V2',
    has => [
	    assembly_directory => {
            is => 'Text',
            doc => 'Path to allpaths assembly.',
        },
        first_tier => {
            is => 'Number',
            doc => 'First tier value',
            is_optional => 1,
        },
        second_tier => {
            is => 'Number',
            doc => 'Second tier value',
            is_optional => 1,
        },
        major_contig_length => {
            is => 'Number',
            is_optional => 1,
            default_value => 500,
            doc => 'Cutoff value for major contig length',
        },
        output_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'Stats output file',
        },
        output_clean_assembly_fasta_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'Clean scaffolds fasta file',
        },
    ],
    has_optional => [
        _metrics => { is_transient => 1, },
        _add_reads_method => { is_transient => 1, },
        _reads_files => { is_transient => 1, },
    ],
};

sub help_brief {
    return 'Produce metrics for allpaths assemblies'
}

sub help_detail {
    return;
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__(@_);
    return @errors if @errors;

    if ( $self->assembly_directory ) {
        if ( not -d $self->assembly_directory ) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ assembly_directory /],
                desc => 'The assembly_directory is not a directory!',
            );
            return @errors;
        }
        if ( not defined $self->output_file ) {
            $self->output_file( $self->assembly_directory."/metrics.out" );
        }
        if ( not defined $self->output_clean_assembly_fasta_file ) {
            $self->output_clean_assembly_fasta_file( $self->assembly_directory.'/clean.final.assembly.fasta' );
        }
    }
    elsif ( not $self->output_file ) { 
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ output_file /],
            desc => 'No output file given and no assembly_directory given to determine the output file!',
        );
    }

    my $scaffolds_efasta = $self->resolve_scaffolds_efasta;
    if ( not $scaffolds_efasta ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ assembly_directory /],
            desc => $self->error_message,
        );
    }

    my @files;
    if ( @files = $self->resolve_reads_metrics_files ) {
        $self->_add_reads_method('add_reads_metrics_file')
    }
    elsif ( @files = $self->resolve_reads_files ) {
        $self->_add_reads_method('add_reads_file')
    }
    else {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ assembly_directory /],
            desc => 'No reads metrics files or reads files in assembly directory: '.$self->assembly_directory,
        );
    }
    $self->_reads_files(\@files);

    return @errors;
}

sub resolve_scaffolds_efasta {
    my $self = shift;

    my @files = glob($self->assembly_directory.'/*/data/*/ASSEMBLIES/*/final.assembly.efasta');
    if ( not @files ) {
        $self->error_message('No scaffold efasta file found in '.$self->assembly_directory);
        return;
    }
    elsif ( @files > 1 ) {
        $self->error_message("More than one scaffold assembly file found!\n".join("\n", @files));
        return;
    }

    return $files[0];
}

sub resolve_reads_files {
    my $self = shift;
    return glob($self->assembly_directory.'/*.fastq');
}

sub resolve_reads_metrics_files {
    my $self = shift;
    return glob($self->assembly_directory.'/*.output_metrics');
}

sub execute {
    my $self = shift;
    $self->debug_message('Allpaths metrics...');

    my $scaffolds_efasta = $self->resolve_scaffolds_efasta;

    # tier values
    my ($t1, $t2);
    if ($self->first_tier and $self->second_tier) {
        $t1 = $self->first_tier;
        $t2 = $self->second_tier;
    }
    else {
        my $est_genome_size = -s $scaffolds_efasta;
        $t1 = int ($est_genome_size * 0.2);
        $t2 = int ($est_genome_size * 0.2);
    }
    $self->debug_message('Tier one: 1 to '.$t1);
    $self->debug_message('Tier two: '.$t1.' to '.($t1 + $t2));

    # metrics
    my $metrics = Genome::Model::Tools::Sx::Metrics::Assembly->create(
        major_contig_threshold => $self->major_contig_length,
        tier_one => $t1,
        tier_two => $t2,
    );
    $self->_metrics($metrics);

    # add reads metrics/reads files
    my $add_reads_methods = $self->_add_reads_method;
    for my $file ( @{$self->_reads_files} ) {
        $self->debug_message('Add reads'.( $file =~ /metrics$/ ? ' metrics ' : ' ').'file: '.$file);
        my $add_reads = $metrics->$add_reads_methods($file);
        return if not $add_reads;
    }

    # output enhanced seq to file
    my $writer = Genome::Model::Tools::Sx::PhredWriter->create(
        file => $self->output_clean_assembly_fasta_file,
    );

    # add scaffolds
    $self->debug_message('Add scaffolds file: '.$scaffolds_efasta);
    my $reader = Genome::Model::Tools::Sx::PhredEnhancedSeqReader->create(file => $scaffolds_efasta);
    return if not $reader;
    while ( my $seq = $reader->read ) {
        $writer->write( $seq );
        my $add_scaffold = $metrics->add_scaffold($seq);
        return if not $add_scaffold;
    }
    
    # transform metrics
    $self->debug_message('Transform metrics to text...');
    my $text = $metrics->transform_xml_to('txt');
    if ( not $text ) {
        $self->error_message('Failed to transform metrics to text!');
        return;
    }
    $self->debug_message('Transform metrics to text...OK');

    # write file
    my $output_file = $self->output_file;
    unlink $output_file if -e $output_file;
    $self->debug_message('Write output file: '.$output_file);
    my $fh = eval{ Genome::Sys->open_file_for_writing($output_file); };
    if ( not $fh ) {
        $self->error_message('Failed to open metrics output file!');
        return;
    }
    $fh->print($text);
    $fh->close;

    $self->debug_message('Allpaths metrics...OK');
    return 1;
}

1;

