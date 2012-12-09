package Genome::Model::Build::AmpliconAssembly;

use strict;
use warnings;

use Genome;

use Carp 'confess';

class Genome::Model::Build::AmpliconAssembly {
    is => 'Genome::Model::Build',
    has => [
    map( { $_ => { via => 'amplicon_assembly' } } Genome::Model::Tools::AmpliconAssembly::Set->helpful_methods ),
    ],
};

sub post_allocation_initialization {
    my $self = shift;
    $self->amplicon_assembly or return 0;
}

#< AA >#
sub amplicon_assembly {
    my $self = shift;

    unless ( $self->{_amplicon_assembly} ) {
        # get
        my $amplicon_assembly = Genome::Model::Tools::AmpliconAssembly::Set->get(
            directory => $self->data_directory,
        ); 
        # create
        unless ( $amplicon_assembly ) {
            $amplicon_assembly = Genome::Model::Tools::AmpliconAssembly::Set->create(
                directory => $self->data_directory,
                description => sprintf(
                    'Model Name: %s Id: %s Build Id: %s', 
                    $self->model->name,
                    $self->model->id,
                    $self->id,
                ),
                (
                    map { $_ => $self->model->$_ } (qw/ 
                        assembly_size sequencing_center sequencing_platform subject_name 
                        /),
                ),
                exclude_contaminated_amplicons => 1,
            );
        }
        # validate
        unless ( $amplicon_assembly ) {
            $self->error_message("Can't get/create amplicon assembly.");
            return;
        }
        $self->{_amplicon_assembly} = $amplicon_assembly;
    }

    return $self->{_amplicon_assembly};
}

#< INTR DATA >#
sub link_instrument_data {
    my ($self, $instrument_data) = @_;

    unless ( $instrument_data ) {
        $self->error_message("No instument data to link");
        return;
    }

    my $chromat_dir = $self->chromat_dir;
    my $instrument_data_dir = $instrument_data->full_path;
    my $dh = Genome::Sys->open_directory($instrument_data_dir)
        or return;

    my $cnt = 0;
    while ( my $trace = $dh->read ) {
        next if $trace =~ m#^\.#;
        $cnt++;
        my $target = sprintf('%s/%s', $instrument_data_dir, $trace);
        my $link = sprintf('%s/%s', $chromat_dir, $trace);
        next if -e $link; # link points to a target that exists
        unlink $link if -l $link; # remove - link exists, but points to something that does not exist
        Genome::Sys->create_symlink($target, $link)
            or return;
    }

    unless ( $cnt ) {
        $self->error_message("No traces found in instrument data directory ($instrument_data_dir)");
    }

    return $cnt;
}

#< Reports >#
sub get_stats_dataset_from_stats_report {
    my $self = shift;

    my $report = $self->get_report('Stats');
    if ( not $report ) { 
        $report = $self->get_report('Assembly Stats'); #old name
        return if not $report;
    }

    my $dataset = $report->get_dataset('stats');
    unless ( $dataset ) {
        $self->error_message('No stats dataset found in stats report.');
        return;
    }

    return $dataset;
}

sub get_value_from_stats_report {
    my ($self, $attr) = @_;

    unless ( defined $attr ) {
        confess "No attribute given to get from stats report.";
    }
    
    my $dataset = $self->get_stats_dataset_from_stats_report;
    return 'NA' unless $dataset;
    
    my ($value) = $dataset->get_row_values_for_header($attr);
    unless ( defined $value ) {
        $self->error_message("No value for attrbute ($attr) found in stats report.");
        return;
    }
    
    return $value;
}

sub amplicons_attempted {
    return $_[0]->get_value_from_stats_report('attempted');
}

sub amplicons_assembled {
    return $_[0]->get_value_from_stats_report('assembled');
}

sub percent_assembled {
    return $_[0]->get_value_from_stats_report('assembly-success');
}

#< Files >#
sub fasta_dir {
    return $_[0]->amplicon_assembly->fasta_dir;
}

sub oriented_fasta {
    return $_[0]->amplicon_assembly->fasta_file_for_type('oriented');
}

sub oriented_qual {
    return $_[0]->amplicon_assembly->qual_file_for_type('oriented');
}

#< These map to MC16s models tmp until all aa go to mgc >#
sub amplicons_processed_success {
    return percent_assembled(@_);
}

sub oriented_fasta_file {
    return oriented_fasta(@_);
}

sub oriented_qual_file {
    return oriented_qual(@_);
}

sub processed_fasta_file {
    return $_[0]->amplicon_assembly->fasta_file_for_type('assembly');
}

sub processed_qual_file {
    return $_[0]->amplicon_assembly->qual_file_for_type('assembly');
}

sub files_ignored_by_diff {
    return qw(
        build.xml
        properties.stor
    );
}

sub dirs_ignored_by_diff {
    return qw(
        logs/
        phd_dir/
        reports/
        edit_dir/
        chromat_dir/
        abi_dir/
    );
}

1;

