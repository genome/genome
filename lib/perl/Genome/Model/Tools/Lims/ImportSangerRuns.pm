package Genome::Model::Tools::Lims::ImportSangerRuns;

use strict;
use warnings;

use Genome;

require Compress::Zlib;
use Regexp::Common;

class Genome::Model::Tools::Lims::ImportSangerRuns {
    is => 'Command::V2',
    has => [
        run_names => {
            is => 'Text',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Run names to import.',
        },
    ],
};

sub help_brief {
    return 'Import sanger runs into genome';
}

sub execute {
    my $self = shift;
    $self->status_message('Import sanger instrument data...');

    my ($att, $succ) = (qw/ 0 0 /);
    for my $run_name ( $self->run_names ) {
        $att++;
        $succ++ if $self->_import_run($run_name);
    }

    $self->status_message('Attempted: '.$att);
    $self->status_message('Success: '.$succ);
    return $succ;
}

sub _import_run {
    my ($self, $run_name) = @_;
    $self->status_message("Run name: $run_name");

    my $sanger = Genome::InstrumentData::Sanger->get($run_name);
    my $created = 0;
    if ( not $sanger ) {
        $sanger = Genome::InstrumentData::Sanger->create(
            id => $run_name,
            run_name => $run_name,
        );
        if ( not $sanger ) {
            $self->error_message('Failed to create sanger instrument data for run name: '.$run_name);
            return;
        }
        $created = 1;
    }

    my $reads = $self->_dump_to_file_system($sanger);
    if ( not $reads ) {
        $self->error_message('Failed to dump lims sanger reads to file system!');
        $sanger->delete if $created;
        return;
    }

    my $library = $self->_set_unique_read_library($sanger, $reads);
    if ( not $library ) {
        $sanger->delete if $created;
        return;
    };

    return 1;
}

sub _dump_to_file_system {
    my ($self, $sanger) = @_;
    $self->status_message('Dump reads to file system...');

    # Allocation get/create
    my $disk_allocation = $sanger->disk_allocation;
    unless ( $disk_allocation ) {
        $disk_allocation = Genome::Disk::Allocation->allocate(
            disk_group_name => $ENV{GENOME_DISK_GROUP_ALIGNMENTS},
            allocation_path => '/instrument_data/sanger'.$sanger->id,
            kilobytes_requested => 10240, # 10 Mb
            owner_class_name => $sanger->class,
            owner_id => $sanger->id
        );
        unless ($disk_allocation) {
            $self->error_message('Failed to create disk allocation for sanger instrument data '.$sanger->id);
            return;
        }
    }
    $self->status_message('Allocation: '.$disk_allocation->id);

    my $data_dir = $disk_allocation->absolute_path;
    unless ( Genome::Sys->validate_existing_directory($data_dir) ) {
        $self->error_message('Absolute path from disk allocation does not exist for sanger instrument data '.$sanger->id);
        return;
    }
    $self->status_message('Directory: '.$data_dir);

    # Read iterator
    my $reads = App::DB::TableRow::Iterator->new(
        class => 'GSC::Sequence::Read',
        params => {
            prep_group_id => $sanger->run_name,
        },
    );
    unless ( $reads ) {
        $self->error_message( sprintf('Could not make read iterartor for run name (%s)', $sanger->run_name) );
        return;
    }

    $self->status_message('Go through read iterator...');
    my %reads; # read_name => library_id
    while ( my $read = $reads->next ) {
        $reads{$read->trace_name} = $read->library_id || 'NA';
        my $scf_name = $read->default_file_name('scf');
        my $scf_file = sprintf('%s/%s.gz', $data_dir, $scf_name);
        my $size = -s $scf_file;
        next if $size and $size > 1000; # if small retry dump
        unlink $scf_file if -e $scf_file; 
        my $scf_fh = IO::File->new($scf_file, 'w');
        unless ( $scf_fh ) {
            $self->error_message("Can't open scf ($scf_file)\n$!");
            return;
        }
        $scf_fh->print( Compress::Zlib::memGzip($read->scf_content) );
        $scf_fh->close;
        $self->error_message("No scf content for $scf_name") unless -s $scf_file;
    }

    unless ( %reads ) {
        $self->error_message("No reads found for run ".$sanger->run_name);
        return;
    }

    # Remove any other files from the directory
    my $dh = eval{ Genome::Sys->open_directory($data_dir); };
    if ( not $dh ) {
        $self->error_message('Failed to open directory! '.$data_dir);
        return;
    }
    for (1..2) { $dh->read; } # . and ..
    while ( my $file = $dh->read ) {
        my $read_name = $file;
        $read_name =~ s/\.gz$//;
        unlink $data_dir.'/'.$file if not exists $reads{$read_name};
    }

    $self->status_message("Read count: ".scalar(keys %reads) );
    $self->status_message("Dump reads to file system...OK");
    return \%reads;
}

sub _set_unique_read_library {
    my ($self, $sanger, $reads) = @_;

    my @library_ids = grep { $_ =~ /^$RE{num}{int}$/ }  values %$reads;
    $self->status_message('Found '.@library_ids.' libraries');
    $self->status_message('Using library id: '.$library_ids[0]);
    my $library = Genome::Library->get(id => $library_ids[0]);
    if ( not $library ) {
        $self->error_message('Failed to get library for id: '.$library_ids[0]);
        return;
    }
    $sanger->library($library) if not $sanger->library;

    return $library;
}

1;

