package Genome::Site::TGI::Synchronize::Classes::SangerRun;

use strict;
use warnings;

use Genome;

require Compress::Zlib;
use Regexp::Common;
require List::MoreUtils;

class Genome::Site::TGI::Synchronize::Classes::SangerRun {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
    table_name => 'select run_name from gsc_run',
    id_by => {
        run_name => { is => 'Text', },
    },
    has_transient => {
        library_ids => { is => 'Array', },
    },
};

sub genome_class_for_create { 'Genome::InstrumentData::Sanger'; }

sub create_in_genome {
    my $self = shift;

    my %params = $self->params_for_create_in_genome;
    return if not %params;

    my $genome_class = $self->genome_class_for_create;
    my $genome_object = eval { $genome_class->create(%params); };
    if ( not $genome_object ) {
        Carp::confess("$@\nFailed to create $genome_class with parmas: ".Data::Dumper::Dumper(\%params));
    }

    my $dump_to_file_system = $genome_object->dump_to_file_system;
    if ( not $dump_to_file_system ) {
        $self->error_message('Failed to dump lims sanger reads to file system!');
        $genome_object->delete;
        return;
    }

    $genome_object->library_id( $self->get_library_id );

    return $genome_object;
}

sub properties_to_copy {
    return ( 'id', 'library_id', 'run_name', );
}

sub _get_read_iterator {
    my $self = shift;

    my $reads = App::DB::TableRow::Iterator->new(
        class => 'GSC::Sequence::Read',
        params => {
            prep_group_id => $self->run_name,
        },
    );

    unless ( $reads ) {
        $self->error_message( sprintf('Could not make read iterartor for run name (%s)', $self->run_name) );
        return;
    }

    return $reads;
}

sub dump_to_file_system {
    my $self = shift;
    $self->debug_message('Dump reads to file system...');

    # Allocation get/create
    my $disk_allocation = Genome::Disk::Allocation->get(owner_id => $self->id);
    unless ( $disk_allocation ) {
        $disk_allocation = Genome::Disk::Allocation->allocate(
            disk_group_name => $ENV{GENOME_DISK_GROUP_ALIGNMENTS},
            allocation_path => '/instrument_data/sanger'.$self->id,
            kilobytes_requested => 10240, # 10 Mb
            owner_class_name => $self->class,
            owner_id => $self->id
        );
        unless ($disk_allocation) {
            $self->error_message('Failed to create disk allocation for sanger instrument data '.$self->id);
            return;
        }
    }
    $self->debug_message('Allocation: '.$disk_allocation->id);

    my $data_dir = $disk_allocation->absolute_path;
    unless ( Genome::Sys->validate_existing_directory($data_dir) ) {
        $self->error_message('Absolute path from disk allocation does not exist for sanger instrument data '.$self->id);
        return;
    }
    $self->debug_message('Directory: '.$data_dir);

    my $reads = $self->_get_read_iterator;
    return unless $reads;

    $self->debug_message('Go through read iterator...');
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
        $self->error_message("No reads found for run ".$self->run_name);
        return;
    }

    $self->library_ids([ values %reads ]);

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

    $self->debug_message("Read count: ".scalar(keys %reads) );
    $self->debug_message("Dump reads to file system...OK");
    return 1;
}

sub library_id {
    my $self = shift;

    if ( not $self->library_ids ) {
        my $reads = $self->_get_read_iterator;
        return if not $reads;

        my %library_ids;
        while ( my $read = $reads->next ) {
            next if not $read->library_id;
            $library_ids{ $read->library_id }++;
        }

        $self->library_ids([ grep { $_ =~ /^$RE{num}{int}$/ } List::MoreUtils::uniq( values %library_ids ) ]);
    }

    my @library_ids = @{$self->library_ids};
    return if not @library_ids;

    my $library = Genome::Library->get(id => $library_ids[0]);
    if ( not $library ) {
        $self->error_message('Failed to get library for id: '.$library_ids[0]);
        return;
    }

    return $library->id;
}

1;

