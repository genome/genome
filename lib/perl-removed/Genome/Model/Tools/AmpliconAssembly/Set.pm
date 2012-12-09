package Genome::Model::Tools::AmpliconAssembly::Set;

use strict;
use warnings;

use Genome;

use Carp 'confess';
use Cwd 'abs_path';
use Data::Dumper 'Dumper';
require File::Copy;
require Genome::Model::Tools::Consed::Directory;
require Genome::Model::Tools::AmpliconAssembly::Amplicon;

use Bio::SeqIO;

my %ATTRIBUTES = (
    directory => {
        is => 'Text',
        is_optional => 1,
        doc => 'Base directory for the amplicon assembly.',
    },
    description => {
        is => 'Text',
        is_optional => 1,
        doc => 'A brief description of the amplicon assembly. Helps identify in messaging. Default is the directory.',
    },
    sequencing_center => {
        is => 'Text',
        is_optional => 1,
        default_value => __PACKAGE__->default_sequencing_center,
        doc => 'Sequencing Center that the amplicons were sequenced. Currently supported centers: '.join(', ', __PACKAGE__->valid_sequencing_centers).' (default: '.__PACKAGE__->default_sequencing_center.').',
    },
    sequencing_platform => {
        is => 'Text',
        is_optional => 1,
        default_value => __PACKAGE__->default_sequencing_platform,
        doc => 'Platform upon whence the amplicons were sequenced. Currently supported platforms: '.join(', ', __PACKAGE__->valid_sequencing_platforms).' (default: '.__PACKAGE__->default_sequencing_platform.').',
    },
    assembly_size => { # amplicon_size
        is => 'Integer',
        is_optional => 1,
        default_value => 1364,
        doc => 'Expected assembly size for an assembled amplicon.',
    },
    subject_name => {
        is => 'Text',
        is_optional => 1,
        doc => 'The subject name of the underlying data. Used primarily in naming files (etc) that represent the entire amplicon assembly, like an assembled fasta of each amplicon.',
    },
    exclude_contaminated_amplicons => {
        is => 'Boolean',
        is_optional => 1,
        default_value => 0,
        doc => 'DEPRECATED! DOES NOT DO ANYTHING!',
    },
    only_use_latest_iteration_of_reads => {
        is => 'Boolean',
        is_optional => 1,
        default_value => 0,
        doc => 'DEPRECATED! DOES NOT DO ANYTHING!',
    },
);

#< UR >#
class Genome::Model::Tools::AmpliconAssembly::Set {
    is => 'UR::Object',
    id_by => 'directory',
    has => [
    %ATTRIBUTES,
    _additional_properties => {
        is => 'HASH',
        default_value => {},
        doc => 'Properties added from operations performed on an amplicon assembly.',
    },
    ],
};

sub get { 
    my ($class, %params) = @_;

    # Do not support getting with any other params...
    my $directory = delete $params{directory};
    if ( %params ) { 
        $class->error_message('Can\'t get with additional parameters, only directory.');
        return;
    }
    
    # Validate directory
    unless ( Genome::Sys->validate_existing_directory($directory) ) {
        $class->error_message("Can't validate amplicon assembly directory. See above error.");
        return;
    }

    # UR Cache
    my ($self) = $class->SUPER::get(directory => $directory);
    return $self if $self;
    
    # Properties file
    return $class->_create_from_saved_properties($directory);
}

sub create {
    my ($class, %params) = @_;

    $params{directory} = Cwd::abs_path($params{directory});
    unless ( Genome::Sys->validate_existing_directory($params{directory}) ) {
        $class->error_message("Can't validate amplicon assembly directory. See above error.");
        return;
    }

    # Get UR cache/properties file
    #my ($self) = $class->get(directory => $params{directory});
    #return $self if $self;

    my $self = $class->SUPER::create(%params)
        or return;
    
    unless ( $self->description ) {
        $self->description( $self->directory );
    }

    $self->_validate_properties
        or return;

    $self->create_directory_structure
        or return;
    
    $self->_save_properties
        or return;

    return $self;
}

#< Properties >#
sub _properties_file {
    my ($self, $directory) = @_;

    if ( not $directory and ref $self ) {
        $directory = $self->directory;
    }

    confess "No directory given to get properties file." unless $directory;

    return $directory.'/properties.stor';
}

sub _create_from_saved_properties {
    my ($class, $directory) = @_;

    my $properties_file = $class->_properties_file($directory);
    return unless -e $properties_file;
    
    my $properties = Storable::retrieve($properties_file);
    unless ( $properties ) {
        confess "No properties in amplicon assembly properties file.";
    }

    return $class->SUPER::create(%$properties)
}

sub _validate_properties {
    my $self = shift;

    for my $attribute (qw/ sequencing_platform sequencing_center /) {
        unless ( $self->validate_attribute_value($attribute) ) {
            $self->delete;
            return;
        }
    }

    return 1;
}

sub _save_properties {
    my $self = shift;

    my %properties;
    for my $attr ( $self->attribute_names ) {
        $properties{$attr} = $self->$attr;
    }

    $properties{_additional_properties} = $self->_additional_properties;

    my $properties_file = $self->_properties_file;
    unlink $properties_file if -e $properties_file;
    eval {
        Storable::nstore(\%properties, $properties_file);
    };
    
    if ( $@ ) {
        confess "Can't save properties to file ($properties_file): $@";
    }
    
    return 1;
}

sub update_additional_properties {
    my ($self, %properties) = @_;

    confess "No properties given to update additional properties." unless %properties;

    for my $property ( keys %properties ) {
        $self->_additional_properties->{$property} = $properties{$property};
    }
    
    return $self->_save_properties;
}

#< Atributes >#
sub attributes {
    return %ATTRIBUTES;
}

sub attribute_names {
    return keys %ATTRIBUTES;
}

sub attributes_without_default_values {
    my $class = shift;

    my %attributes;
    for my $name ( attribute_names() ) {
        $attributes{$name} = $class->get_attribute_for_name($name);
        delete $attributes{$name}->{default_value};
    }

    return %attributes;
}

sub get_attribute_for_name {
    my ($class, $name) = @_;

    confess "No attribute name given to get attribute" unless $name;
    confess "No attribute found for name ($name)" unless exists $ATTRIBUTES{$name};

    return $ATTRIBUTES{$name};
}

sub helpful_methods {
    return (qw/ 
        chromat_dir phd_dir edit_dir
        consed_directory create_directory_structure
        get_amplicons 
        amplicon_fasta_types amplicon_bioseq_method_for_type
        fasta_file_for_type qual_file_for_type
        assembly_fasta reads_fasta processed_fasta 
        /);
}

sub validate_attribute_value {
    my ($self, $attribute) = @_;

    my $value = $self->$attribute;
    unless ( defined $value ) { # all attrs have defaults, make sure it is set
        $self->error_message("No value for $attribute given");
        return;
    }
    
    my $valid_values_method = 'valid_'.$attribute.'s';
    unless ( grep { $value eq $_ } $self->$valid_values_method ) {
        $self->error_message(
            sprintf(
                'Invalid %s %s. Valid %ss: %s',
                $attribute,
                $value,
                $attribute,
                join(', ',valid_sequencing_centers()),
            )
        );
        return;
    }

    return 1;
}

# Sequencing Centers
sub valid_sequencing_centers {
    return (qw/ gsc broad /);
}

sub default_sequencing_center {
    return (valid_sequencing_centers)[0];
}

# Sequencing Platforms
sub valid_sequencing_platforms {
    return (qw/ sanger /);
}

sub default_sequencing_platform {
    return (valid_sequencing_platforms)[0];
}

#< DEPRACATED FIXME #>
sub assembler { return 'phredphrap' };

#< DIRS >#
sub consed_directory {
    my $self = shift;

    unless ( $self->{_consed_directory} ) {
        $self->{_consed_directory} = Genome::Model::Tools::Consed::Directory->create(directory => $self->directory);
    }

    return $self->{_consed_directory};
}

sub create_directory_structure {
    my $self = shift;

    $self->consed_directory->create_extended_directory_structure
        or return;

    return 1;
}

sub edit_dir {
    return $_[0]->consed_directory->edit_dir;
}
    
sub phd_dir {
    return $_[0]->consed_directory->phd_dir;
}
    
sub chromat_dir {
    return $_[0]->consed_directory->chromat_dir;
}

sub fasta_dir {
    return $_[0]->consed_directory->fasta_dir;
}

#< FASTA >#
my %_fastas_and_amplicon_bioseq_methods = (
    reads => 'get_bioseqs_for_raw_reads',
    processed => 'get_bioseqs_for_processed_reads',
    assembly => 'get_assembly_bioseq',
    oriented => 'get_oriented_bioseq',
);

sub amplicon_fasta_types {
    return keys %_fastas_and_amplicon_bioseq_methods;
}

sub amplicon_bioseq_method_for_type {
    return $_fastas_and_amplicon_bioseq_methods{$_[1]};
}

sub fasta_file_for_type {
    my ($self, $type) = @_;

    return sprintf(
        '%s/%s%s.fasta',
        $self->fasta_dir,
        ( defined $self->subject_name ? $self->subject_name.'.' : '' ),
        $type,
    );
}

sub qual_file_for_type {
    return $_[0]->fasta_file_for_type($_[1]).'.qual';
}

#< Amplicons >#
sub get_amplicons {
    my $self = shift;

    my $method = sprintf(
        '_get_amplicons_and_read_names_for_%s_%s', 
        $self->sequencing_center,
        $self->sequencing_platform,
    );

    my $iterator = $self->$method;
    
    my $amplicons = $self->$method;
    unless ( $amplicons and %$amplicons ) {
        $self->error_message(
            sprintf('No amplicons found in chromat_dir of directory (%s)', $self->directory) 
        );
        return;
    }

    my @amplicons;
    my $edit_dir = $self->edit_dir;
    for my $name ( sort { $a cmp $b } keys %$amplicons ) {
        push @amplicons, Genome::Model::Tools::AmpliconAssembly::Amplicon->create(
            name => $name,
            reads => $amplicons->{$name},
            directory => $edit_dir,
        );
    }

    return \@amplicons;
}

sub _get_amplicons_and_read_names_for_gsc_sanger {
    my $self = shift;

    my $read_iterator = $self->_get_read_name_iterator
        or return;

    my %amplicons;
    while ( my $read_name = $read_iterator->() ) {
        my $amplicon_name = $self->_get_amplicon_name_for_gsc_sanger_read_name($read_name)
            or next; # Assume conversion logic is good, and returns undef cuz it is not a read
        push @{$amplicons{$amplicon_name}}, $read_name;
    }

    return \%amplicons;
}

sub _get_amplicons_and_read_names_for_broad_sanger {
    my $self = shift;

    my $read_iterator = $self->_get_read_name_iterator
        or return;

    my %amplicons;
    while ( my $read_name = $read_iterator->() ) {
        my $amplicon = $self->_get_amplicon_name_for_broad_sanger_read_name($read_name)
            or return; # Assume conversion logic is good, and returns undef cuz it is not a read
        push @{$amplicons{$amplicon}}, $read_name;
    }

    return  \%amplicons;
}

sub _get_read_name_iterator {
    my $self = shift;

    my $dh = Genome::Sys->open_directory( $self->chromat_dir );
    unless ( $dh ) {
        $self->error_message("Can't open chromat dir to get reads. See above error.");
        return;
    }

    return sub{
        while ( my $name = $dh->read ) {
            next if $name =~ m#^\.#;
            $name =~ s#\.gz##;
            return $name;
        }
        $dh->close;
        return;
    }
}

#< Amplicon Reads >#
sub get_method_for_get_amplicon_name_for_read_name {
    my $self = shift;

    return sprintf(
        '_get_amplicon_name_for_%s_%s_read_name', 
        $self->sequencing_center,
        $self->sequencing_platform,
    );
}

sub _get_amplicon_name_for_gsc_sanger_read_name {
    my ($self, $read_name) = @_;

    $read_name =~ /^(.+)\.[bg]\d+$/
        or return;

    return $1;
}

sub _get_amplicon_name_for_broad_sanger_read_name {
    my ($self, $read_name) = @_;

    $read_name =~ s#\.T\d+$##;
    $read_name =~ s#[FR](\w\d\d?)$#\_$1#; # or next;

    return $read_name;
}

sub get_all_amplicons_reads_for_read_name {
    my ($self, $read_name) = @_;

    unless ( $read_name ) {
        $self->error_message("No read name given to get all reads for amplicon.");
        return;
    }

    my $amp_method = $self->get_method_for_get_amplicon_name_for_read_name;
    my $amplicon_name = $self->$amp_method($read_name);
    unless ( $amplicon_name ) {
        $self->error_message("Can't get amplicon name for read name ($read_name)");
        return;
    }

    my $reads_method = $self->get_method_for_get_all_amplicons_reads_for_read_name;
    my @read_names = $self->$reads_method($amplicon_name);
    unless ( @read_names ) {
        $self->error_message("No reads found for amplicon name ($amplicon_name)");
        return;
    }


    return @read_names;
}

sub get_method_for_get_all_amplicons_reads_for_read_name {
    sprintf(
        '_get_all_reads_for_%s_%s_amplicon',
        $_[0]->sequencing_center,
        $_[0]->sequencing_platform,
    );
}

sub _get_all_reads_for_gsc_sanger_amplicon {
    my ($self, $amplicon_name) = @_;
    
    my $chromat_dir = $self->chromat_dir;
    my @read_names;
    for my $read_name ( glob("$chromat_dir/$amplicon_name.*") ) {
        $read_name =~ s#$chromat_dir/##;
        $read_name =~ s#\.gz##;
        push @read_names, $read_name;
    }

    return @read_names;
}

sub _get_all_reads_for_broad_sanger_amplicon {
    die "Not implemented\n";
}

#< Contamination Screening >#
sub contamination_dir {
    return $_[0]->directory.'/contamination';
}

sub contamination_reads_dir {
    return $_[0]->contamination_dir.'/reads';
}
    
sub amplicon_fasta_file_for_contamination_screening {
    return $_[0]->contamination_dir.'/amplicon_reads.fasta';
}

sub create_contamination_dir_structure {
    my $self = shift;

    my $contamination_dir = $self->contamination_dir;
    unless ( -d $contamination_dir ) {
        return unless Genome::Sys->create_directory($contamination_dir);
    }

    my $reads_dir = $self->contamination_reads_dir;
    unless ( -d $reads_dir ) {
        return unless Genome::Sys->create_directory($reads_dir);
    }

    return 1;
}

sub create_amplicon_fasta_file_for_contamination_screening {
    my $self = shift;

    my $amplicons = $self->get_amplicons
        or return;

    $self->create_contamination_dir_structure
        or return;
    
    my $fasta_file = $self->amplicon_fasta_file_for_contamination_screening;
    unlink $fasta_file if -e $fasta_file;
    my $fasta_writer = Bio::SeqIO->new(
        '-file' => '>'.$fasta_file,
        '-fomat' => 'fasta',
    );
    for my $amplicon ( @$amplicons ) {
        for my $bioseq ( $amplicon->get_bioseqs_for_processed_reads ) {
            next unless $bioseq->length >= 11;
            $bioseq->seq( uc $bioseq->seq );
            $fasta_writer->write_seq($bioseq);
        }
    }

    return $fasta_file;
}

sub remove_contaminated_amplicons_by_reads_in_file {
    my ($self, $file) = @_;

    $self->create_contamination_dir_structure
        or return;

    my $fh = Genome::Sys->open_file_for_reading($file)
        or return;

    my %amplicons_seen;
    my $chromat_dir = $self->chromat_dir;
    my $contamination_reads_dir = $self->contamination_reads_dir;
    my $amp_method = $self->get_method_for_get_amplicon_name_for_read_name;
    my $reads_method = $self->get_method_for_get_all_amplicons_reads_for_read_name;

    while ( my $read_name = $fh->getline ) {
        chomp $read_name;
        my $amplicon_name = $self->$amp_method($read_name);
        unless ( $amplicon_name ) {
            $self->error_message("Can't get amplicon name for read name ($read_name)");
            return;
        }
        my @read_names = $self->$reads_method($amplicon_name);
        unless ( @read_names ) {
            $self->error_message("Can't get reads for amplicon name ($amplicon_name)");
            return;
        }
        # Move chromats - this effectively removes the amplicon
        # TODO create the amplicon interface here, to support sanger and 454
        #  maybe move this logic to amplicon?
        for my $read_name ( @read_names ) {
            my $from = "$chromat_dir/$read_name.gz"; 
            my $to = "$contamination_reads_dir/$read_name.gz";
            unless ( -e $from ) { # ok, i guess
                $self->error_message("Can't find trace $read_name");
                next;
            }
            unless ( File::Copy::move($from, $to) ) { # not ok
                $self->error_message("Can't move $from to $to\: $!");
                return;
            }
        }
        # TODO move phd_dir edit_dir files?
    }

    return 1;
}

1;

=pod

=head1 Name

Genome::Model::Tools::AmpliconAssembly::Set

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 valid_sequencing_centers

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 default_sequencing_center

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 consed_directory

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 create_directory_structure

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 edit_dir

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 phd_dir

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 chromat_dir

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 fasta_dir

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 amplicon_fasta_types

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 amplicon_bioseq_method_for_type

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 fasta_file_for_type

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 qual_file_for_type

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 get_amplicons

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 contamination_dir

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 contamination_reads_dir

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 amplicon_fasta_file_for_contamination_screening

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 create_contamination_dir

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 create_amplicon_fasta_files_for_contamination_screening

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head2 read_is_contaminated

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2009 Genome Center at Washington University in St. Louis

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@genome.wustl.edu>

=cut

#$HeadURL$
#$Id$
