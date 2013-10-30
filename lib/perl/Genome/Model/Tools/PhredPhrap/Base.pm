package Genome::Model::Tools::PhredPhrap::Base;

use strict;
use warnings;

use Genome;

require Cwd;
use Data::Dumper;
require Genome::Model::Tools::Consed::Directory;
require IO::File;

#- PROPERTIES -#
my %properties = (
    assembly_name  => {
        type => 'String',
        is_optional => 0,
        doc =>'Name of assembly (all files created will have this as a base)',
    },
    directory => {
        type => 'String', 
        is_optional => 0,
        doc =>'Base directory to work in, with (or will create) chromat_dir, phd_dir and edit_dir',
    },
    user  => {
        type => 'String',
        is_optional => 1,
        default => qw(Genome::Sys->username),
        doc =>'User (default is current user)',
    }, 
);

#- PROCESSOR CLASSES - ADD PROPS TO THIS CLASS -#
for my $processor ( pre_assembly_processors(), post_assembly_processors() ) {
    my $class = class_for_pre_assembly_processor($processor);
    eval{ $class->class; };

    $properties{"no_$processor"} = {
        type => 'Boolean',
        is_optional => 1,
        doc => "Do not run $processor",
    };

    my $acc_class_meta = $class->__meta__;
    for my $property ( $acc_class_meta->direct_property_metas ) {
        next if $property->property_name eq 'fasta_file'
            or exists $properties{ $property->property_name };
        $properties{ $property->property_name } = {
            type => $property->property_name,
            is_optional => 1,
            doc => "(for $processor) " . $property->doc,
        };
    }
}

class Genome::Model::Tools::PhredPhrap::Base {
    is => 'Genome::Model::Tools::PhredPhrap::BaseBase',
    is_abstract => 1,
    has => [ %properties ],
};

sub _directory {
    return $_[0]->{_directory};
}

sub _cwd {
    return $_[0]->{_cwd};
}

sub create { 
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    $self->{_cwd} = Cwd::getcwd();

    my $directory = Genome::Model::Tools::Consed::Directory->create(directory => $self->directory)
        or return;
    $directory->create_consed_directory_structure;
    $self->{_directory} = $directory;

    my $create_busy_file_ok = $self->_create_busy_file;
    return if not $create_busy_file_ok;
    
    for my $file_method ( (qw/ acefile singlets_file /), $self->_files_to_remove ) {
       my $file_name = $self->$file_method;
       unlink $file_name if -e $file_name;
    }
    
    return $self;
}

sub DESTROY {
    my $self = shift;

    $self->SUPER::DESTROY;

    unlink $self->busy_file if -e $self->busy_file;

    return chdir $self->_cwd;
}

#- BUSY FILE -#
sub busy_file {
    my $self = shift;
    
    return sprintf('%s/%s.PHREDPHRAP.BUSY', $self->directory, $self->assembly_name);
}

sub _create_busy_file {
    my $self = shift;

    my $busy_file = $self->busy_file;
    $self->error_message('Phred phrap busy file: '.$busy_file);
    if ( -e $self->busy_file ) {
        $self->error_message('Phred phrap busy file exists! There may be another phred phrap process running for this assembly! '.$self->assembly_name);
        $self->error_message("If a phred phrap process is not running, remove the busy file ($busy_file) and retry.");
        return;
    }
    return IO::File->new('>' . $self->busy_file)->close;
}

#- EXECUTE -#
sub execute {
    my $self = shift;

    $self->_handle_input;

    # Save fasta and qual before processing
    my $fasta_file = $self->_get_value_or_default_for_property_name('fasta_file');
    my $bak_fasta_file = sprintf('%s/%s.reads.fasta', $self->_directory->edit_dir, $self->assembly_name);
    unless ( File::Copy::copy($fasta_file, $bak_fasta_file) ) {
        $self->error_message("Can't copy $fasta_file to $bak_fasta_file: $!");
        return;
    }

    my $qual_file = $self->_get_value_or_default_for_property_name('qual_file');
    my $bak_qual_file = $bak_fasta_file.'.qual';
    unless ( File::Copy::copy($qual_file, $bak_qual_file) ) {
        $self->error_message("Can't copy $qual_file to $bak_qual_file: $!");
        return;
    }

    # PRE ASSEMBLY PROCESSING
    $self->status_message("Pre Assembly Processing");
    for my $processor_name ( $self->pre_assembly_processors ) {
        my $no_processor = 'no_' . $processor_name ;
        next if $self->$no_processor;
        my $class = class_for_pre_assembly_processor($processor_name);
        my %params = $self->_params_for_class($class);
        my $processor = $class->create(%params);
        if ( $processor->__errors__ ) { 
            print "\n\n",$processor->error_message,"\n\n";
            #or $self->error_message("Can't create class ($class)");
            return;
        }
        $processor->status_message('Running');
        unless ( $processor->execute ) {
            $self->error_message("Pre Assembly Processing failed, cannot assemble");
            return;
        }
    }

    my %params = $self->_params_for_class('Genome::Model::Tools::PhredPhrap::Fasta');
    my $phrap = Genome::Model::Tools::PhredPhrap::Fasta->create(%params);
    $phrap->status_message("Assembling");
    $phrap->execute
        or return;

    #POST ASEMBLY PROCESS
    $self->status_message("Post Assembly Processing");

    unlink $self->busy_file if -e $self->busy_file;
    chdir $self->_cwd;

    $self->status_message("Finished");

    return 1;
}

my @options_with_defaults = (qw/ scf_file exclude_file phd_file fasta_file qual_file /);
sub _params_for_class {
    my ($self, $class) = @_;

    my %params_for_class;
    for my $property ( grep {
            $_->class_name ne 'UR::Object' and $_->class_name ne 'Command' 
        }$class->__meta__->all_property_metas ) {
        my $property_name = $property->property_name;
        my $value;# = $self->$property_name;
        unless ( $self->can($property_name) and defined($value = $self->$property_name) ) {
            next unless grep { $property_name eq $_ } @options_with_defaults;
            my $default_property_name = sprintf('default_%s', $property_name);
            $value = $self->$default_property_name;
            next unless defined $value;
        }
        $params_for_class{$property_name} = $value;
    }

    return %params_for_class;
}

sub _get_value_or_default_for_property_name {
    my ($self, $property_name) = @_;

    my $value;
    unless ( $self->can($property_name) and defined($value = $self->$property_name) ) {
        return unless grep { $property_name eq $_ } @options_with_defaults;
        my $default_property_name = sprintf('default_%s', $property_name);
        $value = $self->$default_property_name;
    }

    return $value;
}

#- INPUT PROCESSORS -#
sub assembly_processors_and_classes {
    return pre_assembly_processors_and_classes(), post_assembly_processors_and_classes();
}

sub pre_assembly_processors_and_classes { 
    return (
        A_screen_vector => 'Genome::Model::Tools::Fasta::ScreenVector',
        B_trim_quality => 'Genome::Model::Tools::Fasta::TrimQuality',
    );
}

sub pre_assembly_processors {
    my %itnc = pre_assembly_processors_and_classes();
    return grep { s/^\w\_// } sort keys %itnc;
}

sub pre_assembly_processor_classes {
    my %itnc = pre_assembly_processors_and_classes();
    return values %itnc;
}

sub class_for_pre_assembly_processor {
    my ($processor) = @_;
    my %itnc = pre_assembly_processors_and_classes();
    my ($found) = grep { $_ =~ /$processor/ } keys %itnc;
    return $itnc{ $found };
}

#- POST ASSEMBLY PROCESSORS -#
sub post_assembly_processors_and_classes { 
    return (
    );
}

sub post_assembly_processors {
    my %itnc = post_assembly_processors_and_classes();
    return sort keys %itnc;
}

sub post_assembly_processor_classes {
    my %itnc = post_assembly_processors_and_classes();
    return values %itnc;
}

sub class_for_post_assembly_processor {
    my ($processor) = @_;
    my %itnc = post_assembly_processors_and_classes();
    return $itnc{ $processor };
}

#- DIRS, FILES, DEFAULTS, ETC -#
sub default_scf_file {
    my $self = shift;

    return sprintf('%s/%s.include', $self->_directory->edit_dir, $self->assembly_name);
}

sub default_exclude_file {
    my $self = shift;

    return sprintf('%s/%s.exclude', $self->_directory->edit_dir, $self->assembly_name);
}

sub default_phd_file {
    my $self = shift;

    return sprintf('%s/%s.phds', $self->_directory->edit_dir, $self->assembly_name);
}

sub default_fasta_file {
    my $self = shift;

    return sprintf('%s/%s.fasta', $self->_directory->edit_dir, $self->assembly_name);
}

sub default_qual_file {
    my $self = shift;

    return sprintf('%s.qual', $self->default_fasta_file);
}

sub acefile {
    my $self = shift;

    return sprintf('%s.ace', $self->default_fasta_file);
}

sub singlets_file {
    my $self = shift;

    return sprintf('%s.singlets', $self->default_fasta_file);
}

1;

=pod

=head1 Name

ModuleTemplate

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
