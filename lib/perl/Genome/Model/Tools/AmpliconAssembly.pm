package Genome::Model::Tools::AmpliconAssembly;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::AmpliconAssembly {
    is => 'Command',
    is_abstract => 1,
    has => [
        directory => {
            is => 'Text',
            doc => 'Base directory for the amplicon assembly.  It is required that the ' .
                'amplicon assembly have been previously created, saving it\'s properties. ' .
                'See the "create" command.',
        },
        amplicon_assembly => {
            is => 'Genome::Model::Tools::AmpliconAssembly::Set',
            calculate_from => 'directory',
            calculate => q{ return Genome::Model::Tools::AmpliconAssembly::Set->get(directory => $directory); },
        },
        chromat_dir => { via => 'amplicon_assembly' },
        phd_dir => { via => 'amplicon_assembly' },
        edit_dir => { via => 'amplicon_assembly' },
        consed_directory => { via => 'amplicon_assembly' },
        create_directory_structure => { via => 'amplicon_assembly' },
        get_amplicons => { via => 'amplicon_assembly' },
        amplicon_fasta_types => { via => 'amplicon_assembly' },
        amplicon_bioseq_method_for_type => { via => 'amplicon_assembly' },
        fasta_file_for_type => { via => 'amplicon_assembly' },
        qual_file_for_type => { via => 'amplicon_assembly' },
        assembly_fasta => { via => 'amplicon_assembly' },
        reads_fasta => { via => 'amplicon_assembly' },
        processed_fasta => { via => 'amplicon_assembly' },
    ],
};

#< Helps >#
sub help_brief {
    return ucfirst(join(' ', split('-', $_[0]->command_name_brief))).' amplicon assemblies';
}

sub help_synopsis {
};

#< UR >#
sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    unless ( $self->amplicon_assembly ) {
        $self->error_message("Can't get amplicon assembly with given parameters. There maybe an error, or it might just needs to be created first.");
        return;
    }

    return $self;
}

1;

#$HeadURL$
#$Id$
