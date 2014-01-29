package Genome::Model::Tools::Soap::FastaToAgp;

use strict;
use warnings;

use Genome;
use File::Basename;

class  Genome::Model::Tools::Soap::FastaToAgp {
    is => 'Genome::Model::Tools::Soap::Base',
    has => [
        scaffold_size_cutoff => {
            is => 'Integer',
            doc => 'Minimum scaffold size cutoff',
        },
	version => {
	    is => 'String',
	    doc => 'Version of fasta2agp script to run',
	    valid_values => ['9.27.10'],
	},
	assembly_directory => {
	    is => 'Text',
	    doc => 'Assembly directory',
	},
    ],
    has_optional => [
        output_dir => {
            is => 'Text',
            doc => 'Directory to put output files',
        },
	scaffold_sequence_file => {
            is => 'Text',
            doc => 'Soap generated scaffold fasta file, if not specified, tool will derive it',
	},
	file_prefix => {
            is => 'Text',
            doc => 'Output file prefix name, if not specified, tool will derive it from soap output file prefixes',
        },
    ],
    has_optional_transient => [
#	_file_prefix         => { is => 'Text', },
    ],
};

sub help_brief {
    'Tool to run fasta2agp script for soap PGA assemblies';
}

sub help_detail {
    return <<"EOS"
gmt soap fasta-to-agp --version 9.27.10 --scaffold_sequence_file /gscmnt/111/soap_assembly/SRA111_WUGC.scafSeq --scaffold-size-cutoff 100 --file-prefix SRA111_WUGC --output-dir  /gscmnt/111/soap_assembly/PGA
EOS
}

sub execute {
    my $self = shift;
    
    #get version of script
    my $script = $self->_full_path_to_version_script;

    #input files
    my $scaf_seq_file = ( $self->scaffold_sequence_file ) ? $self->scaffold_sequence_file : $self->assembly_scaffold_sequence_file;

    #output directory
    my $output_dir = ($self->output_dir) ? $self->output_dir : $self->assembly_directory;
    Genome::Sys->create_directory( $output_dir ) unless -d $output_dir;

    #output file prefix
    my $file_prefix = ( $self->file_prefix ) ? $self->file_prefix : $self->assembly_file_prefix;

    #create script command string
    my $command = 'perl '.$script.' -i '.$scaf_seq_file.' -o '.$output_dir;
    $command .= ' -size '.$self->scaffold_size_cutoff if $self->scaffold_size_cutoff;
    $command .= ' -name '.$file_prefix;

    $self->debug_message("Running fasta2agp with command: $command");

    system("$command"); #script has no return value

    #check for expected output files
    for my $file_ext ( qw/ contigs.fa agp scaffolds.fa / ) {
	my $file = $output_dir.'/'.$file_prefix.'.'.$file_ext;
	unless ( -e $file ) {
	    $self->error_message("Failed to find output file: $file");
	    return;
	}
    }

    return 1;
}

sub _full_path_to_version_script {
    my $self = shift;

    my $module_path = $self->class;
    $module_path =~ s/::/\//g;

    my $genome_dir = Genome->get_base_directory_name();
    my $inc_dir = substr($genome_dir, 0, -7);
    my $script = $inc_dir.'/'.$module_path.'/'.$self->version.'/'.$self->_script_name;

    unless ( -x $script ) {
        $self->error_message("Failed to find script: $script");
        return;
    }

    return $script;
}

sub _script_name {
    return 'fasta2agp.pl';
}

1;
