package Genome::Model::Tools::Picard::BamToFastqSimple;

use strict;
use warnings;

use Genome;
use Genome::InlineConfig;
use File::Copy;
use File::Basename;

class Genome::Model::Tools::Picard::BamToFastqSimple {
    is  => 'Genome::Model::Tools::Picard',
    has => [ 
        input_sam_file    => { 
            is  => 'String',      
            doc => 'name of sam input file',
        },
    ],
    has_optional => [
        _java_cmd => {
            is_transient => 1,
        }
    ]
};


sub help_brief {
    "Prints out paired data from bams regardless of sorting.";
}


sub help_detail {
    return <<EOS 
Prints out paired data from bams regardless of sorting.  Prints these inline rather than in separate paired files.
EOS
}


sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    my $picard_version = $self->use_version;
    if ($picard_version < 1.31) {
        $self->error_message("bam-to-fastq-simple requires picard 1.31 or better!");
        die $self->error_message;
    }

    my $picard_dir = $self->picard_path;
    my $picard_jar_path = $picard_dir . "/sam-".$picard_version.".jar";
    my $sam_jar_path = $picard_dir . "/picard-".$picard_version.".jar";
    my $tool_jar_path = $class->base_dir . "/BamToFastqSimple.jar";

    my $cp = join ":", ($picard_jar_path, $sam_jar_path, $tool_jar_path);

    my $jvm_options = $self->additional_jvm_options || '';
    my $java_vm_cmd = 'java -Xmx'. $self->maximum_memory .'g -XX:MaxPermSize=' . $self->maximum_permgen_memory . 'm ' . $jvm_options . ' -cp '. $cp . ' edu.wustl.genome.samtools.BamToFastqSimple ' . $self->input_sam_file;
    
    $self->_java_cmd($java_vm_cmd);

    return $self;
}

sub execute {
   my $self = shift;
   Genome::Sys->shellcmd(cmd=>$self->_java_cmd);
}

sub open_stream {
    my $self = shift;

    my $fh = IO::File->new($self->_java_cmd . "|");
}
1;
