package Genome::Model::Tools::Sam::SamHeaderEditor;

use strict;
use warnings;

use Genome;
use Genome::InlineConfig;
use File::Copy;
use File::Basename;

class Genome::Model::Tools::Sam::SamHeaderEditor {
    is  => 'Genome::Model::Tools::Sam',
    has => [ 
        input_sam_file    => { 
            is  => 'String',      
            doc => 'name of sam input file',
        },
        output_sam_file    => { 
            is  => 'String',      
            doc => 'name of output same file',
        },
    ],
    has_optional => [
        bam_version_field => {
            is => 'String',
            doc => 'bam version used to generate the file',
        },
        seq_id_field => {
            is => 'String',
            doc => 'The seq_id off of the Solexa lane summary table for this alignment.  Will be set to the Read Group (RG tag) name.',
        },
        description_field => {
            is => 'String',
            doc => 'The seq_id off of the Solexa lane summary table for this alignment.  Will be set to the Read Group (RG tag) name.',
        },
        insert_size_field => {
            is => 'String',
            doc => 'The median_insert_size field off of the Solexa lane summary table for this alignment.  Will be set to the Paired Insert (PI tag) name.',
        },
        sort_order_field => {
            is => 'String',
            doc => 'sort order of the bam file',
        },
        assembly_field => {
            is => 'String',
            doc => 'the assembly used to generate the bam file e.g. HG18',
        },
        species_field => {
            is => 'String',
            doc => 'species',
        },
        refseq_path_field => {
            is => 'String',
            doc => 'the path to the reference sequence',
        },
        sample_name_field => {
            is => 'String',
            doc => 'the sample name',
        },
        platform_field => {
            is => 'String',
            doc => 'sequencing platform e.g. illumina',
        },
        platform_unit_field => {
            is => 'String',
            doc => 'a run identifier',
        },
        library_field => {
            is => 'String',
            doc => 'the library',
        },
        date_run_field => {
            is => 'String',
            doc => 'a run identifier',
        },
        genome_center_field => {
            is => 'String',
            doc => 'the genome center identifier e.g. gc-wustl',
        },
        aligner_command_field => {
            is => 'String',
            doc => 'the aligner command',
        },
        aligner_version_field => {
            is => 'String',
            doc => 'the aligner version',
        },

    ],
};


sub help_brief {
    "Edits the given sam file header to contain the provided name value pairs.";
}


sub help_detail {
    return <<EOS 
Edits the given sam file header to contain the provided name value pairs.
EOS
}


sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    #find the location of the required jar files based on path of the class            
    my $base_dir = $class->base_dir; 
    my $cp = "$ENV{GENOME_SW_LEGACY_JAVA}/samtools/sam-current.jar".":".$base_dir."/"."GCSam.jar";
            
    #eval the 'use inline' command to inject the correct classpath
    #note:  the "PACKAGE=>'main'" entry allows you to replace "Genome::Model::Tools::Sam::SamHeaderEditor" with "main" when subsequently calling the Java class.  See in execute. 
    $ENV{PERL_INLINE_JAVA_JNI} = 1;
    my $inline_string =  "use Inline ( Java => 'STUDY', CLASSPATH => qw($cp), STUDY => ['edu.wustl.genome.samtools.SamHeaderEditor'], EXTRA_JAVA_ARGS => '-Xmx2g', AUTOSTUDY => 1, PACKAGE => 'main', DIRECTORY => Genome::InlineConfig::DIRECTORY(), JNI=>1 );";
    eval "$inline_string";
    unless ( defined $@ ) {
        $self->error_message("There was a problem evaluating the 'use Inline' code for the required Java libraries: $inline_string.  The classpath was set to: $cp");
        die;
    } 
        
    #test to see if the input sam file exists. 
    $self->error_message('Sam input file does not exist.') and return unless -s $self->input_sam_file;

    return $self;
}

sub execute {

   my $self = shift;
 
   #see note in create() for explanation of 'main' package definition.
   my $she = new main::edu::wustl::genome::samtools::SamHeaderEditor($self->input_sam_file, $self->output_sam_file);
   $self->debug_message("Created SamHeaderEditor");
   my $rg_field_result = $she->setReadGroupName($self->seq_id_field);  #this should also be the id of the read group
   $self->debug_message("Called setReadGroupName.  Return value: $rg_field_result");
 
   my $rg_result = $she->addReadGroup( 
                         $self->seq_id_field, #the id of the read group is the first field 
                         $self->description_field, #the description of the read group, e.g. paired end or fragment 
                         $self->insert_size_field, #setting the PI tag. 
                         $self->sample_name_field, 
                         $self->platform_field, 
                         $self->platform_unit_field, 
                         $self->library_field, 
                         $self->date_run_field,
                         $self->genome_center_field, 
                    );
   $self->debug_message("Called addReadGroup.  Return value: $rg_result");
   
   #The order of the following parameters is important.  Do not change unless you change the corresponding Java class!
   my $aa_result = $she->addAttributes(  
                         $self->bam_version_field, 
                         $self->sort_order_field, 
                         $self->assembly_field, 
                         $self->species_field,
                         $self->refseq_path_field, 
                         $self->aligner_command_field,
                         $self->aligner_version_field
                      );
   
   $self->debug_message("Called addAttributes.  Return value: $aa_result");
   
   $she->execute();
   $self->debug_message("SamHeaderEditor execute() complete.");
   
   my $result = 1; 
   return $result;
}

1;
