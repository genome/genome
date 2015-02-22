package Genome::Model::Tools::Joinx::VcfAnnotate;

use strict;
use warnings;

use Genome;
use Data::Dumper;

our $MINIMUM_JOINX_VERSION = 1.4; #available but horribly broken in 1.3

class Genome::Model::Tools::Joinx::VcfAnnotate {
    is => 'Genome::Model::Tools::Joinx',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'Vcf File to filter',
            shell_args_position => 1,
        },
        annotation_file => {
            is => 'Text',
            doc => 'Vcf File containing annotation',
            shell_args_position => 2,
        },
        info_fields => {
            is => 'Text',
            doc => 'Field ids to embed from the annotation VCF. Use colons to separate multiple field descriptors.',
            #doing the above because UR autosplits on commas with is_many, but joinx uses commas in its field descriptors
        },
    ],
    has_optional_input => [
        output_file => {
            is => 'Text',
            is_output => 1,
            doc => 'The output file (defaults to stdout)',
        },
        use_bgzip => {
            is => 'Boolean',
            doc => 'zcats the input files into stdin, and bgzips the output',
            default => 0,
        },
        identifiers => {
            is => 'Boolean',
            default => 1,
            doc => 'copy identifiers from the annotation file',
        },
        info => {
            is => 'Boolean',
            default => 1,
            doc => 'copy information from info fields from the annotation file',
        },
    ],
};

sub help_brief {
    "Annotate information from one VCF file's identifiers and info fields into another."
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  This tool wraps the joinx vcf-annotate command which places information from one VCF into another. For dbSNP you would commonly want to call it like so: 

  gmt joinx vcf-annotate --input-file a.vcf --annotation-file dbsnp.vcf --info-fields GMAF:dbSNPBuildID=dbSNPBuildID,per-alt --output-file annotated.vcf

  By default the id field of the annotating VCF is used to populate the id field of the input VCF.
  The --info-fields option allows you to specify how items are embedded from the annotation VCF to the input VCF.
    * The name of the tag along (e.g. GMAF) moves it wholesale into the annotating file.
    * You can rename a tag using the syntax OLD=NEW. e.g. GMAF=dbSNPMAF would embed the dbSNP file's GMAF as dbSNPMAF
    * You can also request that the field be populated on a per-alternate basis. This allows you to distinguish which alleles have been annotated. e.g. dbSNPBuildID=dbSNPBuildID,per-alt

  An example output of the above command, illustrating these points is below. Irrelevant aspects of the output VCF are left out for clarity:
    1	196341842	rs185457472	G	C,T	.	PASS	dbSNPBuildID=135,.
    1	196351766	rs2147469	A	T,G	.	PASS	GMAF=0.367459;dbSNPBuildID=96,.
EOS
}

sub execute {
    my $self = shift;

    $self->check_minimum_version($MINIMUM_JOINX_VERSION);

    if($self->use_bgzip && not defined($self->output_file)){
       die $self->error_message("If use_bgzip is set, output_file must also be set, otherwise binary nonsense will spew forth."); 
    }
    my $output = "-";
    $output = $self->output_file if (defined $self->output_file);

    my $input_file = $self->input_file;

    unless(-s $input_file) {
        die $self->error_message("$input_file does not exist");
    }

    my $annotation_file = $self->annotation_file;
    unless(-s $annotation_file) {
        die $self->error_message("$annotation_file does not exist");
    }

    if($self->use_bgzip){
        $input_file = "<(zcat $input_file)";
    }

    my $cmd = $self->joinx_path . " vcf-annotate" . " --input-file $input_file" . " --annotation-file $annotation_file";
    if($self->info_fields) {
        my $info_fields = " --info-fields " . join(" --info-fields ", split /:/, $self->info_fields);
        $cmd .= $info_fields;
    }
    unless($self->identifiers) {
        $cmd .= " --no-identifiers";
    }
    unless($self->info) {
        $cmd .= " --no-info";
    }

    if(defined($self->output_file) && not defined($self->use_bgzip)){
        $cmd .= " --output-file $output" if defined($self->output_file);
    } elsif ( defined($self->use_bgzip) && defined($self->output_file) ){
        $cmd .= " | bgzip -c > $output";
        $cmd = "bash -c \"$cmd\"";
    }
    $self->status_message($cmd);
        
    my %params = (
        cmd => $cmd,
    );
    $params{output_files} = [$output] if $output ne "-";
    $params{skip_if_output_is_present} = 0;
    Genome::Sys->shellcmd(%params);

    return 1;
}


1;
