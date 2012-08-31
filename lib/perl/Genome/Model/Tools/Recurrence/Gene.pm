package Genome::Model::Tools::Recurrence::Gene;

use warnings;
use strict;

use Genome;
use Carp;
use IO::File;
use Genome::Info::IUB;

class Genome::Model::Tools::Recurrence::Gene{
    is => 'Command',
    has => [
    somatic_anno_file => {
        is  => 'String',
        is_input  => 1,
        doc => 'The list of annotated tumor SNPs',
    },
    output_file => {
        is => 'Text',
        is_input => 1,
        is_output => 1,
        doc => "horribly formatted output for possible recurrences"
    },
    source_somatic_model_id=> {
        is => 'Text',
        is_input => 1,
        is_optional => 1,
        doc => "the model id that produced these variants. prevents finding yourself as a recurrence.",
    },
    skip_if_output_present => {
        is => 'Boolean',
        is_optional => 1,
        is_input => 1,
        default => 0,
        doc => 'enable this flag to shortcut through annotation if the output_file is already present. Useful for pipelines.',
    },
    putative => {
        is => 'Boolean',
        default=> 0,
        is_optional=> 1,
    },
    skip_MT => {
        is => 'Boolean',
        default=> 1,
        is_optional=> 1,
    },
]
};

sub help_brief {
    "Teach Dave what REAL recurrence looks like",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
    not filling this part out
EOS
}

sub help_detail {                           
    return <<EOS 
    This tool is one of the greatest tools ever imaginated.
EOS
}

sub execute {
    my $self = shift;

    unless(-f $self->somatic_anno_file) {
        $self->error_message($self->tumor_snp_file . " is not a file");
        die;
    }
    my $fh = IO::File->new($self->somatic_anno_file);
    my $output_fh = IO::File->new($self->output_file, ">");

    my @builds = Genome::Model::Build->get(model_id=>$self->source_somatic_model_id);
    my @build_ids = map { $_->build_id } @builds;
    while (my $line = $fh->getline) {

        my ($chr, $start, $stop, $ref, $var, $type, $gene,) = split "\t", $line;;
        if($chr eq 'MT' && $self->skip_MT) {
            next;
        }
        
        $DB::single=1 if ($gene eq 'KRAS');
        my @EXCITING_VARIANTS = Genome::Model::Variant->get(gene_name=>$gene);
        if(@EXCITING_VARIANTS) {
            my $found_in_other_model=0;
            my $build_id=0;

            for my $variant (@EXCITING_VARIANTS) {
              my @potential_recurrences = Genome::Model::BuildVariant->get(variant_id=>$variant->id);
                for my $link (@potential_recurrences) {
                    my @found = grep /$link->build_id/,  @build_ids;
                    if (scalar(@found) == 1) { 
                        next;
                    }
                    $found_in_other_model =1;
                    $build_id = $link->build_id;
                }
             

            if($found_in_other_model == 1) {    
                my $official_validation= $variant->get_official_validation_for_build($build_id);
                my $model_name=  Genome::Model::Build->get($build_id)->model->name;
                my $validation_status = $official_validation->validation_result;
                my $temp = $variant->chromosome . "\t" . $variant->start_pos . "\t";
                $temp .= $variant->stop_pos . "\t" . $variant->reference_allele . "\t";
                $temp .= $variant->variant_allele . "\t" . $variant->gene_name . "\t" . $validation_status . "\t" . $model_name;
                if($validation_status ne 'P' || $self->putative) {
                    $output_fh->print($line);
                    $output_fh->print( $temp . "\n");
                }
            }
        }
    }
}
    return 1;
}

1;
