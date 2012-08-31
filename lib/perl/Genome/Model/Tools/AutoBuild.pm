package Genome::Model::Tools::AutoBuild;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::AutoBuild { is => 'Command', };

sub sub_command_sort_position { 24 }

sub help_brief {
"Tool to build solexa reference-alignment genome models that are configured for automatic processing"
      ,;
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt auto-build ...
EOS
}

sub help_detail {
    return <<EOS
EOS
}

sub execute {
    my $self = shift;

    #needed to get around a UR bug that tony is working on - je
    Genome::Model::ReferenceAlignment->get();


    foreach my $build_model (
        Genome::Model::ReferenceAlignment->get( auto_build_alignments => 1 ) )
    {
        Genome::Model::Build::Command::ScheduleStage->execute(stage_name => "alignment", model_id => $build_model->genome_model_id());
    }
    return 1;
}

1;
