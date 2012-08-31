package Genome::Model::Tools::ManualReview::ReviewVariants;
use strict;
use warnings;
use Genome::Model::Tools::ManualReview::MRGui;
use Command;
use Data::Dumper;
use IO::File;
use PP::LSF;
use File::Temp;
use File::Basename;
class Genome::Model::Tools::ManualReview::ReviewVariants
{
    is => 'Command',                       
    has => 
    [ 
        project_file => 
        {
            type => 'String',
            is_optional => 1,
            doc => "Manual review csv file",
        },

    ], 
};
############################################################
sub help_brief {   
    return;
}
sub help_synopsis { 
    return;
}
sub help_detail {
    return <<EOS 
    Launches the variant review editing application.
EOS
}
############################################################
sub execute { 
    my $self = shift;
    my $project_file = $self->project_file;
    

    my $mr = Genome::Model::Tools::ManualReview::MRGui->new(project_file => $project_file);    

    Gtk2->main();
    return 1;
}
############################################################
1;
