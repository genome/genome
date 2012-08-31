package Genome::PopulationGroup::Command;

use strict;
use warnings;

use Genome::Command::Crud;
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::PopulationGroup',
    target_name => 'population group',
    create => { 
        exclude => [qw/ member_ids member_hash /], 
        add => {
            models => { is => 'Genome::Model', is_many => 1, is_optional => 1, doc => 'Models to get individuals from. Resolved via string.', require_user_verify => 0, },
            member_ids => { is => 'Number', is_many => 1, is_optional => 1 , doc => 'Member ids', },
        },
        before => sub{
            my ($create) = @_;
            my %seen;
            my @models = $create->models;
            print "\n\n***Found ".@models." models\n\n";
            my @attrs = Genome::SubjectAttribute->get(
                'subject_id in' => [ map { $_->subject_id } @models ],
                attribute_label => 'source_id',
            );
            print "\n\n***Found ".@attrs." attrs\n\n";
            $create->member_ids([ map { $_->attribute_value } @attrs ]);
            $create->models([]);
        },
    },
    list => { show => 'id,name,description', },
    update => { only_if_null => 1, exclude => [qw/ member_hash /], },
    delete => { do_not_init => 1, },
);

1;

