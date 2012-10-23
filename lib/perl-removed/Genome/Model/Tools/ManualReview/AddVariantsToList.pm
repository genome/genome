package Genome::Model::Tools::ManualReview::AddVariantsToList;

use strict;
use warnings;
use Genome::Utility::VariantReviewListReader;
use Genome;

UR::Object::Type->define(
    class_name => __PACKAGE__, is => 'Command',
    has => [
        list => {
            is => 'String',
            doc => 'Master csv to update from', 
        },
    ],
    has_optional => [
        db_list_name =>{
            is => 'String',
            doc => 'Name of list to be updated',
        },
        db_list_id => { 
            is => 'String', 
            doc => 'ID of list to be updated',
        },
        separation_character => {
            is => 'String',
            doc => 'separation character in list of variants to be added',
        },
    ]
);

sub help_brief{
    return "Add new variants to an existing list";
}

sub help_synopsis{
    return "gmt manual-review add-variants-to-list --list <list> --db_list_id <id>";
}

sub help_detail{
"Adds variants to existing lists given the list name or list id. The list must be a character delimited list.  The default delimiter is the '|' char, but this can be specified on the command line.  The list columns must be in the following format:\n".join('|', Genome::Utility::VariantReviewListReader->list_columns)."\n Once uploaded, the variants can be viewed and edited online at this address: https://gscweb.gsc.wustl.edu/view/variant_review_list.html";
}

BEGIN {
    my @det_property_names;
    my @props =
        sort { $a->property_name cmp $b->property_name }
        grep { $_->column_name ne '' }
        grep { defined($_->column_name) }
        Genome::VariantReviewDetail->get_class_object->all_property_metas;
    @det_property_names = map { $_->property_name } @props;
    sub fix_hash_data
    {
        my ($line_hash) = @_;
        my %temp_line_hash;
        foreach my $col (@det_property_names)
        {
            $temp_line_hash{$col} = $line_hash->{$col} if exists $line_hash->{$col};
        }
        $temp_line_hash{start_position} = $line_hash->{begin_position} if exists $line_hash->{begin_position};
        $temp_line_hash{stop_position} = $line_hash->{end_position} if exists $line_hash->{end_position};
        $temp_line_hash{stop_position} = $temp_line_hash{start_position} if !defined $temp_line_hash{stop_position};
        $line_hash = \%temp_line_hash;

        return $line_hash;

    }
}

sub execute{
    my $self = shift;
    my $list = Genome::Utility::VariantReviewListReader->new($self->list, $self->separation_character);
    my $db_list = Genome::VRList->get($self->db_list_name ? (name=>$self->db_list_name) : (id=>$self->db_list_id)); 
    unless ($db_list){
        $self->error_message("List doesn't exist");
        return 0;
    }
    my $list_id = $db_list->id;
    my ($detail) = $db_list->details;
    my $subject_name = $detail->subject_name;
    while (my $line_hash = $list->next_line_data()){
        last unless $line_hash;
        next if $line_hash->{header};

        $line_hash = fix_hash_data($line_hash);
        my $current_member = Genome::VariantReviewDetail->get( start_position => $line_hash->{start_position}, chromosome => $line_hash->{chromosome}, subject_name => $subject_name );
        if ($current_member){
            $self->status_message("Member found, skipping");
        }else{
            $current_member = Genome::VariantReviewDetail->create(
                    %$line_hash   
                    );
        }
        my $member_id = $current_member->id;
        my $review_list_member = Genome::VRListMember->get_or_create(
                list_id => $list_id,
                member_id => $member_id,
                );

    }  #while ( my $line = getline);
    return 1;
}

1;
