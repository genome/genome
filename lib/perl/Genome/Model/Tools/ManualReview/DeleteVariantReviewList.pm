package Genome::Model::Tools::ManualReview::DeleteVariantReviewList;

use strict;
use warnings;

use Genome;

UR::Object::Type->define(
    class_name => __PACKAGE__, is => 'Command',
    has => [
        db_list_name =>{
            is          => 'String',
            doc         => 'Name of list to be deleted',
            is_optional => 1,
        },
        db_list_id => { 
            is          => 'String', 
            doc         => 'ID of list to be deleted',
            is_optional => 1,
        },
    ]
);

sub help_brief{
    return "Delete a variant review list and all of the variant review details that are only in that list";
}

sub help_synopsis{
    return "gmt manual-review delete-variant-review-list --db-list-name <list-name>";
}

sub help_detail{
"Deletes a variant review list from the database. Also deletes all variant review details that have no other links to other lists.  Use with caution.  Takes in either the list id, or the list name";
}


sub execute{
    my ($self) = @_;
    unless($self->db_list_name||$self->db_list_id)
    {
        $self->error_message("Need either a db_list_name or db_list_id.\n");
        return 0;
    }
    my $db_list = Genome::VRList->get($self->db_list_name ? (name=>$self->db_list_name) : (id=>$self->db_list_id)); 
    
    unless ($db_list){
        $self->error_message("List doesn't exist");
        return 0;
    }
    
    my @db_list_members = Genome::VRListMember->get(list_id=>$db_list->id);    
    foreach (@db_list_members){
        my $detail = Genome::VariantReviewDetail->get(id=>$_->member_id);
        my @test_members = Genome::VRListMember->get(member_id=>$detail->id);
        if (@test_members > 1){
           my %list_id_hash; 
           foreach (@test_members){
                $list_id_hash{$_->list_id}++;
            }
            if (scalar keys %list_id_hash == 1){
                $detail->delete;
            }
        }else {
            $detail->delete;
        }

        $_->delete;
        
    }
    $db_list->delete;
    return 1;    
}

1;
