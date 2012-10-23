package Genome::Model::Tools::ManualReview::AssignVariantReviewers; 

use strict;
use warnings;
use Command;
use Genome;

UR::Object::Type->define(
    class_name => __PACKAGE__, is => 'Command',
    has => [
        finishers => {
            is => 'String',
            doc => 'Comma separated list of finishers to assign to variant list',
        },
        review_type => {
            is => 'String',
            doc => '3730 or manual',
        },
    ],
    has_optional => [
        list_name => {
            is => 'String',
            doc => 'Name of list', 
        },
        list_id => { 
            is => 'String', 
            doc => 'Id for list',
        },
    ]
);

sub help_brief{
    return "Assigns a finisher or list of finishers to an aml variant list for review";
}

sub help_synopsis{
    return "gmt manual-review assign-variant-reviewers --finishers 'adukes, pminx, jweible' --review-type 'manual' --list_name 'AML_ROUND_6.list'";
}

sub help_detail{
    return "assigns a finisher or group of finishers to reviewer status for a list in the AML variant db.  If multiple finishers are given (comma-separated list), variants are assigned in a rotating order to the finishers.  There are two review types to assign a finisher to: manual or 3730";
}

sub execute {
    my $self = shift;
    
    $DB::single = $DB::stopper;
    
    my $list_id     = $self->list_id;
    my $list_name   = $self->list_name;
    my @finishers   = split( /,/, $self->finishers);
    my $review_type = $self->review_type;
    my $review_list;
    unless ($list_id || $list_name){
        $self->error_message("Need to supply a list_id or list_name");
        return 0;
    }
    if ($list_id) {
       $review_list = Genome::VariantReviewList->get(id => $list_id); 
       unless ($review_list){
            $self->error_message("Couldn't get list from list_id: $list_id");
            return 0;
        }
        if ($list_name and $list_name ne $review_list->name){ 
            $self->error_message("list_name does not match name derived from list id");
            return 0;
        }
    }else {
        $review_list = Genome::VariantReviewList->get(name => $list_name);
        unless ($review_list){
         $self->error_message("Couldn't get list from list_name: $list_name"); 
         return 0;   
        }
        $list_id = $review_list->id;
    }
    
    my @bridges = Genome::VariantReviewListMember->get(list_id => $list_id); 
    my @details;
    foreach (@bridges){
        my $detail = Genome::VariantReviewDetail->get(id => $_->member_id);
        push @details, $detail;
    }
     
    my $col;
    if ($review_type =~ /3730/){
        $col = "finisher_3730_review";
    }
    elsif ($review_type =~ /manual/i){
        $col = "finisher_manual_review";
    }
    else{
        $self->error_message("Invalid review_type: $review_type");
        return 0;
    }

    foreach (@details){
        if (defined $_->$col){
            $self->status_message("finisher already assigned for detail at ".$_->chromosome." ".$_->begin_position." ".$_->end_position.". skipping");
        }else{
            my $assign = shift @finishers;
            $_->$col($assign);
            push @finishers, $assign;
        }
    }
    return 1;
}

1;
