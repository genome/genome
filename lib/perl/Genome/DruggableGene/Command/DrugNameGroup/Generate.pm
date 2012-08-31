package Genome::DruggableGene::Command::DrugNameGroup::Generate;

use strict;
use warnings;
use Genome;
use List::MoreUtils qw/ uniq /;

=cut
Concept Design Notes:

(  drug 1 )  (  drug 2  )  (drug 3)
  |  |   |    |   |   |      |   |
 a1  a2 a3    a2  a3  a4    a4  a5   <== alternate names, now hash keys to drugs

hash{a1} is drug1
hash{a3} is d1, d2
hash{a4} is d2, d3

for pubchem keys in hash (like a1)
is hash{a1} already in group?, add to group
multiple groups? merge them
else
Create group named after a1's hugo name
add all drugs from hash{a1} to this new group

Next, cycle through all drugs, and add to hugo groups
if name or alt names match exactly 1 group name,
or if alternate names map to exactly 1 group

Pseudo Code:

cycle trusted drug alt names (pubchem_primary_name)
  if group exists
    add this drug (referenced by this alt name)
  else
    create group (name=alternate_name)
    add this drug

cycle drugs
  skip if already in group

  find groups named identical to one of this drug's alt names or name
  if exactly 1 group is found, add drug to that group
      else
  find all groups with alt names identical to this drug's alt names or name
  if exactly 1 group is found, add drug to that group
      This step is made easy with a hash of alt name strings to alt name objects

=cut

class Genome::DruggableGene::Command::DrugNameGroup::Generate{
    is => 'Genome::Command::Base',
    doc => 'Generate a ton of groups to bundle drugs with similar alternate names',
};

sub help_brief { 'Generate a ton of groups to bundle drugs with similar alternate names' }

sub help_synopsis { help_brief() }

sub help_detail { help_brief() }

sub load {
    my $self = shift;

    print "Preloading all drugs\n";
    Genome::DruggableGene::DrugNameReport->get;

    print "Preloading all groups\n";
    Genome::DruggableGene::DrugNameGroup->get;
    Genome::DruggableGene::DrugNameGroupBridge->get;

    print "Loading alternate names and creating hash\n";
    my %alt_to_pubchem;
    my %alt_to_other;

    for (Genome::DruggableGene::DrugAlternateNameReport->get) { #operate on all alternate names
        my $alt = $_->alternate_name;
        print "Skipping $alt\n" and next if $alt =~ /^.$/;    #ignore single character names
        print "Skipping $alt\n" and next if $alt =~ /^\d\d$/; #ignore 2 digit names

        #Save genes with the same alternate name in an array in a hash with key being the alt-name
        if($_->nomenclature eq 'pubchem_primary_name'){
            push @{ $alt_to_pubchem{$alt} }, $_;
        } else {
            push @{ $alt_to_other{$alt} }, $_;
        }
    }

    return \%alt_to_pubchem, \%alt_to_other;
}

sub create_groups {
    my $self = shift;
    my $alt_to_pubchem = shift;
    my $progress_counter = 0;

    print "Putting " . scalar(keys(%{$alt_to_pubchem})) . " pubchem drug names into drug groups\n";
    for my $alt (keys %{$alt_to_pubchem}) {
        $progress_counter++;
        my @drugs = map{$_->drug} @{$alt_to_pubchem->{$alt}};

        my $group = Genome::DruggableGene::DrugNameGroup->get(name => $alt);
        if($group){ #hugo name group already exists
            for my $drug (@drugs){#make sure each drug is already in this group
                Genome::DruggableGene::DrugNameGroupBridge->create(
                    drug_id => $drug->id,
                    group_id => $group->id
                ) if not Genome::DruggableGene::DrugNameGroupBridge->get(drug_id => $drug->id);
            }
        }else{
            $group = Genome::DruggableGene::DrugNameGroup->create(name => $alt);
            Genome::DruggableGene::DrugNameGroupBridge->create(drug_id => $_->id, group_id => $group->id) for @drugs;
        }
        print "$progress_counter : created/refreshed group for $alt\n" if rand() < .001;
    }

    print "\n****\nFinished $progress_counter.\n****\n\n";
}

sub add_members {
    my $self = shift;
    my $alt_to_other = shift;
    my $progress_counter = 0;
    print "Now processing all " . scalar @{[Genome::DruggableGene::DrugNameReport->get]} . " to add members to groups\n";

    for my $drug (Genome::DruggableGene::DrugNameReport->get){
        $progress_counter++;
        next if Genome::DruggableGene::DrugNameGroupBridge->get(drug_id => $drug->id); #if already in a group


        my %indirect_groups;#groups found through alternate names
        my %direct_groups;#groups found through pubchem name
        my $drug_name = $drug->name;

        $direct_groups{$drug_name}++ if Genome::DruggableGene::DrugNameGroup->get(name=>$drug_name); #go drugs for instance have pubchem names

        for my $alt($drug->alternate_names){
            $direct_groups{$alt}++ if Genome::DruggableGene::DrugNameGroup->get(name=>$alt);

            my @alt_drugs = map{$_->drug} @{$alt_to_other->{$alt}};
            for my $alt_drug (@alt_drugs){
                #Get the group if it exists and add it to our list of indirectly found groups
                my $bridge = Genome::DruggableGene::DrugNameGroupBridge->get(drug_id => $alt_drug->id);
                $indirect_groups{$bridge->group->name}++ if $bridge;
            }
        }

        #If we found exactly one group, add this drug to it
        if(scalar keys %direct_groups == 1){
            my ($group_name) = keys %direct_groups;
            my $group_id = Genome::DruggableGene::DrugNameGroup->get(name=>$group_name)->id;
            Genome::DruggableGene::DrugNameGroupBridge->create(drug_id => $drug->id, group_id => $group_id);
            print "$progress_counter : added $drug_name to $group_name directly\n" if rand() < .01;
            next;
        }

        if(scalar keys %direct_groups == 0 and scalar keys %indirect_groups == 1){
            my ($group_name) = keys %indirect_groups;
            my $group_id = Genome::DruggableGene::DrugNameGroup->get(name=>$group_name)->id;
            Genome::DruggableGene::DrugNameGroupBridge->create(drug_id => $drug->id, group_id => $group_id);
            print "$progress_counter : added $drug_name to $group_name indirectly\n" if rand() < .01;
            next;
        }
        print "$progress_counter : failed to add $drug_name to any group. Direct: " .
        join(' ',keys %direct_groups) . "   Indirect: " . join(' ',keys %indirect_groups) . "\n" if rand() < .01;

    }
}

sub execute {
    my $self = shift;

    my ($alt_to_pubchem, $alt_to_other) = $self->load();
    $self->create_groups($alt_to_pubchem);
    $self->add_members($alt_to_other);

    return 1;
}
1;
