#----------------------------------
# $Author: bshore $ 
# $Date: 2008-07-01 09:11:56 -0500 (Tue, 01 Jul 2008) $
# $Revision: 36100 $
# $URL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/MG/MutationDiagram/LayoutManager.pm $
#----------------------------------
package Genome::Model::Tools::Graph::MutationDiagram::MutationDiagram::LayoutManager;
#------------------------------------------------
our $VERSION = '1.0';
#------------------------------------------------
use strict;
use warnings;
use FileHandle;
use Carp;

# If you know about and want to use inheritance:
#use Base::Class;
#our @ISA = qw( Base::Class );

sub new {
    my ($class, %arg) = @_;
    my $self = {
                _iterations     => $arg{iterations} || 100,
                _max_distance   => $arg{max_distance},
                _spring_constant => $arg{spring_constant},
                _spring_force   => $arg{spring_force},
                _attractive_weight => $arg{attractive_weight},
    };
    bless($self, ref($class) || $class);
    return $self;
}
#-------------------------------------------------
sub layout {
  my ($self, @mutations) = @_;
  @mutations = sort {$a->get_label_position <=> $b->get_label_position} @mutations;
  for(my $t = 0; $t < $self->{_iterations}; $t++) {
      $self->repulse_labels(@mutations);
  }
}
#----------------------------------
sub repulse_labels {
    my ($self, @mutations) = @_;
    #mutations should be sorted in the x-axis
    #in ascending order
    my @forces = (0) x scalar(@mutations);
    for(my $i = 0; $i < scalar(@mutations); $i++) {
        for(my $j = 1; $j < scalar(@mutations); $j++) {
            $self->calculate_repulsive_force(mutation1 => {mutation => $mutations[$i],
                    force => \$forces[$i]},
                mutation2 => {mutation => $mutations[$j], force =>
                    \$forces[$j]});
        }
    }
    # for(my $i = 0; $i < scalar(@mutations);  $i++) {
    #     #move based on force
    #     #did not set a maximum movement (mistake?)
    #     $self->calculate_attractive_force(mutation1 => {mutation => $mutations[$i],
    #                 force => \$forces[$i]}),

    #     }

    for(my $i = 0; $i < scalar(@mutations);  $i++) {
        #move based on force
        #did not set a maximum movement (mistake?)
        my $xmove =  $self->{_spring_force} * $forces[$i];
		my $max = 1;
		$xmove = $max if $xmove > $max;
		$xmove = -$max if $xmove < -$max;
        $mutations[$i]->translate_label_position($xmove);
        }
}
#----------------------------------
sub calculate_repulsive_force {
    my ($self, %mutations) = @_;
    #hash should contain the two mutations and their forces
    my $mutation1 = $mutations{mutation1}{mutation};
    my $mutation2 = $mutations{mutation2}{mutation};
    my $force1 = $mutations{mutation1}{force};
    my $force2 = $mutations{mutation2}{force};
    
	my $dx = $mutation2->get_label_position -
	         $mutation1->get_label_position;

	my $d2 = $dx * $dx;
	if ($d2 < 0.01) {
		$dx = rand (0.1) + 0.1;
		$d2 = $dx * $dx;
	}

	my $d = sqrt $d2;

	if ($d < $self->{_max_distance}) {
		my $repulsive_force = $self->{_spring_constant} * $self->{_spring_constant} / $d;
        #changing the forces by reference. Yay!
		$$force2 += $repulsive_force * $dx / $d;
        $$force1 -= $repulsive_force * $dx / $d;
	}
    return;

}
#----------------------------------
sub calculate_attractive_force {
    my ($self, %mutations) = @_;
    #hash should contain the two mutations and their forces
    my $mutation = $mutations{mutation1}{mutation};
    my $force = $mutations{mutation1}{force};
    
	my $dx = $mutation->amino_acid - $mutation->get_label_position;
	         

	my $d2 = $dx * $dx;
    #if ($d2 < 0.01) {
	#	$dx = rand (0.1) + 0.1;
	#	$d2 = $dx * $dx;
	#}

	my $d = sqrt $d2;

	if ($d > $self->{_max_distance}) {
        $d = $self->{_max_distance};
        $d2 = $d * $d;
    }
    if($d != 0) {
		my $attractive_force = ($d2 - $self->{_spring_constant} * $self->{_spring_constant}) / $self->{_spring_constant};
        my $weight = $self->{_attractive_weight};
        $weight = 1 if not $weight or $weight < 1;
        $attractive_force *= log($weight) * 0.5 + 1;

        #changing the forces by reference. Yay!
        $$force += $attractive_force * $dx / $d;
    }
    return;

}

=head1 AUTHOR

David Larson, E<lt>dlarson@watson.wustl.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 by David Larson

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8 or,
at your option, any later version of Perl 5 you may have available.


=cut
1;

