package Genome::Utility::MetagenomicClassifier::Weka;

use strict;
use warnings;

require Bio::Taxon;
use Data::Dumper 'Dumper';
use Genome::InlineConfig;
require Genome::Sys;
require Genome::Utility::MetagenomicClassifier;
require Genome::Utility::MetagenomicClassifier::SequenceClassification;

#$ENV{PERL_INLINE_JAVA_JNI} = 1;

use Inline(
    Java => <<'END', 
// Generated with Weka 3.6.0
//
// This code is public domain and comes with no warranty.
//
// Timestamp: Tue May 19 08:41:36 CDT 2009

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.Capabilities.Capability;
import weka.core.Instance;
import weka.core.Instances;
import weka.classifiers.Classifier;
import java.io.*;
import java.util.*;

public class WekaClassifier {

    public WekaClassifier() {
    }
    public static double classify(Object[] i)
        throws Exception {

        for (int ndx = 0; ndx < i.length;ndx++) {
            i[ndx] = Double.valueOf((String)i[ndx]);
        }

        double p = Double.NaN;
        p = WekaClassifier.N51dee10(i);
        return p;
    }
    static double N51dee10(Object []i) {
        double p = Double.NaN;
        if (i[1] == null) {
            p = 1;
        } else if (((Double) i[1]).doubleValue() <= 0.0) {
            p = WekaClassifier.Nff8f2c1(i);
        } else if (((Double) i[1]).doubleValue() > 0.0) {
            p = WekaClassifier.N10b0fd2(i);
        } 
        return p;
    }
    static double Nff8f2c1(Object []i) {
        double p = Double.NaN;
        if (i[5] == null) {
            p = 0;
        } else if (((Double) i[5]).doubleValue() <= -0.16) {
            p = 0;
        } else if (((Double) i[5]).doubleValue() > -0.16) {
            p = 1;
        } 
        return p;
    }
    static double N10b0fd2(Object []i) {
        double p = Double.NaN;
        if (i[1] == null) {
            p = 2;
        } else if (((Double) i[1]).doubleValue() <= 1.0) {
            p = WekaClassifier.N1d3c9ab3(i);
        } else if (((Double) i[1]).doubleValue() > 1.0) {
            p = WekaClassifier.N15d73d39(i);
        } 
        return p;
    }
    static double N1d3c9ab3(Object []i) {
        double p = Double.NaN;
        if (i[4] == null) {
            p = 0;
        } else if (((Double) i[4]).doubleValue() <= -0.07) {
            p = WekaClassifier.N15d1f6a4(i);
        } else if (((Double) i[4]).doubleValue() > -0.07) {
            p = 2;
        } 
        return p;
    }
    static double N15d1f6a4(Object []i) {
        double p = Double.NaN;
        if (i[5] == null) {
            p = 0;
        } else if (((Double) i[5]).doubleValue() <= -0.01) {
            p = WekaClassifier.Nb9a1445(i);
        } else if (((Double) i[5]).doubleValue() > -0.01) {
            p = WekaClassifier.N1c5d5e48(i);
        } 
        return p;
    }
    static double Nb9a1445(Object []i) {
        double p = Double.NaN;
        if (i[0] == null) {
            p = 2;
        } else if (((Double) i[0]).doubleValue() <= 1.0) {
            p = 2;
        } else if (((Double) i[0]).doubleValue() > 1.0) {
            p = WekaClassifier.N10a9daf6(i);
        } 
        return p;
    }
    static double N10a9daf6(Object []i) {
        double p = Double.NaN;
        if (i[0] == null) {
            p = 0;
        } else if (((Double) i[0]).doubleValue() <= 3.0) {
            p = 0;
        } else if (((Double) i[0]).doubleValue() > 3.0) {
            p = WekaClassifier.N1f67dee7(i);
        } 
        return p;
    }
    static double N1f67dee7(Object []i) {
        double p = Double.NaN;
        if (i[5] == null) {
            p = 0;
        } else if (((Double) i[5]).doubleValue() <= -0.42) {
            p = 0;
        } else if (((Double) i[5]).doubleValue() > -0.42) {
            p = 1;
        } 
        return p;
    }
    static double N1c5d5e48(Object []i) {
        double p = Double.NaN;
        if (i[3] == null) {
            p = 0;
        } else if (((Double) i[3]).doubleValue() <= 0.98) {
            p = 0;
        } else if (((Double) i[3]).doubleValue() > 0.98) {
            p = 2;
        } 
        return p;
    }
    static double N15d73d39(Object []i) {
        double p = Double.NaN;
        if (i[3] == null) {
            p = 0;
        } else if (((Double) i[3]).doubleValue() <= 0.79) {
            p = WekaClassifier.N5f246510(i);
        } else if (((Double) i[3]).doubleValue() > 0.79) {
            p = WekaClassifier.N1ebbb9312(i);
        } 
        return p;
    }
    static double N5f246510(Object []i) {
        double p = Double.NaN;
        if (i[0] == null) {
            p = 0;
        } else if (((Double) i[0]).doubleValue() <= 1.0) {
            p = WekaClassifier.N4e80d311(i);
        } else if (((Double) i[0]).doubleValue() > 1.0) {
            p = 0;
        } 
        return p;
    }
    static double N4e80d311(Object []i) {
        double p = Double.NaN;
        if (i[3] == null) {
            p = 2;
        } else if (((Double) i[3]).doubleValue() <= 0.28) {
            p = 2;
        } else if (((Double) i[3]).doubleValue() > 0.28) {
            p = 0;
        } 
        return p;
    }
    static double N1ebbb9312(Object []i) {
        double p = Double.NaN;
        if (i[4] == null) {
            p = 0;
        } else if (((Double) i[4]).doubleValue() <= -0.17) {
            p = WekaClassifier.N167d3c113(i);
        } else if (((Double) i[4]).doubleValue() > -0.17) {
            p = 2;
        } 
        return p;
    }
    static double N167d3c113(Object []i) {
        double p = Double.NaN;
        if (i[0] == null) {
            p = 2;
        } else if (((Double) i[0]).doubleValue() <= 1.0) {
            p = 2;
        } else if (((Double) i[0]).doubleValue() > 1.0) {
            p = WekaClassifier.N17f612514(i);
        } 
        return p;
    }
    static double N17f612514(Object []i) {
        double p = Double.NaN;
        if (i[3] == null) {
            p = 0;
        } else if (((Double) i[3]).doubleValue() <= 0.99) {
            p = 0;
        } else if (((Double) i[3]).doubleValue() > 0.99) {
            p = 2;
        } 
        return p;
    }
}
      
END

    AUTOSTUDY => 1,
    CLASSPATH => '/gsc/scripts/lib/java/weka.jar',
    STUDY => [
    'java.lang.Double',
    'java.util.ArrayList',
    ],
    PACKAGE => 'main',
    DIRECTORY => Genome::InlineConfig::DIRECTORY(),
    EXTRA_JAVA_ARGS => '-Xmx1000m',
#    JNI => 1,
) ;


sub create {
    my ($class, %params) = @_;
    my $self = bless \%params, $class;
    my $classifier = new WekaClassifier();
    $self->{classifier} = $classifier;
    return $self;
}

our @classes = ('unclassified','clean','chimera');

sub classify {
    my ($self, @args) = @_;

    my $name = shift @args;
    #my $class = pop @args;
    
    my $classifier = $self->{classifier};

    my $result = $classifier->classify(\@args);
    return $classes[$result];
}

1;

#$HeadURL: $
#$Id: $
