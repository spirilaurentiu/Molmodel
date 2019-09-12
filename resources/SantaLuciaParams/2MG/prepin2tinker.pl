#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my %CommandLineOptions;
getopt('tpfrR', \%CommandLineOptions);

my $residueName = $CommandLineOptions{r};
my $pdbResidueName = $CommandLineOptions{R};
my $prepinFile = $CommandLineOptions{p};
my $frcmodFile = $CommandLineOptions{f};
my $olderParamFile = $CommandLineOptions{t};

my $usage = 
"Usage: $0 -r <residueName> -R <pdbResidueName> -p <file.prepin> -f <file.frcmod> -t <tinker.param> [-c]\n
  -c : Generate C++ code for SimTK molmodel, instead of generating tinker parameters\n"
;

die $usage unless defined $residueName;
die $usage unless defined $prepinFile;
die $usage unless defined $frcmodFile;
die $usage unless defined $olderParamFile;
die $usage unless defined $pdbResidueName;

### INPUT ###
#     0    0    2
# 
# This is a remark line
# 2MG.res
# 2MG   INT  0
# CORRECT     OMIT DU   BEG
#   0.0000
#    1  DUMM  DU    M    0  -1  -2     0.000      .0        .0      .00000
#    2  DUMM  DU    M    1   0  -1     1.449      .0        .0      .00000
#    3  DUMM  DU    M    2   1   0     1.522   111.1        .0      .00000
#    4  P     P     M    3   2   1     1.540   111.208   180.000   1.08783
#    5  O1P   O     E    4   3   2     1.476    87.087  -175.953  -0.76666
#    6  O2P   O     E    4   3   2     1.495    34.379   -14.053  -0.76666
#    7  O5'   OS    M    4   3   2     1.599   135.795   -64.492  -0.47131
#    8  C5'   CT    M    7   4   3     1.447   120.649    80.070   0.06353
#    9  H5'1  H1    E    8   7   4     1.070   106.582   -71.398   0.06892
#   10  H5'2  H1    E    8   7   4     1.070   106.568    57.898   0.06892

### OUTPUT ###
# atom      1    14    N       "Glycine N"                 7     14.010     3
# atom      2     1    CT      "Glycine CA"                6     12.010     4

   #####################################################
   ##                                                 ##
   ##  TINKER Atom Class Numbers to Amber Atom Types  ##
   ##                                                 ##
   ##    1  CT      11  CN      21  OW      31  HO    ##
   ##    2  C       12  CK      22  OH      32  HS    ##
   ##    3  CA      13  CQ      23  OS      33  HA    ##
   ##    4  CM      14  N       24  O       34  HC    ##
   ##    5  CC      15  NA      25  O2      35  H1    ##
   ##    6  CV      16  NB      26  S       36  H2    ##
   ##    7  CW      17  NC      27  SH      37  H3    ##
   ##    8  CR      18  N*      28  P       38  HP    ##
   ##    9  CB      19  N2      29  H       39  H4    ##
   ##   10  C*      20  N3      30  HW      40  H5    ##
   ##                                                 ##
   #####################################################

# Create atom records in tinker parameter file format

sub writeTinkerParameters
{
    my $atoms = shift;
    my $atomTypes = shift;
    my $residueName = shift;
    my $atom;
    
    
    
    # atom records
    
    foreach $atom (@$atoms) 
    {
        next unless ref $atom;
        
        my $atomType = $atomTypes->{$atom->{TYPE}};
        
        printf "atom %6d %5d    %-7s %-25s %3d %10.3f %2d\n", 
            $atom->{INDEX}, 
            $atomType->{TINKER_INDEX}, 
            $atomType->{NAME}, 
            "\"".$residueName." ".$atom->{NAME}."\"",
            $atomType->{ATOMIC_NUMBER},
            $atomType->{MASS},
            $atom->{NBONDS} + 0;
    }
            
    # Blank line
    print "\n";
        
    # vdw radius
    
    foreach my $atomTypeName (%$atomTypes) {
        my $atomType = $atomTypes->{$atomTypeName};
        next unless ref $atomType;
        next unless exists $atomType->{VDW_RADIUS};
        
        printf "vdw       %4d          %10.4f %10.4f\n",
            $atomType->{TINKER_INDEX},
            $atomType->{VDW_RADIUS},
            $atomType->{VDW_DEPTH};
    }
            
    # Blank line
    print "\n";
        
    # bond stretch
    
    foreach my $atomTypeName (%$atomTypes) {
        my $atomType = $atomTypes->{$atomTypeName};
        next unless ref $atomType->{STRETCH};
        foreach my $atomTypeName2 ( %{$atomType->{STRETCH}} ) {
            my $atomType2 = $atomTypes->{$atomTypeName2};
            next unless ref $atomType2;
                
                printf "bond      %4d %4d     %10.1f %10.4f\n", 
                    $atomType->{TINKER_INDEX}, 
                    $atomType2->{TINKER_INDEX}, 
                    $atomType->{STRETCH}{$atomTypeName2}{STIFFNESS}, 
                    $atomType->{STRETCH}{$atomTypeName2}{DISTANCE};
        }
    }
    
    # Blank line
    print "\n";
        
    # bond bend
    
    foreach my $atomTypeName (%$atomTypes) {
        my $atomType = $atomTypes->{$atomTypeName};
        next unless ref $atomType->{ANGLE};
        foreach my $atomTypeName2 ( %{$atomType->{ANGLE}} ) {
            my $atomType2 = $atomTypes->{$atomTypeName2};
            next unless ref $atomType2;
            foreach my $atomTypeName3 ( %{$atomType->{ANGLE}{$atomTypeName2}} ) {
                my $atomType3 = $atomTypes->{$atomTypeName3};
                next unless ref $atomType3;
                
                my $angle = $atomType->{ANGLE}{$atomTypeName2}{$atomTypeName3}{ANGLE};
                my $sigma = $atomType->{ANGLE}{$atomTypeName2}{$atomTypeName3}{SIGMA};
                
                printf "angle     %4d %4d %4d%10.2f %10.2f\n", $atomType->{TINKER_INDEX}, $atomType2->{TINKER_INDEX}, $atomType3->{TINKER_INDEX}, $sigma, $angle;
            }
        }
    }
    
    print "\n";
    
    # torsion
    
    foreach my $atomTypeName (%$atomTypes) {
        my $atomType = $atomTypes->{$atomTypeName};
        next unless ref $atomType->{TORSION};
        foreach my $atomTypeName2 ( %{$atomType->{TORSION}} ) {
            my $atomType2 = $atomTypes->{$atomTypeName2};
            next unless ref $atomType2;
            foreach my $atomTypeName3 ( %{$atomType->{TORSION}{$atomTypeName2}} ) {
                my $atomType3 = $atomTypes->{$atomTypeName3};
                next unless ref $atomType3;
                
                foreach my $atomTypeName4 ( %{$atomType->{TORSION}{$atomTypeName2}{$atomTypeName3}} ) {
                    my $atomType4 = $atomTypes->{$atomTypeName4};
                    next unless ref $atomType4;
                
                    my $torsion = $atomType->{TORSION}{$atomTypeName2}{$atomTypeName3}{$atomTypeName4};
                    next unless ref $torsion;
                    
                    printf "torsion   %4d %4d %4d %4d", 
                        $atomType->{TINKER_INDEX}, 
                        $atomType2->{TINKER_INDEX}, 
                        $atomType3->{TINKER_INDEX}, 
                        $atomType4->{TINKER_INDEX};
                    
                    if ( 1 == (scalar @{$torsion}) ) {
                        printf "     %10.3f %6.1f %2d", 
                            $torsion->[0]{SIGMA},
                            $torsion->[0]{ANGLE},
                            $torsion->[0]{PERIOD}
                    }
                    else {
                        foreach my $trs (@{$torsion}) {
                            printf "     %10.3f %6.1f %2d", 
                                $trs->{SIGMA},
                                $trs->{ANGLE},
                                $trs->{PERIOD}
                            }
                    }
                        
                    print "\n";
                }
            }
        }
    }
    
    print "\n";
    
    # improper torsion
    
    foreach my $atomTypeName (%$atomTypes) {
        my $atomType = $atomTypes->{$atomTypeName};
        next unless ref $atomType->{IMPTORS};
        foreach my $atomTypeName2 ( %{$atomType->{IMPTORS}} ) {
            my $atomType2 = $atomTypes->{$atomTypeName2};
            next unless ref $atomType2;
            foreach my $atomTypeName3 ( %{$atomType->{IMPTORS}{$atomTypeName2}} ) {
                my $atomType3 = $atomTypes->{$atomTypeName3};
                next unless ref $atomType3;
                
                foreach my $atomTypeName4 ( %{$atomType->{IMPTORS}{$atomTypeName2}{$atomTypeName3}} ) {
                    my $atomType4 = $atomTypes->{$atomTypeName4};
                    next unless ref $atomType4;
                
                    my $torsion = $atomType->{IMPTORS}{$atomTypeName2}{$atomTypeName3}{$atomTypeName4};
                    next unless ref $torsion;
                    
                    printf "imptors   %4d %4d %4d %4d", 
                        $atomType->{TINKER_INDEX}, 
                        $atomType2->{TINKER_INDEX}, 
                        $atomType3->{TINKER_INDEX}, 
                        $atomType4->{TINKER_INDEX};
                    
                        printf "     %10.3f %6.1f %2d", 
                            $torsion->{SIGMA},
                            $torsion->{ANGLE},
                            $torsion->{PERIOD};
                        
                    print "\n";
                }
            }
        }
    }
    
    print "\n";
    
    # charge records
    # Create charge records in tinker parameter file format

    foreach $atom (@$atoms) 
    {
        next unless ref $atom;
        
        printf "charge %7d %13.4f\n", $atom->{INDEX}, $atom->{CHARGE};
    }

    # Blank line
    print "\n";
        
    # biotype records
    # Create biotype records in tinker parameter file format

    foreach $atom (@$atoms) 
    {
        next unless ref $atom;
        
        printf "biotype %6d    %-7s %-21s %7d\n", 
            $atom->{INDEX}, 
            $atom->{NAME},
            "\"".$residueName."\"",
            $atom->{INDEX};
    }
}

my %atomTypes = ();
my $maxTinkerAtomIndex = 0;

open OLDERPARAMS, "<$olderParamFile" or die;
while (<OLDERPARAMS>) {
    if ($_ =~ m/^atom /) 
    {
        # atom      1    14    N       "Glycine N"                 7     14.010     3
        # atom      2     1    CT      "Glycine CA"                6     12.010     4
        # atom      3     2    C       "Glycine C"                 6     12.010     3
        # atom      4    29    H       "Glycine HN"                1      1.008     1
        die $_ unless $_ =~ m/
            ^atom
            \s+(\d+) # Tinker atom number
            \s+(\d+) # Atom type number
            \s+(\S+) # Atom type name
            \s+\"[^\"]*\" # residue atom name string
            \s+(\d+) # atomic number
            \s+([0-9]+\.[0-9]+) # atomic mass
            \s+(\d+) # number of bonds
            /x;
        
        $atomTypes{$3}{TINKER_INDEX} = $2;
        $atomTypes{$3}{ATOMIC_NUMBER} = $4;
        $atomTypes{$3}{MASS} = $5;
        $atomTypes{$3}{NAME} = $3;
        
        if ($1 > $maxTinkerAtomIndex) {$maxTinkerAtomIndex = $1;}
    }
}
close OLDERPARAMS;

open FRCMOD, "<$frcmodFile" or die;
my $current_record_type = "None";
while (<FRCMOD>) 
{
    if ($_ =~ m/^MASS/) {
        $current_record_type = "MASS";
        next;
    }
    
    if ($_ =~ m/^BOND/) {
        $current_record_type = "BOND";
        next;
    }
    
    if ($_ =~ m/^ANGLE/) {
        $current_record_type = "ANGLE";
        next;
    }
    
    if ($_ =~ m/^DIHE/) {
        $current_record_type = "DIHE";
        next;
    }
    
    if ($_ =~ m/^IMPROPER/) {
        $current_record_type = "IMPROPER";
        next;
    }
    
    if ($_ =~ m/^NONBON/) {
        $current_record_type = "NONBON";
        next;
    }
    
    if ($current_record_type eq "MASS") {
        next if $_ =~ m/^\s*$/; # empty line
        die $_ unless $_ =~ m/^(\S+)\s+([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+).*$/;
        $atomTypes{$1}{MASS} = $2;
        $atomTypes{$1}{POLARIZABILITY} = $2;
        $atomTypes{$1}{NAME} = $1;
    }
    
    elsif ($current_record_type eq "BOND") {
        next if $_ =~ m/^\s*$/; # empty line
        
        # P -O   456.40   1.485       same as o -p4 Autocorrect

        die $_ unless $_ =~ m/^(\S+)\s?-(\S+)\s?\s+([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+).*$/;
        $atomTypes{$1}{STRETCH}{$2}{STIFFNESS} = $3;
        $atomTypes{$1}{STRETCH}{$2}{DISTANCE} = $4;
    }
    
    elsif ($current_record_type eq "ANGLE") {
        next if $_ =~ m/^\s*$/; # empty line
        
        # P -OS-CT   77.600     120.649   same as c3-os-p4 Autocorrect

        die $_ unless $_ =~ m/^(\S+)\s?-(\S+)\s?-(\S+)\s+([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+).*$/;
        $atomTypes{$1}{ANGLE}{$2}{$3}{SIGMA} = $4;
        $atomTypes{$1}{ANGLE}{$2}{$3}{ANGLE} = $5;
    }
    
    elsif ($current_record_type eq "DIHE") {
        next if $_ =~ m/^\s*$/; # empty line
        
        # P -OS-CT   77.600     120.649   same as c3-os-p4 Autocorrect

        die $_ unless $_ =~ m/^(\S+)\s?-(\S+)\s?-(\S+)\s?-(\S+)\s+(\d+)\s+([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+).*$/;
        $atomTypes{$1}{TORSION}{$2}{$3}{$4}[$5-1]{SIGMA} = $6;
        $atomTypes{$1}{TORSION}{$2}{$3}{$4}[$5-1]{ANGLE} = $7;
        $atomTypes{$1}{TORSION}{$2}{$3}{$4}[$5-1]{PERIOD} = $8;
    }
    
    elsif ($current_record_type eq "IMPROPER") {
        next if $_ =~ m/^\s*$/; # empty line
        
        # CB-CK-N*-CT         1.1          180.0         2.0          Using default value

        die $_ unless $_ =~ m/^(\S+)\s?-(\S+)\s?-(\S+)\s?-(\S+)\s+([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+).*$/;
        $atomTypes{$1}{IMPTORS}{$2}{$3}{$4}{SIGMA} = $5;
        $atomTypes{$1}{IMPTORS}{$2}{$3}{$4}{ANGLE} = $6;
        $atomTypes{$1}{IMPTORS}{$2}{$3}{$4}{PERIOD} = $7;
    }
    
    elsif ($current_record_type eq "NONBON") {
        next if $_ =~ m/^\s*$/; # empty line
        
        #  P           2.1000  0.2000             same as p4 

        die $_ unless $_ =~ m/^\s+(\S+)\s+([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+).*$/;
        $atomTypes{$1}{VDW_RADIUS} = $2;
        $atomTypes{$1}{VDW_DEPTH} = $3;
    }
}
close FRCMOD;



open PREPIN, "<$prepinFile" or die;

my $next_atom_index = $maxTinkerAtomIndex + 1;
my @atoms = ();
my %atomsByName = ();

$current_record_type = "ATOM";
while (<PREPIN>) 
{
    if ($_ =~ m/^LOOP/) {
        $current_record_type = "LOOP";
        next;
    }
    
    if ($_ =~ m/^IMPROPER/) {
        $current_record_type = "IMPROPER";
        next;
    }
   
    # Increase bond count for atoms in LOOP records
    
    if ($current_record_type eq "LOOP") {
        next if $_=~ m/^\s*$/; # blank line
        
        die unless $_ =~ m/^\s*(\S+)\s+(\S+)\s*$/;
        
        my $atom1 = $atomsByName{$1};
        my $atom2 = $atomsByName{$2};
        
        if (ref $atom1) {
            ++$atom1->{NBONDS};
            push @{$atom1->{CHILD_ATOMS}}, $atom2 if ref $atom2;
            push @{$atom1->{RING_BONDS}}, $atom2 if ref $atom2;
        }
        if (ref $atom2) {
            ++$atom2->{NBONDS};
            push @{$atom2->{CHILD_ATOMS}}, $atom1 if ref $atom1;
        }
    }
    
    if ($current_record_type eq "ATOM") {
        if ($_ =~ m/^          # anchor at the beginning of the line
                    \s*(\d+)   # atom number
                    \s+(\S+)   # atom name
                    \s+(\S+)   # atom type
                    \s+[MES3B] # main vs. side
                    \s+(\d+)   # parent atom 1
                    \s+(\d+)   # parent atom 2
                    \s+(\d+)   # parent atom 3
                    \s+([012]\.[0-9]{3})+   # bond length
                    \s+([0-9]*\.[0-9]+)+    # bond angle
                    \s+(-?[0-9]*\.[0-9]+)+  # dihedral angle
                    \s+(-?[0-9]*\.[0-9]+)+  # atomic partial charge
                    \s*$       # anchor end of line
                    /x)        # allow whitespace
        {
            my %atom = ();
            $atom{INDEX} = $next_atom_index;
            $atom{PREPIN_INDEX} = $1;
            $atom{NAME} = $2;
            $atom{TYPE} = $3;
            $atom{PARENT} = $4;
            $atom{PARENT2} = $5;
            $atom{PARENT3} = $6;
            $atom{BOND_LENGTH} = $7;
            $atom{BOND_ANGLE} = $8;
            $atom{DIHEDRAL_ANGLE} = $9;
            $atom{CHARGE} = $10;
            
            next if $atom{TYPE} eq "DU";

            die $atom{TYPE} unless exists $atomTypes{$atom{TYPE}};
            
            my $atomType = $atomTypes{$atom{TYPE}};
            
            die $atom{TYPE} unless exists $atomType->{TINKER_INDEX};

            # Increase bond counts
            
            if ($atom{PARENT} > 0) {
                ++$atom{NBONDS};
            }
            my $parentAtom = $atoms[$atom{PARENT}];
            if (ref $parentAtom) {
                ++$parentAtom->{NBONDS};
                push @{$parentAtom->{CHILD_ATOMS}}, \%atom;
            }

            ++$next_atom_index;
            
            $atoms[$atom{PREPIN_INDEX}] = \%atom;
            $atomsByName{$atom{NAME}} = \%atom;
        }
    }
}
close PREPIN;

sub elementName 
{
    my $atomicNumber = shift;
    return "Hydrogen"   if  1 == $atomicNumber;
    return "Carbon"     if  6 == $atomicNumber;
    return "Nitrogen"   if  7 == $atomicNumber;
    return "Oxygen"     if  8 == $atomicNumber;
    return "Phosphorus" if 15 == $atomicNumber;
    return "Sulfur"     if 16 == $atomicNumber;
    die "Unknown element number $atomicNumber";
}

# print out code that initializes a new SingleAtomCompound() in molmodel
sub atomCode
{
    my $atom = shift;
    my $atomTypes = shift;
        
    my $valence = $atom->{NBONDS};
    die unless defined $valence;
    die unless $valence >= 1;
    die unless $valence <= 4;
    
    print "UnivalentAtom"    if 1 == $valence;
    print "BivalentAtom"     if 2 == $valence;
    print "TrivalentAtom"    if 3 == $valence;
    print "QuadrivalentAtom" if 4 == $valence;
    
    print "(\"" . $atom->{NAME} . "\", Element::";
    
    my $atomicNumber = $atomTypes->{$atom->{TYPE}}{ATOMIC_NUMBER};
    print elementName($atomicNumber);
    
    print "()";

    # Set bond angle
    if (2 == $valence) {
        my $angle = $atom->{CHILD_ATOMS}[0]{BOND_ANGLE};
        print ", $angle*Deg2Rad";
    }
    elsif (3 == $valence) {
        my $angle1 = $atom->{CHILD_ATOMS}[0]{BOND_ANGLE};
        print ", $angle1*Deg2Rad";
        # This might be a ring-closing bond...
        my $angle2 = $atom->{CHILD_ATOMS}[1]{BOND_ANGLE};
        print ", $angle2*Deg2Rad";
    }
    
    print ")";
}

sub writeMolmodelCode 
{
    my $atoms = shift;
    my $atomTypes = shift;
    my $residueName = shift;
    my $pdbResidueName = shift;

    print "class $residueName : public Compound {\n";
    print "public:\n";

    # Constructor
    print "    $residueName()\n";
    print "    {\n";
    
    print "        instantiateBiotypes();\n\n";
    print "        setPdbResidueName(\"$pdbResidueName\");\n\n";
    print "        setCompoundName(\"$residueName\");\n\n";

    my $atom;
    my $atomCount = 0;
    foreach $atom (@$atoms) {
        next unless ref $atom;
        $atom->{NEXT_CHILD_BOND} = 2;
        
        # First atom is not bonded
        if ($atomCount == 0) {
            print "        setBaseAtom( ";
            atomCode($atom, $atomTypes);
            print " );\n\n";
        }
        
        else {
            print "        bondAtom( ";
            atomCode($atom, $atomTypes);
            
            my $parentAtom = $atoms->[$atom->{PARENT}];
            
            # Parent bond center
            # TODO - consider chirality for tetrahedral centers
            my $bondNumber;
            if ( ($parentAtom->{NBONDS} == 4) && ($parentAtom->{NEXT_CHILD_BOND} > 2) ) 
            {
                my $angleChange = $atom->{DIHEDRAL_ANGLE} - $parentAtom->{FIRST_DIHEDRAL};
                $angleChange += 360.0 while $angleChange < -180.0;
                $angleChange -= 360.0 while $angleChange > 180.0;
                
                if ( exists $parentAtom->{FILLED_BOND3} ) {
                    $bondNumber = "bond4";
                    $parentAtom->{NEXT_CHILD_BOND} = 3;
                }
                elsif ( exists $parentAtom->{FILLED_BOND4} ) {
                    $bondNumber = "bond3";
                    $parentAtom->{NEXT_CHILD_BOND} = 4;
                }
                elsif ($angleChange > 0) {
                    $bondNumber = "bond3";
                    $parentAtom->{FILLED_BOND3} = 1;
                    $parentAtom->{NEXT_CHILD_BOND} = 4;
                }
                else {
                    $bondNumber = "bond4";
                    $parentAtom->{FILLED_BOND4} = 1;
                    $parentAtom->{NEXT_CHILD_BOND} = 3;
                }
            }
            else {
                $bondNumber = "bond" . $parentAtom->{NEXT_CHILD_BOND};
                ++$parentAtom->{NEXT_CHILD_BOND};
                $parentAtom->{FIRST_DIHEDRAL} = $atom->{DIHEDRAL_ANGLE};
            }
            print ", \"" . $parentAtom->{NAME} . "/$bondNumber\"";
            
            # Bond length
            print ", " . sprintf "%.5f", 0.10 * $atom->{BOND_LENGTH};
            
            # Set dihedral and mobility only if NBONDS > 1
            if ($atom->{NBONDS} > 1) {
                my $childAtom = $atom->{CHILD_ATOMS}[0];

                # Default mobility is rigid if:
                #    a) atom on either side of bond has only 1 bond
                # OR b) atoms on both sides of bond have exactly 3 bonds

                my $rigidBond = undef;
                
                # Case a)
                if ( (1 == $parentAtom->{NBONDS}) || (1 == $atom->{NBONDS}) ) {
                    $rigidBond = 1;
                }
                elsif ( (4 <= $parentAtom->{NBONDS}) || (4 <= $atom->{NBONDS}) ) {
                    $rigidBond = 0;
                }
                elsif ( (2 == $parentAtom->{NBONDS}) && (2 == $atom->{NBONDS}) ) {
                    $rigidBond = 0;
                }
                # 2 or 3 bonds on both atoms, but not both with 2 bonds
                else {
                    $rigidBond = 1;
                }
                
                # Dihedral angle of new bond
                
                my $dihedral_angle = $childAtom->{DIHEDRAL_ANGLE};
                
                print ", " . sprintf "%.2f*Deg2Rad", $dihedral_angle;
                
                # Mobility
                
                if ($rigidBond) {
                    print ", BondMobility::Rigid";
                }
                else {
                    print ", BondMobility::Torsion";
                }
            }
            
            print " );\n";
        }
        
        ++$atomCount;
    }
    
    print "\n";
    
    # Ring closing bonds
    foreach $atom (@$atoms) {
        foreach my $atom2 (@{$atom->{RING_BONDS}}) {
            print "        addRingClosingBond( ";
            print "\"". $atom->{NAME} . "/bond" . $atom->{NEXT_CHILD_BOND} . "\"";
            print ", \"". $atom2->{NAME} . "/bond" . $atom2->{NEXT_CHILD_BOND} . "\"";
            print ", 0.15"; # dummy length for ring closing bonds -- too hard to figure out explicitly
            print ");\n";
            
            $atom->{NEXT_CHILD_BOND} ++;
            $atom2->{NEXT_CHILD_BOND} ++;
        }
    }

    print "\n";

    # Biotypes
    foreach $atom (@$atoms) {
        next unless ref $atom;
        next unless exists $atom->{NAME};
        my $atomName = $atom->{NAME};
        print "        setBiotypeIndex( \"$atomName\", Biotype::get(\"$residueName\", \"$atomName\").getIndex());\n";
    }
    
    print "\n";
    
    # Explicitly set dihedral angles, because the bondAtom arguments are not always right
    
    foreach $atom (@$atoms) 
    {
        next unless ref $atom;
        next unless exists $atom->{PARENT};
        next unless exists $atom->{PARENT2};
        next unless exists $atom->{PARENT3};
        
        next unless $atom->{PARENT} > 0;
        next unless $atom->{PARENT2} > 0;
        next unless $atom->{PARENT3} > 0;
        
        next unless exists $atoms->[$atom->{PARENT}]{NAME};
        next unless exists $atoms->[$atom->{PARENT2}]{NAME};
        next unless exists $atoms->[$atom->{PARENT3}]{NAME};

        my $a1 = $atom->{NAME};
        my $a2 = $atoms->[$atom->{PARENT}]{NAME};
        my $a3 = $atoms->[$atom->{PARENT2}]{NAME};
        my $a4 = $atoms->[$atom->{PARENT3}]{NAME};
        
        next unless defined $a1;
        next unless defined $a2;
        next unless defined $a3;
        next unless defined $a4;        
        
        # set default dihedral angle
        print "        setDefaultDihedralAngle(" . $atom->{DIHEDRAL_ANGLE} . "*Deg2Rad";
        print ", \"$a1\", \"$a2\", \"$a3\", \"$a4\");\n";
        
    }
    
    print "    } // end constructor $residueName \n\n";

    print "    static void instantiateBiotypes() {\n";
    foreach $atom (@$atoms) {
        next unless ref $atom;
        next unless exists $atom->{NAME};
        my $atomName = $atom->{NAME};
        die unless defined $atomName;
        my $nbonds = $atom->{NBONDS};
        my $atomType = $atomTypes->{$atom->{TYPE}};
        my $atomicNumber = $atomType->{ATOMIC_NUMBER};
        my $elementName = elementName($atomicNumber);
        print "        if (! Biotype::exists(\"$residueName\", \"$atomName\") )\n";
        print "            Biotype::defineBiotype(Element::$elementName(), $nbonds, \"$residueName\", \"$atomName\"); \n";
    }
    print "    } // end instantiateBiotypes\n\n";

    print "}; // end class $residueName \n";

}

if (defined $CommandLineOptions{c}) 
{
    writeMolmodelCode(\@atoms, \%atomTypes, $residueName, $pdbResidueName);
}
else 
{
    writeTinkerParameters(\@atoms, \%atomTypes, $residueName);
}

