#!/usr/bin/perl

use strict;
use warnings;

use project;
use Scalar::Util qw(looks_like_number);
use Cwd;

# Read path for positive and negative sequences.
print "Give fasta file path:\n";
my $fastafile=<STDIN>;
print "Give fasta.tab file path:\n"
my $fastatabfile=<STDIN>;

# Read the choice of the user for a dictionary.
my $bool = 0;
my $Q;
my $Qnum;
while ($bool == 0){
	print "Give number Q: ";
	$Q = <STDIN>;
	$Qnum = $Q;
	if(!looks_like_number($Q)){
		print "Q must be a number. ";
	} elsif ($Q == 1) {
		$Q = 'Dictionaries/DNA.txt';
		$bool = 1;
	}elsif ($Q == 2) {
		$Q = 'Dictionaries/RNA.txt';
		$bool = 1;
	}elsif ($Q == 3) {
		$Q = 'Dictionaries/PROTEIN.txt';
		$bool = 1;
	}else {
		print "Give a number between 1-3. ";
	}
}

# Read the choice of the user for markov chain window size
my $K;
$bool=0;
while ($bool == 0){
	print "Give number K: ";
	$K = <STDIN>;
	if(!looks_like_number($K)){
		print "K must be a number. ";
	} elsif(index($K,'.') != -1){
		print "K must be an integer. ";
	}elsif ($K>0) {
		$bool =1;
	} else {
		print "K must be a positive number. ";
	}
}

# Read dictionary from file.
my @Q = @{project::readQ($Q)};

# Read AAtype character from user. The character must be one from the dictionary
$bool = 0;
my $AAt;
do{
	print "Give AAtype: ";
	$AAt = <STDIN>;
	chomp($AAt);
	foreach (@Q){
		my $tmp =  $_;
		if ($tmp eq $AAt){
			$bool=1;
		}
	}
}while ($bool!=1);

# Read user's choice for which tests to execute. The choice must be a number
# between 1-7.
my $test;
do{
	$bool = 0;
	print "\nChoose a test: \n";
	print "1. Self Consistency\n2. Cross Validation\n3. External Data\n4. Self Consistency & Cross Validation\n5. Cross Validation & External Data\n";
	print "6. Self Consistency & External Data\n7. All three tests\n\nAnswer: ";
	$test = <STDIN>;
	if(!looks_like_number($test)){
		print "Answer must be a number.\n";
	} elsif(index($test,'.') != -1){
		print "Answer must be an integer.\n";
	}
	else{
		$bool = 1;
	}
}while( (($test<1) || ($test>7)) || ($bool == 0) );

# Write users choices in a file.
my $cwd = getcwd();
open(OUT, ">>$cwd/RESULTS.txt") || die "Could not open file for writing";
print OUT "\n\n\nK:\t\t",$K,;
print OUT "Q:\t\t",$Q,"\n\n";
print OUT "AAtype:\t\t",$AAt,"\n\n";
print OUT "User's choice number:\t\t",$test,"\n\n";
close(OUT); # closing file handle

# Read external data files path from user.
my $efastafile;
my $efastatabfile;
# A flag for an easier external data selection check.
my $externalflag = 0;
if($test == 3 || $test == 5 || $test == 6 || $test == 7){
	# Read external positive and negative sequences.
	print "Give external fasta file path: ";
	$efastafile = <STDIN>;
	print "Give external fasta.tab file path: ";
	$efastatabfile = <STDIN>;
}

# Read N from the user for cross Validation.
my $N;
if($test == 2 || $test == 4 || $test == 5 || $test == 7){
	$bool = 0;
	while ($bool == 0){
		print "Give number N for cross Validation: ";
		$N = <STDIN>;
		if(!looks_like_number($N)){
			print "N must be a number. ";
		} elsif(index($N,'.') != -1){
			print "N must be an integer. ";
		}elsif ($N>0) {
			$bool =1;
		} else {
			print "N must be a positive number. ";
		}
	}
}

# Read data from fasta files and separate in positive and negative datasets.
my %sequences = %{&project::readFASTAThema2($fastafile)};
my @titles = @{&project::readFASTATabFile($fastatabfile)};

# Create positive and negative datasets.
my @posSeq = ();
my @negSeq = ();
foreach (@titles){
		my @tmp = split /\t/, $_;
		my $SID = shift @tmp;
		my $AAtype = shift @tmp;
		my $start = shift @tmp;
		my $stop = shift @tmp;
		my $score = shift @tmp;
		# If the AAtype of the file is the same with the user's option the part of
		# the sequence will be added in positive sequences. Else the entire sequence
		# will be added in the negative sequences.
		if($AAt eq $AAtype){
			if(exists $sequences{$SID}){
				if (length($sequences{$SID}) < $stop){
					print "Length = ",length($sequences{$SID})," Start = $start Stop = $stop\n";
				}
				my $tmpseq = substr $sequences{$SID}, $start, $stop;
				push @posSeq, $tmpseq;
			} else {
				print "$SID does not exist!\n";
			}
		}else{
			push @negSeq, $sequences{$SID};
		}
}

my @eposSeq = ();
my @enegSeq = ();
# If external data test is selected create positive and negative test.
if ($externalflag == 1){
	# Read data from fasta files and separate in positive and negative datasets.
	my %externalsequences = %{&project::readFASTAThema2($efastafile)};
	my @externaltitles = @{&project::readFASTATabFile($efastatabfile)};

	foreach (@externaltitles){
			my @tmp = split /\t/, $_;
			my $SID = shift @tmp;
			my $AAtype = shift @tmp;
			my $start = shift @tmp;
			my $stop = shift @tmp;
			my $score = shift @tmp;
			# If the AAtype of the file is the same with the user's option the part of
			# the sequence will be added in positive sequences. Else the entire sequence
			# will be added in the negative sequences.
			if($AAt eq $AAtype){
				# Chech if the sequence exists in the fasta file.
				if(exists $externalsequences{$SID}){
					# Check if the length of the sequence is at least the size of number
					# stop.
					if (length($externalsequences{$SID}) < $stop){
						print "Length = ",length($externalsequences{$SID})," Start = $start Stop = $stop\n";
					}
					# Add the part of the sequence we want in a temporary scalar variable.
					my $tmpseq = substr $externalsequences{$SID}, $start, $stop;
					push @eposSeq, $tmpseq;
				} else {
					print "$SID does not exist!\n";
				}
			}else{
				push @enegSeq, $externalsequences{$SID};
			}
	}

}

my %B = ();
# A flag for external dataset test.
my $flag =0;
# Run the program user chooced.
if ($test == 1) {
	&selfCon();
} elsif ($test == 2) {
	&crossVal();
} elsif ($test == 3) {
	&externalTest();
} elsif ($test == 4) {
	&selfCon();
	&crossVal();
} elsif ($test == 5) {
   &crossVal();
   &externalTest();
} elsif ($test == 6) {
	$flag =1;
	&selfCon();
	&externalTest();
} elsif ($test == 7) {
	$flag=1;
	&selfCon();
	&externalTest();
	&crossVal();
}

=item selfCon()
	This function calculates B and then call function selfConsistency. From this
	function the accurancy is returned and is printed on the terminal and writen
	on a file.
=cut
sub selfCon {
# Input: void
# Returns: void
	open(OUT, ">>$cwd/RESULTS.txt") || die "Could not open file for writing";
	print OUT "Self consistency Test:\n";
	close(OUT); # closing file handle
	# Calculate B and the accurancy.
	%B = %{&project::calculateB(\@posSeq, \@negSeq, \@Q, $K, $Qnum)};
	my $acc = &project::selfConsistency(\%B, \@posSeq, \@negSeq, $K);

	# Write results on the terminal and in a file.
	print "Self consistency accurancy = ",$acc,"\n";
	open(OUT, ">>$cwd/RESULTS.txt") || die "Could not open file for writing";
	print OUT "Accurancy:\t", $acc,"\n\n";
	close(OUT); # closing file handle
	return;
}

=item crossVal()
	This function calls the function that calculates the cross Validation from
	the project package. From this function the average accuracy is retuned from
	all cross validations and is printed on the terminal.
=cut
sub crossVal() {
# Input: void
# Returns: void
	open(OUT, ">>$cwd/RESULTS.txt") || die "Could not open file for writing";
	print OUT "Cross Validation:\n";
	print OUT "N:\t\t",$N,"\n";
	close(OUT); # closing file handle
	# Execute cross Validation and get the Overall accuracy.
	my $averageacc = &project::crossValidation($N, \@posSeq, \@negSeq, \@Q, $K, $Qnum);
	# Print the average accuracy on the terminal.
	print ("Average accuracy =",$averageacc,"\n\n" );
	return;
}

=item externalTest()
	This function checks if B is calculated. If not calculateB function is
	executed with the teaching sequences. With the B calculated now our program
	can find the accurancy of the external sequences. The accurancy is returned to
	this function and is printed on the terminal and writen in a file.
=cut
sub externalTest() {
# Input: void
# Returns: void
	open(OUT, ">>$cwd/RESULTS.txt") || die "Could not open file for writing";
	print OUT "External Data Test:\n";
	close(OUT); # closing file handle

	# If a self consistency test is not done calculate B. Else the B is already
	# calculated.
	if($flag==0){
		print "Calculating B.\n";
		%B = %{&project::calculateB(\@posSeq, \@negSeq, \@Q, $K, $Qnum)};
	}
	# Results for the external dataset.
	my $ExAcc = &project::selfConsistency(\%B, \@eposSeq, \@enegSeq, $K);

	# Write results on terminal and in file.
	print "External accuracy = ",$ExAcc,"\n";
	open(OUT, ">>$cwd/RESULTS.txt") || die "Could not open file for writing";
	print OUT "Accuracy:\t", $ExAcc,"\n\n";
	close(OUT); # closing file handle
	return;
}
