#!/usr/bin/perl

package project;

use strict;
use warnings;
use Cwd;
=pod

=item concatbase()
	This function creates all posible combinations of the bases from the
	array bases. Depending on the choice of the user the program select the
  character which have to be excluded.
=cut
sub concatbase{
# Input: a reference to an array (Library - bases), a scalar K for window size, a scalar flag
#        for the bases the user choosed.
# Return: an array of all the combinations of the bases.
	my @bases = @{shift @_};
	my $K = shift @_;
	my $flag = shift @_;
	# Select the character which will be excluded depending on user's choice.
	if ($flag != 0) {
		my $excludeChar;
		if ($flag == 3){
			$excludeChar = "X";
		}else {
			$excludeChar = "N";
		}
		for (my $i=0;$i<scalar(@bases);$i++) {
			if ($bases[$i] eq $excludeChar){
				splice(@bases,$i, 1);
			}
		}
	}

	# Concat K bases, all posible combinations will be created.
	my @concatbases;
	my @concatbases2;
	# Add bases in an array.
	push @concatbases, @bases;
	for(my $i=0; $i<$K; $i++){
		foreach(@concatbases){
			my $conbase = $_;
			for (my $j=0; $j<scalar(@bases);$j++){
				# On every existing combination of bases concatanate all bases from
				# array bases.
				my $base = $conbase.$bases[$j];
				# Add the new base in an array.
				push @concatbases2, $base;
			}
		}
		# Remove everything from this array.
		while(scalar(@concatbases)!= 0){
			pop @concatbases;
		}
		# Add everything from the second array to the first one.
		push @concatbases, @concatbases2;
		# Remove everything from this array.
		while(scalar(@concatbases2)!= 0){
			pop @concatbases2;
		}
	}
	return 	@concatbases;
}

=item estimateA()
   This function calculates the probability of each transicion from all the
  sequences.
=cut
sub estimateA{
# Input: A reference to an array of sequences, a reference to an array of the
#        dictionary, a scalar K for window size, a scalar number for the
#				 balander character.
# Returns: A hashtable: key-> transicion(AA, AT, AC,...), value-> probability
#					 of the transicion
	my $refseq = shift @_;
	my @bases = @{shift @_};
	my $K = shift @_;
	my $flag = shift @_;

	my %P = ();
	my %count = ();
	my %basecount = ();
	# Get all posible transicions.
	my @concatbases = &concatbase(\@bases, $K, 0);

	# Set count of each transicion to a very small number close to 0.
	foreach(@concatbases){
		$count{$_} = 0.00000000001;
	}
	# Count transicions. If an illegal transicion is found exit program.
	for(my $i=0;$i<scalar(@{$refseq});$i++){
		my $seq = $refseq->[$i];
		for(my $j=0;$j<length($seq)-$K;$j++){
			# Get the transition.
			my $transition = substr($seq,$j,$K+1);
			# Increase the count if a transition is found. Else stop the program.
			if(exists($count{$transition})){
				$count{$transition}++;
			}else {
				die "Illegal transition $transition found in sequence $seq";
			}
		}
	}
	# Set the balander character depending on user's choice of dictionary.
	my $tmpchar;
	if ($flag == 3) {
		$tmpchar = 'X';
	}
	else {
		$tmpchar = 'N';
	}

	# Remove all the keys which include balander character
	foreach my $key (sort keys %count) {
		for (my $i = 0; $i<($K+1);$i++){
			my $letter = substr($key, $i, 1);
			if ($letter eq $tmpchar){
				delete $count{$key};
				# print "deleted key: ",$key,"\n";
			}
		}
	}
	# Add each probability in a hash table with key -> transicion and value the
	# probability.
	foreach my $tr (keys %count){
		my $start = substr($tr,0,1);
		my $denom = 0;
		my @newbases = &concatbase(\@bases, $K-1, $flag);
		foreach my $end (@newbases){
			my $tr2 = $start.$end;
			$denom += $count{$tr2};
		}
		if($denom >0){
			$P{$tr} = $count{$tr}/$denom;
		}else{
			die "Zero denominator for transition $tr";
		}
	}

	return (%P);
}

=item estimateP()
	This function calculates the probability of the appearance of each character at the begining
	of the sequence.
=cut
sub estimateP{
# Input: A reference to an array of sequences, a reference to an array of the
#        dictionary.
# Returns: A hash table: key->base(A,T,C,G), value->probability of a sequence starting with each base

	my $refseq = shift @_;
	my @bases = @{shift @_};

	my %P = ();
	my %count = ();

	# Initialize the hash with the counts.
	foreach (@bases){
		$count{$_} = 0;
	}

	# Count bases. If an illegal character is found exit the program.
	for(my $i=0;$i<scalar(@{$refseq});$i++){
		my $base = substr($refseq->[$i], 0, 1);
		if (exists($count{$base})){
			$count{$base}++;
		}else{
			die "Illegal character $base found in sequence ",$refseq->[$i],"\n";
		}
	}

	# Add probability of each base to a hashtable.
	foreach(keys %count){
		$P{$_} = $count{$_}/scalar @{$refseq};
	}
	return(%P);
}


=item readQ()
 This function read the dictionary of bases from a file.
=cut
sub readQ{
# Input: a file name
# Returns: a reference to an array of qenes
	my $Qpath = shift @_;
	open(IN, $Qpath) || die "Could not open file name $Qpath";

	my @seq = ();
	my $i=1;
	while(<IN>){
		chomp($_);
		# Check if q.txt has correct data input
		# if(length($_)>2){
			# die "$_ is not a gen";
		# }
		push @seq, $_;
	}
	close(IN);

	return(\@seq);
}


=item readFASTATabFile()
 This function read and return the data from a fasta.tab file.
=cut
sub readFASTATabFile{
# Input: a scalar with the filename.
# Returns: a reference to an array with the read data.

	my $infile = shift @_;
	open(IN, $infile) || die "Could not open file name $infile";
	# Read data line of file. (First line.)
	my $data = <IN>;
	# Read all the other lines and push them to an array.
	my @seq = ();
	while(<IN>){
		chomp($_);
		push @seq, $_;
	}
	close(IN);
	return(\@seq);
}

=item readFASTAThema2()
 This function is a more specific reading function from a fasta file for thema 2.
  This function read from a fasta file and remove the stop codon (*) at the end of
  each sequence.
=cut
sub readFASTAThema2{
# Input: a file name
# Returns: a reference to a list of sequences

	my $infile = shift @_;
	$/="\n>"; #Change input record separator to read multiple FASTA sequences
	open(IN, $infile) || die "Could not open file name $infile";

	my %seq = ();
	my $count = 0;
	# Read data from file.
	while(<IN>){
		chomp($_);
		my @tmp = split /\n/, $_;
		my $key = "";
		# Add ">" which was removed except the first line.
		if($count > 0 ){
			$key = ">";
		}
		$key .= shift @tmp;
		my $seq = join '', @tmp;
    my $flag=0;
    for(my $i=0;$i<length($seq)-1;$i++){
        my $checkChar=substr($seq,$i,1);
        if ($checkChar eq "*"){
          $flag=1;
        }
    }
    if ($flag==0){
      # remove stop codon * at the end of a sequence.
      if(substr ($seq, -1) eq '*'){
  		 	$seq = substr($seq, 0, length($seq)-1);
  		}
      $seq{$key} = uc($seq);
    }

		$count++;
	}
	close(IN);
	$/="\n"; #Change input record separator to read multiple FASTA sequences
	return(\%seq);
}

=item readFASTA()
  This function read and return the data from a fasta file.
=cut
sub readFASTA{
# Input: a file name
# Returns: a reference to a list of sequences
	my $infile = shift @_;
	$/="\n>"; #Change input record separator to read multiple FASTA sequences
	open(IN, $infile) || die "Could not open file name $infile.";

	my @seq = ();
	# Read data from file.
	while(<IN>){
		chomp($_);
		my @tmp = split /\n/, $_;
		shift @tmp;
		my $seq = join '', @tmp;
    # Add sequen to an array
		push @seq, uc($seq);
	}
	close(IN);
	$/="\n"; #Change input record separator to read multiple FASTA sequences
	return(\@seq);
}

=item createrandomHash()
  This function separate randomly in N (parameter given by user) parts the sequences.
=cut
sub createrandomHash{
# Input: An array with sequences.
# Returns: A hashtable with the sequences in random sets.
	my $N = shift @_;
	my @Seq = @{shift @_};
	my %Hash=();
	foreach (@Seq){
		my $gen = $_;
		my $r = rand($N);
		# Find the correct key of hash and push the sequence.
		for(my $i=1;$i<=$N;$i++){
			if($r<=$i && $r>$i-1){
				push @{$Hash{$i}},$gen;
			}
		}
	}
	return(\%Hash);
}

=item calculateB()
  This function calculates the B values and returns a hash table with this values.
=cut
sub calculateB{
# Input: a reference to an array with the possitive teaching sequences, a
#        reference to an array with the negative teaching sequences, a reference
#        to an array with the dictionary, a scalar number of the markov chain
#        window size and a scalar number for the choice of the dictionary.
# Return: a reference to a hashtable with keys -> bases, values -> score.
	my @posteachSeq = @{shift @_};
	my @negteachSeq = @{shift @_};
	my @Q = @{shift @_};
	my $K = shift @_;
	my $Qnum = shift @_;
	my %B = ();
	my %Apos = &estimateA(\@posteachSeq, \@Q, $K, $Qnum);
	my %Aneg = &estimateA(\@negteachSeq, \@Q, $K, $Qnum);
	foreach (sort keys %Apos){
		$B{$_} = log($Apos{$_}/$Aneg{$_});
		#print OUT $_,"\t",$Apos{$_},"\t",$Aneg{$_},"\t",$B{$_},"\n";
	}
	return (\%B);
}

=item selfConsistency()
  This function calculates the True Positive, False Negative, False Positive
  and True Negative for the positive and negative sequences .
=cut
sub selfConsistency{
# Input: A hash reference with the values of B, a reference to a positive sequence array,
# a reference to an negative sequence array and a scalar number of the markov chain window
#
# Returns: A scalar number of the accurancy.
	my %B = %{shift @_};
	my @posSeq = @{shift @_};
	my @negSeq = @{shift @_};
	my $K = shift @_;
	# TP, FN, FP, TN
	my @selfCon = ('0','0','0','0');

	my $tmpChar;

print "\nPositive Sequence!\n\n";
	for(my $j=0;$j<scalar(@posSeq);$j++){
		my $query = $posSeq[$j];
		my $score = 0;

		my $flag = 0;
		for(my $c=0;$c<length($query)-$K;$c++){

			$tmpChar = substr($query, $c,$K+1);

			foreach my $key (sort keys %B){
				if ($key eq $tmpChar){
					$flag=1;
					last; #break the loop
					print "Broken: $tmpChar\n\n";
				}

			}
			if ($flag == 1) {
				$score= $score+$B{$tmpChar};
				$flag = 0;
			}
			# }else {
			# 	print "Key not found: $tmpChar!!!\n";
			# }
			# $score += $B{substr($query, $c,$K+1)};
		}
		my $normilizedscore = $score /length($query);
		if($score>0){
			$selfCon[0]++;
		}else{
			$selfCon[1]++;
		}
	}

print "Negative Sequence!\n\n";

	for(my $j=0;$j<scalar(@negSeq);$j++){
		my $query = $negSeq[$j];
		my $score = 0;
		my $flag = 0;
		for(my $c=0;$c<length($query)-$K;$c++){

			$tmpChar = substr($query, $c,$K+1);

			foreach my $key (sort keys %B){
				if ($key eq $tmpChar){
					$flag=1;
					last; #break the loop
					print "Broken: $tmpChar\n\n";
				}

			}
			if ($flag == 1) {
				$score= $score+$B{$tmpChar};
				$flag = 0;
			}
			# }else {
			# 	print "Key not found: $tmpChar!!!\n";
			# }

		}
		my $normilizedscore = $score / length($query);
		if($score<0){
			$selfCon[3]++;
		}else{
			$selfCon[2]++;
		}
	}
	print "\nTP\tFN\tFP\tTN\n";
	foreach(@selfCon){
		print $_,"\t";
	}
	print "\n";

	my @ROCheadings = ("TP","FN","FP","TN");
	my $cwd = getcwd();
	open(OUT, ">>$cwd/RESULTS.txt") || die "Could not open file for writing";
	for(my $i = 0; $i<scalar(@selfCon);$i++){
		print OUT $ROCheadings[$i],":\t\t", $selfCon[$i],"\n";
	}
	close(OUT); # closing file handle
	# print scalar (@posSeq), "   ", scalar (@negSeq),"\n";
	my $acc = (($selfCon[0]+$selfCon[3])/($selfCon[0]+$selfCon[3]+$selfCon[1]+$selfCon[2]))*100;
	return ($acc);
}

=item selfConsistency()
 This function calculates the self  consistenct accurancy for cross validation.
=cut
sub calculateAccuracy{
# Input: A reference to a Hashtable of possitive sequences, A reference to a Hashtable of negative sequences,
#        A reference to an array with the dictionary, a scalar number K and a scalar number Qnum.
# Returns: A hashtable with the accurancy of each cross validation.
	my %posHash = %{shift @_};
	my %negHash = %{shift @_};
	my @Q = @{shift @_};
	my $K = shift @_;
	my $Qnum = shift @_;
	my %ret = ();
	my $N = scalar( keys %posHash);
	my $cwd = getcwd();
	for(my $i=1;$i<=$N;$i++){

		open(OUT, ">>$cwd/RESULTS.txt") || die "Could not open file for writing";
		print OUT "Ni:\t\t", $i,"\n";
		close(OUT); # closing file handle

		# Create positive and negatice data sets.
		my @posteachSeq =();
		my @negteachSeq =();
		for(my $j=1;$j<=$N;$j++){
			if($i!=$j){
				push @posteachSeq,@{$posHash{$j}};
				push @negteachSeq,@{$negHash{$j}};
			}
		}

		my $Bref = &calculateB(\@posteachSeq, \@negteachSeq, \@Q, $K, $Qnum);
		#my $acc = &selfConsistency($Bref, \@posteachSeq, \@negteachSeq, $K);
    my $acc = &selfConsistency($Bref, \@{$posHash{$i}}, \@{$negHash{$i}}, $K);
		$ret{$i} = $acc;

		# Write results on terminal and in file.
		print "Accurancy ",$i," = ", $acc,"\n";
		open(OUT, ">>$cwd/RESULTS.txt") || die "Could not open file for writing";
		print OUT "Accuracy:\t", $acc,"\n\n";
		close(OUT); # closing file handle
	}

	return (%ret);
}

=item crossValidation()
  This function randomly separates the dataset in N parts and calculates the
  accuracy of each one and the average accuracy of all dataset.
=cut
sub crossValidation{
# Input: The number of groups, Two references to arrays (positive sequences and negatice sequences),
#		 a reference to an array of the dictionary with the bases, a scalar number of the window size
#		 for markov chain and a scalar number of the choice of dictionary.
# Returns: A scalar with the average accurancy.
	my $N = shift @_;
	my @posSeq = @{shift @_};
	my @negSeq = @{shift @_};
	my @Q = @{shift @_};
	my $K = shift @_;
	my $Qnum = shift @_;
  # Randomly separate the positive and negative sequences.
	my %posHash = %{&createrandomHash($N,\@posSeq)};
	my %negHash = %{&createrandomHash($N,\@negSeq)};
  # Calculate the accurancy of each dataset.
	my %accurancy = &calculateAccuracy( \%posHash, \%negHash, \@Q, $K, $Qnum);
  # Calculate the average accurancy.
  my $sumacc = 0;
	foreach(keys %accurancy){
		$sumacc += $accurancy{$_};
	}
	my $averageacc = $sumacc/$N;
  # Write the results in a file.
	my $cwd = getcwd();
	open(OUT, ">>$cwd/RESULTS.txt") || die "Could not open file for writing";
	print OUT "Overall Accuracy:\t", $averageacc,"\n\n";
	close(OUT); # closing file handle
  # Return average accurancy.
	return ($averageacc);

}

1;
