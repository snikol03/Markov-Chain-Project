#!/usr/bin/perl

package project;

use strict;
use warnings;

sub concatbase{
	my @bases = @{shift @_};
	my $K = shift @_;
	my @concatbases;
	my @concatbases2;
	push @concatbases, @bases;
	for(my $i=0; $i<$K; $i++){
		foreach(@concatbases){
			my $conbase = $_;
			for (my $j=0; $j<scalar(@bases);$j++){
				my $base = $conbase.$bases[$j];
				push @concatbases2, $base;
			}
		}
		while(scalar(@concatbases)!= 0){
			pop @concatbases;
		}
		push @concatbases, @concatbases2;
		while(scalar(@concatbases2)!= 0){
			pop @concatbases2;
		}
	}
	return 	@concatbases;
}

sub estimateA{
# Input: A reference to an array of sequences, a reference to an array (library - Q), a scalar K for window size
# Returns: A hashtable: key-> transmition(AA, AT, AC,...), value-> probability of the transmition
	my $refseq = shift @_;
	my %P = ();
	my %count = ();
	my %basecount = ();
	my @bases = @{shift @_};
	my $K = shift @_;
	
	my @concatbases = &concatbase(\@bases, $K);
	#my @bases = ('A','T','C','G');
	# for (my $i=0; $i<scalar(@bases);$i++){
		# for (my $j=0; $j<scalar(@bases);$j++){
			# $count{$bases[$i].$bases[$j]}=0.000001;
			# #$P{$bases[$i].$bases[$j]}=0;
		# }
	# }
	# push @concatbases, @bases;
	# for(my $i=0; $i<$K; $i++){
		# foreach(@concatbases){
			# my $conbase = $_;
			# for (my $j=0; $j<scalar(@bases);$j++){
				# my $base = $conbase.$bases[$j];
				# push @concatbases2, $base;
			# }
		# }
		# while(scalar(@concatbases)!= 0){
			# pop @concatbases;
		# }
		# push @concatbases, @concatbases2;
		# while(scalar(@concatbases2)!= 0){
			# pop @concatbases2;
		# }
	# }
	
	foreach(@concatbases){
		print $_,"\n";
		$count{$_} = 0.00000000001;
	}
	for(my $i=0;$i<scalar (@{$refseq});$i++){
		my $seq = $refseq->[$i];
		for(my $j=0;$j<length($seq)-$K;$j++){
			my $transition = substr($seq,$j,$K+1);
			if(exists($count{$transition})){
				$count{$transition}++;
			}else{
				die "Illegal transition $transition fount in sequence ", $seq,"\n";
			}
		}

	}

	foreach my $tr (keys %count){
		my $start = substr($tr,0,1);
		my $denom = 0;
		my @newbases = &concatbase(\@bases, $K-1);
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

sub estimateP{
# Input: A reference to an array of sequences
# Returns: A hash table: key->base(A,T,C,G), value->probability of a sequence starting with each base

	my $refseq = shift @_;
	my @bases = @{shift @_};
	
	my %P = ();
	my %count = ();
	foreach (@bases){
		$count{$_} = 0;
	}
	#$count{'A'}=$count{'T'}=$count{'C'}=$count{'G'}=0;

	for(my $i=0;$i<scalar(@{$refseq});$i++){
		my $base = substr($refseq->[$i], 0, 1);
		if (exists($count{$base})){
			$count{$base}++;
		}else{
			die "Illegal character $base found in sequence ",$refseq->[$i],"\n";
		}
	}
	foreach(keys %count){
		$P{$_} = $count{$_}/scalar @{$refseq};
	}
	return(%P);
}

sub readQ{
	# Input: a file name
	# Returns: a reference to a list of qens
	my $Qpath = shift @_;
	open(IN, $Qpath) || die "Could not open file name $Qpath";
	# Read data line of file. (First line.)
	# my $data = <IN>;
	# print "Data: ",$data,"\n";
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
	# while (my $row = <$fh>) {
		# chomp $row;
		# print "$row\n";
	# }
	close(IN);
	return(\@seq);
}

sub readFASTA{
# Input: a file name
# Returns: a reference to a list of sequences
	my $infile = shift @_;
	$/="\n>"; #Change input record separator to read multiple FASTA sequences
	open(IN, $infile) || die "Could not open file name $infile";
	# Read data line of file. (First line.)
	#my $data = <IN>;
	my @seq = ();
	while(<IN>){
		chomp($_);
		my @tmp = split /\n/, $_;
		shift @tmp;
		my $seq = join '', @tmp;
		#print $seq, "\n\n";
		push @seq, $seq;
		#$seq .= $_;
	}
	close(IN);
	return(\@seq);
}

sub createrandomHash{
	#Input: An array with sequences
	#Output: A hashtable with the sequences in random sets
	my $N = shift @_;
	my @Seq = @{shift @_};
	my %Hash=();
	foreach (@Seq){
		my $gen = $_;
		my $r = rand($N);
		for(my $i=1;$i<=$N;$i++){
			if($r<=$i && $r>$i-1){
				push @{$Hash{$i}},$gen;
			}
		}
	}
	return(\%Hash);
}

sub calculateB{
# Input: an array with the possitive teaching sequences, an array with the negative teaching sequences
# Return: a reference to a hash with the values of B
	my @posteachSeq = @{shift @_};
	my @negteachSeq = @{shift @_};
	my @Q = @{shift @_};
	my $K = shift @_;
	my %B = ();
	my %Apos = &estimateA(\@posteachSeq, \@Q, $K);
	my %Aneg = &estimateA(\@negteachSeq, \@Q, $K);
	foreach (sort keys %Apos){
		$B{$_} = log($Apos{$_}/$Aneg{$_});
		#print OUT $_,"\t",$Apos{$_},"\t",$Aneg{$_},"\t",$B{$_},"\n";
	}
	return (\%B);
}
sub selfConsistency{
# Input: A hash reference with the values of B, a reference to a positive sequence array, a reference to an negative sequence array
# Output: 	
	my %B = %{shift @_};
	my @posSeq = @{shift @_};
	my @negSeq = @{shift @_};
	my $K = shift @_;
	my @selfCon = ('0','0','0','0');
	for(my $j=0;$j<scalar(@posSeq);$j++){
		my $query = $posSeq[$j];
		my $score = 0;
		for(my $c=0;$c<length($query)-$K;$c++){
			$score += $B{substr($query, $c,$K+1)};
		}
		if($score>0){
			$selfCon[0]++;
		}else{
			$selfCon[1]++;
		}
	}
	for(my $j=0;$j<scalar(@negSeq);$j++){
		my $query = $negSeq[$j];
		my $score = 0;
		for(my $c=0;$c<length($query)-$K;$c++){
			$score += $B{substr($query, $c,$K+1)};
		}
		if($score<0){
			$selfCon[3]++;
		}else{
			$selfCon[2]++;
		}
	}
	foreach(@selfCon){
		print $_,"\t";
	}
	print "\n";
	# print scalar (@posSeq), "   ", scalar (@negSeq),"\n";
	my $acc = (($selfCon[0]+$selfCon[3])/($selfCon[0]+$selfCon[3]+$selfCon[1]+$selfCon[2]))*100;
	return ($acc);
}

sub calculateAccurancy{
	#Input: A reference to a Hashtable of possitive sequences, A reference to a Hashtable of negative sequences 
	#Output: Accurancy
	my %posHash = %{shift @_};
	my %negHash = %{shift @_};
	my @Q = @{shift @_};
	my $K = shift @_;
	my %ret = ();
	my $N = scalar( keys %posHash);
	my $sumacc =0;
	for(my $i=1;$i<=$N;$i++){
		my @posteachSeq =();
		my @negteachSeq =();
		for(my $j=1;$j<=$N;$j++){
			if($i!=$j){
				push @posteachSeq,@{$posHash{$j}};
				push @negteachSeq,@{$negHash{$j}};
			}
		}
		my $Bref = &calculateB(\@posteachSeq, \@negteachSeq, \@Q, $K);
		my $acc = &selfConsistency($Bref, \@posteachSeq, \@negteachSeq, $K);
		$ret{$i} = $acc;
		$sumacc += $acc;
		print "Overall Accurancy ",$i," = ", $acc,"\n";
	}
	return (%ret);
}

sub crossValidation{
#Input: the number of groups, Two references to arrays (positive sequences and negatice sequences)
#Output: Average accurancy (scalar)
	my $N = shift @_;
	my @posSeq = @{shift @_};
	my @negSeq = @{shift @_};
	my @Q = @{shift @_};
	my $K = shift @_;
	
	my %posHash = %{&createrandomHash($N,\@posSeq)};
	my %negHash = %{&createrandomHash($N,\@negSeq)};
	
	my %accurancy = &calculateAccurancy( \%posHash, \%negHash, \@Q, $K);
	my $sumacc = 0;
	foreach(keys %accurancy){
		$sumacc += $accurancy{$_};
	}
	my $averageacc = $sumacc/$N;
	return ($averageacc);

}


1;
