#!/usr/bin/perl -w

use strict;

sub runCmd {
	my ($cmd) = @_;
	print STDERR sprintf("Running cmd:\t%s\n", $cmd);
	open(IN, sprintf("%s | ", $cmd)) or die sprintf("FAILED:\t%s\n", $cmd) or return 1;
	while(defined(my $line = <IN>)) {
		chomp $line;
	}
	return close(IN);
}

sub runCmDie {
	my ($cmd) = @_;
	
	unless(&runCmd($cmd)) {
		die sprintf("Failed command:\t%s\n", $cmd);
	}
}

sub trimAsm {
	
	my ($infilename, $outfilename) = @_;
	my $print = 0;
	
	open(IN, "<", $infilename) or die "Filename:\t$infilename\n$!\n";
	open(OUT, ">", $outfilename) or die "Filename:\t$outfilename\n$!\n";
	
	while(defined(my $line = <IN>)) {
		if($line =~	m/^>NODE_\d+_length_(\d+)_cov_(\d+\.\d+)/o) {
			$print = 512 < $1 && 1.0 < $2;
		}
		if($print) {
			print OUT $line;
		}
	}
	
	close IN;
	close OUT;
}

sub matchAsm {
	
	my ($infilename, $outfilename) = @_;
	my %uniqMatch = ();
	
	open(IN, sprintf("gunzip -c %s | cut -f 1,7 | uniq | ", $infilename)) or die "Filename:\t$infilename\n$!\n";
	open(OUT, ">", $outfilename) or die "Filename:\t$outfilename\n$!\n";
	
	while(defined(my $line = <IN>)) {
		chomp $line;
		my @fields = split('\t', $line);
		unless(exists($uniqMatch{$fields[1]})) {
			$uniqMatch{$fields[1]} = 1;
			print OUT sprintf(">%s\n%s\n", $fields[1], $fields[0]);
		}
	}
	%uniqMatch = ();
	close IN;
	close OUT;
}

sub checkTarget {
	
	my ($filename) = @_;
	my $empty = 1;
	open(IN, "<", $filename) or return $empty;#die sprintf("Filename: %s\n%s\n", $filename, $!);
	if(defined(my $line = <IN>)) {
		$empty = 0;
	}
	close IN;
	
	return $empty;
}

sub runIteration {
	
	my ($output, $target, $db, $kma_out, $target_reads, $SPAdes_out, @reads) = @_;
	my $se = scalar(@reads);
	
	######################
	## get target reads ##
	######################
	
	# Index target
	&runCmDie(sprintf("kma index -i %s -o %s.target -m 14", $target, $db));
	
	# align reads to target
	&runCmDie(sprintf("kma -i %s -o %s -tmp -t_db %s.target -mem_mode -1t1 -mrc 0.25 -ID 50 -na -nc -ef -t 4", join(" ", @reads), $kma_out, $db));
	
	# extract aligning reads
	&runCmDie(sprintf("gunzip -c %s.frag.gz | cut -f 7 | cut -f 1 -d \" \" | sort | uniq > %s", $kma_out, $target_reads));
	my $filename = sprintf("%s.target.reads", $output);
	my $fqgrep_cmd = sprintf("fqgrep -f %s -o %s ", $target_reads, $filename);
	if($se == 1) {
		$fqgrep_cmd .= sprintf(" -i %s", $reads[0]);
	} elsif($se == 2) {
		$fqgrep_cmd .= sprintf(" -p %s %s", $reads[0], $reads[1]);
	} else {
		$fqgrep_cmd .= sprintf(" -p %s %s -i %s", $reads[0], $reads[1], $reads[2]);
	}
	if(&checkTarget($target_reads) == 0) {
		&runCmDie($fqgrep_cmd);
	} else {
		$filename = sprintf("%s.fasta", $db);
		open(IN, ">", $filename) or die sprintf("Filename: %s\n%s\n", $filename, $!);
		close(IN);
		return $filename;
	}
	
	##################################
	## SPAdes assemble target reads ##
	##################################
	
	# run SPAdes
	#my $SPAdes_cmd = sprintf("spades.py -o %s ", $SPAdes_out);
	#my $SPAdes_cmd = sprintf("spades.py --only-assembler -o %s ", $SPAdes_out);
	my $SPAdes_cmd = sprintf("spades.py --only-assembler -k 21,33,55,87 -o %s ", $SPAdes_out);
	#my $SPAdes_cmd = sprintf("spades.py --only-assembler -k 27,47,67,87,107,127 -t 39 -m 188 -o %s ", $SPAdes_out);
	my $valid_files = 0;
	$valid_files++ if(&checkTarget(sprintf("%s.fq", $filename)) == 0);
	$valid_files += 2 if(&checkTarget(sprintf("%s%s.fq", $filename, "_1")) == 0);
	if($valid_files == 1) {
		$SPAdes_cmd .= sprintf("--sc -s %s.fq", $filename);
	} elsif($valid_files == 2) {
		$SPAdes_cmd .= sprintf("--meta -1 %s%s.fq -2 %s%s.fq", $filename, "_1", $filename, "_2");
	} elsif($valid_files == 3) {
		$SPAdes_cmd .= sprintf("--meta -1 %s%s.fq -2 %s%s.fq -s %s.fq", $filename, "_1", $filename, "_2", $filename);
	}
	if($valid_files != 0) {
		&runCmDie($SPAdes_cmd);
	}
	
	# Get contigs with original matches
	if(&checkTarget(sprintf("%s/scaffolds.fasta", $SPAdes_out)) == 0) {
		&runCmDie(sprintf("kma -i %s -o %s -tmp -t_db %s -ID 50 -na -nc -a -proxi -0.7", sprintf("%s/scaffolds.fasta", $SPAdes_out), $kma_out, $db));
		&matchAsm(sprintf("%s.frag.gz", $kma_out), sprintf("%s.fasta", $db));
	} else {
		$filename = sprintf("%s.fasta", $db);
		open(IN, ">", $filename) or die sprintf("Filename: %s\n%s\n", $filename, $!);
		close(IN);
	}
	
	# Set assembly as new target
	return $target = sprintf("%s.fasta", $db);
}

sub countMatches {
	
	my ($mapstat) = @_;
	
	open(IN, "<", $mapstat) or die sprintf("Filename: %s\n%s\n", $mapstat, $!);
	my $matches = 0;
	while(defined(my $line = <IN>)) {
		if(!($line =~ m/^#/o)) {
			my @fields = split('\t', $line);
			$matches += $fields[13] + 0;
		}
	}
	close IN;
	
	return $matches;
}

sub main {
	
	my ($argc, @argv) = @_;
	
	# check input
	if($argc < 4 || 6 < $argc) {
		print STDERR "targetAsm.pl creates assemblies around target sequences.\n";
		print STDERR "Requirements:\tKMA\tSPAdes\n";
		print STDERR "Usage:\n";
		print STDERR "targetAsm.pl iter output target read(s)\n";
		print STDERR "iter:\tNumber of target assembly iterations\n";
		print STDERR "output:\tOutput destination\n";
		print STDERR "target:\tTarget fasta\n";
		print STDERR "read(s):\tInput read(s), paired- or single-end, ordered as: 1,2,s (only 2nd generation)\n";
		return 1;
	}
	
	# get input
	my ($iter, $output, $target, @reads) = @ARGV;
	$iter = $iter + 0;
	
	# init variables
	my ($matches, $prev_matches, $reAsm, $iternumber) = (0, 0, 0, 0);
	my ($db, $kma_out, $tmpoutput, $target_reads, $SPAdes_out, $SPAdes_cmd, $mapstat) = ("", "", "", "", "", "", "");
	my @targetReads = ();
	
	# set variables
	&runCmDie(sprintf("mkdir -p %s.tmp/", $output));
	$tmpoutput = sprintf("%s.tmp/tmp", $output);
	$db = sprintf("%s.target.index", $tmpoutput);
	$kma_out = sprintf("%s.target.kma", $tmpoutput);
	$target_reads = sprintf("%s.target.reads", $tmpoutput);
	$SPAdes_out = sprintf("%s.target.SPAdes", $tmpoutput);
	$mapstat = sprintf("%s.mapstat", $kma_out);
	
	# Index target
	&runCmDie(sprintf("kma index -i %s -o %s -m 14", $target, $db));
	
	# Get reads hitting the target
	$target = &runIteration($tmpoutput, $target, $db, $kma_out, $target_reads, $SPAdes_out, @reads);
	
	# Run iterative assembly
	$reAsm = $iter < 0 ? 1 : $iter;
	# check if target is empty
	if(&checkTarget($target)) {
		$reAsm = 0;
	}
	while($reAsm) {
		# Run iteration assembly
		$target = &runIteration($tmpoutput, $target, $db, $kma_out, $target_reads, $SPAdes_out, @reads);
		
		#######################
		## validate assembly ##
		#######################
		
		# check for improvement
		$matches = &countMatches($mapstat);
		if($prev_matches < $matches) { # Improvement
			$prev_matches = $matches;
		} else { # no improvement
			$reAsm = 0;
		}
		
		# update iter
		if(0 < $iter) {
			--$reAsm;
		}
		
		# check if target is empty
		if(&checkTarget($target)) {
			$reAsm = 0;
		}
		++$iternumber;
	}
	
	# Save wanted files and clean up
	&runCmDie(sprintf("mv %s %s.fasta", $target, $output));
	$target = sprintf("%s/assembly_graph_with_scaffolds.gfa", $SPAdes_out);
	if(&checkTarget($target) == 0) {
		&runCmDie(sprintf("mv %s %s.gfa", $target, $output));
		&runCmDie(sprintf("mv %s.tmp/tmp.target.kma.frag.gz %s.frag.gz", $output, $output));
		&runCmDie(sprintf("mv %s.tmp/tmp.target.kma.frag_raw.gz %s.frag_raw.gz", $output, $output));
	} else {
		&runCmDie(sprintf("touch %s.gfa", $output));
		&runCmDie(sprintf("touch %s.frag.gz", $output));
		&runCmDie(sprintf("touch %s.frag_raw.gz", $output));
	}
	&runCmDie(sprintf("rm -rf %s.tmp/", $output));
	
	print STDERR sprintf("Total number of iterations:\t%d\n", $iternumber);
	
	return 0;
}

exit(&main(scalar(@ARGV), @ARGV));
