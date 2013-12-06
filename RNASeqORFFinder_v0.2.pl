#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw(time);
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw (sum shuffle min max);
use threads;
use threads::shared;
use Statistics::Descriptive;
use URI::Escape;
use Cwd 'abs_path';
use Bio::Range;

#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a perl script to find the ORF and define coding and nonCoding transcripts based on RNA-Seq data.
#
#	Input
#		--gffPath=						file path[compulsory]; path of the reference GFF for gene annotation;
#		--fastaPath=					file path [compulsory]; the path fasta file contains the genome sequence, for generating blank perl storables;
#		--getorfFastaPath=				file path [compulsory]; the path of a fasta file that contains the getorf result;
#		--fullPileupStorableIndexPath=	file path[compulsory]; path of the index file of the pileupStorable to be counted, in full mode;
#		--end5PileupStorableIndexPath=	file path[compulsory]; path of the index file of the pileupStorable to be counted, in end 5 mode;
#		--maxThread=					integer [4]; max number of threads to be used;
#		--outDir=						directory path ['./BAMToReadEndPerlStorable/']; output directory;
#
#	Usage
#		/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/RNASeqORFFinder/v0.1/RNASeqORFFinder_v0.1.pl --getorfFastaPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/trypanosome/inUse/927/getorf/getorf.fa --gffPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/trypanosome/inUse/927/TbruceiTreu927_TriTrypDB-4.2.gff --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/trypanosome/inUse/927/TbruceiTreu927Genomic_TriTrypDB-4.2.fasta --maxThread=7
#
#	v0.1
#	[Mon 26 Aug 2013 16:44:23 CEST] debut;
#
#	v0.2
#	[Fri 27 Sep 2013 15:23:53 CEST] will generate percentage of read 5 end at 1st codon position
#
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-09-27 17:27]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/RNASeqORFFinder/v0.2/RNASeqORFFinder_v0.2.pl --getorfFastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/trypanosome/inUse/927/v4.2/getorf_9nt/getorf.fa --gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/trypanosome/inUse/927/v4.2/TbruceiTreu927_TriTrypDB-4.2.gff --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/trypanosome/inUse/927/v4.2/TbruceiTreu927Genomic_TriTrypDB-4.2.fasta --maxThread=4 --end5PileupStorableIndexPath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS026_27_29_30_pooled/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.no/cntgCovPls/index.hsh.pls --fullPileupStorableIndexPath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS026_27_29_30_pooled/BAMToReadEndPerlStorable/countMode.full.offset.0.baseComp.no/cntgCovPls/index.hsh.pls --outDir=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS026_27_29_30_pooled/RNASeqORFFinder_with_codonPos/
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/RNASeqORFFinder/v0.2/RNASeqORFFinder_v0.2.pl
#	--getorfFastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/trypanosome/inUse/927/v4.2/getorf_9nt/getorf.fa
#	--gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/trypanosome/inUse/927/v4.2/TbruceiTreu927_TriTrypDB-4.2.gff
#	--fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/trypanosome/inUse/927/v4.2/TbruceiTreu927Genomic_TriTrypDB-4.2.fasta
#	--maxThread=4
#	--end5PileupStorableIndexPath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS026_27_29_30_pooled/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.no/cntgCovPls/index.hsh.pls
#	--fullPileupStorableIndexPath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS026_27_29_30_pooled/BAMToReadEndPerlStorable/countMode.full.offset.0.baseComp.no/cntgCovPls/index.hsh.pls
#	--outDir=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS026_27_29_30_pooled/RNASeqORFFinder_final/
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $globalScriptDirPath = dirname(rel2abs($0));
open DEBUGLOG, ">", "$globalScriptDirPath/debug.log.txt";
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|1161, readParameters|1526
#	secondaryDependOnSub: currentTime|636
#
#<section ID="startingTasks" num="0">
########################################################################## 
&printCMDLogOrFinishMessage("CMDLog");#->1161

my ($gffPath, $fastaPath, $getorfFastaPath, $fullPileupStorableIndexPath, $maxThread, $end5PileupStorableIndexPath, $outDir) = &readParameters();#->1526

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="1">
my $minValidCov = 1;
my $minCovPerNt = 2;
my $minCovPosPct = 70;
my $uORFMinNtLen = 9; #---minLen of uORF (i.e. the ORF that is overlapping the annotated 5’UTR)
my $nonAnnoORFMinNtLen = 30; #---minLen of nonAnno ORF (i.e. the ORF that is not overlapping with any existing annotations)
my $keepOvrlpNovelORF = 'no';
my $margin = 20; #---the margin to define nonAnno ORF, i.e. the nonAnno ORF have to be at least 20 nt away from any annotated features
my $paramTag = "VC$minValidCov.CN$minCovPerNt.CP$minCovPosPct.UL$uORFMinNtLen.NL.$nonAnnoORFMinNtLen.KO_$keepOvrlpNovelORF.MG.$margin";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="2">
my @mkDirAry;
my $resultDir = "$outDir/$paramTag/"; push @mkDirAry, $resultDir;
my $resultStorableDir = "$resultDir/storable/"; push @mkDirAry, $resultStorableDir;
my $resultLogDir = "$resultDir/log/"; push @mkDirAry, $resultLogDir;
my $resultGFFDir = "$resultDir/GFF/"; push @mkDirAry, $resultGFFDir;
foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_defineOutFilePath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutFilePath" num="3">
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_processInputData
#	primaryDependOnSub: addCDSRngGeneWithoutCDS|248, checkGeneInfo|269, getIndivCntgCovPlsPath|833, readGFF_oneRNAPerGene|1373, readMultiFasta|1472, zipUnzipCntgCovInPlsPathHsh|1653
#	secondaryDependOnSub: currentTime|636, reportStatus|1564
#
#<section ID="processInputData" num="4">
########################################################################## 
#----------Read Fasta
#my ($fastaHsh_ref) = &readMultiFasta($fastaPath);
my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);#->1373
&checkGeneInfo($geneInfoHsh_ref);#->269
&addCDSRngGeneWithoutCDS($geneInfoHsh_ref);#->248

my ($fullPileupStorablePathHsh_ref) = &getIndivCntgCovPlsPath($fullPileupStorableIndexPath);#->833
my ($end5PileupStorablePathHsh_ref) = &getIndivCntgCovPlsPath($end5PileupStorableIndexPath);#->833
&zipUnzipCntgCovInPlsPathHsh('unzip', $fullPileupStorablePathHsh_ref);#->1653
&zipUnzipCntgCovInPlsPathHsh('unzip', $end5PileupStorablePathHsh_ref);#->1653
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_getAllPossibleORFs
#	primaryDependOnSub: parseGetorfResults|1092
#	secondaryDependOnSub: readMultiFasta|1472, reportStatus|1564
#
#<section ID="getAllPossibleORFs" num="5">
my ($getorfInfoHsh_ref) = &parseGetorfResults($getorfFastaPath, $resultStorableDir);#->1092
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_checkOverlappingAndCoverageOfNovelORF
#	primaryDependOnSub: countNovelORFCoverage|601, getNonOverlappingORF|866, getNovelORFCoverage|944, getORFCodonPosCount|982, tagOverlappingShorterTranscribedORF|1585
#	secondaryDependOnSub: checkOverlapAndProximity_withMargin|297, checkOvrlpTwoRangePair|534, generateGeneByCntgHsh|654, getCoverageOfItemRngType|702, getCtgryGeneInfo|796, reportStatus|1564
#
#<section ID="checkOverlappingAndCoverageOfNovelORF" num="6">
my ($novelORFInfoHsh_ref, $novelORFByCntgHsh_ref) = &getNonOverlappingORF($getorfInfoHsh_ref, $geneInfoHsh_ref, $maxThread, $resultStorableDir, $resultGFFDir, $keepOvrlpNovelORF, $margin, $nonAnnoORFMinNtLen, $uORFMinNtLen);#->866

my ($novelORFCovHsh_ref) = &getNovelORFCoverage($fullPileupStorablePathHsh_ref, $novelORFInfoHsh_ref, $novelORFByCntgHsh_ref, $resultStorableDir, $maxThread);#->944

&countNovelORFCoverage($novelORFCovHsh_ref, $novelORFInfoHsh_ref, $minValidCov, $minCovPerNt, $minCovPosPct, $resultGFFDir);#->601

my ($codonPosCountHsh_ref) = &getORFCodonPosCount($end5PileupStorablePathHsh_ref, $novelORFInfoHsh_ref, $novelORFByCntgHsh_ref, $geneInfoHsh_ref, $maxThread, $resultStorableDir);#->982

&tagOverlappingShorterTranscribedORF($keepOvrlpNovelORF, $novelORFInfoHsh_ref, 1);#->1585
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 7_printTheTranscribeNovelORFLog
#	primaryDependOnSub: outputAllGff|1039, printGeneBaseduORFInfo|1265, printTranscribedNovelORFInfo|1331, readMultiFasta|1472
#	secondaryDependOnSub: printGFF_oneRNAPerGene_chooseStrnd_filterAry|1194, reportStatus|1564
#
#<section ID="printTheTranscribeNovelORFLog" num="7">
&printTranscribedNovelORFInfo($novelORFInfoHsh_ref, $codonPosCountHsh_ref, $resultLogDir);#->1331

#----just to get the name of cntig
my ($fastaHsh_ref) = &readMultiFasta($fastaPath);#->1472
&printGeneBaseduORFInfo($novelORFInfoHsh_ref, $resultLogDir, $geneInfoHsh_ref, $fastaHsh_ref);#->1265
&outputAllGff($novelORFInfoHsh_ref, $getorfInfoHsh_ref, $resultGFFDir);#->1039
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 8_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|1161
#	secondaryDependOnSub: currentTime|636
#
#<section ID="finishingTasks" num="8">
&printCMDLogOrFinishMessage("finishMessage");#->1161
close DEBUGLOG;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	coverage [n=1]:
#		getCoverageOfItemRngType
#
#	fasta [n=1]:
#		readMultiFasta
#
#	general [n=5]:
#		currentTime, getCtgryGeneInfo, printCMDLogOrFinishMessage
#		readParameters, reportStatus
#
#	gff [n=5]:
#		addCDSRngGeneWithoutCDS, checkGeneInfo, generateGeneByCntgHsh
#		printGFF_oneRNAPerGene_chooseStrnd_filterAry, readGFF_oneRNAPerGene
#
#	multithread [n=1]:
#		generateThreadHshWithRandomCntg
#
#	range [n=2]:
#		checkOverlapAndProximity_withMargin, checkOvrlpTwoRangePair
#
#	specific [n=8]:
#		countNovelORFCoverage, getNonOverlappingORF, getNovelORFCoverage
#		getORFCodonPosCount, outputAllGff, printGeneBaseduORFInfo
#		printTranscribedNovelORFInfo, tagOverlappingShorterTranscribedORF
#
#	storable [n=2]:
#		getIndivCntgCovPlsPath, zipUnzipCntgCovInPlsPathHsh
#
#	thridPartyApp [n=1]:
#		parseGetorfResults
#
#====================================================================================================================================================#

sub addCDSRngGeneWithoutCDS {
#....................................................................................................................................................#
#	subroutineCategory: gff
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|134
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref
#	output: 
#	toCall: &addCDSRngGeneWithoutCDS($geneInfoHsh_ref);
#	calledInLine: 144
#....................................................................................................................................................#
	my ($geneInfoHsh_ref) = @_;
	
	#---the purpose if to prevent the getORF hit on non-mRNA region like rRNA
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		$geneInfoHsh_ref->{$geneID}{'CDSRng'} = $geneInfoHsh_ref->{$geneID}{'geneRng'} if not $geneInfoHsh_ref->{$geneID}{'CDSRng'};
	}

	return ();
}
sub checkGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: gff
#	dependOnSub: reportStatus|1564
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|134
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref
#	output: 
#	toCall: &checkGeneInfo($geneInfoHsh_ref);
#	calledInLine: 143
#....................................................................................................................................................#
	
	my ($geneInfoHsh_ref) = @_;
	
	&reportStatus("Checking gene categories", 0, "\n");#->1564
	my $ctrgyCountHsh_ref = {};
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		$ctrgyCountHsh_ref->{$ctgry}++;
	}
	
	foreach my $ctgry (sort keys %{$ctrgyCountHsh_ref}) {
		&reportStatus("Item in $ctgry = $ctrgyCountHsh_ref->{$ctgry}", 0, "\n");#->1564
	}
	
	return ();
}
sub checkOverlapAndProximity_withMargin {
#....................................................................................................................................................#
#	subroutineCategory: range
#	dependOnSub: currentTime|636, generateThreadHshWithRandomCntg|675, reportStatus|1564
#	appearInSub: getNonOverlappingORF|866, tagOverlappingShorterTranscribedORF|1585
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_checkOverlappingAndCoverageOfNovelORF|164
#	input: $checkPrxmty, $maxThread, $qryInfoHsh_ref, $qryMargin, $qryRngType, $refInfoHsh_ref, $refMargin, $refRngType, $reportExactMatch
#	output: $hitAndPrxmtyByQryHsh_ref, $hitAndPrxmtyByRefHsh_ref
#	toCall: my ($hitAndPrxmtyByRefHsh_ref, $hitAndPrxmtyByQryHsh_ref) = &checkOverlapAndProximity_withMargin($refInfoHsh_ref, $qryInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType, $refMargin, $qryMargin);
#	calledInLine: 894, 1614
#....................................................................................................................................................#
	#---incoming variables
	my ($refInfoHsh_ref, $qryInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType, $refMargin, $qryMargin) = @_;

	#---outgoing variables
	my $refGeneNumTotal = 0;
	
	#---make a tmpHsh to contain all cntgs the have either ref and qry
	my $tmpCntgHsh_ref = {};
	my $refCntgHsh_ref = {};
	my $qryCntgHsh_ref = {};
	
	foreach my $refGeneID (keys %{$refInfoHsh_ref}) {
		if ($refInfoHsh_ref->{$refGeneID}{$refRngType}) {
			@{$refInfoHsh_ref->{$refGeneID}{$refRngType}} = sort {$a <=> $b} @{$refInfoHsh_ref->{$refGeneID}{$refRngType}};
			$refGeneNumTotal++;
			$tmpCntgHsh_ref->{$refInfoHsh_ref->{$refGeneID}{'cntg'}}++;
			$refCntgHsh_ref->{$refInfoHsh_ref->{$refGeneID}{'cntg'}}{$refGeneID}++;
		}
	}
	
	foreach my $qryGeneID (keys %{$qryInfoHsh_ref}) {
		if ($qryInfoHsh_ref->{$qryGeneID}{$qryRngType}) {
			@{$qryInfoHsh_ref->{$qryGeneID}{$qryRngType}} = sort {$a <=> $b} @{$qryInfoHsh_ref->{$qryGeneID}{$qryRngType}};
			$tmpCntgHsh_ref->{$qryInfoHsh_ref->{$qryGeneID}{'cntg'}}++;
			$qryCntgHsh_ref->{$qryInfoHsh_ref->{$qryGeneID}{'cntg'}}{$qryGeneID}++;
		}
	}
	
	my $totalCntgNum = keys %{$tmpCntgHsh_ref};
	my @cntgAry = keys %{$tmpCntgHsh_ref};

	my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($maxThread, \@cntgAry);#->675
	my $refGeneNumProc :shared = 0;
	my %threadHsh = ();

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
		my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
		&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 20, "\n");#->1564

		#---spawn a new thread
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
		
			sub {
				my ($cntgAry_ref) = @_;
	
				my $hitAndPrxmtyByRefHsh_InThr_ref = {};
				my $hitAndPrxmtyByQryHsh_InThr_ref = {};

				foreach my $cntg (@{$cntgAry_ref}) {

					#---update the on-screen progress
					#&reportStatus("Finding overlaping on $cntg", 20, "\r");#->1564
		
					my $tmpPrxmtyByRefHsh_ref = {};
					my $tmpPrxmtyByQryHsh_ref = {};

					if ((exists $qryCntgHsh_ref->{$cntg}) and (exists $refCntgHsh_ref->{$cntg})) {#---if there are both ref and qry can both be found on cntg
						foreach my $refGeneID (keys %{$refCntgHsh_ref->{$cntg}}) {#--- all ftur on the $strnd of $cntg of refGff
							my ($refStart, $refEnd) = ($refInfoHsh_ref->{$refGeneID}{$refRngType}->[0]-$refMargin, $refInfoHsh_ref->{$refGeneID}{$refRngType}->[-1]+$refMargin);
				
							$refGeneNumProc++;
							
							&reportStatus("$refGeneNumProc of $refGeneNumTotal reference genes checked", 20, "\r");#->1564

							foreach my $qryGeneID (keys %{$qryCntgHsh_ref->{$cntg}}) {#--- all ftur on the $strnd of $cntg of QryGtf
	
								my $samestrnd = "no";
								$samestrnd = "yes" if ($refInfoHsh_ref->{$refGeneID}{'strnd'} eq $qryInfoHsh_ref->{$qryGeneID}{'strnd'});
								my ($qryStart, $qryEnd) = ($qryInfoHsh_ref->{$qryGeneID}{$qryRngType}->[0]-$qryMargin, $qryInfoHsh_ref->{$qryGeneID}{$qryRngType}->[-1]+$qryMargin);

								my $scene;
								my $ovrlpSize;

								if (($refStart == $qryStart) && ($refEnd == $qryEnd)) {#---scene 0
									$scene = 'exactMatch';
									$ovrlpSize = $qryEnd - $qryStart;
						
								} elsif (($refStart<=$qryStart)&&($refEnd>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 1
									$scene = 'overlapTail';
									$ovrlpSize = $refEnd - $qryStart;

								} elsif (($refStart>=$qryStart)&&($refStart<=$qryEnd)&&($refEnd>=$qryEnd)) {#---scene 2
									$scene = 'overlapHead';
									$ovrlpSize = $qryEnd - $refStart;

								} elsif (($refStart<=$qryStart)&&($refEnd>=$qryEnd)) {#---scene 3
									$scene = 'cover';
									$ovrlpSize = $qryEnd - $qryStart;

								} elsif (($refStart>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 4
									$scene = 'within';
									$ovrlpSize = $refEnd - $refStart;

								#------Proximity with ref's tail proximal to qry's head
								} elsif (($refEnd<=$qryStart)&&($refEnd<$qryEnd)) {#---scene 5 ---> ref Tail, qry Head

									$scene = 'prxmtyTail';

									if ($checkPrxmty eq "yes") {
										my $tmpPrmxty = $qryStart - $refEnd;
										$tmpPrxmtyByRefHsh_ref->{'XS'}{$refGeneID}{"T"}{$qryGeneID} = $tmpPrmxty;
										$tmpPrxmtyByQryHsh_ref->{'XS'}{$qryGeneID}{"H"}{$refGeneID} = $tmpPrmxty;

										if ($samestrnd eq "yes") {
											$tmpPrxmtyByRefHsh_ref->{'SS'}{$refGeneID}{"T"}{$qryGeneID} = $tmpPrmxty;
											$tmpPrxmtyByQryHsh_ref->{'SS'}{$qryGeneID}{"H"}{$refGeneID} = $tmpPrmxty;
										}
									}

								#------Proximity with ref's head proximal to qry's tail
								} elsif (($refStart>=$qryEnd)&&($refStart>$qryStart)) {#---scene 6 ---> ref Head, qry Tail

									$scene = 'prxmtyHead';

									if ($checkPrxmty eq "yes") {
										my $tmpPrmxty = $refStart - $qryEnd;
										$tmpPrxmtyByRefHsh_ref->{'XS'}{$refGeneID}{"H"}{$qryGeneID} = $tmpPrmxty;
										$tmpPrxmtyByQryHsh_ref->{'XS'}{$qryGeneID}{"T"}{$refGeneID} = $tmpPrmxty;

										if ($samestrnd eq "yes") {
											$tmpPrxmtyByRefHsh_ref->{'SS'}{$refGeneID}{"H"}{$qryGeneID} = $tmpPrmxty;
											$tmpPrxmtyByQryHsh_ref->{'SS'}{$qryGeneID}{"T"}{$refGeneID} = $tmpPrmxty;
										}
									}

								} else {#---BUG! possibly other scene?
									#print "[".&currentTime()."] refStart=$refStart; refEnd=$refEnd; qryStart=$qryStart; qryEnd=$qryEnd\n";#->636
									die "Unexpected overlapping scene between $refGeneID and $qryGeneID. It's a Bug. Program qutting.\n";
								}
								
								if ($scene ne 'prxmtyTail' and $scene ne 'prxmtyHead' and not ($reportExactMatch eq 'no' and $scene eq 'exactMatch')) {

									@{$hitAndPrxmtyByRefHsh_InThr_ref->{'XS'}{'hit'}{$refGeneID}{$qryGeneID}} = ($scene, $ovrlpSize);
									@{$hitAndPrxmtyByQryHsh_InThr_ref->{'XS'}{'hit'}{$qryGeneID}{$refGeneID}} = ($scene, $ovrlpSize);
	
									#print DEBUGLOG join "", (join "\t", ("XS", $qryGeneID, $scene, $ovrlpSize)), "\n" if ($refGeneID eq "getorf_Tb927_02_v4_+_162");

									if ($samestrnd eq "yes") {
	
										#print DEBUGLOG join "", (join "\t", ("SS", $qryGeneID, $scene, $ovrlpSize)), "\n" if ($refGeneID eq "getorf_Tb927_02_v4_+_162");
										
										@{$hitAndPrxmtyByRefHsh_InThr_ref->{'SS'}{'hit'}{$refGeneID}{$qryGeneID}} = ($scene, $ovrlpSize);
										@{$hitAndPrxmtyByQryHsh_InThr_ref->{'SS'}{'hit'}{$qryGeneID}{$refGeneID}} = ($scene, $ovrlpSize);
									}
								}
							}
						}
					}

					#---find the closest proximity for all refs
					if ($checkPrxmty eq "yes") {
						my $refQryRefHsh_ref = {};

						$refQryRefHsh_ref->{'ref'}{'tmpPrxmtyHsh_ref'} = $tmpPrxmtyByRefHsh_ref;
						$refQryRefHsh_ref->{'ref'}{'cntgHsh_ref'} = $refCntgHsh_ref;
						$refQryRefHsh_ref->{'ref'}{'hitAndPrxmtyHsh_ref'} = $hitAndPrxmtyByRefHsh_InThr_ref;

						$refQryRefHsh_ref->{'qry'}{'tmpPrxmtyHsh_ref'} = $tmpPrxmtyByQryHsh_ref;
						$refQryRefHsh_ref->{'qry'}{'cntgHsh_ref'} = $qryCntgHsh_ref;
						$refQryRefHsh_ref->{'qry'}{'hitAndPrxmtyHsh_ref'} = $hitAndPrxmtyByQryHsh_InThr_ref;
			
						foreach my $refOrQry ('ref', 'qry') {

							my $cntgHsh_ref = $refQryRefHsh_ref->{$refOrQry}{'cntgHsh_ref'};
							my $tmpPrxmtyHsh_ref = $refQryRefHsh_ref->{$refOrQry}{'tmpPrxmtyHsh_ref'};
							my $hitAndPrxmtyHsh_ref = $refQryRefHsh_ref->{$refOrQry}{'hitAndPrxmtyHsh_ref'};

							foreach my $ftur (keys %{$cntgHsh_ref->{$cntg}}) {
								foreach my $XSOrSS ('XS', 'SS') {
									foreach my $HOrT ('H', 'T') {
										$tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{"edge"} = -999 if (not exists $tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT});
										foreach my $otherFtur (sort {$tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{$a} <=> $tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{$b}} keys %{$tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}}) {
											@{$hitAndPrxmtyHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT}} = ($tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{$otherFtur}, $otherFtur);
											last; #---sample the smallest only
										}
									}
								}
							}
						}
					}
				}
				
				return ($hitAndPrxmtyByRefHsh_InThr_ref, $hitAndPrxmtyByQryHsh_InThr_ref);
			}
			,($cntgAry_ref)
		);
	}
	
	
	my %tmpTransferThrDataHsh = ();
	$tmpTransferThrDataHsh{'ref'}{'all'} = {};
	$tmpTransferThrDataHsh{'qry'}{'all'} = {};
	
	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				($tmpTransferThrDataHsh{'ref'}{'thr'}, $tmpTransferThrDataHsh{'qry'}{'thr'}) = $thr->join;
				foreach my $refOrQry (keys %tmpTransferThrDataHsh) {
					my ($allHsh_ref, $thrHsh_ref) = ($tmpTransferThrDataHsh{$refOrQry}{'all'}, $tmpTransferThrDataHsh{$refOrQry}{'thr'});
					foreach my $XSOrSS ('XS', 'SS') {
						foreach my $ftur (keys %{$thrHsh_ref->{$XSOrSS}{'hit'}}) {
							foreach my $hitftur (keys %{$thrHsh_ref->{$XSOrSS}{'hit'}{$ftur}}) {
								@{$allHsh_ref->{$XSOrSS}{'hit'}{$ftur}{$hitftur}} = @{$thrHsh_ref->{$XSOrSS}{'hit'}{$ftur}{$hitftur}};
							}
						}
						
						if ($thrHsh_ref->{$XSOrSS}{'prxmty'}) {
							foreach my $ftur (keys %{$thrHsh_ref->{$XSOrSS}{'prxmty'}}) {
								foreach my $HOrT (keys %{$thrHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}}) {
									@{$allHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT}} = @{$thrHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT}} if $thrHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT};
								}
							}
						}
					}
				}
				delete $threadHsh{$threadNum};
			}
		}
		sleep 1;
	}

	print "\n";

	my $hitAndPrxmtyByRefHsh_ref = $tmpTransferThrDataHsh{'ref'}{'all'};
	my $hitAndPrxmtyByQryHsh_ref = $tmpTransferThrDataHsh{'qry'}{'all'};

	return ($hitAndPrxmtyByRefHsh_ref, $hitAndPrxmtyByQryHsh_ref);
}
sub checkOvrlpTwoRangePair {
#....................................................................................................................................................#
#	subroutineCategory: range
#	dependOnSub: >none
#	appearInSub: getNonOverlappingORF|866
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_checkOverlappingAndCoverageOfNovelORF|164
#	input: $qryRngAry_ref, $refRngAry_ref
#	output: $ovrlpSize, $prmxty, $scene
#	toCall: my ($scene, $ovrlpSize, $prmxty) = &checkOvrlpTwoRangePair($refRngAry_ref, $qryRngAry_ref);
#	calledInLine: 916, 918
#....................................................................................................................................................#
	my ($refRngAry_ref, $qryRngAry_ref) = @_;
	
	my %tmpHsh = ();
	
	$tmpHsh{'ref'}{'rngAry_ref'} = $refRngAry_ref;
	$tmpHsh{'qry'}{'rngAry_ref'} = $qryRngAry_ref;
	
	#---get the start and end
	foreach (keys %tmpHsh) {
		@{$tmpHsh{$_}{'rngAry_ref'}} = sort {$a <=> $b} @{$tmpHsh{$_}{'rngAry_ref'}};
		$tmpHsh{$_}{'start'} = $tmpHsh{$_}{'rngAry_ref'}->[0];
		$tmpHsh{$_}{'end'} = $tmpHsh{$_}{'rngAry_ref'}->[-1];
	}
	
	my $scene = undef;
	my $ovrlpSize = my $prmxty = 0;

	if (($tmpHsh{'ref'}{'start'} == $tmpHsh{'qry'}{'start'}) && ($tmpHsh{'ref'}{'end'} == $tmpHsh{'qry'}{'end'})) {#---scene 0
		$scene = 'exactMatch';
		$ovrlpSize = $tmpHsh{'qry'}{'end'} - $tmpHsh{'qry'}{'start'};
						
	} elsif (($tmpHsh{'ref'}{'start'}<=$tmpHsh{'qry'}{'start'})&&($tmpHsh{'ref'}{'end'}>=$tmpHsh{'qry'}{'start'})&&($tmpHsh{'ref'}{'end'}<=$tmpHsh{'qry'}{'end'})) {#---scene 1
		$scene = 'overlapTail';
		$ovrlpSize = $tmpHsh{'ref'}{'end'} - $tmpHsh{'qry'}{'start'};

	} elsif (($tmpHsh{'ref'}{'start'}>=$tmpHsh{'qry'}{'start'})&&($tmpHsh{'ref'}{'start'}<=$tmpHsh{'qry'}{'end'})&&($tmpHsh{'ref'}{'end'}>=$tmpHsh{'qry'}{'end'})) {#---scene 2
		$scene = 'overlapHead';
		$ovrlpSize = $tmpHsh{'qry'}{'end'} - $tmpHsh{'ref'}{'start'};

	} elsif (($tmpHsh{'ref'}{'start'}<=$tmpHsh{'qry'}{'start'})&&($tmpHsh{'ref'}{'end'}>=$tmpHsh{'qry'}{'end'})) {#---scene 3
		$scene = 'cover'; #---reference cover query
		$ovrlpSize = $tmpHsh{'qry'}{'end'} - $tmpHsh{'qry'}{'start'};

	} elsif (($tmpHsh{'ref'}{'start'}>=$tmpHsh{'qry'}{'start'})&&($tmpHsh{'ref'}{'end'}<=$tmpHsh{'qry'}{'end'})) {#---scene 4
		$scene = 'within'; #---reference within query
		$ovrlpSize = $tmpHsh{'ref'}{'end'} - $tmpHsh{'ref'}{'start'};

	#------Proximity with ref's tail proximal to qry's head
	} elsif (($tmpHsh{'ref'}{'end'}<=$tmpHsh{'qry'}{'start'})&&($tmpHsh{'ref'}{'end'}<$tmpHsh{'qry'}{'end'})) {#---scene 5 ---> ref Tail, qry Head

		$scene = 'prxmtyTail';
		$prmxty = $tmpHsh{'qry'}{'start'} - $tmpHsh{'ref'}{'end'};

	#------Proximity with ref's head proximal to qry's tail
	} elsif (($tmpHsh{'ref'}{'start'}>=$tmpHsh{'qry'}{'end'})&&($tmpHsh{'ref'}{'start'}>$tmpHsh{'qry'}{'start'})) {#---scene 6 ---> ref Head, qry Tail

		$scene = 'prxmtyHead';
		$prmxty = $tmpHsh{'ref'}{'start'} - $tmpHsh{'qry'}{'end'};

	} else {#---BUG! possibly other scene?
		die "Unexpected overlapping scene. qutting.\n";
	}

	return ($scene, $ovrlpSize, $prmxty);
}
sub countNovelORFCoverage {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 6_checkOverlappingAndCoverageOfNovelORF|164
#	secondaryAppearInSection: >none
#	input: $minCovPerNt, $minCovPosPct, $minValidCov, $novelORFCovHsh_ref, $novelORFInfoHsh_ref, $resultGFFDir
#	output: 
#	toCall: &countNovelORFCoverage($novelORFCovHsh_ref, $novelORFInfoHsh_ref, $minValidCov, $minCovPerNt, $minCovPosPct, $resultGFFDir);
#	calledInLine: 173
#....................................................................................................................................................#
	my ($novelORFCovHsh_ref, $novelORFInfoHsh_ref, $minValidCov, $minCovPerNt, $minCovPosPct, $resultGFFDir) = @_;

	foreach my $ORFID (keys %{$novelORFInfoHsh_ref}) {
		$novelORFInfoHsh_ref->{$ORFID}{'covPosPct'} = 0;
		$novelORFInfoHsh_ref->{$ORFID}{'covPerNt'} = 0;
		$novelORFInfoHsh_ref->{$ORFID}{'transcribed'} = 'no';

		if ($novelORFCovHsh_ref->{$ORFID}) {
			my $length = @{$novelORFCovHsh_ref->{$ORFID}};
			my $covPos = 0;
			my $covSum = 0;
			foreach my $cov (@{$novelORFCovHsh_ref->{$ORFID}}) {
				$covSum += $cov;
				$covPos++ if $cov >= $minValidCov;
			}
			$novelORFInfoHsh_ref->{$ORFID}{'covPosPct'} = sprintf "%.2f", 100*$covPos/$length;
			$novelORFInfoHsh_ref->{$ORFID}{'covPerNt'} = sprintf "%.2f", $covSum/$length;
			$novelORFInfoHsh_ref->{$ORFID}{'transcribed'} = 'yes' if $novelORFInfoHsh_ref->{$ORFID}{'covPosPct'} >= $minCovPosPct and $novelORFInfoHsh_ref->{$ORFID}{'covPerNt'} >= $minCovPerNt;
		}
	}

	return ();
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: checkOverlapAndProximity_withMargin|297, printCMDLogOrFinishMessage|1161, readGFF_oneRNAPerGene|1373, reportStatus|1564
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|79, 4_processInputData|134, 8_finishingTasks|197
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 435, 1181, 1184, 1189, 1393, 1580
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub generateGeneByCntgHsh {
#....................................................................................................................................................#
#	subroutineCategory: gff
#	dependOnSub: reportStatus|1564
#	appearInSub: getNonOverlappingORF|866
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_checkOverlappingAndCoverageOfNovelORF|164
#	input: $geneInfoHsh_ref
#	output: $geneByCntgHsh_ref
#	toCall: my ($geneByCntgHsh_ref) = &generateGeneByCntgHsh($geneInfoHsh_ref);
#	calledInLine: 939
#....................................................................................................................................................#
	my ($geneInfoHsh_ref) = @_;
	
	&reportStatus("Generating geneByCntgHsh", 0, "\n");#->1564

	my $geneByCntgHsh_ref = {};
	$geneByCntgHsh_ref->{$geneInfoHsh_ref->{$_}{'cntg'}}{$_}++ foreach (keys %{$geneInfoHsh_ref});

	return ($geneByCntgHsh_ref);
}
sub generateThreadHshWithRandomCntg {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: >none
#	appearInSub: checkOverlapAndProximity_withMargin|297, getCoverageOfItemRngType|702
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $cntgAry_ref, $threadToSpawn
#	output: $randCntgInThreadHsh_ref
#	toCall: my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($threadToSpawn, $cntgAry_ref);
#	calledInLine: 339, 716
#....................................................................................................................................................#

	my ($threadToSpawn, $cntgAry_ref) = @_;

	my @shuffleCntgAry = shuffle(@{$cntgAry_ref});
	my $threadNum = 1;
	my $randCntgInThreadHsh_ref = {};
	foreach my $cntg (@{$cntgAry_ref}) {
		$threadNum = 1 if $threadNum > $threadToSpawn;
		push @{$randCntgInThreadHsh_ref->{$threadNum}}, $cntg;
		$threadNum++;
	}
	
	return $randCntgInThreadHsh_ref;

}
sub getCoverageOfItemRngType {
#....................................................................................................................................................#
#	subroutineCategory: coverage
#	dependOnSub: generateThreadHshWithRandomCntg|675, reportStatus|1564
#	appearInSub: getNovelORFCoverage|944, getORFCodonPosCount|982
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_checkOverlappingAndCoverageOfNovelORF|164
#	input: $dirtn, $fullPileupStorablePathHsh_ref, $itemByCntgHsh_ref, $itemInfoHsh_ref, $margin, $maxThread, $rngType
#	output: $covHsh_ref
#	toCall: my ($covHsh_ref) = &getCoverageOfItemRngType($fullPileupStorablePathHsh_ref, $itemInfoHsh_ref, $itemByCntgHsh_ref, $maxThread, $rngType, $margin, $dirtn);
#	calledInLine: 974, 1020
#....................................................................................................................................................#
	my ($fullPileupStorablePathHsh_ref, $itemInfoHsh_ref, $itemByCntgHsh_ref, $maxThread, $rngType, $margin, $dirtn) = @_;
	
	my @cntgAry = (keys %{$fullPileupStorablePathHsh_ref});
	my $randCntgInThreadHsh_ref = &generateThreadHshWithRandomCntg($maxThread, \@cntgAry);#->675
	my $cntgProc :shared = 0;
	my %threadHsh = ();

	my %strndCovToGetHsh = ();

	$strndCovToGetHsh{'+'}{'s'} = '+';
	$strndCovToGetHsh{'+'}{'a'} = '-';
	$strndCovToGetHsh{'-'}{'s'} = '-';
	$strndCovToGetHsh{'-'}{'a'} = '+';

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
		my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
		&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 20, "\n");#->1564

		#---spawn a new thread
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
		
			sub {
				my ($cntgAry_ref) = @_;
				my $covInThrHsh_ref = {};
				foreach my $cntg (@{$cntgAry_ref}) {
					$cntgProc++;
					my $cntgCovAry_ref = retrieve($fullPileupStorablePathHsh_ref->{$cntg});
					next if not exists $itemByCntgHsh_ref->{$cntg};
					foreach my $itemID (keys %{$itemByCntgHsh_ref->{$cntg}}) {
						next if not $itemInfoHsh_ref->{$itemID}{$rngType};

						my @rng = sort {$a <=> $b} @{$itemInfoHsh_ref->{$itemID}{$rngType}};
						unshift @rng, ($rng[0]-1-$margin, $rng[0]-1);
						push @rng, ($rng[-1]+1, $rng[-1]+1+$margin);
	
						my $posCovAry_ref = ();
						for (my $i=0; $i < $#rng; $i += 2) {
							foreach my $j ($rng[$i]-1..$rng[$i+1]-1) {
								my %tmpCovHsh = ('+'=>0, '-'=>0);
								($tmpCovHsh{'+'}, $tmpCovHsh{'-'}) = split /,/, $cntgCovAry_ref->[$j] if ($cntgCovAry_ref->[$j]);
								push @{$posCovAry_ref}, $tmpCovHsh{$strndCovToGetHsh{$itemInfoHsh_ref->{$itemID}{'strnd'}}{$dirtn}};
							}
						}
	
						if (($itemInfoHsh_ref->{$itemID}{'strnd'} eq '-' and $dirtn eq 's') or ($itemInfoHsh_ref->{$itemID}{'strnd'} eq '+' and $dirtn eq 'a')) {
							@{$covInThrHsh_ref->{$itemID}} = reverse @{$posCovAry_ref};
						} else {
							@{$covInThrHsh_ref->{$itemID}} = @{$posCovAry_ref};
						}
					}
					
					&reportStatus("Finished counting items on $cntgProc cntg", 20, "\r");#->1564
				}
				return ($covInThrHsh_ref);
			}
			,($cntgAry_ref)
		);
	}
	
	my $covHsh_ref = {};
	
	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				my ($covInThrHsh_ref) = $thr->join;
				foreach my $itemID (keys %{$covInThrHsh_ref}) {
					$covHsh_ref->{$itemID} = $covInThrHsh_ref->{$itemID};
				}
				delete $threadHsh{$threadNum};
			}
		}
		sleep 1;
	}

	my $numItem = scalar(keys %{$covHsh_ref});
	&reportStatus("$rngType coverage of $numItem items in $dirtn direction were stored", 10, "\n");#->1564
	
	return ($covHsh_ref);

}
sub getCtgryGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: reportStatus|1564
#	appearInSub: getORFCodonPosCount|982
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_checkOverlappingAndCoverageOfNovelORF|164
#	input: $ctgryAry_ref, $geneInfoHsh_ref
#	output: $geneCtgryByCntgHsh_ref, $geneCtgryInfoHsh_ref
#	toCall: my ($geneCtgryInfoHsh_ref, $geneCtgryByCntgHsh_ref) = &getCtgryGeneInfo($geneInfoHsh_ref, $ctgryAry_ref);
#	calledInLine: 1006, 1007
#....................................................................................................................................................#

	my ($geneInfoHsh_ref, $ctgryAry_ref) = @_;
	
	my $geneCtgryInfoHsh_ref = {};
	my $geneCtgryByCntgHsh_ref = {};
	
	my $ctgryStr = join ",", @{$ctgryAry_ref};

	&reportStatus("Filtering GFF on cgtry $ctgryStr.\n", 10, "\n");#->1564
	
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		my $cntg = $geneInfoHsh_ref->{$geneID}{'cntg'};
		if (grep /^$ctgry$/, @{$ctgryAry_ref}) {
			%{$geneCtgryInfoHsh_ref->{$geneID}} = %{$geneInfoHsh_ref->{$geneID}};
			$geneCtgryByCntgHsh_ref->{$cntg}{$geneID}++;
		}
	}
	
	my $numGene = keys %{$geneCtgryInfoHsh_ref};
	
	&reportStatus("$numGene gene filtered on cgtry $ctgryStr", 10, "\n");#->1564
	
	return ($geneCtgryInfoHsh_ref, $geneCtgryByCntgHsh_ref);
}
sub getIndivCntgCovPlsPath {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: reportStatus|1564
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|134
#	secondaryAppearInSection: >none
#	input: $cntgCovPlsIndexPath
#	output: $cntgCovInPlsPathHsh_ref
#	toCall: my ($cntgCovInPlsPathHsh_ref) = &getIndivCntgCovPlsPath($cntgCovPlsIndexPath);
#	calledInLine: 146, 147
#....................................................................................................................................................#
	
	my ($cntgCovPlsIndexPath) = @_;

	$cntgCovPlsIndexPath =~ s/\.gz$//;

	my $cntgCovInPlsPathHsh_ref = {};
	
	system ("gzip -fd $cntgCovPlsIndexPath.gz") if -s "$cntgCovPlsIndexPath.gz";
	my %plsIndexHsh = %{retrieve($cntgCovPlsIndexPath)};
	
	my (undef, $cntgCovStroableDir, undef) = fileparse($cntgCovPlsIndexPath, qr/\.[^.]*/);
	foreach my $cntg (keys %plsIndexHsh) {
		my $cntgCovPlsPath = "$cntgCovStroableDir/$plsIndexHsh{$cntg}";
		die "cntgCovPlsPath $cntgCovPlsPath is invalid\n" if ((not -s $cntgCovPlsPath) and (not -s $cntgCovPlsPath.".gz"));
		$cntgCovInPlsPathHsh_ref->{$cntg} = $cntgCovPlsPath;
	}
	my $numCntg = keys %{$cntgCovInPlsPathHsh_ref};
	&reportStatus("pls path of $numCntg contig stored.", 0, "\n");#->1564
	
	return $cntgCovInPlsPathHsh_ref;
}
sub getNonOverlappingORF {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: checkOverlapAndProximity_withMargin|297, checkOvrlpTwoRangePair|534, generateGeneByCntgHsh|654, reportStatus|1564
#	appearInSub: >none
#	primaryAppearInSection: 6_checkOverlappingAndCoverageOfNovelORF|164
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref, $getorfInfoHsh_ref, $keepOvrlpNovelORF, $margin, $maxThread, $nonAnnoORFMinNtLen, $resultGFFDir, $resultStorableDir, $uORFMinNtLen
#	output: $novelORFByCntgHsh_ref, $novelORFInfoHsh_ref
#	toCall: my ($novelORFInfoHsh_ref, $novelORFByCntgHsh_ref) = &getNonOverlappingORF($getorfInfoHsh_ref, $geneInfoHsh_ref, $maxThread, $resultStorableDir, $resultGFFDir, $keepOvrlpNovelORF, $margin, $nonAnnoORFMinNtLen, $uORFMinNtLen);
#	calledInLine: 169
#....................................................................................................................................................#
	my ($getorfInfoHsh_ref, $geneInfoHsh_ref, $maxThread, $resultStorableDir, $resultGFFDir, $keepOvrlpNovelORF, $margin, $nonAnnoORFMinNtLen, $uORFMinNtLen) = @_;
	
	my $hitAndPrxmtyByGetORFHsh_ref = {};
	my $novelORFInfoHsh_ref = {};
	
	my $hitAndPrxmtyByGetORFHshPlsPath = "$resultStorableDir/hitAndPrxmtyByGetORFHsh.pls";
	if (-s $hitAndPrxmtyByGetORFHshPlsPath) {
		&reportStatus("Retrieving hitAndPrxmtyByGetORFHsh", 0,"\n");#->1564
		$hitAndPrxmtyByGetORFHsh_ref = retrieve($hitAndPrxmtyByGetORFHshPlsPath);
	} else {
		&reportStatus("Finding overlapping ORFs with genes", 0,"\n");#->1564
		my $checkPrxmty = 'no';
		my $reportExactMatch = 'yes';
		my $refRngType = 'geneRng';
		my $qryRngType = 'geneRng';
		my $refMargin = 0;
		my $qryMargin = $margin;
		($hitAndPrxmtyByGetORFHsh_ref, undef) = &checkOverlapAndProximity_withMargin($getorfInfoHsh_ref, $geneInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType, $refMargin, $qryMargin);#->297
		store($hitAndPrxmtyByGetORFHsh_ref, $hitAndPrxmtyByGetORFHshPlsPath);
	}
	
	&reportStatus("Recording overlapping ORFs with genes", 0,"\n");#->1564
	foreach my $ORFID (keys %{$getorfInfoHsh_ref}) {
		
		if (not $hitAndPrxmtyByGetORFHsh_ref->{'XS'}{'hit'}{$ORFID}) {#---not overlapping with any genes in any direction

			if ($getorfInfoHsh_ref->{$ORFID}{'CDSRng'}[1]-$getorfInfoHsh_ref->{$ORFID}{'CDSRng'}[0] >= $nonAnnoORFMinNtLen) {
				$novelORFInfoHsh_ref->{$ORFID} = $getorfInfoHsh_ref->{$ORFID};
				$novelORFInfoHsh_ref->{$ORFID}{'nonAnno'} = 'yes';
				$novelORFInfoHsh_ref->{$ORFID}{'uORF'} = 'no';
				$novelORFInfoHsh_ref->{$ORFID}{'geneWithuORF'} = 'no'; #---assume one ORFID will hit only 1 gene
				$novelORFInfoHsh_ref->{$ORFID}{'longestOvrlping'} = 'yes'; #---will be changed to no later if keepOvrlpNovelORF is no and it not the longest
			}

		} elsif ($hitAndPrxmtyByGetORFHsh_ref->{'SS'}{'hit'}{$ORFID}) {#---hit some genes in sense dirtn
			
			if ($getorfInfoHsh_ref->{$ORFID}{'CDSRng'}[1]-$getorfInfoHsh_ref->{$ORFID}{'CDSRng'}[0] >= $uORFMinNtLen) {
				foreach my $geneID (keys %{$hitAndPrxmtyByGetORFHsh_ref->{'SS'}{'hit'}{$ORFID}}) {
					if ($geneInfoHsh_ref->{$geneID}{'UTR5Rng'} and $geneInfoHsh_ref->{$geneID}{'CDSRng'}) {#---the genes has UTR5 annotated
						my ($scene, $UTR5OvrlpSize, undef) = &checkOvrlpTwoRangePair($getorfInfoHsh_ref->{$ORFID}{'CDSRng'}, $geneInfoHsh_ref->{$geneID}{'UTR5Rng'});#->534
						if ($scene eq 'within') {#---overlap with UTR5
							my (undef, $CDSOvrlpSize, undef) = &checkOvrlpTwoRangePair($getorfInfoHsh_ref->{$ORFID}{'CDSRng'}, $geneInfoHsh_ref->{$geneID}{'CDSRng'});#->534
							if ($CDSOvrlpSize == 0) {#---not overlaping with UTR5
								$novelORFInfoHsh_ref->{$ORFID} = $getorfInfoHsh_ref->{$ORFID};
								$novelORFInfoHsh_ref->{$ORFID}{'nonAnno'} = 'no';
								$novelORFInfoHsh_ref->{$ORFID}{'uORF'} = 'yes';
								$novelORFInfoHsh_ref->{$ORFID}{'geneWithuORF'} = $geneID; #---assume one ORFID will hit only 1 gene
								$novelORFInfoHsh_ref->{$ORFID}{'longestOvrlping'} = 'yes'; #---will be changed to no later if keepOvrlpNovelORF is no and it not the longest
							}
						}
					}
				}
			}
		}
	}

	my ($novelORFByCntgHsh_ref) = &generateGeneByCntgHsh($novelORFInfoHsh_ref);#->654

	return ($novelORFInfoHsh_ref, $novelORFByCntgHsh_ref);
}
sub getNovelORFCoverage {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: getCoverageOfItemRngType|702, reportStatus|1564
#	appearInSub: >none
#	primaryAppearInSection: 6_checkOverlappingAndCoverageOfNovelORF|164
#	secondaryAppearInSection: >none
#	input: $fullPileupStorablePathHsh_ref, $maxThread, $novelORFByCntgHsh_ref, $novelORFInfoHsh_ref, $resultStorableDir
#	output: $novelORFCovHsh_ref
#	toCall: my ($novelORFCovHsh_ref) = &getNovelORFCoverage($fullPileupStorablePathHsh_ref, $novelORFInfoHsh_ref, $novelORFByCntgHsh_ref, $resultStorableDir, $maxThread);
#	calledInLine: 171
#....................................................................................................................................................#
	my ($fullPileupStorablePathHsh_ref, $novelORFInfoHsh_ref, $novelORFByCntgHsh_ref, $resultStorableDir, $maxThread) = @_;
	
	my $novelORFCovHsh_ref = {};
	my $novelORFCovHshPlsPath = "$resultStorableDir/novelORFCovHsh.pls";
	
	if (-s $novelORFCovHshPlsPath) {
		
		&reportStatus("Retrieving novelORFCovHshPlsPath", 10, "\n");#->1564
		$novelORFCovHsh_ref = retrieve($novelORFCovHshPlsPath);
		
	} else {

		&reportStatus("Getting coverage on new ORF", 10, "\n");#->1564
		my $itemInfoHsh_ref = $novelORFInfoHsh_ref;
		my $itemByCntgHsh_ref = $novelORFByCntgHsh_ref;
		my $rngType = 'CDSRng';
		my $margin = 0;
		my $dirtn = 's';
		
		($novelORFCovHsh_ref) = &getCoverageOfItemRngType($fullPileupStorablePathHsh_ref, $itemInfoHsh_ref, $itemByCntgHsh_ref, $maxThread, $rngType, $margin, $dirtn);#->702
		
		store($novelORFCovHsh_ref, $novelORFCovHshPlsPath);
	}

	return ($novelORFCovHsh_ref);
}
sub getORFCodonPosCount {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: getCoverageOfItemRngType|702, getCtgryGeneInfo|796, reportStatus|1564
#	appearInSub: >none
#	primaryAppearInSection: 6_checkOverlappingAndCoverageOfNovelORF|164
#	secondaryAppearInSection: >none
#	input: $end5PileupStorablePathHsh_ref, $geneInfoHsh_ref, $maxThread, $novelORFByCntgHsh_ref, $novelORFInfoHsh_ref, $resultStorableDir
#	output: $codonPosCountHsh_ref
#	toCall: my ($codonPosCountHsh_ref) = &getORFCodonPosCount($end5PileupStorablePathHsh_ref, $novelORFInfoHsh_ref, $novelORFByCntgHsh_ref, $geneInfoHsh_ref, $maxThread, $resultStorableDir);
#	calledInLine: 175
#....................................................................................................................................................#
	my ($end5PileupStorablePathHsh_ref, $novelORFInfoHsh_ref, $novelORFByCntgHsh_ref, $geneInfoHsh_ref, $maxThread, $resultStorableDir) = @_;

	my $codonPosCountHsh_ref = {};
	my $codonPosCountHshPlsPath = "$resultStorableDir/codonPosCountHsh.pls";
	
	if (-s $codonPosCountHshPlsPath) {
		
		$codonPosCountHsh_ref = retrieve($codonPosCountHshPlsPath);
		
	} else {

		#----will process three types of item, novel ORF, mRNA and snRNA
		my %tmpItemHsh = ();
		($tmpItemHsh{'control'}{'itemInfoHsh_ref'}, $tmpItemHsh{'control'}{'itemByCntgHsh_ref'}) = &getCtgryGeneInfo($geneInfoHsh_ref, [qw/snRNA rRNA tRNA/]);#->796
		($tmpItemHsh{'mRNA'}{'itemInfoHsh_ref'}, $tmpItemHsh{'mRNA'}{'itemByCntgHsh_ref'}) = &getCtgryGeneInfo($geneInfoHsh_ref, [qw/mRNA/]);#->796
		$tmpItemHsh{'novelORF'}{'itemInfoHsh_ref'} = $novelORFInfoHsh_ref;
		$tmpItemHsh{'novelORF'}{'itemByCntgHsh_ref'} = $novelORFByCntgHsh_ref;
	
		foreach my $item (keys %tmpItemHsh) {

			&reportStatus("Counting codon position ratio of $item", 10, "\n");#->1564

			my $itemInfoHsh_ref = $tmpItemHsh{$item}{'itemInfoHsh_ref'};
			my $itemByCntgHsh_ref = $tmpItemHsh{$item}{'itemByCntgHsh_ref'};
			my $rngType = 'CDSRng';
			my $margin = 0;
			my $dirtn = 's';
			my ($covHsh_ref) = &getCoverageOfItemRngType($end5PileupStorablePathHsh_ref, $itemInfoHsh_ref, $itemByCntgHsh_ref, $maxThread, $rngType, $margin, $dirtn);#->702
			foreach my $geneID (keys %{$covHsh_ref}) {
				my $posFirst = 0;
				for (my $i = 0; $i < $#{$covHsh_ref->{$geneID}}; $i += 3) {
					$posFirst += $covHsh_ref->{$geneID}->[$i];
				}
				$codonPosCountHsh_ref->{$item}{$geneID}{'posFirst'} = $posFirst;
				$codonPosCountHsh_ref->{$item}{$geneID}{'posAll'} = sum(@{$covHsh_ref->{$geneID}});
				$codonPosCountHsh_ref->{$item}{$geneID}{'posFirst'} = $posFirst;
				$codonPosCountHsh_ref->{$item}{$geneID}{'firstPct'} = -1;
				$codonPosCountHsh_ref->{$item}{$geneID}{'firstPct'} = sprintf "%.2f", 100*$codonPosCountHsh_ref->{$item}{$geneID}{'posFirst'}/$codonPosCountHsh_ref->{$item}{$geneID}{'posAll'} if $codonPosCountHsh_ref->{$item}{$geneID}{'posAll'} > 0;
			}
		}
		store($codonPosCountHsh_ref, $codonPosCountHshPlsPath);
	}
	
	return ($codonPosCountHsh_ref);
}
sub outputAllGff {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: printGFF_oneRNAPerGene_chooseStrnd_filterAry|1194, reportStatus|1564
#	appearInSub: >none
#	primaryAppearInSection: 7_printTheTranscribeNovelORFLog|182
#	secondaryAppearInSection: >none
#	input: $getorfInfoHsh_ref, $novelORFInfoHsh_ref, $resultGFFDir
#	output: 
#	toCall: &outputAllGff($novelORFInfoHsh_ref, $getorfInfoHsh_ref, $resultGFFDir);
#	calledInLine: 192
#....................................................................................................................................................#
	my ($novelORFInfoHsh_ref, $getorfInfoHsh_ref, $resultGFFDir) = @_;
	
	foreach my $transcribed (qw/yes no/) {
		foreach my $longestOvrlping (qw/yes no/) {
			foreach my $ORFType (qw/nonAnno uORF both/) {
				my $outGFFPath = "$resultGFFDir/novelORF.transcribed_$transcribed.longest_$longestOvrlping.$ORFType.gff";
				&reportStatus("Printing gff $outGFFPath", 10, "\n");#->1564
				my $gffGeneLineOnly = 'no';
				my $strndInWord = 'both';
				my $filterAry_ref = [];
				if ($transcribed eq 'yes') {
					push @{$filterAry_ref}, 'transcribed';
				} else {
					$gffGeneLineOnly = 'yes';
				}
				push @{$filterAry_ref}, 'longestOvrlping' if $longestOvrlping eq 'yes';
				push @{$filterAry_ref}, $ORFType if $ORFType eq 'nonAnno' or $ORFType eq 'uORF';
				&printGFF_oneRNAPerGene_chooseStrnd_filterAry($novelORFInfoHsh_ref, $outGFFPath, $gffGeneLineOnly, $strndInWord, $filterAry_ref);#->1194
			}
		}
	}

	{
		my $outGFFPath = "$resultGFFDir/novelORF.all.gff";
		my $gffGeneLineOnly = 'yes';
		my $strndInWord = 'both';
		my $filterAry_ref = [];
		&printGFF_oneRNAPerGene_chooseStrnd_filterAry($novelORFInfoHsh_ref, $outGFFPath, $gffGeneLineOnly, $strndInWord, $filterAry_ref);#->1194
	}

	{
		my $outGFFPath = "$resultGFFDir/getorf.all.gff";
		my $gffGeneLineOnly = 'yes';
		my $strndInWord = 'both';
		my $filterAry_ref = [];
		&printGFF_oneRNAPerGene_chooseStrnd_filterAry($getorfInfoHsh_ref, $outGFFPath, $gffGeneLineOnly, $strndInWord, $filterAry_ref);#->1194
	}


	return ();
}
sub parseGetorfResults {
#....................................................................................................................................................#
#	subroutineCategory: thridPartyApp
#	dependOnSub: readMultiFasta|1472, reportStatus|1564
#	appearInSub: >none
#	primaryAppearInSection: 5_getAllPossibleORFs|154
#	secondaryAppearInSection: >none
#	input: $getorfFastaPath, $resultStorableDir
#	output: $getorfInfoHsh_ref
#	toCall: my ($getorfInfoHsh_ref) = &parseGetorfResults($getorfFastaPath, $resultStorableDir);
#	calledInLine: 159
#....................................................................................................................................................#
	my ($getorfFastaPath, $resultStorableDir) = @_;

	my $getorfInfoHsh_ref = {};
	my $getorfInfoHshPlsPath = "$resultStorableDir/getorfInfoHsh.pls";

	if (-s $getorfInfoHshPlsPath) {
	
		&reportStatus("Retrieving getorfInfoHsh", 0,"\n");#->1564
		$getorfInfoHsh_ref = retrieve($getorfInfoHshPlsPath);
		
	} else {

		&reportStatus("Parsing getorf results", 0,"\n");#->1564
		my ($getorfFastaHsh_ref) = &readMultiFasta($getorfFastaPath);#->1472

		foreach my $fastaSeqID (keys %{$getorfFastaHsh_ref}) {
			next if $getorfFastaHsh_ref->{$fastaSeqID} =~ m/X/ or $getorfFastaHsh_ref->{$fastaSeqID} !~ m/M/;
			
			#---trim aa before start codon
			my $aaSeq = $getorfFastaHsh_ref->{$fastaSeqID};
			
			if ($fastaSeqID =~ m/^(\S+)_(\d+)_\[(\d+)_-_(\d+)\]_([\(\)_\w]*)$/) {
				my ($start, $end, $strnd);
				if ($5 eq '(REVERSE_SENSE)_') {
					$strnd = '-';
					$start = $4 - 3;
					$end = $3;
				} else {
					$strnd = '+';
					$start = $3;
					$end = $4 + 3;
				}

				my $ORFID = "getorf_".$1."_".$strnd."_".$2;
				$getorfInfoHsh_ref->{$ORFID}{'aaSeq'} = $aaSeq;
				$getorfInfoHsh_ref->{$ORFID}{'cntg'} = $1;
				$getorfInfoHsh_ref->{$ORFID}{'strnd'} = $strnd;

				@{$getorfInfoHsh_ref->{$ORFID}{'geneRng'}} = sort {$a <=> $b} ($start, $end);
				@{$getorfInfoHsh_ref->{$ORFID}{'CDSRng'}} = @{$getorfInfoHsh_ref->{$ORFID}{'geneRng'}};
				@{$getorfInfoHsh_ref->{$ORFID}{'RNARng'}} = @{$getorfInfoHsh_ref->{$ORFID}{'geneRng'}};
				@{$getorfInfoHsh_ref->{$ORFID}{'exonRng'}} = @{$getorfInfoHsh_ref->{$ORFID}{'geneRng'}};
				$getorfInfoHsh_ref->{$ORFID}{"ctgry"} = "mRNA";
				$getorfInfoHsh_ref->{$ORFID}{"description"} = "getORF";
				$getorfInfoHsh_ref->{$ORFID}{"RNAID"} = "RNA_$ORFID";
				#print DEBUGLOG $aaSeq."\n";
			} else {
				die "unexpected getorf formatm quitting\n";
			}
		}
	
		store($getorfInfoHsh_ref, $getorfInfoHshPlsPath);
	}
	
	return ($getorfInfoHsh_ref);
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|636
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|79, 8_finishingTasks|197
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 85, 202
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->636
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->636
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->636
		print "=========================================================================\n\n";
	}
}
sub printGFF_oneRNAPerGene_chooseStrnd_filterAry {
#....................................................................................................................................................#
#	subroutineCategory: gff
#	dependOnSub: >none
#	appearInSub: outputAllGff|1039
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 7_printTheTranscribeNovelORFLog|182
#	input: $filterAry_ref, $geneInfoHsh_ref, $gffGeneLineOnly, $outGFFPath, $strndInWord
#	output: none
#	toCall: &printGFF_oneRNAPerGene_chooseStrnd_filterAry($geneInfoHsh_ref, $outGFFPath, $gffGeneLineOnly, $strndInWord, $filterAry_ref);
#	calledInLine: 1067, 1077, 1085
#....................................................................................................................................................#

	my ($geneInfoHsh_ref, $outGFFPath, $gffGeneLineOnly, $strndInWord, $filterAry_ref) = @_;

	#------gffGeneLineOnly: print only the gene line, which is perfect for trnsfrgs
	#------strndInWord: plus, minus or both
		
	$outGFFPath =~ s/\/+/\//g;
	
	my $strndOut = '+';
	$strndOut = '-' if $strndInWord eq 'minus';
	
	open (GFFOUT, ">", "$outGFFPath");
	print GFFOUT "##gff-version\t3\n";

	foreach my $geneID (sort {$a cmp $b} keys %{$geneInfoHsh_ref}) {
		
		my $toPrint = 'yes';
		foreach my $filter (@{$filterAry_ref}) {
			$toPrint = 'no' if ($geneInfoHsh_ref->{$geneID}{$filter}) eq 'no';
		}
		
		if ($toPrint eq 'yes') {
		
			my $strnd = $geneInfoHsh_ref->{$geneID}{"strnd"};

			next if $strndOut ne $strnd and $strndInWord ne 'both';
		
			my $cntg = $geneInfoHsh_ref->{$geneID}{"cntg"};
			my $ctgry = $geneInfoHsh_ref->{$geneID}{"ctgry"};
			my $description = $geneInfoHsh_ref->{$geneID}{"description"};
		
			my ($geneStart, $geneEnd) = sort {$a <=> $b} @{$geneInfoHsh_ref->{$geneID}{'geneRng'}};
			print GFFOUT join "", (join "\t", ($cntg, 'BCP', 'gene', $geneStart, $geneEnd, ".", $strnd, ".", "ID=$geneID;Name=$geneID;description=$description")), "\n";

			if ($gffGeneLineOnly eq 'no') {#-----will print the RNA and exon lines also, aim to avoid display annoyance on IGV
				my $RNAID = $geneInfoHsh_ref->{$geneID}{"RNAID"};
				my ($RNAStart, $RNAEnd) = sort {$a <=> $b} @{$geneInfoHsh_ref->{$geneID}{'RNARng'}};
				print GFFOUT join "", (join "\t", ($cntg, 'BCP', $ctgry, $RNAStart, $RNAEnd, ".", $strnd, ".", "ID=$RNAID;Name=$geneID;Parent=$geneID;description=$description;")), "\n";
			
				foreach my $rngType ('exonRng', 'CDSRng') {
					my $gffRng = 'exon';
					$gffRng = 'CDS' if $rngType eq 'CDSRng';
				
					if ($geneInfoHsh_ref->{$geneID}{$rngType}) {
						my @rngAry = sort {$a <=> $b} @{$geneInfoHsh_ref->{$geneID}{$rngType}};
						my $num = 1;
						for (my $i=0; $i < $#rngAry; $i += 2) {
							$num++;
							my ($start, $end) = ($rngAry[$i], $rngAry[$i+1]);
							my $ID = "$rngType\_$num\_$RNAID";
							print GFFOUT join "", (join "\t", ($cntg, "BCP", $gffRng, $start, $end, ".", $strnd, ".", "ID=$ID;Name=$ID;Parent=$RNAID;description=.;")), "\n";
						}
					}
				} #---end of foreach my $rngType ('exonRng', 'CDSRng') 
			} #---end of if ($gffGeneLineOnly eq 'no') 
		} #---end of if ($toPrint eq 'yes')
	}
	close GFFOUT;
}
sub printGeneBaseduORFInfo {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 7_printTheTranscribeNovelORFLog|182
#	secondaryAppearInSection: >none
#	input: $fastaHsh_ref, $geneInfoHsh_ref, $novelORFInfoHsh_ref, $resultLogDir
#	output: 
#	toCall: &printGeneBaseduORFInfo($novelORFInfoHsh_ref, $resultLogDir, $geneInfoHsh_ref, $fastaHsh_ref);
#	calledInLine: 191
#....................................................................................................................................................#
	my ($novelORFInfoHsh_ref, $resultLogDir, $geneInfoHsh_ref, $fastaHsh_ref) = @_;
	
	
	#---set set uORFNum and uORFTranscribed to 0
	foreach my $geneID (sort {$a cmp $b} keys %{$geneInfoHsh_ref}) {
		if ($geneInfoHsh_ref->{$geneID}{'UTR5Rng'}) {
			$geneInfoHsh_ref->{$geneID}{'uORFNum'} = 0;
			$geneInfoHsh_ref->{$geneID}{'uORFTranscribed'} = 0;
			$geneInfoHsh_ref->{$geneID}{'maxuORFLen'} = 0;
			$geneInfoHsh_ref->{$geneID}{'avguORFLen'} = 0;
		}
	}
	
	foreach my $ORFID (sort {$a cmp $b} keys %{$novelORFInfoHsh_ref}) {
		if ($novelORFInfoHsh_ref->{$ORFID}{'geneWithuORF'} ne 'no') {
			my $geneID = $novelORFInfoHsh_ref->{$ORFID}{'geneWithuORF'};
			$geneInfoHsh_ref->{$geneID}{'uORFNum'}++;
			$geneInfoHsh_ref->{$geneID}{'uORFTranscribed'}++ if ($novelORFInfoHsh_ref->{$ORFID}{'transcribed'} eq 'yes');
			my $uORFLength = $novelORFInfoHsh_ref->{$ORFID}{'CDSRng'}[1] - $novelORFInfoHsh_ref->{$ORFID}{'CDSRng'}[0];
			push @{$geneInfoHsh_ref->{$geneID}{'uORFLength'}}, $uORFLength;
			$geneInfoHsh_ref->{$geneID}{'maxuORFLen'} = max @{$geneInfoHsh_ref->{$geneID}{'uORFLength'}};
			$geneInfoHsh_ref->{$geneID}{'avguORFLen'} = sum (@{$geneInfoHsh_ref->{$geneID}{'uORFLength'}})/@{$geneInfoHsh_ref->{$geneID}{'uORFLength'}};
		}
	}

	open UORFINFO, ">", "$resultLogDir/gene.based.uORF.info.xls";
	my @outputAry = ('geneID', 'locationTag', 'cntg', 'strnd', 'geneLen', 'UTR5Len', 'uORFNum', 'uORFTranscribed', 'maxuORFLen', 'avguORFLen');
	print UORFINFO join "", ((join "\t", @outputAry), "\n");
	foreach my $geneID (sort {$a cmp $b} keys %{$geneInfoHsh_ref}) {
		my $cntg = $geneInfoHsh_ref->{$geneID}{'cntg'};
		next if not $fastaHsh_ref->{$cntg};
		my ($geneStart, $geneEnd) = ($geneInfoHsh_ref->{$geneID}{'geneRng'}[0], $geneInfoHsh_ref->{$geneID}{'geneRng'}[1]);
		my $locationTag = $cntg.":".$geneStart."-".$geneEnd;
		my $geneLen = $geneEnd - $geneStart;
		my $strnd = $geneInfoHsh_ref->{$geneID}{'strnd'};
		my $UTR5Len = -1;
		my $uORFNum = -1;
		my $uORFTranscribed = -1;
		my $maxuORFLen = -1;
		my $avguORFLen = -1;
		if ($geneInfoHsh_ref->{$geneID}{'UTR5Rng'}) {
			$UTR5Len = $geneInfoHsh_ref->{$geneID}{'UTR5Rng'}[1] - $geneInfoHsh_ref->{$geneID}{'UTR5Rng'}[0];
			$uORFNum = $geneInfoHsh_ref->{$geneID}{'uORFNum'};
			$uORFTranscribed = $geneInfoHsh_ref->{$geneID}{'uORFTranscribed'};
			$maxuORFLen = $geneInfoHsh_ref->{$geneID}{'maxuORFLen'};
			$avguORFLen = $geneInfoHsh_ref->{$geneID}{'avguORFLen'};
		}
		my @outputAry = ($geneID, $locationTag, $cntg, $strnd, $geneLen, $UTR5Len, $uORFNum, $uORFTranscribed, $maxuORFLen, $avguORFLen);
		print UORFINFO join "", ((join "\t", @outputAry), "\n");
	}
	close UORFINFO;
	
	return ();
}
sub printTranscribedNovelORFInfo {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 7_printTheTranscribeNovelORFLog|182
#	secondaryAppearInSection: >none
#	input: $codonPosCountHsh_ref, $novelORFInfoHsh_ref, $resultLogDir
#	output: 
#	toCall: &printTranscribedNovelORFInfo($novelORFInfoHsh_ref, $codonPosCountHsh_ref, $resultLogDir);
#	calledInLine: 187
#....................................................................................................................................................#
	my ($novelORFInfoHsh_ref, $codonPosCountHsh_ref, $resultLogDir) = @_;
	
	open NOVELORFIFNO, ">", "$resultLogDir/transcribed.novel.ORF.info.xls";
	my @outputAry = (qw/ORFID locationTag cntg CDSStart CDSEnd strnd nonAnno uORF geneWithuORF longestOvrlping aaLength covPerNt covPosPct aaSeq posFirst posAll firstPct/);
	print NOVELORFIFNO join "", ((join "\t", @outputAry), "\n");
	foreach my $ORFID (sort {$a cmp $b} keys %{$novelORFInfoHsh_ref}) {
		next if $novelORFInfoHsh_ref->{$ORFID}{'transcribed'} eq 'no';
		my $cntg = $novelORFInfoHsh_ref->{$ORFID}{'cntg'};
		my ($CDSStart, $CDSEnd) = @{$novelORFInfoHsh_ref->{$ORFID}{'CDSRng'}};
		my $locationTag = $cntg.":".$CDSStart."-".$CDSEnd;
		my $strnd = $novelORFInfoHsh_ref->{$ORFID}{'strnd'};
		my $nonAnno = $novelORFInfoHsh_ref->{$ORFID}{'nonAnno'};
		my $uORF = $novelORFInfoHsh_ref->{$ORFID}{'uORF'};
		my $geneWithuORF = $novelORFInfoHsh_ref->{$ORFID}{'geneWithuORF'};
		my $longestOvrlping = $novelORFInfoHsh_ref->{$ORFID}{'longestOvrlping'};
		my $covPerNt = $novelORFInfoHsh_ref->{$ORFID}{'covPerNt'};
		my $covPosPct = $novelORFInfoHsh_ref->{$ORFID}{'covPosPct'};
		my $aaSeq = $novelORFInfoHsh_ref->{$ORFID}{'aaSeq'};
		my $aaLength = length $aaSeq;
		my $posFirst = $codonPosCountHsh_ref->{'novelORF'}{$ORFID}{'posFirst'};
		my $posAll = $codonPosCountHsh_ref->{'novelORF'}{$ORFID}{'posAll'};
		my $firstPct = $codonPosCountHsh_ref->{'novelORF'}{$ORFID}{'firstPct'};

		my @outputAry = ($ORFID, $locationTag, $cntg, $CDSStart, $CDSEnd, $strnd, $nonAnno, $uORF, $geneWithuORF, $longestOvrlping, $aaLength, $covPerNt, $covPosPct, $aaSeq, $posFirst, $posAll, $firstPct);
		print NOVELORFIFNO join "", ((join "\t", @outputAry), "\n");
	}
	close NOVELORFIFNO;
	
	return ();
}
sub readGFF_oneRNAPerGene {
#....................................................................................................................................................#
#	subroutineCategory: gff
#	dependOnSub: currentTime|636
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|134
#	secondaryAppearInSection: >none
#	input: $gffPath
#	output: $geneInfoHsh_ref
#	toCall: my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);
#	calledInLine: 142
#....................................................................................................................................................#

	my ($gffPath) = @_;

	my $geneInfoHsh_ref = {};
	
	#---read the gff
	my $geneByRNAHsh_ref = {};

	open (GFF, $gffPath);
	print "[".&currentTime()."] Reading: $gffPath\n";#->636
	while (my $theLine = <GFF>) {

		chomp $theLine;
		
		last if $theLine =~ m/^##FASTA/;
		
		if ($theLine !~ m/^\#|^\@/ and $theLine !~ m/\tsupercontig\t/) {

			my ($seq, undef, $geneCategory, $featureStart, $featureEnd, undef, $geneStrd, undef, $dscrptns) = split (/\t/, $theLine);
			
			#----assigne all non -/+ will be treated as plus
			$geneStrd = "+" if (($geneStrd ne "-") and ($geneStrd ne "+"));
			
			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent);
			my $geneName = "unknown";
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^description=/) {$geneName = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}

			if ($geneCategory eq "gene") {#---gene
				
				my $geneID = $unqID;
				
				$geneInfoHsh_ref->{$geneID}{'strnd'} = $geneStrd;
				$geneInfoHsh_ref->{$geneID}{'cntg'} = $seq;
				$geneInfoHsh_ref->{$geneID}{'description'} = uri_unescape($geneName);
				$geneInfoHsh_ref->{$geneID}{'description'} =~ s/\+/ /g;
				@{$geneInfoHsh_ref->{$geneID}{'geneRng'}} = ($featureStart, $featureEnd);

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes
				
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'CDSRng'}}, ($featureStart, $featureEnd);
				
			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'exonRng'}}, ($featureStart, $featureEnd);
				
			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh_ref->{$RNAID} = $geneID;
				$geneInfoHsh_ref->{$geneID}{'ctgry'} = $geneCategory;
				@{$geneInfoHsh_ref->{$geneID}{'RNARng'}} = ($featureStart, $featureEnd);
				$geneInfoHsh_ref->{$geneID}{'RNAID'} = $RNAID;
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close GFF;
	
	#---get the UTR if any
	my $minUTRLength = 10;
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {

		if (exists $geneInfoHsh_ref->{$geneID}{'CDSRng'}) {
			my $exonMin = min(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $exonMax = max(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $CDSMin = min(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});
			my $CDSMax = max(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});

			if ($geneInfoHsh_ref->{$geneID}{'strnd'} eq '+') {
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			} else {
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			}
		}
	}
	
	return ($geneInfoHsh_ref);
}
sub readMultiFasta {
#....................................................................................................................................................#
#	subroutineCategory: fasta
#	dependOnSub: reportStatus|1564
#	appearInSub: parseGetorfResults|1092
#	primaryAppearInSection: 4_processInputData|134, 7_printTheTranscribeNovelORFLog|182
#	secondaryAppearInSection: 5_getAllPossibleORFs|154
#	input: $fastaPath
#	output: $fastaHsh_ref
#	toCall: my ($fastaHsh_ref) = &readMultiFasta($fastaPath);
#	calledInLine: 190, 1116
#....................................................................................................................................................#

	my ($fastaPath) = @_;

	my ($seq, $seqName);
	my $fastaHsh_ref = {};
	my $i = 0;

	&reportStatus("Reading: $fastaPath", 0, "\n");#->1564
	
	open (INFILE, $fastaPath);
	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/ *\|/, $curntLine);
			$seqName = $theLineSplt[0]; #---get the first tag
			$seqName =~ s/ /_/g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$curntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$seq =~ tr/a-z/A-Z/;
			$fastaHsh_ref->{$seqName} = $seq;
			$seq = "";
		} elsif (eof(INFILE)) {#---this is the last line
			$seq =~ tr/a-z/A-Z/;
			$seq = $seq.$nextLine;
			$fastaHsh_ref->{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}

	close INFILE;
	return ($fastaHsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|79
#	secondaryAppearInSection: >none
#	input: none
#	output: $end5PileupStorableIndexPath, $fastaPath, $fullPileupStorableIndexPath, $getorfFastaPath, $gffPath, $maxThread, $outDir
#	toCall: my ($gffPath, $fastaPath, $getorfFastaPath, $fullPileupStorableIndexPath, $maxThread, $end5PileupStorableIndexPath, $outDir) = &readParameters();
#	calledInLine: 87
#....................................................................................................................................................#
	
	my ($gffPath, $fastaPath, $getorfFastaPath, $fullPileupStorableIndexPath, $maxThread, $end5PileupStorableIndexPath, $outDir);
	
	$maxThread = 4;
	my $dirPath = dirname(rel2abs($0));
	$outDir = "$dirPath/RNASeqORFFinder/";
	
	GetOptions 	("gffPath=s"  => \$gffPath,
				 "fastaPath=s"  => \$fastaPath,
				 "getorfFastaPath=s"  => \$getorfFastaPath,
				 "fullPileupStorableIndexPath=s"  => \$fullPileupStorableIndexPath,
				 "end5PileupStorableIndexPath=s"  => \$end5PileupStorableIndexPath,
				 "maxThread:i"  => \$maxThread,
				 "outDir:s"  => \$outDir) 

	or die		("Error in command line arguments\n");
	
	#---check file
	foreach my $fileToCheck ($gffPath, $fastaPath, $getorfFastaPath, $fullPileupStorableIndexPath, $end5PileupStorableIndexPath) {
		die "Can't read $fileToCheck" if not -s $fileToCheck;
	}

	system "mkdir -p -m 777 $outDir/";
	
	return($gffPath, $fastaPath, $getorfFastaPath, $fullPileupStorableIndexPath, $maxThread, $end5PileupStorableIndexPath, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|636
#	appearInSub: checkGeneInfo|269, checkOverlapAndProximity_withMargin|297, generateGeneByCntgHsh|654, getCoverageOfItemRngType|702, getCtgryGeneInfo|796, getIndivCntgCovPlsPath|833, getNonOverlappingORF|866, getNovelORFCoverage|944, getORFCodonPosCount|982, outputAllGff|1039, parseGetorfResults|1092, readMultiFasta|1472, tagOverlappingShorterTranscribedORF|1585, zipUnzipCntgCovInPlsPathHsh|1653
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_processInputData|134, 5_getAllPossibleORFs|154, 6_checkOverlappingAndCoverageOfNovelORF|164, 7_printTheTranscribeNovelORFLog|182
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 283, 291, 346, 360, 371, 667, 730, 765, 790, 815, 828, 861, 884, 887, 898, 962, 967, 1013, 1056, 1110, 1115, 1490, 1607, 1668
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->636

	return ();
}
sub tagOverlappingShorterTranscribedORF {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: checkOverlapAndProximity_withMargin|297, reportStatus|1564
#	appearInSub: >none
#	primaryAppearInSection: 6_checkOverlappingAndCoverageOfNovelORF|164
#	secondaryAppearInSection: >none
#	input: $keepOvrlpNovelORF, $maxThread, $novelORFInfoHsh_ref
#	output: 
#	toCall: &tagOverlappingShorterTranscribedORF($keepOvrlpNovelORF, $novelORFInfoHsh_ref, $maxThread);
#	calledInLine: 177
#....................................................................................................................................................#
	my ($keepOvrlpNovelORF, $novelORFInfoHsh_ref, $maxThread) = @_;
	
	#---remove overlapping novel ORF if keepOvrlpNovelORF eq no
	
	if ($keepOvrlpNovelORF eq 'no') {

		my $transcribedORFInfoHsh_ref = {};
		foreach (keys %{$novelORFInfoHsh_ref}) {#---check only transcribed, otherwise would take ages
			$transcribedORFInfoHsh_ref->{$_} = $novelORFInfoHsh_ref->{$_} if ($novelORFInfoHsh_ref->{$_}{'transcribed'} eq 'yes');
		}
		
		&reportStatus("Removing overlapping novel ORFs", 0,"\n");#->1564
		my $checkPrxmty = 'no';
		my $reportExactMatch = 'no';
		my $refRngType = 'CDSRng';
		my $qryRngType = 'CDSRng';
		my $refMargin = 0;
		my $qryMargin = 0;
		my ($hitAndPrxmtySelfHitGetORFHsh_ref, undef) = &checkOverlapAndProximity_withMargin($transcribedORFInfoHsh_ref, $transcribedORFInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType, $refMargin, $qryMargin);#->297
		
		my %checkedIDHsh = ();
		foreach my $refORFID (keys %{$novelORFInfoHsh_ref}) {
			next if $checkedIDHsh{$refORFID};
			if ($hitAndPrxmtySelfHitGetORFHsh_ref->{'SS'}{'hit'}{$refORFID}) {
				my $longestLength = 0;
				my %tmpORFLenHsh = ();
				my $refLen = $novelORFInfoHsh_ref->{$refORFID}{'CDSRng'}[1]-$novelORFInfoHsh_ref->{$refORFID}{'CDSRng'}[0];
				$tmpORFLenHsh{$refORFID} = $refLen;
				$longestLength = $refLen if $refLen > $longestLength;
				foreach my $qryORFID (keys %{$hitAndPrxmtySelfHitGetORFHsh_ref->{'SS'}{'hit'}{$refORFID}}) {
					my $qryLen = $novelORFInfoHsh_ref->{$qryORFID}{'CDSRng'}[1]-$novelORFInfoHsh_ref->{$qryORFID}{'CDSRng'}[0];
					$tmpORFLenHsh{$qryORFID} = $qryLen;
					$longestLength = $qryLen if $qryLen > $longestLength;
				}
				my @longestOvrlpingIDAry = (); #---to avoid ORF with same longest length
				foreach (keys %tmpORFLenHsh) {
					$checkedIDHsh{$_}++;
					if ($tmpORFLenHsh{$_} < $longestLength) {
						$novelORFInfoHsh_ref->{$_}{'longestOvrlping'} = 'no';
					} else {
						push @longestOvrlpingIDAry, $_;
					}
				}
				
				#---i.e. more than one ORF with same longest length
				if (@longestOvrlpingIDAry > 1) {
					for my $i (1..$#longestOvrlpingIDAry) {#---only the ID in index 0 will keep the yes
						$novelORFInfoHsh_ref->{$longestOvrlpingIDAry[$i]}{'longestOvrlping'} = 'no';
					}
				}
			}
		}
	}

	return ();
}
sub zipUnzipCntgCovInPlsPathHsh {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: reportStatus|1564
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|134
#	secondaryAppearInSection: >none
#	input: $cntgCovInPlsPathHsh_ref, $zipUnzip
#	output: none
#	toCall: &zipUnzipCntgCovInPlsPathHsh($zipUnzip, $cntgCovInPlsPathHsh_ref);
#	calledInLine: 148, 149
#....................................................................................................................................................#

	my ($zipUnzip, $cntgCovInPlsPathHsh_ref) = @_;
	
	foreach my $cntg (sort keys %{$cntgCovInPlsPathHsh_ref}) {
		&reportStatus("Trying to $zipUnzip cntg ary", 20, "\r");#->1564

		my $cntgCovPlsPath = "$cntgCovInPlsPathHsh_ref->{$cntg}";

		if ($zipUnzip eq 'unzip') {
			system ("gzip -df $cntgCovPlsPath.gz") if (-s "$cntgCovPlsPath.gz");
		} else {
			system ("gzip -f $cntgCovPlsPath") if (-s "$cntgCovPlsPath");
		}
	}
	print "\n";
}

exit;
