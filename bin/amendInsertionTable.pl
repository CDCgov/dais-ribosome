#!/usr/bin/env perl
# Sam Shepard

use Getopt::Long;
GetOptions(	'use-std-out|U' => \$useSTDOUT, 'split-id-field|T' => \$splitID 
	);
if ( scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <insertion_table.txt> [-U|--use-std-out] [-T|--split-id-field]\n";
	die($message."\n");
}

$insertionTable = $ARGV[0]; $/ = "\n"; %inserts = ();
open(INS,'<',$insertionTable) or die("Cannot open $insertionTable for reading.\n");
@lines = <INS>; chomp(@lines);
foreach $line ( @lines ) {
	@fields = split("\t",$line);
	if ( scalar(@fields) < 3 || scalar(@fields) > 4 ) {
		die("Expected 3 fields in this format:\n\tID<tab>POSITION<tab>INSERT\n");
	} elsif ( scalar(@fields) == 4 ) {
		die("Likely already translated.\n");
	} else {
		($id,$pos,$insert) = @fields;
		$inserts{$id}{$pos} = lc($insert);
	}
}
close(INS);

# augmented translation table
%gc = (
	'TAA'=>'*','TAG'=>'*','TAR'=>'*','TGA'=>'*','TRA'=>'*','GCA'=>'A','GCB'=>'A','GCC'=>'A','GCD'=>'A','GCG'=>'A','GCH'=>'A',
	'GCK'=>'A','GCM'=>'A','GCN'=>'A','GCR'=>'A','GCS'=>'A','GCT'=>'A','GCV'=>'A','GCW'=>'A','GCY'=>'A','TGC'=>'C','TGT'=>'C',
	'TGY'=>'C','GAC'=>'D','GAT'=>'D','GAY'=>'D','GAA'=>'E','GAG'=>'E','GAR'=>'E','TTC'=>'F','TTT'=>'F','TTY'=>'F','GGA'=>'G',
	'GGB'=>'G','GGC'=>'G','GGD'=>'G','GGG'=>'G','GGH'=>'G','GGK'=>'G','GGM'=>'G','GGN'=>'G','GGR'=>'G','GGS'=>'G','GGT'=>'G',
	'GGV'=>'G','GGW'=>'G','GGY'=>'G','CAC'=>'H','CAT'=>'H','CAY'=>'H','ATA'=>'I','ATC'=>'I','ATH'=>'I','ATM'=>'I','ATT'=>'I',
	'ATW'=>'I','ATY'=>'I','AAA'=>'K','AAG'=>'K','AAR'=>'K','CTA'=>'L','CTB'=>'L','CTC'=>'L','CTD'=>'L','CTG'=>'L','CTH'=>'L',
	'CTK'=>'L','CTM'=>'L','CTN'=>'L','CTR'=>'L','CTS'=>'L','CTT'=>'L','CTV'=>'L','CTW'=>'L','CTY'=>'L','TTA'=>'L','TTG'=>'L',
	'TTR'=>'L','YTA'=>'L','YTG'=>'L','YTR'=>'L','ATG'=>'M','AAC'=>'N','AAT'=>'N','AAY'=>'N','CCA'=>'P','CCB'=>'P','CCC'=>'P',
	'CCD'=>'P','CCG'=>'P','CCH'=>'P','CCK'=>'P','CCM'=>'P','CCN'=>'P','CCR'=>'P','CCS'=>'P','CCT'=>'P','CCV'=>'P','CCW'=>'P',
	'CCY'=>'P','CAA'=>'Q','CAG'=>'Q','CAR'=>'Q','AGA'=>'R','AGG'=>'R','AGR'=>'R','CGA'=>'R','CGB'=>'R','CGC'=>'R','CGD'=>'R',
	'CGG'=>'R','CGH'=>'R','CGK'=>'R','CGM'=>'R','CGN'=>'R','CGR'=>'R','CGS'=>'R','CGT'=>'R','CGV'=>'R','CGW'=>'R','CGY'=>'R',
	'MGA'=>'R','MGG'=>'R','MGR'=>'R','AGC'=>'S','AGT'=>'S','AGY'=>'S','TCA'=>'S','TCB'=>'S','TCC'=>'S','TCD'=>'S','TCG'=>'S',
	'TCH'=>'S','TCK'=>'S','TCM'=>'S','TCN'=>'S','TCR'=>'S','TCS'=>'S','TCT'=>'S','TCV'=>'S','TCW'=>'S','TCY'=>'S','ACA'=>'T',
	'ACB'=>'T','ACC'=>'T','ACD'=>'T','ACG'=>'T','ACH'=>'T','ACK'=>'T','ACM'=>'T','ACN'=>'T','ACR'=>'T','ACS'=>'T','ACT'=>'T',
	'ACV'=>'T','ACW'=>'T','ACY'=>'T','GTA'=>'V','GTB'=>'V','GTC'=>'V','GTD'=>'V','GTG'=>'V','GTH'=>'V','GTK'=>'V','GTM'=>'V',
	'GTN'=>'V','GTR'=>'V','GTS'=>'V','GTT'=>'V','GTV'=>'V','GTW'=>'V','GTY'=>'V','TGG'=>'W','TAC'=>'Y','TAT'=>'Y','TAY'=>'Y'
);

if ( $useSTDOUT ) {
	*INS = *STDOUT;
} else {
	open(INS,'>',$insertionTable) or die("Cannot open $insertionTable for writing.\n");
}
foreach $id ( sort { $a cmp $b } keys(%inserts) ) {
	foreach $pos ( sort { $a <=> $b } keys(%{$inserts{$id}}) ) {
		$insert = $inserts{$id}{$pos};
		$L = length($insert);
		if ( $L % 3 == 0 ) {
			$aa = '';
			for($i=0;$i<$L;$i+=3) {
				$codon = uc(substr($insert,$i,3));
				$aa .= defined($gc{$codon}) ? $gc{$codon} : 'X';
			}
		} else {	
			$aa = "?" x ($L%3);
		}

		if ( $pos % 3 != 0 ) {
			$pos = int($pos/3) + ($pos%3) * 0.3;
		} else {
			$pos = int($pos/3);
		}

		if ( defined($splitID) ) { $id = join("\t",split('\|',$id)); }
		print INS $id,"\t",$pos,"\t",$insert,"\t$aa\n";
	}
}
close(INS);
