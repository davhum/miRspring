###################################################################################################################################
# Intermediate_to_miRspring.pl
#      This program converts the intermediate file into a miRspring document.  
# 	For more information please visit http://mirspring.victorchang.edu.au.
#
#      The following command will list the options available with this script:
#               perl Intermediate_to_miRspring.pl
#
# Version 1.2
#		- Added non-template addition feature (up to 2nt)
#		- Script settings are now saved and displayed in final miRspring document.
#			
#	
# Author: David Thomas Humphreys
# Copyright (C) 2013, Victor Chang Cardiac Research Institute
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    see <http://www.gnu.org/licenses/> for GNU license.

####################################################################################################################################


#!/usr/bin/perl

use strict;


my $parameters = {};

####------SETUP A-------------------------------------------------------------------------------------------------------------------------
#### Before running the scripts it is recommended to set the following FOUR global variables.

## The three variables below should contain the DIRECTORY of the miRbase files.
my $Precursor_GFF_Dir	= "ENTER DIRECTORY where gff file resides";
my $Precursor_Seq_Dir	= "ENTER DIRECTORY where hairpin.fa (or equivalent) resides";
my $Mature_Seq_Dir		= "ENTER DIRECTORY where mature.fa resides";
my $JavaScriptTemplateFile      = "./javascriptTEMPLATE.txt";	# Set the full path so that you do not have to execute this script within the installed directory.

####-------------------------------------------------------------------------------------------------------------------------------


my $HeaderInfo = "var HeaderInfo = new Object();\n";	# Will store header information PLUS script settings for miRspring document.

sub load_default_reference_files()
{
####------SETUP B -------------------------------------------------------------------------------------------------------------------------
### Before running the scripts recommend to set the default mirbase file names

#$parameters->{gff} = "$Precursor_GFF_Dir/ ENTER FILENAME HERE, if using default miRbase files use: $parameters->{species}.gff2", if (! defined $parameters->{gff});
#$parameters->{precursor_seq_file} = "$Precursor_Seq_Dir/ ENTER FILENAME HERE, if using miRBase files use: hairpin.fa", if (! defined $parameters->{precursor_seq_file});
#$parameters->{mature_seq_file} = "$Mature_Seq_Dir/ ENTER FILENAME HERE, if using miRBase files use: mature.fa", if (! defined $parameters->{mature_seq_file});
####-------------------------------------------------------------------------------------------------------------------------------

	#    $miRspringHeader   .= "GFF coordinate file (-gff): ".$parameters->{gff}."\t <br>";
	#    $miRspringHeader   .= "Mature sequence file (-mat): ".$parameters->{mature_seq_file}."\t<br>";
	#    $miRspringHeader   .= "Hairpin precursor file (-pre): ".$parameters->{precursor_seq_file}.".<br>";
	if ((! defined  $parameters->{gff}) || (! defined  $parameters->{mature_seq_file}) || (! defined  $parameters->{percursor_seq_file}))
	{	# Search through input file to see if the required files were previously declared and saved
		my @Check;
		$Check[0] = 1, if (! defined  $parameters->{gff});
		$Check[1] = 1, if (! defined $parameters->{mature_seq_file});
		$Check[2] = 1, if (! defined $parameters->{precursor_seq_file});
		open MIRSPRING_DATA, $parameters->{input} or  die "Cannot open the raw miRspring data file ($parameters->{input})";
		while (<MIRSPRING_DATA>)
		{
			if ($_ =~ m/^\@MIRSPRING/)
			{	chomp();
				$HeaderInfo .= "\nHeaderInfo.MIRSPRING = \"$_\";";
				if ($_ =~ m/-gff\): (.*?)\t </)
				{	$parameters->{gff} = $1, if (! defined  $parameters->{gff});	
					$Check[0]--;
				}
				if ($_ =~ m/-mat\): (.*?)\t<br>/)
				{	$parameters->{mature_seq_file} = $1, if (! defined  $parameters->{mature_seq_file});	
					$Check[1]--;
				}
				if ($_ =~ m/-pre\): (.*?)\.</)
				{	$parameters->{precursor_seq_file} = $1, (! defined  $parameters->{percursor_seq_file});	
					$Check[2]--;
				}
				last, if (($Check[0] < 1) && ($Check[1] < 1) && ($Check[2] < 1)); # Save reading through the whole file. Also means first entry gets saved over later entries
			}
		}
		close MIRSPRING_DATA;
	}


	$parameters->{mintag} = 0, if (! defined $parameters->{mintag});
	$parameters->{tagfreq} = 0, if (! defined $parameters->{tagfreq});	
	$parameters->{starttagfreq} = 0, if (! defined $parameters->{starttagfreq});	
	$parameters->{flank} = 0, if (! defined $parameters->{flank});
	$parameters->{output} =  "$parameters->{input}".".html", if (! defined $parameters->{output});
	$parameters->{mbase} = 3, if (! defined $parameters->{mbase});
	$parameters->{mismatch} = 1, if (! defined $parameters->{mismatch});
	$parameters->{showflank} = 0, if (! defined $parameters->{showflank});
	$parameters->{seedclass} = 0, if (! defined $parameters->{seedclass});
	$parameters->{normalisation} = 0, if (! defined $parameters->{normalisation});
}



sub usage {
	print "\nUsage: $0 \n\n\t ";

	print "REQUIRED \n\t";
	print "-in <input file> \n\t";
	print "-s <three letter species code>\n\t";

	print "OPTIONAL \n\t";
	print "-out <output file name with full path> \n\n\t";
	print "-mm <mismatches: '0' or '1', default is 1> \n\t";
	print "-gff <Gene features file (gff), default is defined in script> \n\t";
	print "-mat <Mature sequence file, default is defined in script> \n\t";
	print "-pre <Precursor file, default is defined in script> \n\t";
	print "-flank <length of flanking sequence>\n\t";
	print "-ref <format of reference: \'0\' (genome) or \'1\' (custom), 0 (genome) is default\n\t";
	print "-mintag <only include precursors with at least this number of tags, default = 0>\n\t";
	print "-tagfreq <minimum sequence count for all unique sequence tags, default = 0>\n\t";
	print "-starttagfreq <minimum sequence count for all unique sequence tags, default = 0>\n\t";
	print "-comment <message for top of miRspring document>\n\t";
	print "-mbase <number of nucleotides to extend up and downstream to accept tag as a processed miRNA, default = 3>\n\t";
	print "-mismatch <number of nucleotide mismatches to display, default = 1>\n\t";
	print "-showflank <set the miRspring document to display the flanking sequence as default. No = 0, Yes = 1. Default = 0>\n\t";
	print "-seedclass <0 = seeds only from miRs; 1 = seeds from all classes, default = 0>\n\t";
	print "-normalisation < 0 = raw counts; 1 = RPM miRs; 2 = RPM all RNA classes in document, default = 0>\n\n\n";

	exit(1);
}

if(scalar(@ARGV) == 0)
{    usage();	}

# Parse the Command Line
&parse_command_line($parameters, @ARGV);

sub parse_command_line
{

    my($parameters, @ARGV) = @_;
    my $CommentFlag = 0;
    my $next_arg;

    while(scalar @ARGV > 0)
    {
        $next_arg = shift(@ARGV);

	if ($next_arg =~ m/^-.*?/)
	{	$CommentFlag = 0;	}

	if ($next_arg eq "-s"){ $parameters->{species} = shift(@ARGV);  }
	elsif($next_arg eq "-in"){ $parameters->{input} = shift(@ARGV); }
	elsif($next_arg eq "-ml"){ $parameters->{min_length} = shift(@ARGV); }
	elsif($next_arg eq "-mm"){ $parameters->{mismatch} = shift(@ARGV); }
	elsif($next_arg eq "-gff"){ $parameters->{gff} = shift(@ARGV); }
	elsif($next_arg eq "-mat"){ $parameters->{mature_seq_file} = shift(@ARGV); }
	elsif($next_arg eq "-pre"){ $parameters->{precursor_seq_file} = shift(@ARGV); }
	elsif($next_arg eq "-out"){ $parameters->{output} = shift(@ARGV); }
	elsif($next_arg eq "-flank"){ $parameters->{flank} = shift(@ARGV); }
	elsif($next_arg eq "-ref"){ $parameters->{reference_format} = shift(@ARGV); }
	elsif($next_arg eq "-mintag"){ $parameters->{mintag} = shift(@ARGV); }
	elsif($next_arg eq "-tagfreq"){ $parameters->{tagfreq} = shift(@ARGV); }
	elsif($next_arg eq "-starttagfreq"){ $parameters->{starttagfreq} = shift(@ARGV); }
	elsif($next_arg eq "-mbase") { $parameters->{mbase} = shift(@ARGV); }
	elsif($next_arg eq "-mismatch") { $parameters->{mismatch} = shift(@ARGV); }
	elsif($next_arg eq "-showflank") { $parameters->{showflank} = shift(@ARGV); }
	elsif($next_arg eq "-seedclass") { $parameters->{seedclass} = shift(@ARGV); }
	elsif($next_arg eq "-normalisation") { $parameters->{normalisation} = shift(@ARGV); }
	elsif($next_arg eq "-comment")
	{ 	$parameters->{comment} = shift(@ARGV);
		$CommentFlag = 1;
	}
	elsif ($CommentFlag == 1) { $parameters->{comment} .= " ".$next_arg;	}
        else { print "Invalid argument: $next_arg"; usage(); }
    }

    # Check all essential options have been entered.
    my $Error_Log = '';
    $Error_Log .= "Species undefined (-s)\n\t", if (! defined $parameters->{species});
    $Error_Log .= "No output file defined (-in)\n\t", if (! defined $parameters->{input});
    if ($Error_Log ne '')
    {       print "\nThe following essential option(s) have not been defined:\n\t$Error_Log\n\n";
                usage();
    }

    # Set the sequence files if user did not specify
    $parameters->{reference_format} = 0, if (! defined $parameters->{reference_format});

    &load_default_reference_files($parameters->{species});

    print "\n\n**************    All parameters for Intermediate_to_miRspring script ********************************\n";
    while (my ($key, $value) = each %$parameters)
    {       print "$key\t$value\n"; }
	print "\n********************************************************************************************";
}


my $IsomiR_Range = 3;		# The accepted window upstream or downstream relative to miRBase to accept a miR as genuin


############################################### Load precursor file into memory
open MIRBASE, "$parameters->{precursor_seq_file}" or die "Cannot open $parameters->{precursor_seq_file}";
my $miRs;
my $Multi_loci;	# Hash that will save all multi-loci miRs
my $Order =0;	# The index, miRs are listed in chromosome order
my $TotalMultiLoci = 0;
my $No_Key;	# quick way for me to work out how many mulit-loci miRs.
while(<MIRBASE>)
{
	if ($_ =~ m/>(.*?),/)
    {       my $ID = $1;
		$ID =~ s/miR/mir/;
		$miRs->{$ID} = <MIRBASE>;
        chomp($miRs->{$ID});
        if ($ID =~ m/(\w\w\w-\w\w\w-.*?)-.*?$/)
        {	my $ShortID = $1;
			my $key = "$ShortID TOTAL";
			$Multi_loci->{$key} = 0, if ($Multi_loci->{$key} <= 0);
            $No_Key->{$key} = 0;
            $Multi_loci->{"$ShortID $Multi_loci->{$key}"}= $ID;
            $Multi_loci->{$key}++;
            $TotalMultiLoci++;
        }
        $miRs->{"$ID Middle"} = length($miRs->{$ID})/2;	# Estimate the middle of precursor and save
    }
}
close MIRBASE;
my @Allkeys = keys(%$No_Key);

#print "\nA total of $TotalMultiLoci miRs that are clustered from a total of ".scalar(@Allkeys)." miRs\n";

open GENOME, $parameters->{gff} or die "Cannot open mirbase genome feature file ($parameters->{gff})";
my $Genome_Coords;
while (<GENOME>)
{	# 1	.	miRNA	30366	30503	.	+	.	ACC="MI0006363"; ID="hsa-mir-1302-2";
	my @mir_coord = split('\t',$_);
        if ($mir_coord[8] =~ m/\"(\w\w\w-\w\w\w-.*?)\"/)
        {	my $ID = $1;
        	$Genome_Coords->{$1} = "'$mir_coord[0]', $mir_coord[3], '$mir_coord[6]'";
        	$miRs->{$Order++} = $ID;
        }
}
close GENOME;


#### Load mature file into memory. Need to know coordinates derived from precursor.
open MATURE, "$parameters->{mature_seq_file}" or die "Cannot open mirbase mature file ($parameters->{mature_file})";
my $Mature_Seq;
while(<MATURE>)
{
	if ($_ =~ m/>(.*?) MI/)
    {	my $Mature_ID = $1;
		$Mature_ID =~ tr/R/r/;
        my $Mature_Seq = <MATURE>;
        chomp($Mature_Seq);
        $Mature_Seq =~ tr/Uu/Tt/;
		$Mature_Seq =~ s/\r//;

#     	if ($Mature_ID =~ m/($parameters->{species}-\w\w\w-\w+)(.*?)/)
        if ($Mature_ID =~ m/(\w\w\w-\w\w\w-\w+)(.*?)/)
		{	my $miR_Number = $2; # This will contain a '-' character if multi-loci
			my $Temp = $1;
            $Mature_ID = $1;
			my $Mature_ID_Master = $1;
            my $i = 0;
            my $Max_Loop = 0;
            $Max_Loop = $Multi_loci->{"$Mature_ID TOTAL"}, if (defined $Multi_loci->{"$Mature_ID TOTAL"});
            do
            {
				if ($miRs->{$Mature_ID} =~ m/($Mature_Seq)/i)
                {	my $Match_Pos = length($`);
					my $Side = '5p';
                    $Side = '3p', if ($Match_Pos > $miRs->{"$Mature_ID Middle"});
					# Check to see if already defined
					if (defined $miRs->{"$Mature_ID $Side"})
					{	if ($Match_Pos > $miRs->{"$Mature_ID $Side"})
						{	$miRs->{"$Mature_ID 5p"} = $miRs->{"$Mature_ID $Side"};
							$miRs->{"$Mature_ID 5p END"} = $miRs->{"$Mature_ID $Side END"};
							$Side = '3p';
						}
						elsif ($Match_Pos < $miRs->{"$Mature_ID $Side"})
						{   $miRs->{"$Mature_ID 3p"} = $miRs->{"$Mature_ID $Side"};
							$miRs->{"$Mature_ID 3p END"} = $miRs->{"$Mature_ID $Side END"};
							$Side = '5p';
						}
						#print "\nSwapped $Mature_ID coordinates";
					}
					$miRs->{"$Mature_ID $Side"} = $Match_Pos;	# Record miRbase position
                    $miRs->{"$Mature_ID $Side END"} = $Match_Pos + length($Mature_Seq);
                }
                $Mature_ID = $Multi_loci->{"$Mature_ID_Master $i"};
                $i++;
            } while ($i <= $Max_Loop);
        }
    }
}
close MATURE;


#### Load raw data into memory
open MIRDSPRING_DATA, $parameters->{input} or  die "Cannot open the raw mirdspring data file ($parameters->{input})";
my $Raw_data;
$HeaderInfo .= ScriptSettings();
while(<MIRDSPRING_DATA>)
{
	if ($_ =~ m/^\@RG/) # @RG     ID:Library6_6   PL:SOLID        PU:2_4  LB:Library6     PI:0    DS:50F  DT:2012-06-01T15:42:33+1000     SM:Traude_Nico_Bene_3
	{ 	my $count =0;
		chomp();
		s{(\t)}{((++$count == 6) || ($count == 12)) ? "<br>" : $1 }ige;
		$HeaderInfo .= "\nHeaderInfo.RG = \"$_\";";
	}
	elsif ($_ =~ m/^\@PG/) # @PG     ID:Victor_2012-06-01T15:42:33+1000_2_4  PN:LifeScope    VN:2.5
	{	my $count =0;
		chomp();
		s{(\t)}{((++$count == 6) || ($count == 12)) ? "<br>" : $1 }ige;
		$HeaderInfo .= "\nHeaderInfo.PG = \"$_\";";
	}
	elsif ($_ =~ m/^\@CO/) #@CO     RG:Library6_6   IX:6    II:6    LD:Bene B       LT:Fragment     AT:Small RNA    DE:2012-06-04T16:10:19+1000     CU:18743080     CT:18743080     SP:RNA-Seq SD:      BX:0    TN:25   TX:50   EC:0    ER:0    CO:Administrator        PJ:Hentze Lab   SO:Traude Bene Nico     PN:Victor
	{	my $count =0;
		chomp();
		s{(\t)}{((++$count == 6) || ($count == 12)) ? "<br>" : $1 }ige;
		$HeaderInfo .= "\nHeaderInfo.CO = \"$_\";";
	}
	elsif ($_ =~ m/^\@MIRSPRING/)
	{		}  	# Have alread saved this information. Do nothing as don't wont to double up.
	elsif ($_ =~ m/^\@/)
	{	my $count =0;
		chomp();
		s{(\t)}{((++$count == 6) || ($count == 12)) ? "<br>" : $1 }ige;
		$HeaderInfo .= "\nHeaderInfo.Other = \"$_\";";
	}


	# fileindex	Start pos	miRNA name	Freq	Length	ERROR
	if ($_ =~ m/^(.*?)\t(\d+)\t(\w\w\w-\w\w\w-.*?)\t(\d+)\t(\d+)\t(.*?)$/)
	{	
		my ($Index,$Start, $miR_Name, $Freq, $Length, $Error) = ($1,$2,$3,$4,$5,$6);
		next, if ($Freq < $parameters->{tagfreq});
		if ($Error =~ m/^E(\d+)E\d+(\w)to(\w)/) #old format
		{	$Error = "$1$2$3";	}
		elsif ($Error =~ m/^(\d+\w\w)/) # mismatch or possible editing
		{	$Error = $1;	}
		elsif ($Error =~ m/^(_\w\w)/) # Non template addition
		{	$Error = $1;	}
        else
        {	$Error = 'NONE';	}
		$Raw_data->{$miR_Name}->{$Start}->{$Length}->{$Error} += $Freq;
    }
}
close MIRDSPRING_DATA;

#### Now to export data into javascript format
## STEP 1 Sort the raw data and save total counts
my $New_Index = 0;
my $DataSet_miRs;	# Hash to save relevant miRs

# Javascript array object:
#/* GLOBAL VARIABLES */
#/* [Name][Sequence][Secondary structure][Total counts][5p counts][3p counts][Other Counts]] */
my $JS_premiR_Var;	# Holds all info about precursor
#/* [miR_Index][Start][Length][count][Editing/End modification]*/
my $JS_Data_Var;	# Holds all deep sequencing tag coordinates.
my $JS_Genome_Coords;   # Holds genome coordinates of all miRs

for (my $i = 0; $i < $Order;$i++)
{	my $miR_ID = $miRs->{$i};
	my $Current_miR = $Raw_data->{$miR_ID};
    my $Hairpin_Total_Counts = 0;
    my $Output;	#Debug purposes
    my ($FiveP, $ThreeP, $Non_Canonical, $miR_Avg_Length, $Editing) = (0,0,0,0,0);
    my ($FiveP_5pIsomiR, $FiveP_3pIsomiR, $ThreeP_5pIsomiR,$ThreeP_3pIsomiR) = (0,0,0,0);

    # Now sort keys
    foreach (sort {$Current_miR->{$b} <=> $Current_miR->{$a}} keys %$Current_miR)
    {   my $Start_Pos = $_;
		my $Current_Start = $Current_miR->{$Start_Pos};
		
		# Check to see if start_tagfreq is set. If so then do precount to see if we need to process.
		if ($parameters->{starttagfreq} > 0)
		{	my $CummulativeStartCount = 0;
			foreach (sort {$Current_Start->{$b} <=> $Current_Start->{$a}} keys %$Current_Start)
			{	my $Length = $_;
				my $Mismatch_Data = $Current_Start->{$Length};
				foreach (keys %$Mismatch_Data)
				{	my $Error_Code = $_;
					my $Freq = $Mismatch_Data->{$Error_Code};
					next, if (($parameters->{mismatch} == 0) && ($Error_Code ne ''));
					$CummulativeStartCount += $Freq;
				}
			}
			#print "\nDebug $CummulativeStartCount", if ($miR_ID eq "hsa-mir-34c");
			next, if ($parameters->{starttagfreq} > $CummulativeStartCount);
			#print " .....pass", if ($miR_ID eq "hsa-mir-34c");
		}
		

        foreach (sort {$Current_Start->{$b} <=> $Current_Start->{$a}} keys %$Current_Start)
        {	my $Length = $_;
			# Need to check length and start position are not outside boundaries of precursor + flank

			my $Mismatch_Data = $Current_Start->{$Length};
            foreach (keys %$Mismatch_Data)
            {	my $Error_Code = $_;
				my $Freq = $Mismatch_Data->{$Error_Code};
				$Error_Code = '', if ($Error_Code eq 'NONE');
				$JS_Data_Var .=  "[$New_Index,$Start_Pos,$Length,$Freq,\'$Error_Code\'],";

				next, if (($parameters->{mismatch} == 0) && ($Error_Code ne ''));

				$miR_Avg_Length += ($Length * $Freq); # NEED TO PUT CONDITION TO SAVE LENGTHS ONLY WITHIN PRECURSOR  and not flanking sequences!!
				$Editing += $Freq, if ($Error_Code ne '');
				$Hairpin_Total_Counts += $Freq;

				if (abs($Start_Pos -  $miRs->{"$miR_ID 5p"}) <= $IsomiR_Range)
				{	$FiveP += $Freq;
					# Calculate all miRs that start EXACTLY as described by miRBase. Later I can work out the % of isomiRS
					if ( (abs($Start_Pos + $Length -  $miRs->{"$miR_ID 5p END"})== 0) || (abs($Start_Pos + $Length -  $miRs->{"$miR_ID 3p END"}) == 0) )
	                {	$FiveP_3pIsomiR += $Freq;	}

	                if ( (abs($Start_Pos -  $miRs->{"$miR_ID 5p"})== 0) ||(abs($Start_Pos - $miRs->{"$miR_ID 3p"}) == 0) )
        	        {	$FiveP_5pIsomiR += $Freq;	}
                }
                elsif (abs($Start_Pos - $miRs->{"$miR_ID 3p"}) <= $IsomiR_Range)
                {	$ThreeP += $Freq;
	                if ( (abs($Start_Pos + $Length - $miRs->{"$miR_ID 5p END"})== 0) ||(abs($Start_Pos + $Length - $miRs->{"$miR_ID 3p END"}) == 0) )
	                {	$ThreeP_3pIsomiR += $Freq;	}

                    if ( (abs($Start_Pos -  $miRs->{"$miR_ID 5p"})== 0) ||(abs($Start_Pos - $miRs->{"$miR_ID 3p"}) == 0) )
        	        {	$ThreeP_5pIsomiR += $Freq;	}
                }
                else
                {	$Non_Canonical += $Freq	}

			} # foreach (keys %$Error_Code)
		} # foreach (sort {$Current_Start->{$b} <=> $Current_Start->{$a}} keys %$Current_Start)
	} # foreach (sort {$Current_miR->{$b} <=> $Current_miR->{$a}} keys %$Current_miR)

    if ($Hairpin_Total_Counts > $parameters->{mintag})
    {
		$DataSet_miRs->{$New_Index} = $miRs->{$i};	# Save miR index
		$DataSet_miRs->{"$New_Index TOTAL"} = $Hairpin_Total_Counts;
		$DataSet_miRs->{"$New_Index AVG LENGTH"} = (int($miR_Avg_Length/$Hairpin_Total_Counts*10))/10;
#		my $AverageLength_Fix = (($miR_Avg_Length/$Hairpin_Total_Counts) *10);
		
#		print "\nI see $AverageLength_Fix vs ".$DataSet_miRs->{"$New_Index AVG LENGTH"};
		
		# Save 5p, 3p and non-canonical counts
        $DataSet_miRs->{"$New_Index 5p"} = $FiveP;
        $DataSet_miRs->{"$New_Index 3p"} = $ThreeP;
        $DataSet_miRs->{"$New_Index Non_Canonical"} = $Non_Canonical;
		$DataSet_miRs->{"$New_Index Mismatch"} = $Editing;

        $FiveP_5pIsomiR = $FiveP - $FiveP_5pIsomiR;
		$FiveP_3pIsomiR = $FiveP - $FiveP_3pIsomiR;
        $ThreeP_5pIsomiR = $ThreeP - $ThreeP_5pIsomiR;
        $ThreeP_3pIsomiR =  $ThreeP - $ThreeP_3pIsomiR;

        #  /* GLOBAL VARIABLES */
		#/* [Name][Sequence][Secondary structure]	[Total counts][5p counts][3p counts][Other Counts]] */
        my $Total_Count = $FiveP + $ThreeP + $Non_Canonical;
        my @Start_n_Stop = ($miRs->{"$miRs->{$i} 5p"},$miRs->{"$miRs->{$i} 5p END"},$miRs->{"$miRs->{$i} 3p"},$miRs->{"$miRs->{$i} 3p END"});

		$JS_premiR_Var .= "[$New_Index,'$miRs->{$i}','$miRs->{$miRs->{$i}}','',$Start_n_Stop[0],$Start_n_Stop[1],$Start_n_Stop[2],$Start_n_Stop[3],$Total_Count,$FiveP, $ThreeP, $Non_Canonical,".$DataSet_miRs->{"$New_Index AVG LENGTH"};
        $JS_premiR_Var .= ",$FiveP_5pIsomiR,$FiveP_3pIsomiR,$ThreeP_5pIsomiR,$ThreeP_3pIsomiR,$Editing],\n";
		$JS_premiR_Var =~ s/\r//g;
        $JS_Genome_Coords .= "[$New_Index, ".$Genome_Coords->{$miRs->{$i}}."],";

        $New_Index++;
        #  print $Output;

    }
    else
    {	$JS_Genome_Coords .= "['$miRs->{$i}', ".$Genome_Coords->{$miRs->{$i}}."],";	}

}


open MIRDSPRING_DATA, $JavaScriptTemplateFile or  die "Cannot open javascript template file ($JavaScriptTemplateFile)";
open OUTPUTFILE, ">".$parameters->{output} or die "Cannot create javascript output file: $parameters->{output}";
my $miRspringVersion = "0.9";
while (<MIRDSPRING_DATA>)
{	if ($_ =~ m/^MIRDSPRING DATA HERE/)
	{	print OUTPUTFILE "$HeaderInfo\n";

        print OUTPUTFILE "var Options = ['Length (nt)', 'Editing (no data)', '5\\' IsomiR (%)', '3\\' IsomiRs (%)', 'Arm processing 5p/(5p+3p)','Non-Canonical Processing (%)', 'Polycistronic miR analysis', 'Cummulative Reads', 'Threshold (miRs > minimum abundance)', 'Seed Frequency (data set)','% of tags edited'];var GraphOption = 7;\n";
		print OUTPUTFILE "var WindowSize = [1000,5000,12000,25000,50000];\nvar WindowOption = 3;	// The index for WindowSize.\n";
		print OUTPUTFILE "var miRBaseWindow = $parameters->{mbase};  // the mirbase window used to count tags\n";
		print OUTPUTFILE "var RNAEditing = ";
                my $mm = 1;
                $mm = 0, if ($parameters->{mismatch} == 1);

                print OUTPUTFILE "$mm;		// 0 = show all, 1 = only show perfect, 2 = only show editing\n";
		print OUTPUTFILE "var FlankingSeq = $parameters->{showflank}; 	// 0 = show only when selected; 1 = always show\n";
		print OUTPUTFILE "var SeedClass =$parameters->{seedclass}; 		// 0 = seeds only from miRs; 1 = seeds from all classes\n";
		print OUTPUTFILE "var RPM = $parameters->{normalisation};			// 0 = raw counts; 1 = RPM miRs; 2 = RPM all RNA classes in document\n";
		print OUTPUTFILE "var Colour = ['Blue','DarkOrchid','Green','SteelBlue','Coral','Violet','Indigo','Orange','Sienna','Purple','DarkOrchid','Green','SteelBlue','Coral','Violet','Indigo','Orange','Sienna','Purple'];\n";
		print OUTPUTFILE "var SpacerMax =  $parameters->{flank};\n";

		print OUTPUTFILE "var miRBaseDataSet = [\n$JS_premiR_Var];\n\n";
		print OUTPUTFILE "var data = [\n$JS_Data_Var];\n\n";
		print OUTPUTFILE "var Genome = [\n$JS_Genome_Coords];\n\n";


        }
        elsif ($_ =~ m/^DESCRIPTION HERE/)
        {	print OUTPUTFILE "$parameters->{comment}";	}
	 elsif ($_ =~ m/^Version (\d+\.\d+).*?$/)
	 {	$miRspringVersion = $1;				
		$HeaderInfo .= "\nHeaderInfo.VER = \"$miRspringVersion\";";
		print OUTPUTFILE $_;
	 }
        else
        {	print OUTPUTFILE $_;	}

}
close OUTPUTFILE;
close MIRDSPRING_DATA;

print "\n$parameters->{output} file created\n";
exit(0);

sub ScriptSettings()
{	my $miRspringHeader =  "\nHeaderInfo.MIRSPRING += \"<br><br>\@MIRSPRING Intermediate_to_miRspring.pl settings<br>";
	
	$miRspringHeader   .= "Species (-s): ".$parameters->{species}."<br>";
	$miRspringHeader   .= "Min length (-ml): ".$parameters->{min_length}."<br>";
	$miRspringHeader   .= "Number of mismatches (-mm): ".$parameters->{mismatch}."<br>";
	$miRspringHeader   .= "Number of flanking nucleotides (-flank): ".$parameters->{flank}."<br>";
	$miRspringHeader   .= "Reference format (-ref): ".$parameters->{reference_format}."<br>";
	$miRspringHeader   .= "Minimum number of tags (-mintag): ".$parameters->{mintag}."<br>";
	$miRspringHeader   .= "tag frequency (-tagfreq): ".$parameters->{tagfreq}."<br>";
	$miRspringHeader   .= "Start tag frequency (-starttagfreq): ".$parameters->{starttagfreq}."<br>\";\n";
	
	return $miRspringHeader;
}


