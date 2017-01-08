(* ::Package:: *)
(* Wolfram Language Package *)
(* Created with the Mathematica plugin for IntelliJ IDEA,
 http://wlplugin.halirutan.de/ *)

(* :Title: Insertion Mapping *)
(* :Context: InsertionMapping` *)
(* :Description: Functions for analysis of insertion screen data *)
(* :Author: Frank Buermann *)

(* :Package Version: 0.2.1 *)
(* :Mathematica Version: 11.0.1 *)
(* :Copyright: (c) 2017 Frank Buermann *)

(***********************************************************************)
BeginPackage["InsertionMapping`"];
ClearAll["`*"];

(***********************************************************************)
(* :Documentation: *)

reverseDNA::usage =
    "reverseDNA[\"seq\"] returns the reverse complementary sequence of\
 \"seq\".";

translateCDS::usage =
    "translateCDS[\"seq\"] translates the DNA sequence \"seq\" into an\
 amino acid sequence.";

importFastq::usage =
    "importFastq[file] parses a FastQ file or stream into an FastQ\
 association.
importFastq[file, n] parses n entries.";

insertionMapper::usage =
    "insertionMapper[\"target\", \"ins\", reads]\
 takes reads <|\"ID\" -> \"Sequence\", ... |> and maps those reads\
 containing \"ins\" to \"target\".\
 Output is <|\"ID\" -> <|\"Mappings\" -> coordinates,\
 \"Direction\" -> dir|>|>. 
 The coordinates list contains four elements that have been generated\
 with the following comparisons:\
 {ins vs. target, reverseDNA@ins vs. target, ins vs. reverseDNA@target,\
 reverseDNA@ins vs. reverseDNA@target}. Coordinate pairs report the\
 first and last matching bases, or {0, 0} if unmapped.
 dir is \"Fwd\" or \"Rev\" (for unambiguously mapped reads) or\
 \"NA\" (for unmappable or ambiguous reads).
 All coordinates correspond to native bases flanking the insert.";

findUpstreamCodon::usage =
    "findUpstreamCodon[\"cds\", \"ins\", pos, \"dir\"] gets the\
 closest codon position upstream of the insert that does not lead to\
 an amino acid change. \"cds\" is the coding sequence of the gene,\
 \"ins\" is the insert sequence, pos is the insert position, \"dir\"\
 is the insert direction.";

findCodon::usage =
    "findCodon[base] returns the codon number belonging to a base\
 position.";

findFrame::usage =
    "findFrame[pos, \"dir\"] returns the reading frame of the insert\
 with position pos and orientation \"dir\".";

coiledCoilFraction::usage =
    "coiledCoilFraction[list] computes the fraction of residues in\
 list which are part of the B. subtilis Smc coiled coil.";

coiledCoilEnrichment::usage =
    "coiledCoilEnrichment[s, ref] computes the coiled-coil enrichment\
 of a sample s relative to a reference ref.";

confidenceInterval::usage =
    "confidenceInterval[sample] estimates the 95 % confidence interval\
 of a sample.";

permutationTest::usage =
    "permutationTest[f, s1, s2] estimates the p-value for\
 f[p1, p2] <= f[s1, s2] using samples s1 and s2 and 10^5 permutations\
 p1 and p2.
permutationTest[f, s1, s2, n] uses n permutations.";

permutationTestConfideneceInterval::usage =
    "permutationTestConfidenceInterval[f, s1, s2, n, m] estimates\
 the 95 % confidence interval for the one-sided p value for\
 (f[p1, p2] <= f[s1, s2]) from m tests with n permutations {p1, p2}.";

(***********************************************************************)

(***********************************************************************)
Begin["`Private`"];
ClearAll["`*"];

ccRegions = {{167, 499}, {673, 1027}};

translationRules =
    {
      "UCA" -> "S",
      "UCC" -> "S",
      "UCG" -> "S",
      "UCU" -> "S",
      "UUC" -> "F",
      "UUU" -> "F",
      "UUA" -> "L",
      "UUG" -> "L",
      "UAC" -> "Y",
      "UAU" -> "Y",
      "UAA" -> "Z",
      "UAG" -> "Z",
      "UGC" -> "C",
      "UGU" -> "C",
      "UGA" -> "Z",
      "UGG" -> "W",
      "CUA" -> "L",
      "CUC" -> "L",
      "CUG" -> "L",
      "CUU" -> "L",
      "CCA" -> "P",
      "CCC" -> "P",
      "CCG" -> "P",
      "CCU" -> "P",
      "CAC" -> "H",
      "CAU" -> "H",
      "CAA" -> "Q",
      "CAG" -> "Q",
      "CGA" -> "R",
      "CGC" -> "R",
      "CGG" -> "R",
      "CGU" -> "R",
      "AUA" -> "I",
      "AUU" -> "I",
      "AUC" -> "I",
      "AUG" -> "M",
      "ACA" -> "T",
      "ACC" -> "T",
      "ACG" -> "T",
      "ACU" -> "T",
      "AAC" -> "N",
      "AAU" -> "N",
      "AAA" -> "K",
      "AAG" -> "K",
      "AGC" -> "S",
      "AGU" -> "S",
      "AGA" -> "R",
      "AGG" -> "R",
      "GUA" -> "V",
      "GUC" -> "V",
      "GUG" -> "V",
      "GUU" -> "V",
      "GCA" -> "A",
      "GCC" -> "A",
      "GCG" -> "A",
      "GCU" -> "A",
      "GAC" -> "D",
      "GAU" -> "D",
      "GAA" -> "E",
      "GAG" -> "E",
      "GGA" -> "G",
      "GGC" -> "G",
      "GGG" -> "G",
      "GGU" -> "G"
    };

(***********************************************************************)
(* :Helper functions: *)

reverseDNA[seq_String] :=
    StringReplace[
      StringReverse[seq],
      {
        "A" -> "T",
        "C" -> "G",
        "G" -> "C",
        "T" -> "A",
        "a" -> "t",
        "c" -> "g",
        "g" -> "c",
        "t" -> "a"
      }
    ];

translateCDS[seq_String] :=
    StringJoin[
      StringJoin /@
          Partition[Characters[ToUpperCase@seq] /. "T" -> "U", 3]
          /. translationRules
    ];

importFastq[file_] := importFastq[file, Hold[ReadList[file, Record]]];

importFastq[file_, n_Integer] :=
    importFastq[file, Hold[ReadList[file, Record, n * 4]]];

importFastq[file_, readListExpression_Hold] :=
    (
      ReleaseHold@readListExpression
          // Partition[#, 4]&
          //
          Map[
            First@# -> AssociationThread[
              {"Sequence", "Separator", "QualityString"} -> Rest@#
            ]&
          ]
          // Association
    );

extractFlankingSequences[target_String, seq_String] :=
    (
      StringPosition[target, seq, IgnoreCase -> True]
          /.
          {
            {{__}, __} -> {},
            {a_, b_} :>
                Sequence[
                  StringTake[target, a - 1], StringDrop[target, b]
                ]
          }
    );

trimSequences[{a_String, b_String}, n_Integer] :=
    {
      If[StringLength@a < n, a, StringTake[a, -n]],
      If[StringLength@b < n, b, StringTake[b, n]]
    };

trimSequences[a_, n_Integer] := a;

Options[mapSequencePosition] = {"Homology" -> 10};
mapSequencePosition[target_, seq_, opts : OptionsPattern[]] :=
    If[
      StringLength@seq >= OptionValue["Homology"]
      ,
      StringPosition[target, seq, IgnoreCase -> True]
          /.
          {
            {_, __List} -> {0, 0},
            {} -> {0, 0},
            {{a_, b_}} :> {a, b}
          }
      ,
      {0, 0}
    ];

getInsertionCoordinates[positions_List] :=
    {Last@First@#, First@Last@#}& /@ positions;

Options[mapInsertionPoints] = {"Homology" -> 10};
mapInsertionPoints[
  targetSeq_String,
  flankingSequences_List,
  opts : OptionsPattern[]
] :=
    (
      flankingSequences
          /.
          {
            s : {_, _} :>
                (
                  mapSequencePosition[targetSeq, #,
                    FilterRules[
                      {opts, Options[mapInsertionPoints]},
                      Options[mapSequencePosition]
                    ]
                  ]&
                      /@ s
                ),
            {} -> {{0, 0}, {0, 0}}
          }
          // getInsertionCoordinates
    );

ttInterpreter[tab_] :=
    Switch[
      tab,
      {False, False, False, False},
      "NA"
      ,
      {True, False, False, False},
      "Fwd"
      ,
      {False, True, False, False},
      "Rev"
      ,
      {False, False, True, False},
      "Rev"
      ,
      {False, False, False, True},
      "Fwd"
      ,
      _,
      "NA"
    ];

convertToFwdPosition[pos_, seqLen_] :=
    (
      pos /.
          {
            s : {a_?NumericQ, b_?NumericQ} /; FreeQ[s, 0]
                :> {seqLen - b + 1, seqLen - a + 1}
            ,
            s : {a_?NumericQ, b_?NumericQ} /; a =!= 0
                :> {0, seqLen - a + 1}
            ,
            s : {a_?NumericQ, b_?NumericQ} /; b =!= 0
                :> {seqLen - b + 1, 0}
          }
    );

(***********************************************************************)
(* :Main function: *)

Options[insertionMapper] = {"Homology" -> 10, "Trimming" -> False};

insertionMapper[
  targetSeq_String,
  ins_String,
  reads_Association,
  opts : OptionsPattern[]
] :=
    AssociationThread[
      Keys@reads -> insertionMapper[targetSeq, ins, Values@reads, opts]
    ];

insertionMapper[
  targetSeq_String,
  ins_String,
  reads_List,
  opts : OptionsPattern[]
] :=
    Module[
      {
        len = StringLength@targetSeq,
        trimmingFunction = Identity,
        fwdFlanks, revFlanks, mappings, truthTable, directions
      },

      If[
        OptionValue["Trimming"] =!= False,
        trimmingFunction = trimSequences[#, OptionValue["Trimming"]]&
      ];

      fwdFlanks =
          trimmingFunction[extractFlankingSequences[#, ins]]& /@ reads;
      revFlanks =
          trimmingFunction[
            extractFlankingSequences[#, reverseDNA@ins]]& /@ reads;

      mappings =
          Outer[
            mapInsertionPoints[##,
              FilterRules[
                {opts, Options[insertionMapper]},
                Options[mapInsertionPoints]]
            ]&
            ,
            {targetSeq, reverseDNA@targetSeq},
            {fwdFlanks, revFlanks}
            , 1
          ]
              ~ Flatten ~ 1
              // MapAt[convertToFwdPosition[#, len]&, #, {{3}, {4}}]&
              // Transpose;

      truthTable = mappings /. {{0, 0} -> False, {_, _} -> True};

      directions = ttInterpreter /@ truthTable;

      {
        Thread["Mappings" -> mappings],
        Thread["Direction" -> directions]
      }
          // Transpose
          // Map[Association]

    ];

(***********************************************************************)
(* :Analysis: *)

findCodon[base_Integer] :=
    Ceiling[base / 3];

affectedCodon[cds_String, insPos_Integer] :=
    translateCDS@StringTake[cds, 3 * findCodon[insPos] + {-2, 0}];

findUpstreamCodon[
  cds_String,
  ins_String,
  pos_Integer,
  direction_String
] :=
    Catch[
      Module[
        {insert, refAA, mutAA},

        If[Divisible[pos, 3], Throw[pos / 3, "findInsertionPoint"]];

        insert =
            Which[
              direction === "Fwd",
              ins
              ,
              direction === "Rev",
              reverseDNA@ins
              ,
              True,
              Throw[$Failed, "findInsertionPoint"]
            ];

        refAA = affectedCodon[cds, pos];
        mutAA = affectedCodon[StringInsert[cds, insert, pos + 1], pos];

        If[refAA == mutAA, findCodon[pos], findCodon[pos] - 1]
      ]
      , "findInsertionPoint"
    ];

findFrame[pos_Integer, direction_String] :=
    Which[
      direction === "Fwd",
      Mod[pos, 3] + 1
      ,
      direction === "Rev",
      -(Mod[pos, 3] + 1)
      ,
      True,
      $Failed
    ];

coiledCoilFraction[positions_List] :=
    (
      positions
          // Select[IntervalMemberQ[Interval @@ ccRegions, #]&]
          // Length[#] / Length[positions]&
    );

coiledCoilEnrichment[s_List, ref_List] :=
    coiledCoilFraction[s] / coiledCoilFraction[ref];

permutationStatistics[
  f_,
  sample1_List,
  sample2_List,
  n_Integer
] :=
    With[
      {joined = sample1 ~ Join ~ sample2},

      Table[
        f @@ TakeDrop[
          RandomSample[joined, Length@joined],
          Length@sample1
        ],
        n
      ]
    ];

permutationTest[
  f_,
  sample1_List,
  sample2_List,
  n_Integer : 1000
] :=
    With[
      {
        statistic = f[sample1, sample2],
        permutations = permutationStatistics[f, sample1, sample2, n]
      },

      Select[permutations, # <= statistic&]
          // Length[#] / Length[permutations]&
    ];

confidenceInterval[sample_List, p_Real : 0.95] :=
    {Quantile[sample, (1 - p) / 2], Quantile[sample, 1 - (1 - p) / 2]};

permutationTestConfideneceInterval[
  f_,
  sample1_List,
  sample2_List,
  n_Integer : 10,
  m_Integer : 10
] :=
    (
      Table[permutationTest[f, sample1, sample2, n], m]
          // confidenceInterval
    );

(***********************************************************************)
End[];

EndPackage[];