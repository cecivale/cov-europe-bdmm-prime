#!/bin/sh
WORKDIR="/Users/maceci/code/mt-analysis/"
ANALYSIS="200910_dsEurope0"
N=10

# Results Processing 
# Analysis is run on Euler and entire analysis folder move back to local folder in /mt-analysis
mkdir $WORKDIR/$ANALYSIS/pResults

# 1. Combine log files with Log Combiner 
# Check first convergence in Tracer
LOGCOMBINER="/Applications/BEAST 2.6.3/bin/logcombiner"
$LOGCOMBINER \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.1.log \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.2.log \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.3.log \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.4.log \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.5.log \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.6.log \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.7.log \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.8.log \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.9.log \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.10.log \
-o $WORKDIR/pResults/${ANALYSIS}.log \
-b 10 

# 3. Log summary table
# Visualize combined log file and export data table from Tracer
# Transform table to .md
Rscript $WORKDIR/txttomd.R \
--input $WORKDIR/$ANALYSIS/pResults/${ANALYSIS}_logSummary 


# 4. Combine typed node trees with Log Combiner
$LOGCOMBINER \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.1.typed.node.trees \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.2.typed.node.trees \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.3.typed.node.trees \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.4.typed.node.trees \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.5.typed.node.trees \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.6.typed.node.trees \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.7.typed.node.trees \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.8.typed.node.trees \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.9.typed.node.trees \
-log $WORKDIR/$ANALYSIS/rResults/${ANALYSIS}.10.typed.node.trees \
-o $WORKDIR/pResults/${ANALYSIS}.typed.node.trees \
-b 10 

# 5. Summary tree
TREEANNOTATOR="/Applications/BEAST 2.6.3/bin/treeannotator"
$TREEANNOTATOR $WORKDIR/pResults/${ANALYSIS}.typed.node.trees  \
$WORKDIR/pResults/${ANALYSIS}.typed.node.tree





