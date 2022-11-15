cd manual_fixes_pombase/
python format_manual_changes.py
cd ../
python perform_qc.py --alleles manual_fixes_pombase/manual_changes_applied.tsv --output manual_fixes_pombase/allele_results_manual_changes.tsv
cd manual_fixes_pombase
cat allele_results_manual_changes_errors.tsv|grep -v nucle|grep -v disrupt|grep -v wtf|grep -v canto|grep -v remove|grep -v CTD|grep -v zas1|grep -v cannot|grep -v "lengths don't match"|grep -v asked  > alleles_to_fix.tsv