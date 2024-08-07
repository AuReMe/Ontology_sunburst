
cl_tst1 = 'ontosunburst -i tests/global_test/test1_int.tsv -ir tests/global_test/test1_ref.tsv ' \
          '--onto metacyc --pcut deeper --kwargs bg_color=pink title="C\'est un titre" mx_dpt=3'

cl_tst2 = ' ontosunburst -i Files/test2_int.tsv -ir Files/test1_ref.tsv --onto metacyc ' \
          '-a enrichment -o test2 -sl --pcut deeper'
