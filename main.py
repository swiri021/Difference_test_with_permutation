from difference_test.diff_stats import diff_stats
import pandas as pd

df = pd.read_csv('Some_of_your_data.csv', index_col=0)
group1 = ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5', 'Sample7']
group2 = ['Somegroup']*100

diff_s = diff_stats(df, small_group=group1, large_group=group2)

# Combined result
result = diff_s.combined_pvalues(perm=100)

# T-test result
result = diff_s._permuted_ttest(df, small_group=group1, large_group=group2, perm=100)

# Median-test result
result = diff_s._permuted_median_test(df, small_group=group1, large_group=group2, perm=100)

print result
