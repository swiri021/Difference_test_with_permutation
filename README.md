# Difference test by Permutation method
This code calculates p-value of difference by using permuted sampling T-test and Median test in case of unbalanced sample size between 2 groups. For example, group1 has only 7 samples, and group2 has only 200 samples.

# Example :
```Python
from difference_test.diff_stats import diff_stats
import pandas as pd

# Your Pandas file
df = pd.read_csv('Some_of_your_data.csv', index_col=0)
# Small group
group1 = ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5', 'Sample7']
# Large group
group2 = ['Somegroup']*100

diff_s = diff_stats(df, group1=group1, group2=group2)



# T-test result
result = diff_s._permuted_ttest(df, small_group=group1, large_group=group2, perm=100)

# Median-test result, this function is needed to be fixed later(Wrong algorithm)
# result = diff_s._permuted_median_test(df, small_group=group1, large_group=group2, perm=100)
# Combined result, this function is needed to be fixed later(Wrong algorithm)
# result = diff_s.combined_pvalues(perm=100)


print result
```
