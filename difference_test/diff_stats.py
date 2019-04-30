import pandas as pd
import numpy as np
from scipy import stats
import random
import itertools


class diff_stats:

	def _storey_fdr(self, df, p_name):
		neg_ps = float(df.loc[df[p_name] > 0.5].count()[p_name])
		total_ps = float(df[p_name].count())
		q_cal = lambda x : x/(neg_ps/total_ps)

		df[p_name+'_adj'] = df[p_name]
		df[p_name+'_adj'] = df[p_name+'_adj'].apply(q_cal)
		df.loc[df[p_name+'_adj'] >= 1.0,p_name+'_adj'] = 0.999999

		return df


	def _permutation(self, **kwargs): ##### Permutation samples
		return [random.sample(kwargs['large_group'], len(kwargs['small_group'])) for a in range(0,kwargs['perm'])]

	def _permuted_ttest(self, df, **kwargs): ##### T-test
		small_group = kwargs['small_group']
		large_group = kwargs['large_group']
		perm = kwargs['perm']

		assert perm>2, "Wrong Permutation number"

		test_genes = df.index.tolist()

		permuted_pool = self._permutation(small_group=small_group, large_group=large_group, perm=perm)
		p_arr = []

		for x in test_genes: ##### Each genes
			t_stats = []
			g2_arr = []
			g1 = df[small_group].loc[x]
			for p in permuted_pool: ##### Permutation table
				g2 = df[p].loc[x]
				s, p = stats.ttest_ind(g1, g2)
				t_stats.append(s) # T-statistics
				g2_arr.append(g2.values.tolist()) # Expression of Permuted samples

			permuted_pval = stats.t.sf(abs(np.mean(t_stats)), len(t_stats)-1)*2 # Going to p values from Mean of T- statistics

			g2_arr = list(itertools.chain(*g2_arr))
			fc = g1.mean() - np.mean(g2_arr) # Difference between Mean of small population and Mean of permuted population (Fold Change)
			p_arr.append([x, fc, permuted_pval])
		p_df = pd.DataFrame(data=p_arr, columns=['Gene', 'Ttest_FC', 'Ttest_Permuted_pval']).set_index('Gene')

		return p_df

	def _permuted_median_test(self, df, **kwargs):
		small_group = kwargs['small_group']
		large_group = kwargs['large_group']
		perm = kwargs['perm']

		assert perm>2, "Wrong Permutation number"

		test_genes = df.index.tolist()

		permuted_pool1 = self._permutation(small_group=small_group, large_group=large_group, perm=perm)
		permuted_pool2 = self._permutation(small_group=small_group, large_group=large_group, perm=perm)
		p_arr = []

		for x in test_genes:
			null_median = []
			g1 = df[small_group].loc[x].median()

			for p1,p2 in zip(permuted_pool1, permuted_pool2):  ##### Permutation table
				g2_1 = df[p1].loc[x]
				g2_2 = df[p2].loc[x]
				null_median.append(g2_1.median()-g2_2.median()) # Null distribution

			null_median = np.array(null_median)
			fc = g1.mean() - df[large_group].loc[x].mean() # Difference between Mean of small population and original large population (Fold Change)

			if fc >= 0:
				pval_cut = len(np.where(null_median>g1)[0])/float(len(null_median)) # Position of the median of small group in null dist
			else:
				pval_cut= len(np.where(null_median<g1)[0])/float(len(null_median)) # Position of the median of small group in null dist (opposite side) 2 tails

			p_arr.append([x, fc, pval_cut])

		p_df = pd.DataFrame(data=p_arr, columns=['Gene', 'Mtest_FC', 'Mtest_Permuted_pval']).set_index('Gene')
		return p_df

	def combined_pvalues(self, perm):
		assert perm>2, "Wrong Permutation number"

		t_result = self._permuted_ttest(self._df, small_group=self._small_group, large_group=self._large_group, perm=perm) # Permutated T-test
		m_result = self._permuted_median_test(self._df, small_group=self._small_group, large_group=self._large_group, perm=perm) # Permutated Median test
		combined_df = pd.concat([t_result, m_result], axis=1)

		# Merged 2 different tests in Z value space, but same hypothesis
		combined_df['Combined_pval'] = combined_df['Mtest_Permuted_pval'].apply(lambda x : stats.norm.ppf(1-x)) * combined_df['Ttest_Permuted_pval'].apply(lambda x : stats.norm.ppf(1-x))
		# Merged Z value to P-value
		combined_df['Combined_pval'] = combined_df['Combined_pval'].apply(lambda x : stats.norm.sf(abs(x))*2)
		#Storey FDR
		combined_df = self._storey_fdr(combined_df, p_name='Combined_pval')

		return combined_df


	def __init__(self, df, small_group, large_group):

		#####type checking
		assert df.empty==False, "Data is empty"
		assert type(small_group)==list, "Input group must be list type"
		assert type(large_group)==list, "Input group must be list type"
		assert len(small_group)<len(large_group), "Small group must be put in attribute in front of large group, ex) diff_stats(df, smallgroup, largegroup)"
		#####type checking

		self._df = df
		self._small_group = small_group
		self._large_group = large_group
