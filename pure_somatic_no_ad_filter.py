#!/usr/bin/env python
import sys
import pandas as pd
from pyfaidx import Fasta

usage = '\t --------\n' \
        '\t| usage  : python pure_somatic_no_ad_filter.py f1 f2 f3\n' \
        '\t| input  : f1 = mutect2_variants.tsv\n' \
        '\t| input  : f2 = reference_genome.fasta (and associated .fai in same directory)\n' \
        '\t| input  : f3 = depth_per_pos_original_bam.tsv\n' \
        '\t| output : pure_somatic_no_AD_filter_CTI_spe.tsv\n' \
        '\t| output : pure_somatic_no_AD_filter_TT_spe.tsv\n' \
        '\t| output : pure_somatic_no_AD_filter_in_both_tumors.tsv\n' \
        '\t| output : pure_somatic_no_AD_filter_across_patients_and_tumors_legit.tsv\n' \
        '\t| output : pure_somatic_no_AD_filter_across_patients_and_tumors_all.tsv\n' \
        '\t --------'

if len(sys.argv) != 4:
    print(usage)
    sys.exit()

##################
### Keep_legit ###
##################
# Load variants into a dataframe
#variants_df = pd.read_csv(sys.argv[1], sep='\t')
var_df = pd.read_csv('./all.tsv', sep='\t')
# df_var contains among others flags "alt_allele_in_normal" & "multi_event_alt_allele_in_normal", remove associated variants if Allele_Frequency_Normal > 0.01 (1%)
var1_df = var_df[(var_df.Flag == 'PASS') | (var_df.Allele_Frequency_Normal <= 0.01)]
# Remove if potential contamination for the call
var2_df = var1_df[var1_df.Contamination != 'plausible']
# Need initial BAM depth > 99 for both tumor_sample & normal_sample
var3_df = var2_df[(var2_df.Reads_Tumor > 99) & (var2_df.Reads_Normal > 99)]
# Need AF > 0.1 (10%), otherwise it's problematic to validate in Sanger
var4_df = var3_df[var3_df.Allele_Frequency_Tumor >= 0.1]
# Need alt-allele depth in tumor >= 10 (MuTect2) (with Sanger in mind too)
var5_df = var4_df[var4_df.AD_Tumor_Alt > 9]
# Need coherent SAC field values (don't account for strand biais here, but ask for at least some depth : sum of Alt+/- SAC should be >= 10, which is expected to be coherent with previous filtering on AD_Tumor_Alt)
var6_df = var5_df[(((var5_df['Tumor_Alt+']) + (var5_df['Tumor_Alt-'])) > 9)]
# Remove variants out of BED position
var7_df = var6_df[~var6_df.Amplicon_position.str.contains('position')]

# Get -3/+3 genome context for Ref & Var (pyfaidx module) and report it in new columns
# Then for InDel, remove variants which may originate from Ion Proton homopolymer biais
# Ressources : goo.gl/sbjx4t, goo.gl/uMEMXq
# Nb : No Type = MNP or CLUMPED in our data, thus we don' try to handle them
#hg19_faidx = Fasta(sys.argv[2], sequence_always_upper=True)
hg19_faidx = Fasta('../../Data/Hg19/ucsc.hg19.fasta', sequence_always_upper=True)
var7_df.insert(len(var7_df.columns), 'Context_Var', None)
var7_df.insert(len(var7_df.columns), 'Context_Ref', None)
var7_snp_df = var7_df[var7_df.Type.str.contains('SNP')]
var7_ins_df = var7_df[var7_df.Type.str.contains(' ins')]
var7_del_df = var7_df[var7_df.Type.str.contains(' del')]

for index, row in var7_snp_df.iterrows():
	seq_ref_li = list(hg19_faidx[row["Chrom"]][(int(row["Position"])-4):(int(row["Position"])+3)].seq)
	seq_var_li = seq_ref_li[:]
	seq_var_li[3] = row["Var"]
	var7_df.at[index, 'Context_Ref'] = "".join(seq_ref_li)
	var7_df.at[index, 'Context_Var'] = "".join(seq_var_li)

drop_indel_set = set()
for index, row in var7_ins_df.iterrows():
	seq_ref_li = list(hg19_faidx[row["Chrom"]][(int(row["Position"])-4):(int(row["Position"])+3)].seq)
	seq_var_li = seq_ref_li[:]
	del seq_var_li[3]
	var_li = list(row["Var"])
	seq_var_li[3:3] = var_li
	var7_df.at[index, 'Context_Ref'] = "".join(seq_ref_li)
	var7_df.at[index, 'Context_Var'] = "".join(seq_var_li)
	if seq_var_li[1:4].count(seq_var_li[1]) == 3:
		drop_indel_set.add(index)
	if seq_var_li[(2+len(var_li)):(5+len(var_li))].count(seq_var_li[(2+len(var_li))]) == 3:
		drop_indel_set.add(index)

for index, row in var7_del_df.iterrows():
	seq_ref_li = list(hg19_faidx[row["Chrom"]][(int(row["Position"])-4):(int(row["Position"])+(3+(len(row["Ref"])-1)))].seq)
	seq_var_li = seq_ref_li[:]
	seq_var_li[3:(3+len(row["Ref"]))] = row["Var"]
	var7_df.at[index, 'Context_Ref'] = "".join(seq_ref_li)
	var7_df.at[index, 'Context_Var'] = "".join(seq_var_li)
	if seq_var_li[1:4].count(seq_var_li[1]) == 3:
		drop_indel_set.add(index)
	if seq_var_li[3:6].count(seq_var_li[3]) == 3:
		drop_indel_set.add(index)

var8_df = var7_df.drop(drop_indel_set)

# In initial BAM, we want enough depth (>= 100) in associated CTI/TT (computed with GATK DepthOfCoverage) and we report this depth in a new column named "Reads_Other_Tumor"
#depth_df = pd.read_csv(sys.argv[3], sep='\t')
depth_df = pd.read_csv('../0_Controls/depth_per_position.csv', sep='\t')
drop_depth_set = set()
var8_df.insert(len(var8_df.columns), 'Reads_Other_Tumor', None)
coherence_name_table={'108885':'p1','12A13991':'p2','12A17639':'p3','14A02580':'p4','14A13198':'p5','14A13351':'p6','15A14':'p7','15A04217':'p8','15A12637':'p9','15A18359':'p10','10016297':'p11','14A15780':'p12','14A19434':'p13'}
depth_idx_df = depth_df.set_index('Locus')
for index, row in var8_df.iterrows():
	if row["Tissu(TT|CTI)"] == "TT":
		depth_CTI = depth_idx_df.loc[row["Chrom"]+":"+str(row["Position"])]["Depth_for_"+coherence_name_table[row["Patient"]]+"_CTI"]
		var8_df.at[index, 'Reads_Other_Tumor'] = depth_CTI
		if depth_CTI < 100:
			drop_depth_set.add(index)
	else:
		depth_TT = depth_idx_df.loc[row["Chrom"]+":"+str(row["Position"])]["Depth_for_"+coherence_name_table[row["Patient"]]+"_TT"]
		var8_df.at[index, 'Reads_Other_Tumor'] = depth_TT
		if depth_TT < 100:
			drop_depth_set.add(index)

legit_df = var8_df.drop(drop_depth_set)

###########################
### Legit filters recap ###
###########################
print(var_df.shape, depth_df.shape, len(var1_df), len(var2_df), len(var3_df), len(var4_df), len(var5_df), len(var6_df), len(var7_df), len(var8_df), len(legit_df))

#####################
### CTI_or_TT_spe ###
#####################
## Has to be done in a patient + variant_position specific manner
# Create non-legit dic
non_legit_df = var_df.drop(set(legit_df.index.values.tolist()))
# For non_legit_df, use .groupby() to get some kind of "multiindex" on long_key, then count the number of element for each long_key, store result in a dic with key = long_key (Patient_Chrom_Position_Tissu) and value = count (1->4)
long_key_non_legit_dic = {}
non_legit_series = non_legit_df.groupby(['Patient', 'Chrom', 'Position', 'Tissu(TT|CTI)'])['Version(1|2)'].size()
for index, value in non_legit_series.iteritems():
	long_key = '_'.join(map(str, index))
	if long_key in long_key_non_legit_dic:
		long_key_non_legit_dic[long_key]+=1
	else:
		long_key_non_legit_dic[long_key]=1

# Create long_key_legit_dic for legit_df, with value = count (1->4)
# And long_key_legit_seen_in_non_dic to store number of associated non-legit variants in the other tumor sample
long_key_legit_dic = {}
long_key_legit_seen_in_non_dic = {}
for index, row in legit_df.iterrows():
	long_key = str(row['Patient']) + '_' + str(row['Chrom']) + '_' + str(row['Position']) + '_' + str(row['Tissu(TT|CTI)'])
	long_key_legit_seen_in_non_dic.setdefault(long_key, None)
	if long_key in long_key_legit_dic:
		long_key_legit_dic[long_key]+=1
	else:
		long_key_legit_dic[long_key]=1

# For each legit variant, check if the same event was seen in the other non-legit tumoral tissu variants. If so, add value to dic
for key, value in long_key_legit_seen_in_non_dic.items():
	parts = key.split('_')
	if (row['Tissu(TT|CTI)'] == 'TT'):
		long_key_other = '_'.join(map(str,parts[:-1])) + '_CTI'
		if long_key_other in long_key_non_legit_dic:
			long_key_legit_seen_in_non_dic[key]=(long_key_non_legit_dic[long_key_other])
	else:
		long_key_other = '_'.join(map(str,parts[:-1])) + '_TT'
		if long_key_other in long_key_non_legit_dic:
			long_key_legit_seen_in_non_dic[key]=(long_key_non_legit_dic[long_key_other])

# Change None value in tissu_spe dic to 0
for key, value in long_key_legit_seen_in_non_dic.items():
	if value == None:
		long_key_legit_seen_in_non_dic[key] = 0

# Insert and fill new columns in legit_df (number of the same legit events & number of non-legit event in the other tumoral sample)
legit_df.insert(len(legit_df.columns), 'Occurence_of_this_legit_event(max 4)', None)
legit_df.insert(len(legit_df.columns), 'Seen_in_non-legit_2nd_tissu(max 4)', None)
for index, row in legit_df.iterrows():
	long_key = str(row['Patient']) + '_' + str(row['Chrom']) + '_' + str(row['Position']) + '_' + str(row['Tissu(TT|CTI)'])
	legit_df.at[index, 'Seen_in_non-legit_2nd_tissu(max 4)'] = long_key_legit_seen_in_non_dic[long_key]
	legit_df.at[index, 'Occurence_of_this_legit_event(max 4)'] = long_key_legit_dic[long_key]

# Create a set containing events seen in the other tumoral tissu for the same Patient_Chrom_Position (checked for each event in legit_df)
both_idx_set = set()
for index, row in legit_df.iterrows():
	if (row['Tissu(TT|CTI)'] == 'TT'):
		long_key_other = str(row['Patient']) + '_' + str(row['Chrom']) + '_' + str(row['Position']) + '_CTI'
	else:
		long_key_other = str(row['Patient']) + '_' + str(row['Chrom']) + '_' + str(row['Position']) + '_TT'
	if long_key_other in long_key_legit_dic:
		both_idx_set.add(index)

# Use previous set to obtain a dataset of tissu specific variants
tissu_spe_df = legit_df.drop(both_idx_set)

# Keep only one (best) patient_spe variant if multiple
# Find the best ones
tissu_spe_best_idx_dic = {}
for index, row in tissu_spe_df.iterrows():
	long_key = str(row['Patient']) + '_' + str(row['Chrom']) + '_' + str(row['Position']) + '_' + str(row['Tissu(TT|CTI)'])
	if not long_key in tissu_spe_best_idx_dic:
		tissu_spe_best_idx_dic[long_key] = index
	else:
		worse = True
		while worse:
			if row['Flag'].count(';') < tissu_spe_df.loc[tissu_spe_best_idx_dic[long_key]].Flag.count(';'):
				tissu_spe_best_idx_dic[long_key] = index
				worse = False
				break
			if (row['Flag'] == 'PASS') and (tissu_spe_df.loc[tissu_spe_best_idx_dic[long_key]].Flag != 'PASS'):
				tissu_spe_best_idx_dic[long_key] = index
				worse = False
				break
			break
			if row['Allele_Frequency_Tumor'] > tissu_spe_df.loc[tissu_spe_best_idx_dic[long_key]].Allele_Frequency_Tumor:
				tissu_spe_best_idx_dic[long_key] = index
				worse = False
				break
			if row['AD_Tumor_Alt'] > tissu_spe_df.loc[tissu_spe_best_idx_dic[long_key]].AD_Tumor_Alt:
				tissu_spe_best_idx_dic[long_key] = index
				worse = False
				break
			break

# Create idx set from best ones and use this set to create "tissu_spe_best_df"
tissu_spe_best_idx_set = set()
for key, value in tissu_spe_best_idx_dic.items():
	tissu_spe_best_idx_set.add(value)

tissu_spe_best_df = tissu_spe_df.loc[tissu_spe_best_idx_set]

# Output corresponding CTI and TT specific tables, with columns reordering & sorting
tissu_spe_best_df = tissu_spe_best_df[['Patient', 'Chrom', 'Position', 'Ref', 'Var', 'Allele_Frequency_Tumor', 'Gene_ID', 'Func', 'Exo_func', 'AAChange', 'Cosmic', 'Clinvar (CLINSIG/CLNDBN/CLNACC/CLNDSDB/CLNDSDBID)', 'Flag', 'Occurence_of_this_legit_event(max 4)', 'Seen_in_non-legit_2nd_tissu(max 4)', 'Type', 'Context_Ref', 'Context_Var', 'Amplicon_position', 'AD_Tumor_Ref', 'AD_Tumor_Alt', 'Reads_Tumor', 'Reads_Other_Tumor', 'Tumor_Ref+', 'Tumor_Ref-', 'Tumor_Alt+', 'Tumor_Alt-', 'Allele_Frequency_Normal', 'AD_Normal_Ref', 'AD_Normal_Alt', 'Reads_Normal', 'Normal_Ref+', 'Normal_Ref-', 'Normal_Alt+', 'Normal_Alt-', 'Version(1|2)', 'MinPruning', 'Tissu(TT|CTI)', 'Cytoband', 'Contamination']]
tissu_spe_best_srt_df = tissu_spe_best_df.sort_values(['Patient', 'Chrom', 'Position'], ascending=[True, True, True])
CTI_spe_df = tissu_spe_best_srt_df[tissu_spe_best_srt_df['Tissu(TT|CTI)'] == 'CTI']
TT_spe_df = tissu_spe_best_srt_df[tissu_spe_best_srt_df['Tissu(TT|CTI)'] == 'TT']
CTI_spe_df.to_csv('./pure_somatic_no_AD_filter_CTI_spe.tsv', sep='\t', index=False)
TT_spe_df.to_csv('./pure_somatic_no_AD_filter_TT_spe.tsv', sep='\t', index=False)

####################
### Seen_in_both ###
####################
# Temp otuput of raw data
both_legit_df = legit_df.loc[both_idx_set]
both_legit_df = both_legit_df[['Patient', 'Chrom', 'Position', 'Ref', 'Var', 'Allele_Frequency_Tumor', 'Gene_ID', 'Func', 'Exo_func', 'AAChange', 'Cosmic', 'Clinvar (CLINSIG/CLNDBN/CLNACC/CLNDSDB/CLNDSDBID)', 'Flag', 'Occurence_of_this_legit_event(max 4)', 'Seen_in_non-legit_2nd_tissu(max 4)', 'Type', 'Context_Ref', 'Context_Var', 'Amplicon_position', 'AD_Tumor_Ref', 'AD_Tumor_Alt', 'Reads_Tumor', 'Reads_Other_Tumor', 'Tumor_Ref+', 'Tumor_Ref-', 'Tumor_Alt+', 'Tumor_Alt-', 'Allele_Frequency_Normal', 'AD_Normal_Ref', 'AD_Normal_Alt', 'Reads_Normal', 'Normal_Ref+', 'Normal_Ref-', 'Normal_Alt+', 'Normal_Alt-', 'Version(1|2)', 'MinPruning', 'Tissu(TT|CTI)', 'Cytoband', 'Contamination']]
both_legit_srt_df = both_legit_df.sort_values(['Patient', 'Chrom', 'Position'], ascending=[True, True, True])
both_legit_srt_df.to_csv('./pure_somatic_no_AD_filter_in_both_tumors.tsv', sep='\t', index=False)

# Do better stuff here : multiindex try
#both_legit_idx_df = both_legit_df.set_index(['Patient', 'Chrom', 'Position'])
#both_legit_idx_df = both_legit_idx_df.sort_index()
#for index, row in both_legit_idx_df.iterrows():
#	print(index)
#	print(len(row))

# Do better stuff here : regular try with inc heavy dics
#for index, row in both_legit_df.iterrows():
#	long_key = str(row['Patient']) + '_' + str(row['Chrom']) + '_' + str(row['Position']) + '_' + str(row['Tissu(TT|CTI)'])
#	long_key_other = ''
#	if (row['Tissu(TT|CTI)'] == 'TT'):
#		long_key_other = '_'.join(map(str,parts[:-1])) + '_CTI'
#	else:
#		long_key_other = '_'.join(map(str,parts[:-1])) + '_TT'

###########################
### Seen_across_samples ###
###########################
## Perform this for all_variants & for legit_variants (looking for the most seen non-patient, non-tissu spe variants)
# Use groupby + size on "key" = chr+pos to compute associated occurence
occ_all_series = var_df.groupby(['Chrom', 'Position'])['Version(1|2)'].size()
occ_legit_series = legit_df.groupby(['Chrom', 'Position'])['Reads_Normal'].size()
# Define top quantile and add True to variants greater than (gt) quantile
top_quantile_all_series = occ_all_series.gt(occ_all_series.quantile(q=0.9975, interpolation='higher'))

top_quantile_legit_series = occ_legit_series.gt(occ_legit_series.quantile(q=0.9, interpolation='nearest'))
# Remove variants being False in top_quantile series (cool synthax)
top_quantile_all_series = top_quantile_all_series[top_quantile_all_series]
top_quantile_legit_series = top_quantile_legit_series[top_quantile_legit_series]
# In the occurence series, keep only chr+pos keys being True in top_quantile series
# Use indexes sets intersection to keep good candidate keys in first sets (the one associated with occ series), then reindex accordingly
all_intersec_set = set(occ_all_series.index).intersection(set(top_quantile_all_series.index))
legit_intersec_set = set(occ_legit_series.index).intersection(set(top_quantile_legit_series.index))
occ_all_top_quantile_series = occ_all_series.reindex(all_intersec_set)
occ_legit_top_quantile_series = occ_legit_series.reindex(legit_intersec_set)
# Now, iterate one last time on initial dfs and sample one random row to illustrate each top_quantile chr+pos
across_samples_all_df = pd.DataFrame(columns=var_df.columns)
for key in occ_all_top_quantile_series.keys():
	across_samples_all_df = across_samples_all_df.append((var_df[(var_df.Chrom == key[0]) & (var_df.Position == key[1])]).sample(n=1))

across_samples_legit_df = pd.DataFrame(columns=legit_df.columns)
for key in occ_legit_top_quantile_series.keys():
	across_samples_legit_df = across_samples_legit_df.append((legit_df[(legit_df.Chrom == key[0]) & (legit_df.Position == key[1])]).sample(n=1))

# Add occurence row in output dfs
across_samples_all_df.insert(0, 'Occurence_in_all_patients_and_tumors', None)
for index, row in across_samples_all_df.iterrows():
	across_samples_all_df.at[index, 'Occurence_in_all_patients_and_tumors'] = occ_all_top_quantile_series.loc[row['Chrom'], row['Position']]

across_samples_legit_df.insert(0, 'Occurence_in_all_patients_and_tumors', None)
for index, row in across_samples_legit_df.iterrows():
	across_samples_legit_df.at[index, 'Occurence_in_all_patients_and_tumors'] = occ_legit_top_quantile_series.loc[row['Chrom'], row['Position']]

# Remove patient specific rows
across_samples_legit_df = across_samples_legit_df.drop(['Occurence_of_this_legit_event(max 4)', 'Seen_in_non-legit_2nd_tissu(max 4)'], axis=1)

# Output across_samples variants with columns sorting
across_samples_all_df = across_samples_all_df[['Occurence_in_all_patients_and_tumors', 'Patient', 'Chrom', 'Position', 'Ref', 'Var', 'Allele_Frequency_Tumor', 'Gene_ID', 'Func', 'Exo_func', 'AAChange', 'Cosmic', 'Clinvar (CLINSIG/CLNDBN/CLNACC/CLNDSDB/CLNDSDBID)', 'Flag', 'Occurence_of_this_legit_event(max 4)', 'Seen_in_non-legit_2nd_tissu(max 4)', 'Type', 'Amplicon_position', 'AD_Tumor_Ref', 'AD_Tumor_Alt', 'Reads_Tumor', 'Tumor_Ref+', 'Tumor_Ref-', 'Tumor_Alt+', 'Tumor_Alt-', 'Allele_Frequency_Normal', 'AD_Normal_Ref', 'AD_Normal_Alt', 'Reads_Normal', 'Normal_Ref+', 'Normal_Ref-', 'Normal_Alt+', 'Normal_Alt-', 'Version(1|2)', 'MinPruning', 'Tissu(TT|CTI)', 'Cytoband', 'Contamination']]
across_samples_legit_df = across_samples_legit_df[['Occurence_in_all_patients_and_tumors', 'Patient', 'Chrom', 'Position', 'Ref', 'Var', 'Allele_Frequency_Tumor', 'Gene_ID', 'Func', 'Exo_func', 'AAChange', 'Cosmic', 'Clinvar (CLINSIG/CLNDBN/CLNACC/CLNDSDB/CLNDSDBID)', 'Flag', 'Type', 'Context_Ref', 'Context_Var', 'Amplicon_position', 'AD_Tumor_Ref', 'AD_Tumor_Alt', 'Reads_Tumor', 'Reads_Other_Tumor', 'Tumor_Ref+', 'Tumor_Ref-', 'Tumor_Alt+', 'Tumor_Alt-', 'Allele_Frequency_Normal', 'AD_Normal_Ref', 'AD_Normal_Alt', 'Reads_Normal', 'Normal_Ref+', 'Normal_Ref-', 'Normal_Alt+', 'Normal_Alt-', 'Version(1|2)', 'MinPruning', 'Tissu(TT|CTI)', 'Cytoband', 'Contamination']]

across_samples_all_srt_df = across_samples_all_df.sort_values(['Chrom', 'Position'], ascending=[True, True])
across_samples_legit_srt_df = across_samples_legit_df.sort_values(['Chrom', 'Position'], ascending=[True, True])
across_samples_all_srt_df.to_csv('./pure_somatic_no_AD_filter_across_patients_and_tumors_all.tsv', sep='\t', index=False)
across_samples_legit_srt_df.to_csv('./pure_somatic_no_AD_filter_across_patients_and_tumors_legit.tsv', sep='\t', index=False)



### DA) If filter on allelic disequilibrium
# Discard if strand biais for tumor reads supporting alt-allele is worse than 65/35
# Not good enough, need to account for the number of reads too plz
#df4 = df3[(df3['Tumor_Alt+'] >= 35/100*(df3.AD_Tumor_Alt)) & (df3['Tumor_Alt-'] >= 35/100*(df3.AD_Tumor_Alt))]
