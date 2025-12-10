import argparse
import sys
import os
import pandas as pd
import bisect
import re
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
matplotlib.use('Agg')
def parse_args():
	parser = argparse.ArgumentParser(description='update pericentric position based on other satellite and barplot for each sample/chromosome')
	parser.add_argument('-anno', dest="fcenanno",
						help="bedfile for centromere annotation>",
						required=False)
	parser.add_argument('-pos', dest="fperipos",
						help="bedfile for pericentric position based on alpha satellite>",
						required=False)
	parser.add_argument('-fai', dest="fchromlen",
						help="index file generate by samtools faidx>",
						required=False)
	parser.add_argument('-o', dest="outfperipos",
						help="output bedfile for updated pericentric postion [default='out.bed']", default="out.bed",
						required=False)
	parser.add_argument('-p', dest="prefix",
						help="add sample name to figure title", default="",
						required=False)
	parser.add_argument('-plotOutdir', dest="outdir", help="outDir for plot pdf ", default="plot_output",
						required=False)


	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit()
	else:
		return parser.parse_args()


def stat_data(fcenanno, fchromlen, fperipos, outfperipos):
	data = pd.read_csv(fcenanno, sep='\t', header=None)
	data.columns = ['Chromosome', 'Start', 'End', 'Type', 'Value1', 'Strand', 'Value2', 'Value3', 'Color']

	chromosome_lengths = pd.read_csv(fchromlen, sep='\t', header=None)
	chromosome_lengths.columns = ['Chromosome', 'TotalLength', 'Start', 'End', 'Extra']
	# chrs = chromosome_lengths['Chromosome'].apply(lambda x: x.split('#')[2])
	# chromosome_lengths.drop('Chromosome', axis=1, inplace=True)
	# chromosome_lengths['Chromosome'] = chrs
	# print(chromosome_lengths.head())

	peri_data = pd.read_csv(fperipos, sep = '\t', header=None)
	peri_data.columns = ['Chromosome', 'Start', 'End']

	sortdata = data.groupby('Chromosome').apply(lambda x: x.sort_values(by='Start')).reset_index(drop=True)
	# print(sortdata.head())

	new_peri_start = {}
	new_peri_end = {}
	for chrom in sortdata['Chromosome'].unique():
		positions_list = sortdata.loc[sortdata['Chromosome'] == chrom][['Start', 'End']].values.tolist()
		if chrom in peri_data['Chromosome'].unique():

			peri_start = peri_data.loc[peri_data['Chromosome'] == chrom, 'Start'].values[0]
			peri_end = peri_data.loc[peri_data['Chromosome'] == chrom, 'End'].values[0]
			chrom_len = chromosome_lengths[chromosome_lengths['Chromosome'] == chrom]['TotalLength'].values[0]
			# if chrom in ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']:
			# 	new_peri_start[chrom] = 0
			# else:
			# 	for i in range(len(positions_list)):
			# 		anno_min_start = positions_list[i][0]
			# 		if  (anno_min_start < peri_start) and \
			# 		    (peri_start - anno_min_start <= 2000000) and \
			# 		    (positions_list[i][1] - positions_list[i][0] > 1000):
			# 			if anno_min_start - 2000000 > 0:
			# 				new_peri_start[chrom] = anno_min_start - 2000000
			#
			# 			else:
			# 				new_peri_start[chrom] = 0
			# 			break
			# 		elif anno_min_start > peri_start:
			# 			new_peri_start[chrom] = peri_start
			# 			break
			# 		else:
			# 			pass
			iterflag = True
			while iterflag:
				index = bisect.bisect_left([positions_list[i][0] for i in range(len(positions_list))], peri_start)
				if bool(re.search(r'(chr13|chr14|chr15|chr21|chr22)', chrom)):
					new_peri_start[chrom] = 0
					break
				elif index == 0:
					new_peri_start[chrom] = peri_start
					break
				else:
					for i in range(index, 0, -1):
						anno_start = positions_list[i-1][0]
						if (peri_start-anno_start <= 2000000 ) or (positions_list[i-1][1] - positions_list[i-1][0] > 10000):
							if anno_start > 100000:
								peri_start = anno_start - 100000
								break
							else:
								new_peri_start[chrom] = 0
								iterflag = False
								break
						else:
							iterflag = False
							new_peri_start[chrom] = peri_start - 1900000
							break


			iterflag = True
			while iterflag:
				index = bisect.bisect_right([positions_list[i][1] for i in range(len(positions_list))], peri_end)
				if index == len(positions_list):
					new_peri_end[chrom] = peri_end
					break
				else:
					for i in range(index, len(positions_list)):
						anno_end = positions_list[i][1]
						if (anno_end - peri_end <= 2000000) or (positions_list[i][1] - positions_list[i][0] > 10000):
							if anno_end + 100000 < chrom_len:
								peri_end  = anno_end + 100000
								break
								# new_peri_end[chrom] = anno_max_end + 2000000
							else:
								peri_end = chrom_len
								new_peri_end[chrom] = chrom_len
								iterflag = False
								break
						else:
							iterflag = False
							new_peri_end[chrom] = peri_end + 1900000
							break



	peri_data['Start'] = peri_data['Chromosome'].map(new_peri_start)
	peri_data['End'] = peri_data['Chromosome'].map(new_peri_end)

	peri_data.to_csv(outfperipos, sep = "\t", header = False, index=False)
	return data, chromosome_lengths, peri_data

def sort_chrom(chromosome):
	if re.match(r"(.+)chr(.+)", chromosome):
		try:
			return int(re.search(r"(.+)chr(.+)", chromosome).group(2))
		except ValueError:
			return float('inf')
	else:
		return float('inf')

def chrom_plot(prefix, outdir, annodf, chromlendf, peridf):
	print(annodf.head())

	# formatter = ticker.ScalarFormatter(useMathText=True)
	# formatter.set_scientific(True)
	# formatter.set_powerlimits((-3, 6))
	# plt.gca().xaxis.set_major_formatter(formatter)

	types = annodf['Type'].unique()
	chromosomes = annodf['Chromosome'].unique()

	plt.figure(figsize=(15, int(len(chromosomes)*0.25)))

	color_map = dict(zip(annodf['Type'], annodf['Color']))
	# print(color_map)

	# peridf = peridf.sort_values(by='Chromosome').copy()
	for index, row in peridf.iterrows():
		chrom = row['Chromosome']
		start = row['Start'] / 1000000
		end = (row['End'] - row['Start']) / 1000000
		plt.barh(chrom, end, left=start, color='lightgrey')

	# 绘制每个类型的satellite
	for chrom in chromosomes:
		chrom_data = annodf[annodf['Chromosome'] == chrom]
		for index, row in chrom_data.iterrows():
			start = row['Start'] / 1000000
			end = (row['End'] - row['Start']) / 1000000
			color = color_map[row['Type']]
			plt.barh(chrom, end, left=start, color=color)

	# chromlendf = chromlendf.sort_values(by='Chromosome').copy()
	for index, row in chromlendf.iterrows():
		chrom = row['Chromosome']
		length = row['TotalLength']  / 1000000
		plt.barh(chrom, length,  edgecolor='black', fill=False, linewidth=1.5)

	plt.tick_params(axis='y', labelsize=10)
	newlabels = [re.search(r"(.+)?(chr.+)", i).group(2) if re.match(r"(.+)?(chr.+)", i) else i for i in chromosomes]
	plt.yticks(ticks=list(chromosomes), labels=newlabels)
	plt.xlabel('Chromosome Length(Mb)')
	plt.ylabel('Chromosome')
	plt.title('Centromere Position on Chromosome (%s)' % (prefix))

	legend_patches = []
	new_color_map = {('ASat' if ty.startswith('S') else ty): color for ty, color in color_map.items()}
	for ty, color in new_color_map.items():
		legend_patches.append(mpatches.Patch(color=color, label=ty))
	plt.legend(handles=legend_patches, loc='center right')

	# 显示图
	plt.savefig("%s/%s.genomewide.cent.round1.barplot.png" % (outdir, prefix), dpi=300)
	# plt.show()

def single_plot(outdir, annodf, peridf, chrom, prefix):
	iannodf = annodf[annodf['Chromosome'] == chrom].copy()
	iperidf = peridf[peridf['Chromosome'] == chrom].copy()

	divisor = 1000000
	min_pos = (iperidf['Start'].values[0] // divisor) * divisor
	max_pos = (iperidf['End'].values[0] // divisor) * divisor

	subannodf = iannodf[(iannodf['Start'] >= min_pos) & (iannodf['End'] <= max_pos)].copy()
	# types = subannodf['Type'].unique()
	# color_map = {type: color for type, color in zip(types, subannodf['Color'].unique())}

	scale = 25000
	subannodf.loc[:,'adj_start'] = (subannodf['Start'] - min_pos) / scale
	subannodf.loc[:,'adj_end'] = (subannodf['End'] - min_pos) / scale
	iperidf.loc[:,'adj_start'] = (iperidf['Start']- min_pos) / scale
	iperidf.loc[:,'adj_end'] = (iperidf['End'] - min_pos) / scale
	# print(subannodf.head())

	plt.figure(figsize=(8, 2))
	plt.ylim([298, 302])
	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
	plt.gca().spines['left'].set_visible(False)
	plt.gca().get_yaxis().set_visible(False)
	plt.gca().spines['bottom'].set_position(('data', 299.5))

	for index, row in iperidf.iterrows():
		start = row['adj_start']
		end = row['adj_end'] - row['adj_start']
		plt.barh(300, end, left=start, color='lightgrey', height=0.5)

	for index, row in subannodf.iterrows():
		start = row['adj_start']
		end = row['adj_end'] - row['adj_start']
		plt.barh(300, end, left=start, color=row['Color'], height=0.5)

	original_x_ticks = plt.xticks()
	new_x_ticks = [round(((i * scale + min_pos) / 1000000),2)for i in original_x_ticks[0]]
	plt.xticks(original_x_ticks[0], new_x_ticks)
	plt.title('%s:%s-%s' % ( chrom, iperidf['Start'].values[0], iperidf['End'].values[0]), y=0.6)
	# 添加图例
	# legend_patches = []
	# for type, color in color_map.items():
	# 	legend_patches.append(mpatches.Patch(color=color, label=type))
	# plt.legend(handles=legend_patches, bbox_to_anchor=(1, 0.2), loc=3, borderaxespad=0)
	plt.savefig("%s/%s_%s.round1.barplot.pdf" % (outdir, prefix, chrom))
	plt.cla()
	plt.close("all")
	# plt.show()

if __name__ == "__main__":
	args = parse_args()
	annodf, chromlendf, peridf = stat_data(args.fcenanno, args.fchromlen, args.fperipos, args.outfperipos)
	print("[log]: updated pericentric postion done!")

	if not os.path.exists(args.outdir):
		os.makedirs(args.outdir)

	# filter chrM and another contigs#
	chromlendf = chromlendf[chromlendf['Chromosome'].str.contains('chr') & ~chromlendf['Chromosome'].str.contains('chrM')]
	peridf = peridf[peridf['Chromosome'].str.contains('chr') & ~chromlendf['Chromosome'].str.contains('chrM')]
	peridf = peridf.iloc[peridf['Chromosome'].map(sort_chrom).argsort()].copy()
	chromlendf = chromlendf.iloc[chromlendf['Chromosome'].map(sort_chrom).argsort()].copy()

	chrom_plot(args.prefix, args.outdir, annodf, chromlendf, peridf)
	for chrom in peridf['Chromosome'].unique():
		single_plot( args.outdir, annodf, peridf, chrom, args.prefix)
		print("[log]: %s pericentric distribution plot done!" % (chrom))
	# single_plot(prefix, outdir, annodf, peridf, 'chr6')

