#!/usr/bin/python


# python peak_bedgraph2.py win genome inputfile "wig_file_name"

# align bed file to wig format chr-by-chr using dictionary and list

import sys
from time import time, localtime, strftime
import concurrent.futures
import numpy as np
import math


def peak_seq(inputfile, name):
    # genome chrom sizes
    dm3 = {
        "chr2L": 23011544,
        "chr2LHet": 368872,
        "chr2R": 21146708,
        "chr2RHet": 3288761,
        "chr3L": 24543557,
        "chr3LHet": 2555491,
        "chr3R": 27905053,
        "chr3RHet": 2517507,
        "chr4": 1351857,
        "chrM": 16569,
        "chrX": 22422827,
        "chrXHet": 204112,
        "chrYHet": 347038,
        "chrU": 10049037,
        "chrUextra": 29004656,
    }

    mm8 = {
        "chr1": 197069962,
        "chr2": 181976762,
        "chr3": 159872112,
        "chr4": 155029701,
        "chr5": 152003063,
        "chr6": 149525685,
        "chr7": 145134094,
        "chr8": 132085098,
        "chr9": 124000669,
        "chr10": 129959148,
        "chr11": 121798632,
        "chr12": 120463159,
        "chr13": 120614378,
        "chr14": 123978870,
        "chr15": 103492577,
        "chr16": 98252459,
        "chr17": 95177420,
        "chr18": 90736837,
        "chr19": 61321190,
        "chrM": 16299,
        "chrX": 165556469,
        "chrY": 16029404,
    }

    mm9 = {
        "chr1": 197195432,
        "chr2": 181748087,
        "chr3": 159599783,
        "chr4": 155630120,
        "chr5": 152537259,
        "chr6": 149517037,
        "chr7": 152524553,
        "chr8": 131738871,
        "chr9": 124076172,
        "chr10": 129993255,
        "chr11": 121843856,
        "chr12": 121257530,
        "chr13": 120284312,
        "chr14": 125194864,
        "chr15": 103494974,
        "chr16": 98319150,
        "chr17": 95272651,
        "chr18": 90772031,
        "chr19": 61342430,
        "chrM": 16299,
        "chrX": 166650296,
        "chrY": 15902555,
    }

    hg18 = {
        "chr1": 247249719,
        "chr2": 242951149,
        "chr3": 199501827,
        "chr4": 191273063,
        "chr5": 180857866,
        "chr6": 170899992,
        "chr7": 158821424,
        "chr8": 146274826,
        "chr9": 140273252,
        "chr10": 135374737,
        "chr11": 134452384,
        "chr12": 132349534,
        "chr13": 114142980,
        "chr14": 106368585,
        "chr15": 100338915,
        "chr16": 88827254,
        "chr17": 78774742,
        "chr18": 76117153,
        "chr19": 63811651,
        "chr20": 62435964,
        "chr21": 46944323,
        "chr22": 49691432,
        "chrM": 16569,
        "chrX": 154913754,
        "chrY": 57772954,
    }

    hg19 = {
        "chr1": 249250621,
        "chr2": 243199373,
        "chr3": 198022430,
        "chr4": 191154276,
        "chr5": 180915260,
        "chr6": 171115067,
        "chr7": 159138663,
        "chr8": 146364022,
        "chr9": 141213431,
        "chr10": 135534747,
        "chr11": 135006516,
        "chr12": 133851895,
        "chr13": 115169878,
        "chr14": 107349540,
        "chr15": 102531392,
        "chr16": 90354753,
        "chr17": 81195210,
        "chr18": 78077248,
        "chr19": 59128983,
        "chr20": 63025520,
        "chr21": 48129895,
        "chr22": 51304566,
        "chrM": 16571,
        "chrX": 155270560,
        "chrY": 59373566,
    }

    hg38 = {
        "chr1": 248956422,
        "chr2": 242193529,
        "chr3": 198295559,
        "chr4": 190214555,
        "chr5": 181538259,
        "chr6": 170805979,
        "chr7": 159345973,
        "chrX": 156040895,
        "chr8": 145138636,
        "chr9": 138394717,
        "chr11": 135086622,
        "chr10": 133797422,
        "chr12": 133275309,
        "chr13": 114364328,
        "chr14": 107043718,
        "chr15": 101991189,
        "chr16": 90338345,
        "chr17": 83257441,
        "chr18": 80373285,
        "chr20": 64444167,
        "chr19": 58617616,
        "chrY": 57227415,
        "chr22": 50818468,
        "chr21": 46709983,
        "chrM": 16569,
    }

    win = int(25)
    genome = "hg38"

    if genome == "dm3":
        genome = dm3

    if genome == "mm8":
        genome = mm8

    if genome == "mm9":
        genome = mm9

    if genome == "hg18":
        genome = hg18

    if genome == "hg19":
        genome = hg19

    if genome == "hg38":
        genome = hg38

    currtime = time()

    chrom = []  # chr name
    chroms = []  # chr names
    ext_dict = {}  # inputfile
    ext_refseq = {}  # refseq
    s = {}  #
    peaklist = {}  #

    for i in genome.keys():
        i = i.lstrip("chr")
        if i.isdigit():
            i = int(i)
        chrom.append(i)

    for i in chrom:
        chrom = "chr" + str(i)
        chroms.append(chrom)
        ext_refseq[chrom] = []
        ext_dict[chrom] = []
        ext_dict[chrom].extend([0] * (int(genome.get(chrom)) // win))
        s[chrom] = []
        peaklist[chrom] = []

    # print('time elapsed:', time() - currtime)

    # print('Importing alignment file ...')

    currtime = time()

    for line in inputfile:
        a = line.rstrip().split("\t")
        chrom, start, end, num = a[0], int(a[1]) + 1, int(a[2]) + 1, float(a[3])
        if chrom in chroms:
            ext_refseq[chrom].append([start, end, num])

    # inputfile = inputfile[0]
    # for line in open(inputfile, 'r').readlines():
    #     a = line.rstrip().split('\t')
    #     chrom, start, end, num = a[0], int(a[1])+1, int(a[2])+1, float(a[3])
    #     if chrom in chroms:
    #         ext_refseq[chrom].append([start, end, num])

    # print('Total sequences =', sum(len(ext_refseq[i]) for i in chroms))
    # print('time elapsed:', time() - currtime)

    # print('Coverting alignments ...')

    currtime = time()

    for chrom in chroms:
        if len(ext_refseq[chrom]) != 0:
            for jj in range(len(ext_refseq[chrom])):
                start, end, num = ext_refseq[chrom][jj]
                # 		if start < 0: # start pos less than 0
                # 			ext_dict[chrom][end//win] += (end%win)*num
                # 			for j in range(0, end//win, 1):
                # 				ext_dict[chrom][j] += win*num
                # 			continue
                ##		if start//win == genome.get(chrom)//win: # since taglen at least 26 bases, impossible for win >= 25 bp
                ##			continue
                if (start // win) != (end // win):
                    if (start // win < genome.get(chrom) // win) and (
                        end // win >= genome.get(chrom) // win
                    ):
                        ext_dict[chrom][start // win] += (win - (start % win)) * num
                        for k in range(start // win + 1, genome.get(chrom) // win):
                            ext_dict[chrom][k] += win * num
                        continue
                    if (start // win < genome.get(chrom) // win) and (
                        end // win < genome.get(chrom) // win
                    ):
                        ext_dict[chrom][start // win] += (win - (start % win)) * num
                        ext_dict[chrom][end // win] += (end % win) * num
                        for i in range(1, (end // win - start // win)):
                            ext_dict[chrom][start // win + i] += win * num
                else:
                    ext_dict[chrom][start // win] += (end - start) * num

    # print('time elapsed:', time() - currtime)

    # print('Printing WIG file ...' )

    currtime = time()

    result = open(name + "2.wig", "w")
    result.write("browser hide all\n")
    result.write("browser pack refGene\n")
    result.write(
        'track type=wiggle_0 visibility=full name="'
        + name
        + '" autoScale=on color=0,0,0 windowingFunction=maximum\n'
    )
    for i in chroms:
        result.write("variableStep chrom=" + str(i) + " span=25\n")
        for k in range(len(ext_dict[i])):
            if float(ext_dict[i][k]) / win >= 1:
                # 		if (ext_dict[i][k] > 0) and ((float(ext_dict[i][k])/win) >= 1):
                result.write("%d\t" % (k * win + 1))
                result.write(
                    "%.2f\n" % (float(ext_dict[i][k]) / win)
                )  # forward+reverse

    result.close()

    print("time elapsed:", time() - currtime)
    print

    print("Peak finding ...")

    import scipy
    from scipy.optimize import leastsq  # Levenberg-Marquandt: leastsq
    import warnings

    # Returns (height, center_x, width_x), the gaussian parameters of distribution found by a fit
    def fitgaussian(array):
        height = array[:, 1].max()
        mean = sum(array[:, 0] * array[:, 1]) / sum(array[:, 1])
        width = np.sqrt(
            abs(sum((array[:, 0] - mean) ** 2 * array[:, 1]) / sum(array[:, 1]))
        )
        p0 = scipy.c_[height, mean, width]
        # >>> p0 #array([[  5.02000000e+02,   1.17327342e+03,   6.43156134e-01]])
        warnings.simplefilter("ignore", Warning)  # scipy 0.10
        p1, success = leastsq(
            errfunc, p0.copy()[0], args=(array[:, 0], array[:, 1]), maxfev=0
        )
        warnings.simplefilter("default", Warning)  # scipy 0.10
        # 	p1, success = leastsq(errfunc, p0.copy()[0], args=(array[:,0],array[:,1]), maxfev=2000, warning=False) # scipy 0.7
        # 	>>> p1 #array([  4.74090817e+02,   1.17323869e+03,   2.02731712e-01])
        # 	print height, p1[0], '\n', 'p1 =', p1
        if (
            p1[0] < 2 * height
            and p1[1] > array[:, 0].min()
            and p1[1] < array[:, 0].max()
        ):
            return p1
        return np.array([0, 0, 0])

    def diff(x):
        d = []
        for i in range(len(x) - 1):
            d.append(x[i + 1] - x[i])
        return d

    def find_peaks(peak_region, threshold):
        # extract peaks from enriched regions
        # 	c = [i[0] for i in peak_region] # chr
        p = [i[1] for i in peak_region]  # pos
        s = [i[2] for i in peak_region]  # num
        ds = diff(s)
        l = len(ds)
        maxima = []
        # find all local maxima
        for i in range(l - 1):
            # the sign of the slope changes from + to -, we have a max or inflect
            if ds[i] > 0 and ds[i + 1] <= 0:
                maxima.append([s[i + 1], i + 1])
        # 			maxima.append((s[i+1], i+1))
        # 			maxima.append([c[i+1], p[i+1], s[i+1]])
        # 	print 'maxima =', maxima
        # if enriched region is too small for max finding just return blank (first)
        if maxima == []:
            # 		tiny_peak = peak_region[0][1:]
            # 		tiny_peak.reverse()
            # 		print 'tiny_peak'
            # 		return [tiny_peak]
            return []
        # 	print 'maxima =', maxima
        # sort maxima in descending order
        # 	maxima = sorted(maxima, reverse=True)
        # 	print 'maxima =', maxima
        # 	# the first peak is defined as the highest maximum
        # 	highest_max = maxima[0]
        # 	first_peak = (highest_max[0], p[highest_max[1]])
        # 	# position of first peak within enriched region
        # 	max_local_pos = highest_max[1]
        t = maxima[0]  # t[1]: first max pos
        # 	m = [[t[0], p[t[1]]]]
        n = [peak_region[0][1]]
        for i in maxima[1:]:
            # 		second_max_score = i[0]
            # 		second_max_pos = i[1]
            # find the deepest dip between the two peaks
            dip = min(s[t[1] : i[1]])
            # 		print 'dip =', dip
            for j in range(t[1], i[1] + 1):
                if s[j] == dip:
                    n.append(p[j])
                    break
            # return 2 peaks if the distance between them are >= 200
            # return 2 peaks if dip is deep enough
            # 		print 'j = ', j '\n', 'p[j] = ', p[j], '\n', 't[1] = ', t[1], '\n', 'p[t[1]] = ', p[t[1]], '\n', 'p[i[1]] = ', p[i[1]]
            if dip / min(t[0], i[0]) < threshold:
                # 			print min(t[0], i[0])
                t = i
            # 			m.append([i[0], p[i[1]]])
            else:
                if (i[1] - t[1]) * 25 >= 200:
                    t = i
                else:
                    n.pop()
                    if t[0] <= i[0]:
                        t = i
        # 					print 't =', t, '\n', 'm =', m
        # 					m[-1] = [i[0], p[i[1]]]
        # 	print 'm =', m
        n.append(peak_region[-1][1])
        # 	print 'n =', n
        for i in range(len(n) - 1):
            # 		print 'n[i] =', n[i]
            if n[i + 1] - n[i] <= 50:
                n[i] = n[i + 1]
        # 	print 'n =', n
        n = list(set(n))
        n.sort()
        # 	print 'n =', n
        return n

    # define a gaussian fitting function where
    # p[0] = amplitude/height, p[1] = mean, p[2] = sigma
    # fitfunc = lambda p, x: p[0]*scipy.exp(-(x-p[1])**2/(2.0*p[2]**2))
    def fitfunc(p, x):
        if p[2] != 0:
            return p[0] * np.exp(-((x - p[1]) ** 2) / (2.0 * p[2] ** 2))
        else:
            return 0

    errfunc = lambda p, x, y: fitfunc(p, x) - y

    print("Importing Alignment File ...")

    currtime = time()

    for chrom in chroms:
        for k in range(len(ext_dict[chrom])):
            if float(ext_dict[chrom][k]) / win >= 3:
                start = int(k * win + 1)
                num = float(ext_dict[chrom][k]) / win  # forward+reverse
                # 			num = ('%.2f' %(float(ext_dict[chrom][k])/win)) # forward+reverse
                s[chrom].append([chrom, start, num])

    print("time elapsed:", time() - currtime)
    print("Finding Peaks ...")  # on %s ...' % s[0][0]

    def test(chrom):  # chr-by-chr peak finding
        d = s[chrom]
        if d == []:
            return []
        peaklist[chrom] = []
        threshold = 0.2
        count = d[0][1]
        peak_region = d[0:1]
        # 	print 'count =', count, '\n', 'peak_region =', peak_region
        for i in range(len(d) - 1):
            if d[i + 1][1] - d[i][1] > 25:
                # 			print '25+'
                if d[i][1] - count >= 200:
                    # 				print 'peak finding', '\n', 'peak_region =', peak_region, '\n', 'len(peak_region) =', len(peak_region)
                    peaks = find_peaks(peak_region, threshold)
                    # 				print 'peaks =', peaks
                    u = []
                    for q in range(len(peaks)):
                        # 					print 'q =', q
                        for k in range(len(d) - 1):
                            # 						if d[k][0] == d[i][0] and d[k][1] == peaks[q]:
                            if d[k][1] == peaks[q]:
                                # 							print 'k =', k, '\n', 'd[k] =', d[k]
                                u.append(k)
                                break
                    # 				print 'u =', u
                    for x in range(len(peaks) - 1):
                        # 					print 'x =', x, '\n', 'u[x] =', u[x]
                        w = np.array([[y[1], y[2]] for y in d[u[x] : u[x + 1]]])
                        # 					print w, '\n', 'len(w) =', len(w)
                        p1 = fitgaussian(w)
                        # 					print 'p1 =', p1
                        corrfit = fitfunc(p1, w[:, 0])
                        b = 0.0
                        for e in range(13):
                            # 						print 'e =', e
                            # 						print int(p1[1]) + (e - 6) * 25
                            # 						print fitfunc(p1, int(p1[1]) + (e - 6) * 25)
                            b += fitfunc(p1, int(p1[1]) + (e - 6) * 25)
                        # 					print 'b =', b
                        peak = [d[i + 1][0], int(p1[1]) + 10, b / 8]
                        # 					print 'peak =', peak
                        peaklist[chrom].append(peak)
                # 					print 'peaklist =', peaklist[chrom]
                # 			print 'new assignment'
                peak_region = d[i + 1 : i + 2]
                count = d[i + 1][1]
            # 			print 'peak_region =', peak_region, '\n', 'count =', count
            else:
                # 			print 'else'
                peak_region.append(d[i + 1])

    # 	return peaklist[chrom]

    currtime = time()

    for i in chroms:
        print(i)
        test(i)

    print("time elapsed:", time() - currtime)

    print("Printing Peaks List ...")  # on %s ... \n' % s[0][0]

    currtime = time()

    result = open(name + "_Peaks2.wig", "w")
    result.write("browser hide all\n")
    result.write("browser pack refGene\n")
    result.write(
        'track type=wiggle_0 visibility=full name="'
        + name
        + ' Peaks" autoScale=on color=0,0,0 windowingFunction=maximum alwaysZero=on\n'
    )
    for i in chroms:
        result.write("variableStep chrom=" + str(i) + " span=10\n")
        if len(peaklist[i]) != 0:
            peaklist[i].sort()
            for item in peaklist[i]:
                if item[1] > 10:
                    # 				result.write('%s\t' % item[0])
                    result.write("%d\t" % item[1])
                    result.write("%.2f\n" % item[2])

    result.close()

    print("time elapsed:", time() - currtime)


def chunker(iterable, chunk_size):
    """
    Yield successive chunk_size chunks from iterable.
    """
    it = iter(iterable)
    while True:
        chunk = []
        for _ in range(chunk_size):
            try:
                chunk.append(next(it))
            except StopIteration:
                break
        if not chunk:
            return
        yield chunk


def process_file(file_name, chromosomes):
    with open(file_name, "r") as f:
        lines = f.readlines()

    # Create a dictionary to hold each chromosome's lines
    chr_dict = {chr_name: [] for chr_name in chromosomes}
    for line in lines:
        chr_name = line.rstrip().split("\t")[0]
        if chr_name in chr_dict:
            chr_dict[chr_name].append(line)

    # Process each chromosome's lines
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=len(chromosomes)
    ) as executor:
        futures = {
            executor.submit(peak_seq, chr_lines, chr_name): chr_name
            for chr_name, chr_lines in chr_dict.items()
            if chr_lines
        }

        for future in concurrent.futures.as_completed(futures):
            chr_name = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                print(f"{chr_name} generated an exception: {exc}")
            else:
                print(f"{chr_name} returned results: {result}")


if __name__ == "__main__":
    file_name = "/media/raid/RAWDATA/SUP_B15/Data/MAPPED/BED/tab-dex-minus-nodex-GR-ChIP-Seq-NegtoZero.bg"
    chromosomes = [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrM",
        "chrX",
        "chrY",
    ]
    process_file(file_name, chromosomes)
