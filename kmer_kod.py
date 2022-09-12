# -- coding: utf-8 --

import matplotlib.pylab as plt
fig = plt.figure()
from collections import Counter

#sequence = ""
#ksize = 0

def build_kmers(sequence, ksize):
    kmers = {}
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1

    return kmers


def kmer_freq_top(sequence, ksize):
    kmers = build_kmers(sequence, ksize)
    mostcommon = Counter(kmers).most_common(5)

    return mostcommon


def reverse_freq(freq_list):
    bolu_freq = []

    for freqs in freq_list:
        a = pow(freqs,-1)
        bolu_freq.append(a)

    return bolu_freq

def kmer_graph_1(sequence, ksize):
    kmers = build_kmers(sequence, ksize)
    freqs = []

    for kmer, freq in kmers.items():
        freqs.append(freq)

    reversed_freq = reverse_freq(freqs)

    return reversed_freq, freqs
    #plt.plot(bolu_freq, freqs,'ro')
    #plt.show()

def kmer_graph_2(sequence, ksize):
    kmers = build_kmers(sequence, ksize)
    kmer_list=[]
    freq_list=[]

    for kmer, freq in kmers.items():
        kmer_list.append(kmer)
        freq_list.append(freq)

    count = 0
    kmer_freq = []

    for i in kmer_list:
        count += 1
        kmer_freq.append(count)

    return kmer_freq, freq_list

def kmer_graph_3(sequence, ksize):
    kmers = build_kmers(sequence, ksize)
    mostcommon = Counter(kmers).most_common(5)

    kmer_string = []
    freq = []
    #print("Most Common Kmers: ", *mostcommon, sep="\n")

    for a, b in mostcommon:
        kmer_string.append(a)
        freq.append(b)

    count = 0
    kmer_no = []

    for a in kmer_string:
        count += 1
        kmer_no.append(count)

    return kmer_no, freq

sequence="ACTGTGTAGCTGA"
ksize=3
a=kmer_freq_top(sequence, ksize)
print(a)
for kmer, freq in a:
    txt = "kmer: {} - freq: {}".format(kmer, freq)
    #self.ui.kmer_top5_freq_txt.append("{}".format(txt))
    print("{}".format(txt))
