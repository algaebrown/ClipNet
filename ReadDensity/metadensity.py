import pyBigWig
import pandas as pd
import sys
sys.path.append('/home/hsher/rbp-maps/maps/')
from density.ReadDensity import ReadDensity
import os
from pybedtools import BedTool
import matplotlib.pyplot as plt
import numpy as np
basedir = '/home/hsher/seqdata/eclip_raw/'
from scipy.stats import entropy


transcript = BedTool('/home/hsher/projects/gencode_transcript.gff3')
gencode_feature = BedTool('/home/hsher/projects/gencode_combine_sorted.gff3')

def make_density(series, basedir = basedir):
    ''' Generate 3 ReadDensity Object from pd.Series from encode_data_id.pickle'''
    all_den = []
    for types in ['0','1','control']: # rep1, rep2, control
        neg = basedir + series['minus_'+types].values[0]
        pos = basedir + series['plus_'+types].values[0]
        bam = basedir + series['bam_'+types].values[0]
        
        density = ReadDensity(pos, neg, bam = bam, name = str(series['RBP']))
        all_den.append(density)
    return all_den[0], all_den[1], all_den[2]
        
class eCLIP:
    def __init__(self, density_rep1 = None, density_rep2 = None, density_ctrl = None):
        self.rep1 = density_rep1 # Read Density Objects
        self.rep2 = density_rep2
        self.ctrl = density_ctrl

    def from_pd_series(self, series):
        ''' from the above encode data dataframe'''
        rep1, rep2, ctrl = make_density(series)
        self.__init__(rep1, rep2, ctrl)
        self.uID = series['uID'].values[0]
        self.name = series['RBP'].values[0]
    def add_peaks(self,
                  idr_path = '/home/hsher/seqdata/eclip_bed/sorted/',
                  indv_path = '/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20180205/hg38/'):
        ''' add peaks to eCLIP object'''
        idr_suffix = '.01v02.IDR.out.0102merged.bed.blacklist_removed.bed.narrowPeak.bed'
        self.idr = BedTool(idr_path + self.uID + idr_suffix)
        
        indv_suffix = '{0}_0{1}.basedon_{0}_0{1}.peaks.l2inputnormnew.bed.compressed.bed.blacklist_removed.bed'
        self.peak1 = BedTool(indv_path + indv_suffix.format(self.uID, 1))
        self.peak2 = BedTool(indv_path + indv_suffix.format(self.uID, 2))
    def find_idr_transcript(self, genome_coord = transcript):
        ''' find positive example for generateing data'''
        self.idr_transcript = genome_coord.intersect(self.idr, s=True)
    def find_negative_example(self, genome_coord = transcript):
        ''' find negative regions as no IDR neither individual peaks'''
        self.no_peak = genome_coord.intersect(self.idr, s=True, v=True).intersect(self.peak1, s = True, v =True).intersect(self.peak2, s = True, v =True)
     
    def build_metagene(self, sample_no = 200):
        ''' Use function `Build_many_metagenes()` to get regions in UTR, intron, exon for each transcripts in self.idr_transcript and self.no_peak; store metagene in self.idr_metagene/self.neg_metagene'''
        self.idr_metagene = Build_many_metagene(self.idr_transcript, sample_no = sample_no)
        self.neg_metagene = Build_many_metagene(self.no_peak, sample_no = sample_no)
    def get_metadensity(self):
        ''' calculate metadensity for each metagene in self.idr_metagene or self.neg_metagene.
        store density in each metagene object'''
        _ = [m.metadensity(self) for m in self.idr_metagene.values()]
        _ = [m.metadensity(self) for m in self.neg_metagene.values()]
    def get_feature_density_array(self, feature, target_len, align, pad_value = np.nan, example = 'positive'):
        d1 = []
        d2 = []
        if example == 'positive':
            metagenes = self.idr_metagene.values()
        else:
            metagenes = self.neg_metagene.values()
        for m in metagenes:
            d1.append(trim_or_pad(m.densities[self.uID]['rep1'][feature], target_len, align, pad_value))
            d2.append(trim_or_pad(m.densities[self.uID]['rep2'][feature], target_len, align, pad_value))
        self.density_array[example,feature,align,'rep1']= np.stack(d1)
        self.density_array[example,feature,align,'rep2']= np.stack(d2)
                
        
    def get_density_array(self, five_utr_len=500, three_utr_len=1000, intron_len = 1500, exon_len = 1000):
        ''' extract metadensity from each metagene, zero pad or trim to get np array '''
        self.density_array = {}
        for feature, l in zip(['five_utr', 'exon', 'intron', 'three_utr'], [five_utr_len, exon_len, intron_len, three_utr_len]):
            for example in ['positive', 'negative']:
                for align in ['left', 'right']:
                    self.get_feature_density_array(feature, l, align, example = example)
    
    
    def RBP_centric_approach(self, series, sample_no = 200):
        ''' create eCLIP object, get peaks, get examples, metagene and metadensity with one function. return eCLIP object'''
        self.from_pd_series(series)
        print('adding peaks')
        self.add_peaks()
        print('finding negative/positive examples')
        self.find_idr_transcript()
        self.find_negative_example()
        print('Building metagene and metadensity')
        self.build_metagene(sample_no = sample_no)
        self.get_metadensity()
        
        

########################## Metagene class #####################################

class Metagene:
    def __init__(self, esnt, chro, start, end, strand):
        self.ensembl_id = esnt
        self.chrom = chro
        self.start = start
        self.stop = end
        self.strand = strand
        self.five_utr = set()
        self.three_utr = set()
        self.intron = set()
        self.exon = set()
        self.densities = {} # key: eCLIP uID, values: metadensity
    def multi_feature_avg(self, feature, eCLIP, align = 'left', max_len = None):
        '''
        average and align (zero padding) multiple intron/exon
        return nanmean for both replicates
        '''
        den1 = []
        den2 = []
        
        # feature length
        if max_len == None:
            max_len = max([f[1]-f[0] for f in feature])
        
        for f in feature:
            minus1, minus2 = subtract_input(eCLIP, self.chrom, f[0], f[1], self.strand)
            
        
            den1.append(trim_or_pad(minus1, max_len, pad_value = np.nan))
            den2.append(trim_or_pad(minus2, max_len, pad_value = np.nan))
        return np.nanmean(np.stack(den1), axis = 0), np.nanmean(np.stack(den2), axis = 0)
        
    def metadensity(self, eCLIP):
        ''' get metadensity from eCLIP object (containing density) 
        store. self.densities[eCLIP.uid] as dictionary {'rep1': exon: np.array}'''
        feature_den1 = []
        feature_den2 = []
        for feature in [self.five_utr, self.exon, self.intron, self.three_utr]:
            if len(feature) == 0: # no such feature
                minus1 = np.empty(1)
                minus2 = np.empty(1)
                
            elif len(feature) == 1:
                feature = list(feature)[0]
                minus1, minus2 = subtract_input(eCLIP, self.chrom, feature[0], feature[1], self.strand)
                minus1 = np.array(minus1)
                minus2 = np.array(minus2)
            else:
                minus1, minus2 = self.multi_feature_avg(feature, eCLIP) 
            feature_den1.append(minus1)
            feature_den2.append(minus2)
        
        # transcript level normalization
        if np.nansum(np.concatenate(feature_den1))!= 0:
            n1 = [f/np.nansum(np.concatenate(feature_den1)) for f in feature_den1]
        else: 
            n1 = feature_den1
        if np.nansum(np.concatenate(feature_den2))!= 0:
            n2 = [f/np.nansum(np.concatenate(feature_den2)) for f in feature_den2]
        else:
            n2 = feature_den2
    
    
        fnames = ['five_utr', 'exon', 'intron', 'three_utr']
        
        self.densities[eCLIP.uID] = {
                                    'rep1':dict(zip(fnames, n1)),
                                    'rep2':dict(zip(fnames, n2))
                                    }
        

    


def Build_many_metagene(key_transcript, gencode_feature = gencode_feature, sample_no = 200):
    ''' Create Metagene object for regions for key transcript '''
    if len(key_transcript)< sample_no:
        n = len(key_transcript)
    else: 
        n = sample_no
    
    # put into transcript
    all_metagene = {}
    for i in range(n):
        enst = key_transcript[i].fields[-1].split(';')[3].replace('transcript_id=', '')
        metagene = Metagene(enst, key_transcript[i].chrom, key_transcript[i].start, key_transcript[i].stop, key_transcript[i].strand)
        all_metagene[enst] = metagene
    
    # write where is intron, exon...
    for i in gencode_feature.filter(lambda x: x.fields[-1].split(';')[3].replace('transcript_id=', '') in all_metagene.keys()):
        enst = i.fields[-1].split(';')[3].replace('transcript_id=', '')
        feature_type = i.fields[2]
        if feature_type == 'transcript':
            all_metagene[enst].intron.update([(i.start, i.stop)])
        if feature_type == 'exon':
            all_metagene[enst].exon.update([(i.start, i.stop)])
        if feature_type == 'five_prime_UTR':
            all_metagene[enst].five_utr.update([(i.start, i.stop)])
        if feature_type == 'three_prime_UTR':
            all_metagene[enst].three_utr.update([(i.start, i.stop)])
    return all_metagene

########################################## calculate meta-density ##########################################
def subtract_input(eclip, chrom, start, stop, strand):
    '''
    return IP density minus input density for 2 biological replicates
    '''
    density1 = np.nan_to_num(eclip.rep1.values(chrom, start, stop, strand),0)
    density2 = np.nan_to_num(eclip.rep2.values(chrom, start, stop, strand),0)
    
    
    
    # get input
    density_ctrl = np.nan_to_num(eclip.ctrl.values(chrom, start, stop, strand),0)
    
    if strand == '-':
        density1 = -density1
        density2 = -density2
        density_ctrl = -density_ctrl
        
    minus1 = np.array(density1) - np.array(density_ctrl)
    minus2 = np.array(density2) - np.array(density_ctrl)
    
    
    # no negative value 
    minus1[minus1 < 0] = 0
    minus2[minus2 < 0] = 0
    
    
    return minus1, minus2

def trim_or_pad(density, target_length, align = 'left', pad_value = 0):
    ''' make density your target length by trimming or padding'''
    if len(density) == target_length:
        return density
    elif len(density) > target_length:
        if align == 'left':
            return density[:target_length]
        if align == 'right':
            return density[-target_length:]
    else:
        discrepency = target_length - len(density)
        if align == 'left':
            density = np.array(list(density) + [pad_value]*discrepency)
            
            
        if align == 'right':
            density = np.array([pad_value]*discrepency +list(density))
        return density


######################################## Probability Distribution #########################################
def bind_strength_discretize(density_array, bins = 8, ymax = 0.03):
    ''' given density array(sample * length of feature), count discretization'''
    pseudocount = [np.histogram(d[~np.isnan(d)], range = (0, ymax), bins = bins)[0] + 1 for d in density_array.T] # pseudocount 
    den_prob = np.stack([p/np.sum(p) for p in pseudocount]).T
    
    return den_prob

def pos_spec_bind_strength(eCLIP, peak_max = 0.03, bins = 20):
    ''' return probability distribution on \"binding stregth/normalized peak height\" for each position'''
    bind_strength = {}
    for k in eCLIP.density_array.keys():
        bind_strength[k] = bind_strength_discretize(eCLIP.density_array[k], bins = bins, ymax = peak_max)
    return bind_strength
######################################### Relative entropy ##########################################


def density_array_entropy(f, fn):
    '''
    relative entropy for each position, f = positive binding density; fn = negative binding examples's density
    '''
    return np.array([entropy(pk = f[1:, pos], qk = fn[1:, pos]) for pos in range(f.shape[1])])

def pos_relative_entropy(eCLIP_prob):
    '''
    return relative entropy throughout whole transcripts
    '''
    rel_ent = {}
    for feat in ['five_utr', 'exon', 'intron', 'three_utr']:
        for align in ['left', 'right']:
            pos = np.mean(np.array([eCLIP_prob['positive', feat,align, r] for r in ['rep1', 'rep2']]), axis = 0)
            neg = np.mean(np.array([eCLIP_prob['negative', feat,align, r] for r in ['rep1', 'rep2']]), axis = 0)
            rel_ent[feat, align] = density_array_entropy(pos, neg)
    return rel_ent
