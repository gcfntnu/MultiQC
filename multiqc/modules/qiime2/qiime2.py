#!/usr/bin/env python

""" MultiQC module to parse output from qiime2 """

from __future__ import print_function
from collections import OrderedDict, defaultdict
import logging
import json
import re
import pandas as pd

from multiqc import config
from multiqc.plots import bargraph, scatter
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

default_colors = ['#7cb5ec', '#434348', '#90ed7d', '#f7a35c', '#8085e9',
        '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1']


class MultiqcModule(BaseMultiqcModule):
    """
    qiime2 module class
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='QIIME2', anchor='qiime2',
            href='https://docs.qiime2.org/',
            info="Microbiome workflows for variable region sequencing analysis.")
        
        self.db = 'custom'
        self.qiime2_dada2 = defaultdict(dict)
        
        # Dada2 filtering stats
        self.qiime2_dada2['Total'] = defaultdict(dict)
        for f in self.find_log_files('qiime2/dada2', filehandles=True):
            parsed_data = self.parse_qiime2_dada2_export(f)
            if parsed_data is not None:
                for region, sample_data in parsed_data.items():
                    sample_data = self.ignore_samples(sample_data)
                    for s_name, vals in sample_data.items():
                        self.add_data_source(f, s_name)
                        self.qiime2_dada2[region][s_name] = vals
                        if s_name in self.qiime2_dada2['Total']:
                            for k, v in vals.items():
                                self.qiime2_dada2['Total'][s_name][k] += v
                        else:
                            self.qiime2_dada2['Total'][s_name] = vals.copy()
        if len(self.qiime2_dada2.keys()) == 0:
            raise UserWarning
                
        log.info("Found {} dada2 regions".format(len(self.qiime2_dada2)))

        self.write_data_file(self.qiime2_dada2, 'multiqc_qiime2_dada2')

        self.add_section(
            name = "Dada2 Filtered Reads",
            anchor = "qiime2-dada2-stats",
            description = "Filtering statsitics after dada2 denoising. [DADA2](https://benjjneb.github.io/dada2/index.html): Fast and accurate sample inference from amplicon data with single-nucleotide resolution.",
            helptext = """
            * Passed Filter: Number of reads going into taxonomic classification.  
            * Low Quality: Number of reads filtered out due to too many basepairs with low sequencing quality.  
            * No Overlap: Number of reads filtered out because the paired end reads do not overlap (min_overlap < 12).    
            * Failed denoising: Number of reads mot identified as a true sequence variant.  
            * Chimeric: Number of reads identified as with chimeric origin.
            """,
            plot = self.dada2_filtered_reads_plot()
        )


        # qiime2 taxonomy bargraphs
        self.qiime2_taxa = defaultdict(dict)
        for f in self.find_log_files('qiime2/taxonomy', filehandles=True):
            parsed_data = self.parse_qiime2_taxonomy_export(f)
            if parsed_data is not None:
                for taxa, sample_data in parsed_data.items():
                    sample_data = self.ignore_samples(sample_data)
                    for s_name, vals in sample_data.items():
                        self.add_data_source(f, s_name)
                        self.qiime2_taxa[taxa][s_name] = vals

        if len(self.qiime2_taxa) == 0:
            raise UserWarning
        log.info("Found {} taxonomy reports".format(len(self.qiime2_taxa.keys())))

        self.write_data_file(self.qiime2_taxa, 'multiqc_qiime2_taxonomy')

        self.add_section(
            name = "Taxonomic Classification",
            anchor = "qiime2-taxonomy",
            description = "Predicted taxonomy for qiime2 classification using {} reference database".format(self.db.upper()),
            helptext = "The taxonomy is classified with a naive bayes classifier trained on region specific sequences of the above reference database and the result is joined. <br/>It is important to note that classifications at the lowest level (species) must be considered with care as the variable ribosomal regions may not carry enough information for to classify correctly. ",
            plot = self.taxonomy_plot()
        )

        # qiime2 robust pca
        self.qiime2_rpca = dict()
        for f in self.find_log_files('qiime2/rpca', filehandles=True):
            parsed_data = self.parse_qiime2_rpca(f)
            if parsed_data is not None:
                sample_data = self.ignore_samples(parsed_data)
                for s_name, vals in sample_data.items():
                    self.add_data_source(f, s_name)
                    self.qiime2_rpca[s_name] = vals

        if len(self.qiime2_rpca) == 0:
            raise UserWarning
        log.info("Found {} ordination reports".format(len(self.qiime2_rpca)))

        self.write_data_file(self.qiime2_rpca, 'multiqc_qiime2_taxonomy')

        self.add_section(
            name = "Deicode",
            anchor = "qiime2-deicode",
            description = "DEICODE is a form of Aitchison Distance that is robust to high levels of sparsity. DEICODE utilizes a natural solution to the zero problem formulated in recommendation systems called matrix completion. A simple way to interpret the method is, as a robust compositional PCA (via SVD) where zero values do not influence the resulting ordination.",
            helptext = "The figure shows all samples (each point is a sample) in the two first principal components. By hovering over a point you may see the sample annotation.</br>Points close to each other will have similar bacterial content and point on opposite side of origin will have opposite bacterial content. PCA is usually used to identify clusters of samples and outlying samples.",
            plot = self.rpca_plot()
            )
        
        
        
        # General Stats Table
        self.qiime2_set_table_headers()


    def parse_qiime2_dada2_export(self, f):
        """ Parse the dada2 stats export """
        region_data = defaultdict(dict)
        count_data = defaultdict(dict)
        try:
            txt = f['f'].read().splitlines()
            header = [i.strip() for i in txt.pop(0).split('\t')[1:]]
            for row in txt:
                if row.startswith('#'):
                    continue
                els = row.split('\t')
                s_name = els.pop(0)
                els = map(float, els)
                d = s_name.split('_')
                s_name = '_'.join(d[:-1])
                region = d[-1]
                region_data[region][s_name] = {i:j for i, j in zip(header, els)}
        except:
            log.warn("Could not parse qiime2 metadata: '{}'".format(f['fn']))
            return None
        for region, sample_data in region_data.items():
            for s_name, vals in sample_data.items():
                keep = ['input', 'filtered', 'denoised', 'merged', 'non-chimeric']
                cumulative_counts = [vals[k] for k in keep]
                counts = {}
                for i in range(1, len(cumulative_counts)):
                    counts[keep[i]] = cumulative_counts[i-1] - cumulative_counts[i]
                counts['passed'] = cumulative_counts[i]
                count_data[region][s_name] = counts

        return count_data

    def _extract_scores(self, txt):
        end_row = -1
        rows = {}
        name = None
        xmax, xmin, ymax, ymin = 0,0,0,0
        for i, line in enumerate(txt):
            if name is not None and line == '':
                return rows
            if line.startswith('Site'):
                name, num_rows, num_cols = line.split('\t')
                end_row = i + int(num_rows)
                continue
            if name:
                els = line.split('\t')
                s_name = els.pop(0)
                x, y = [float(i) for i in els[:2]]
                rows[s_name] = {'x': x, 'y': y}
        return None

    def _extract_sample_info(self, f):
        f.seek(0)
        while not f.readline().startswith("Samples:"):
            continue
        return pd.read_csv(f, header=0, sep='\t')
        
    def parse_qiime2_rpca(self, f):
        scores = defaultdict(dict)
        try:
            txt = f['f'].read().splitlines()
            rows = self._extract_scores(txt)
            samples = self._extract_sample_info(f['f'])
            if 'Sample_Group' in samples.columns:
                groups = set(samples['Sample_Group'])
                for i,g in enumerate(groups):
                    g_df = samples[samples["Sample_Group"] == g]
                    for s_id in g_df['Sample_ID']:
                        if s_id in rows.keys():
                            rows[s_id]['color'] = default_colors[i]
                            rows[s_id]['name'] = g
            return rows
        except Exception as e:
            log.warn(e)
            log.warn("Could not parse qiime2 robust pca metadata: '{}'".format(f['fn'])) 
            return None
       
    def parse_qiime2_taxonomy_export(self, f):
        """ Parse the taxonomy bar export """
        count_data = defaultdict(dict)
        taxa_names = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        try:
            txt = f['f'].read().splitlines()
            m = re.match('level-(\d+).csv', f['fn'])
            if m:
                level, = m.groups()
                taxa = taxa_names[int(level)-1]
            else:
                log.warn("Failed to identify correct region and level from filename.")
            if 'D_0__Bacteria' in txt[0]:
                self.db = 'silva'
            elif 'k__Bacteria' in txt[0]:
                self.db = 'greengenes'
            elif 'k__Fungi' in txt[0]:
                self.db = 'unite'
            
            header = [i.strip() for i in txt.pop(0).split(',')]
            if self.db == 'silva':
                db_patt = 'D_0__'
                sep = ';'
            else:
                db_patt = 'k__'
                sep = ';'
            labels = []
            index = []
            for i, name in enumerate(header):
                if db_patt in name:
                    name = name.split(sep)[-1].split('__')[-1].strip()
                    labels.append(name)
                    index.append(i)
            for row in txt:
                if row.startswith('#'):
                    continue
                els = row.split(',')
                s_name = els.pop(0)
                els = [els[i-1] for i in index] 
                els = map(float, els)
                if not taxa == 'Kingdom':
                    count_data[taxa][s_name] = {i:j for i, j in zip(labels, els)}
            return {k: count_data[k] for k in taxa_names[1:]}
        except:
            log.warn("Could not parse qiime2 taxonomy metadata: '{}'".format(f['fn']))

    
    def qiime2_set_table_headers(self):
        """Add dada2 Total data to general stats table.
        """
        headers = OrderedDict()
        headers['passed'] = {'title': 'Dada2 PF PE reads',
                             'description': 'Dada2: number of paired end reads passed filter ({})'.format(config.read_count_desc),
                             'min': 0,
                             'modify': lambda x: x * config.read_count_multiplier,
                             'scale': 'GnBu',
                             'shared_key': 'read_count'}
        self.general_stats_addcols(self.qiime2_dada2['Total'], headers)

    
    def dada2_filtered_reads_plot(self):
        keys = OrderedDict()
        keys['passed'] = { 'name': 'Passed Filter' }
        keys['filtered'] = { 'name': 'Low Quality' }
        keys['merged'] =  { 'name': 'No Overlap' }
        keys['denoised'] = { 'name': 'Failed denoising' }
        keys['non-chimeric'] = { 'name': 'Chimeric' }

        regions = list(self.qiime2_dada2.keys())

        if 'Total' in regions:
            regions.remove('Total')
            regions = sorted(regions)
            if len(regions) > 1:
                regions.insert(0, 'Total')
        else:
            regions = sorted(regions)

        #log.info("Region labels: {} ".format(str(regions))) 
        data = [self.qiime2_dada2[r] for r in regions]
        pconfig = {'id': 'dada2-stats-plot',
                   'title': 'Dada2: Sequence Stats',
                   'xlab': '# Reads',
                   'data_labels': regions}

        return bargraph.plot(data, [keys] * len(data), pconfig)

    def rpca_plot(self):
        """
        """
        pconfig = {'id': 'dada2-rpca-plot',
                   'title': 'Sample Ordination',
                   'xlab': 'PC1',
                   'ylab': 'PC2',
                   'tt_label' :  'PC1: {point.x:.2f}<br/>PC2: {point.y:.2f}'
        }
        return scatter.plot(self.qiime2_rpca, pconfig=pconfig)
    
    def taxonomy_plot(self):
        """Generate barplot of taxonomic classification.
        """
        taxa_names = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        data_labels = [n for n in taxa_names if n in self.qiime2_taxa.keys()]
        data = [self.qiime2_taxa[k] for k in data_labels]

        #set keys to get order wrt. total number of reads
        cats = []
        for taxa in data_labels:
            sample_data = self.qiime2_taxa[taxa]
            tot = defaultdict(float)
            keys = OrderedDict()
            for s_name, vals in sample_data.items():
                for tax, count in vals.items():
                    tot[tax] += count
            sorted_tax = [k[0] for k in sorted(tot.items(), key=lambda item: item[1])[::-1]]
            for n in sorted_tax:
                keys[n] = { 'name': n }
            cats.append(keys)
            
        pconfig = {
            'id': 'taxonomy-classification',
            'title': 'Taxonomy',
            'xlab': '# Reads',
            'data_labels': data_labels
        }
        
        return bargraph.plot(data, cats, pconfig=pconfig)
