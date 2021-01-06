#!/usr/bin/env python

""" MultiQC module to parse output from RNA-SeQC """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph

# Initialise the logger
log = logging.getLogger(__name__)


SEQLEN_MATCH = {'miRNA': lambda f: re.match('^unitas\\.miR\\.(.*)\\.info$',f['fn']),
                'tRNA': lambda f: f['fn']== 'unitas.tRNA.info',
                'rRNA': lambda f: f['fn'] == 'unitas.rRNA.info',
                'snoRNA': lambda f: f['fn'] == 'unitas.snoRNA.info',
                'no annotation': lambda f: f['fn'] == 'unitas.no-annotation.info',
                'protein_coding': lambda f: f['fn'] == 'unitas.protein_coding.info',
                'lincRNA': lambda f: f['fn'] == 'unitas.lincRNA.info',
                'antisense': lambda f: f['fn'] == 'unitas.antisense.info',
                'processed_transcript': lambda f: f['fn'] == 'unitas.processed_transcript.info',
                'bidirectional_promoter_lncRNA': lambda f: f['fn'] == 'unitas.bidirectional_promoter_lncRNA.info'
}

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(name='unitas', anchor='unitas',
        href='https://www.smallrnagroup.uni-mainz.de/software.html',
        info=" is a convenient tool for efficient annotation of small non-coding RNA sequence datasets produced by Next Generation Sequencing.")
        self.biotypes = ['miRNA', 'tRNA', 'rRNA', 'protein_coding', 'lincRNA', 'snoRNA', 'no annotation']
        self.metrics = dict()
        self.annotations = dict()
        self.detected = dict()
        self.seqlen = dict()
        
        # parse annotation information.
        for f in self.find_log_files('unitas/annotation'):
            self.parse_annotation(f)
        self.annotations = self.ignore_samples(self.annotations)
        if len(self.annotations) == 0:
            raise UserWarning

        # parse sequence length distribution from .info files
        for f in self.find_log_files('unitas/seqlen'):
            for biotype in self.biotypes:
                match = SEQLEN_MATCH.get(biotype, lambda x: False)
                if match(f): 
                    self.parse_seqlen(f, biotype)
        self.seqlen = self.ignore_samples(self.seqlen)
        if len(self.seqlen) == 0:
            raise UserWarning
        # parse the simplified mirna counts
        for f in self.find_log_files('unitas/mirna'):
            self.parse_mirna_simplified(f)
        self.metrics = self.ignore_samples(self.metrics)
        if len(self.metrics) == 0:
            raise UserWarning
        
        self.unitas_general_stats()
        self.annotation_plot()
        self.seqlen_lineplot()

    def parse_mirna_simplified(self, f, min_detect=3):
        """parse the simplified mirna counts.

        The requirement for a `present mirna`:  count >= min_detect
        """
        s_name = os.path.basename(f['root'])
        present = 0
        for i, line in enumerate(f['f'].splitlines()):
            if i > 1:
                name, val = line.split()
                if float(val) >= min_detect:
                    present += 1
        if not s_name in self.metrics:
            self.metrics[s_name] = {}
        self.metrics[s_name]['miRNA Detected'] = present

    def parse_seqlen(self, f, key):
        """parse the unitas metrics info of miRNA"""
        s_name = os.path.basename(os.path.dirname(f['root']))
        if not key in self.seqlen:
            self.seqlen[key] = {}
        data = {int(i):0 for i in range(1,51)}
        for line in f['f'].splitlines()[2:]:
            if not line:
                break
            x, y = line.split('\t')
            data[int(x)] = float(y or 0)
        self.seqlen[key][s_name] = data

    def parse_annotation(self, f):
        """parse the unitas.annotation_summary.txt file from unitas
        """
        s_name = os.path.basename(f['root'])
        self.annotations[s_name] = dict()
        total = 0.0
        for line in f['f'].splitlines():
            line = line.strip()
            els = line.split('\t')
            if els[0].startswith(' '):
                continue
            else:
                val = float(line.split()[-1].strip())
                if els[0] in self.biotypes:
                    self.annotations[s_name][els[0]] = val
                else:
                    self.annotations[s_name]['other'] = val
                total += val
                    
        self.biotypes.append('other')
        if not s_name in self.metrics:
            self.metrics[s_name] = {}
        for s_name, anno in self.annotations.items():
            self.metrics[s_name]['miRNA Percentage'] = anno['miRNA'] / total
        
    def unitas_general_stats(self):
        """add mirna alignment rate and number of detected mirna to the general stats table."""
        headers = OrderedDict()
        headers['miRNA Percentage'] = {
            'title': '% miRNA',
            'description': 'Percentage of filtered and trimmed sequences annotated to miRNA type',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn',
            'modify': lambda x: float(x) * 100.0
        }
        headers['miRNA Detected'] = {
            'title': '# Genes',
            'description': 'Number of miRNA genes detected with at least 3 reads.',
            'min': 0,
            'scale': 'Bu',
            'format': '{:,.0f}'
        }
      
        self.general_stats_addcols(self.metrics, headers)

    def annotation_plot(self):
        """ Plot a bargraph showing the annotation type distributuion. """
        keys = OrderedDict()
        for k in self.biotypes:
            keys[k] = { 'name': k}
        pconfig = {
            'id': 'unitas_annotation',
            'title': 'Unitas: Gene annotations',
            'ylab': 'XX',
            'cpswitch': True,
            'tt_decimals': 1,
            'cpswitch_c_active': True
        }
        self.add_section (
            name = 'Annotations',
            anchor = 'annotations',
            helptext = 'All of the above rates are per mapped read. ',
            plot = bargraph.plot(self.annotations, keys, pconfig)
        )

    def seqlen_lineplot (self):
        """ Make HTML for sequence length line plots """
        biotypes = [i for i in self.biotypes if i != 'other']
        data = [self.seqlen[k] for k in biotypes if k in self.seqlen]
        data_labels = []
        for biotype in biotypes:
            data_labels.append({'name': biotype, 'ylab': 'Number of {} Reads'.format(biotype)})
        pconfig = {
            'id': 'unitas_seqlen',
            'title': 'Sequence length distribution',
            'ylab': '# Reads',
            'xlab': 'Nucleotide position',
            'xmin': 0,
            'xmax': 50,
            'tt_label': "<strong>{point.x}% from 5'</strong>: {point.y:.2f}",
            'data_labels': data_labels
        }
        self.add_section (
            name = 'Sequence Length Distribution',
            anchor = 'unitas_seqlen',
            helptext = 'aa',
            plot = linegraph.plot(data, pconfig)
        )
