#!/usr/bin/env python

""" MultiQC module to parse output from cellranger, cellranger-atac, spaceranger """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc import config
from multiqc.plots import bargraph, linegraph, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    cellranger module class
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='cellranger', anchor='cellranger',
            href='https://support.10xgenomics.com/',
            info="10X genomics workflows for single cell sequencing analysis.")

        # Find and load any cellranger reports
        self.cellranger_data = dict()
        self.cellranger_all_data = dict()

        for f in self.find_log_files('cellranger', filehandles=True):
            self.parse_cellranger_log(f)

        # Filter to strip out ignored sample names
        self.cellranger_data = self.ignore_samples(self.cellranger_data)

        if len(self.cellranger_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.cellranger_data)))

        # Write parsed report data to a file
        ## Parse whole JSON to save all its content
        self.write_data_file(self.cellranger_data, 'multiqc_cellranger')



        # General Stats Table
        self.cellranger_set_table_headers()

        self.add_section(
            name = "10X genomics mkfastq",
            anchor = "cellranger-mkfastq-qc",
            description = "QC metrics from 10X genomics mkfastq pipelines (cellranger, spaceranger, cellranger-atac, longranger)",
            plot = table.plot(self.cellranger_data, self.cellranger_qc_headers, {})
        )


    def modify_sample_metrics(sample_metrics):
        for k, v in sample_metrics.items():
            if k in ['number_reads', 'gem_count_estimate']:
                sample_metrics[k] = int(v)
            elif k in ['barcode_exact_match_ratio', 'barcode_q30_base_ratio', 'read1_q30_base_ratio', 'read2_q30_base_ratio', 'bc_on_whitelist']:
                sample_metrics[k] = round(float(v)*100.0,2)
            elif k == 'mean_barcode_qscore':
                sample_metrics[k] = round(float(v),2)
        return sample_metrics


    def parse_cellranger_log(self, f):
        """ Parse the JSON output from cellranger and save the summary statistics """
        try:
            parsed_json = json.load(f['f'])
        except:
            log.warn("Could not parse cellranger JSON: '{}'".format(f['fn']))
            return None

        for k in parsed_json['sample_qc'].keys():
            s_name = k
            self.add_data_source(f, s_name)
            self.cellranger_data[s_name] = parsed_json['sample_qc'][k]['all']

        # Don't delete dicts with subkeys, messes up multi-panel plots


    def cellranger_set_table_headers(self):
        """ Take the parsed stats from the cellranger report and add it to the
        General Statistics table at the top of the report """

        headers = OrderedDict()
        headers['number_reads'] = {
            'title': '{} Reads'.format(config.read_count_prefix),
            'description': 'Total reads before filtering ({})'.format(config.read_count_desc),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'GnBu',
            'shared_key': 'read_count',
        }
        if 'gem_count_estimate' in self.cellranger_data[list(self.cellranger_data.keys())[0]].keys():
            headers['gem_count_estimate'] = {
                'title': '{} GEM estimate'.format(config.read_count_prefix),
                'description': 'GEM count estimate ({})'.format(config.read_count_desc),
                'modify': lambda x: x * config.read_count_multiplier,
                'min': 0,
                'scale': 'GnBu',
                'shared_key': 'read_count',
            }
        headers['barcode_exact_match_ratio'] = {
            'title': 'Ratio perfect BC',
            'description': 'Ratio of perfect match for barcodes',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100.0,
            'format': '{:,.1f}',
            'suffix': '%',
            'scale': 'RdYlGn-rev'
        }
        headers['barcode_q30_base_ratio'] = {
            'title': 'BC % > Q30',
            'description': 'Percentage of barcode reads > Q30',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100.0,
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['mean_barcode_qscore'] =  {
            'title': 'BC qscore',
            'description': 'Mean barcode Q score',
            'format': '{:,.1f}',
            'min': 0,
            'scale': 'GnBu'
        }
        headers['read1_q30_base_ratio'] = {
            'title': 'R1 % > Q30',
            'description': 'Percentage of R1 reads > Q30',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100.0,
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['read2_q30_base_ratio'] = {
            'title': 'R2 % > Q30',
            'description': 'Percentage of R2 reads > Q30',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100.0,
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        self.cellranger_qc_headers = headers
        #self.general_stats_addcols(self.cellranger_data, headers)

