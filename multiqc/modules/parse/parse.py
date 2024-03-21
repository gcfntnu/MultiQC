#!/usr/bin/env python

""" MultiQC module to parse Parse Biosciences summary csv """

from __future__ import print_function
from collections import OrderedDict
import logging
import pandas as pd
from multiqc import config
from multiqc.plots import bargraph, linegraph, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Parse Biosciences module class
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Parse Bioscences', anchor='parse',
            href='https://www.parsebiosciences.com',
            info="Parse Biosciences workflow for single cell sequencing analysis.")

        # Find and load any parse reports
        self.parse_data = dict()
        self.parse_all_data = dict()

        for f in self.find_log_files('parse', filehandles=True):
            self.parse_log(f)

        # Filter to strip out ignored sample names
        self.parse_data = self.ignore_samples(self.parse_data)

        if len(self.parse_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.parse_data)))

        # Write parsed report data to a file
        ## Parse whole JSON to save all its content
        self.write_data_file(self.parse_data, 'multiqc_parse')



        # General Stats Table
        self.parse_set_table_headers()
        self.add_section(
            name = "Summary",
            anchor = "parse-summary",
            description = "QC metrics from Parse Biosciences single cell workflow.",
            plot = table.plot(self.parse_data, self.parse_qc_headers, {})
        )


    def parse_log(self, f):
        """ Parse the csv Summary output from parse and save the summary statistics """
        try:
            parsed_csv = pd.read_csv(f['f'], header=0, index_col=0)
        except:
            log.warn("Could not parse csv: '{}'".format(f['fn']))
            return None
        s_name = f['fn'].replace(".agg_samp_ana_symmary.csv", "")
        s_data = parsed_csv.to_dict(orient='index')
        self.parse_data[s_name] = s_data


    def parse_set_table_headers(self):
        """ Take the parsed stats from the Parse Biosciences summary and add it to the
        General Statistics table at the top of the report """

        headers = OrderedDict()
        headers['number_of_reads'] = {
            'title': '{} Reads'.format(
                config.read_count_prefix,
                ),
            'description': 'Total reads ({})'.format(
                config.read_count_desc
                ),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'GnBu',
            'shared_key': 'read_count',
        }
        headers['valid_barcode_fraction'] = {
            'title': '% Valid BC',
            'description': '% Reads With Valid Barcodes',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100,
            'format': '{:,.1f}',
            'suffix': '%',
            'scale': 'RdYlGn-rev'
        }
        headers['sequencing_saturation'] = {
            'title': 'Saturation',
            'description': 'Sequencing Saturation',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100,
            'format': '{:,.1f}',
            'suffix': '%',
            'scale': 'RdYlGn-rev'
        }
        headers['bc1_Q30'] = {
            'title': 'BC1 % > Q30',
            'description': 'Percentage of Q30 Bases in Barcode Read',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100,
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['bc2_Q30'] = {
            'title': 'BC2 % > Q30',
            'description': 'Percentage of Q30 Bases in Barcode Read',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100,
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['bc3_Q30'] = {
            'title': 'BC3 % > Q30',
            'description': 'Percentage of Q30 Bases in Barcode Read',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100,
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['cDNA_Q30'] = {
            'title': 'RNA % > Q30',
            'description': 'Percentage of Q30 Bases in RNA',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100,
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['transcriptome_map_fraction'] = {
            'title': '% Mapped to Transcriptome',
            'description': 'Percentage Reads Mapped to Transcriptome',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100,
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['number_of_cells'] = {
            'title': '{} # Cells'.format(
                config.read_count_prefix,
                ),
            'description': 'Estimated Number of Cells ({})'.format(
                config.read_count_desc
                ),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'GnBu',
            'shared_key': 'read_count',
        }
        headers['mean_reads_per_cell'] = {
            'title': '{} Mean Reads/Cell'.format(
                config.read_count_prefix,
                ),
            'description': 'Mean reads in cells mapped to unique genes ({})'.format(
                config.read_count_desc
                ),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'GnBu',
            'shared_key': 'read_count',
        }
        self.parse_qc_headers = headers
        gen_stat_cols = ['mean_reads_per_cell', 'number_of_cells', 'sequencing_saturation', 'transcriptome_map_fraction']
        gen_stat_headers = {k: headers[k] for k in gen_stat_cols}
        self.general_stats_addcols(self.parse_data, gen_stat_headers)

