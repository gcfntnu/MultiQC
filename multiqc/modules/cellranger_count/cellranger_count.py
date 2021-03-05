#!/usr/bin/env python

""" MultiQC module to parse STAR solo Summary.csv """

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
    cellranger_count module class
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='cellranger_count', anchor='cellranger_count',
            href='https://support.10xgenomics.com',
            info="10X genomics workflows for single cell sequencing analysis.")

        # Find and load any cellranger_count reports
        self.cellranger_count_data = dict()
        self.cellranger_count_all_data = dict()

        for f in self.find_log_files('cellranger_count', filehandles=True):
            self.parse_cellranger_count_log(f)

        # Filter to strip out ignored sample names
        self.cellranger_count_data = self.ignore_samples(self.cellranger_count_data)

        if len(self.cellranger_count_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.cellranger_count_data)))

        # Write parsed report data to a file
        ## Parse whole JSON to save all its content
        self.write_data_file(self.cellranger_count_data, 'multiqc_cellranger_count')



        # General Stats Table
        self.cellranger_count_set_table_headers()

        self.add_section(
            name = "cellranger count summary",
            anchor = "cellranger_count-summary",
            description = "QC metrics from cellranger single cell gene quantification.",
            plot = table.plot(self.cellranger_count_data, self.cellranger_count_qc_headers, {})
        )


    def parse_cellranger_count_log(self, f):
        """ Parse the csv Summary output from STAR solo and save the summary statistics """
        try:
            parsed_csv = pd.read_csv(f['f'], header=0, index_col=None)
        except:
            log.warn("Could not parse cellranger count csv: '{}'".format(f['fn']))
            return None
        s_name = f['fn'].replace(".metrics_summary.csv", "")
        s_data = parsed_csv.to_dict(orient='index')[0]
        self.cellranger_count_data[s_name] = s_data



    def cellranger_count_set_table_headers(self):
        """ Take the parsed stats from the cellranger_count summary and add it to the
        General Statistics table at the top of the report """

        headers = OrderedDict()
        headers['Number of Reads'] = {
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
        headers['Valid Barcodes'] = {
            'title': '% Valid BC',
            'description': '% Reads With Valid Barcodes',
            'min': 0,
            'max': 100,
            'modify': lambda x: x.replace("%",""),
            'format': '{:,.1f}',
            'suffix': '%',
            'scale': 'RdYlGn-rev'
        }
        headers['Sequencing Saturation'] = {
            'title': 'Saturation',
            'description': 'Sequencing Saturation',
            'min': 0,
            'max': 100,
            'modify': lambda x: x.replace("%",""),
            'format': '{:,.1f}',
            'suffix': '%',
            'scale': 'RdYlGn-rev'
        }
        headers['Q30 Bases in Barcode'] = {
            'title': 'BC % > Q30',
            'description': 'Percentage of Q30 Bases in Barcode Read',
            'min': 0,
            'max': 100,
            'modify': lambda x: x.replace("%",""),
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['Q30 Bases in RNA Read'] = {
            'title': 'RNA % > Q30',
            'description': 'Percentage of Q30 Bases in RNA',
            'min': 0,
            'max': 100,
            'modify': lambda x: x.replace("%",""),
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['Reads Mapped to Genome'] = {
            'title': '% Mapped to Genome',
            'description': 'Percentage Reads Mapped to Genome: Unique+Multiple',
            'min': 0,
            'max': 100,
            'modify': lambda x: x.replace("%",""),
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['Estimated Number of Cells'] = {
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
        headers['Mean Reads per Cell'] = {
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
        headers['Median UMI Counts per Cell'] = {
            'title': '{} Median UMI/Cells'.format(
                config.read_count_prefix,
                ),
            'description': 'Median UMI per Cell ({})'.format(
                config.read_count_desc
                ),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'GnBu',
            'shared_key': 'read_count',
        }
        headers['Median Genes per Cell'] = {
            'title': '{} Median Genes/Cell'.format(
                config.read_count_prefix,
                ),
            'description': 'Median Genes per Cell ({})'.format(
                config.read_count_desc
                ),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'GnBu',
            'shared_key': 'read_count',
        }
        headers['Total Genes Detected'] = {
            'title': '{} Detected Genes'.format(
                config.read_count_prefix,
                ),
            'description': 'Total Genes Detected ({})'.format(
                config.read_count_desc
                ),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'GnBu',
            'shared_key': 'read_count',
        }
        self.cellranger_count_qc_headers = headers
        gen_stat_cols = ['Mean Reads per Cell', 'Estimated Number of Cells', 'Sequencing Saturation']
        gen_stat_headers = {k: headers[k] for k in gen_stat_cols}
        self.general_stats_addcols(self.cellranger_count_data, gen_stat_headers)

