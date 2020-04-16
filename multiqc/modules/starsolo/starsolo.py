#!/usr/bin/env python

""" MultiQC module to parse output from cellranger, cellranger-atac, spaceranger """

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
    cellranger module class
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='starsolo', anchor='starsolo',
            href='https://github.com/alexdobin/STAR',
            info="STAR solo single cell gene quantification.")

        # Find and load any cellranger reports
        self.starsolo_data = dict()
        self.starsolo_all_data = dict()

        for f in self.find_log_files('starsolo', filehandles=True):
            self.parse_starsolo_log(f)

        # Filter to strip out ignored sample names
        self.starsolo_data = self.ignore_samples(self.cellranger_data)

        if len(self.starsolo_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.starsolo_data)))

        # Write parsed report data to a file
        ## Parse whole JSON to save all its content
        self.write_data_file(self.starsolo_data, 'multiqc_cellranger')



        # General Stats Table
        self.starsolo_set_table_headers()

        self.add_section(
            name = "STAR solo summary",
            anchor = "starsolo-summary",
            description = "QC metrics from STAR solo single cell gene quantification",
            plot = table.plot(self.starsolo_data, self.starsolo_qc_headers, {})
        )


    def parse_starsolo_log(self, f):
        """ Parse the csv Summary output from STAR solo and save the summary statistics """
        try:
            parsed_csv = pd.read_csv(f['f'], header=None)
        except:
            log.warn("Could not parse STAR solo csv: '{}'".format(f['fn']))
            return None
        s_name = os.path.basename(f['f']).split('_')[0]
        s_data = dict()
        for i, row in parsed_csv.iterrows():
            s_data[row[0]] = row[1]
            self.cellranger_data[s_name] = s_data



    def cellranger_set_table_headers(self):
        """ Take the parsed stats from the cellranger report and add it to the
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
        headers['Reads With Valid Barcodes'] = {
            'title': '% Valid BC',
            'description': '% Reads With Valid Barcodes',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100.0,
            'format': '{:,.1f}',
            'suffix': '%',
            'scale': 'RdYlGn-rev'
        }
        headers['Sequencing Saturation'] = {
            'title': 'Saturation',
            'description': 'Sequencing Saturation',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100.0,
            'format': '{:,.1f}',
            'suffix': '%',
            'scale': 'RdYlGn-rev'
        }
        headers['Q30 Bases in CB+UMI'] = {
            'title': 'BC % > Q30',
            'description': 'Percentage of Q30 Bases in CB+UMI',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100.0,
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['Q30 Bases in RNA read'] = {
            'title': 'RNA % > Q30',
            'description': 'Percentage of Q30 Bases in RNA',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100.0,
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['Reads Mapped to Genome: Unique+Multiple'] = {
            'title': '% Mapped to Genome',
            'description': 'Percentage Reads Mapped to Genome: Unique+Multiple',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100.0,
            'format': '{:,.1f}',
            'scale': 'GnBu',
            'suffix': '%',
        }
        headers['Reads Mapped to Transcriptome: Unique+Multipe Genes'] = {
            'title': '% Mapped to transcriptome',
            'description': 'Percentage Reads Mapped to Transcriptome: Unique+Multiple',
            'min': 0,
            'max': 100,
            'modify': lambda x: x * 100.0,
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
            'description': 'Mean Reads per Cell ({})'.format(
                config.read_count_desc
                ),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'GnBu',
            'shared_key': 'read_count',
        }
        headers['Mean UMI per Cell'] = {
            'title': '{} Mean UMI/Cells'.format(
                config.read_count_prefix,
                ),
            'description': 'Mean UMI per Cell ({})'.format(
                config.read_count_desc
                ),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'GnBu',
            'shared_key': 'read_count',
        }
        headers['Mean Genes per Cell'] = {
            'title': '{} Mean Genes/Cell'.format(
                config.read_count_prefix,
                ),
            'description': 'Mean Genes per Cell ({})'.format(
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
        self.starsolo_qc_headers = headers
        #self.general_stats_addcols(self.starsolo_data, starsolo_genstats_headers)

