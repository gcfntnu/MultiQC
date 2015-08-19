#!/usr/bin/env python

""" MultiQC module to parse output from Bismark """

from collections import defaultdict, OrderedDict
import json
import logging
import mmap
import os
import re

import multiqc

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self, report):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "Bismark"
        self.anchor = "bismark"
        self.intro = '<p><a href="http://www.bioinformatics.babraham.ac.uk/projects/bismark/" target="_blank">Bismark</a> \
            is a tool to map bisulfite converted sequence reads and determine cytosine methylation states.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any Bismark reports
        self.bismark_raw_data = defaultdict(lambda:dict())
        for root, dirnames, filenames in os.walk(self.analysis_dir):
            for fn in filenames:
                if fn.endswith('_PE_report.txt') or fn.endswith('_SE_report.txt'):
                    with open (os.path.join(root,fn), "r") as f:
                        r_data = f.read()
                        fn_search = re.search("Bismark report for: (\S+)", r_data)
                        if fn_search:
                            s_name = fn_search.group(1)
                            s_name = s_name.split(".gz",1)[0]
                            s_name = s_name.split(".fastq",1)[0]
                            s_name = s_name.split(".fq",1)[0]
                            s_name = s_name.split("_val_1",1)[0]
                            self.bismark_raw_data[s_name]['alignment'] = r_data
                        else:
                            logging.warn("Found bismark alignment report, but couldn't recognise contents: {}".format(fn))

                if fn.endswith('.deduplication_report.txt'):
                    with open (os.path.join(root,fn), "r") as f:
                        r_data = f.read()
                        fn_search = re.search("Total number of alignments analysed in (\S+)", r_data)
                        if fn_search:
                            s_name = fn_search.group(1)
                            s_name = s_name.split(".gz",1)[0]
                            s_name = s_name.split(".fastq",1)[0]
                            s_name = s_name.split(".fq",1)[0]
                            s_name = s_name.split("_val_1",1)[0]
                            self.bismark_raw_data[s_name]['dedup'] = r_data
                        else:
                            logging.warn("Found bismark deduplication report, but couldn't recognise contents: {}".format(fn))

                if fn.endswith('_splitting_report.txt'):
                    with open (os.path.join(root,fn), "r") as f:
                        r_data = f.read()
                        s_name = r_data.splitlines()[0]
                        s_name = s_name.split(".gz",1)[0]
                        s_name = s_name.split(".fastq",1)[0]
                        s_name = s_name.split(".fq",1)[0]
                        s_name = s_name.split("_val_1",1)[0]
                        self.bismark_raw_data[s_name]['methextract'] = r_data

        if len(self.bismark_raw_data) == 0:
            logging.debug("Could not find any Bismark reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} Bismark reports".format(len(self.bismark_raw_data)))

        # Parse the raw reports
        self.parse_bismark_reports()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.bismark_stats_table(report)

        self.sections = list()

        # Section 1 - Column chart of alignment stats
        self.parse_alignment_chart_data()
        self.sections.append({
            'name': 'Alignment Rates',
            'anchor': 'bismark-alignment',
            'content': self.bismark_alignment_chart()
        })

        # Section 2 - Methylation percentages
        self.parse_methylation_chart_data()
        self.sections.append({
            'name': 'Cytosine Methylation',
            'anchor': 'bismark-methylation',
            'content': self.bismark_methlyation_chart()
        })

    def parse_bismark_reports (self):
        """ Search the three types of Bismark report files for
        numbers needed later in the module. """

        regexes = {
            'alignment': {
                'total_reads': r"^Sequence(?:s| pairs) analysed in total:\s+(\d+)$",
                'aligned_reads': r"^Number of(?: paired-end) alignments with a unique best hit:\s+(\d+)$",
                'no_alignments': r"^Sequence(?:s| pairs) with no alignments under any condition:\s+(\d+)$",
                'ambig_reads': r"^Sequence(?:s| pairs) did not map uniquely:\s+(\d+)$",
                'discarded_reads': r"^Sequence(?:s| pairs) which were discarded because genomic sequence could not be extracted:\s+(\d+)$",
                'aln_total_c': r"^Total number of C's analysed:\s+(\d+)$",
                'aln_meth_cpg': r"^Total methylated C's in CpG context:\s+(\d+)",
                'aln_meth_cph': r"^Total methylated C's in CHG context:\s+(\d+)",
                'aln_meth_chh': r"^Total methylated C's in CHH context:\s+(\d+)",
                'aln_unmeth_cpg': r"^Total unmethylated C's in CpG context:\s+(\d+)",
                'aln_unmeth_cph': r"^Total unmethylated C's in CHG context:\s+(\d+)",
                'aln_unmeth_chh': r"^Total unmethylated C's in CHH context:\s+(\d+)",
                'aln_percent_cpg_meth': r"^C methylated in CpG context:\s+([\d\.]+)%",
                'aln_percent_chg_meth': r"^C methylated in CHG context:\s+([\d\.]+)%",
                'aln_percent_chh_meth': r"^C methylated in CHH context:\s+([\d\.]+)%"
            },
            'dedup': {
                # 'aligned_reads' overwrites previous, but I trust this more
                # Leave the number from the alignment report in case deduplication is not run
                'aligned_reads': r"^Total number of alignments analysed in .+:\s+(\d+)$",
                'dup_reads': r"^Total number duplicated alignments removed:\s+(\d+)",
                'dup_reads_percent': r"^Total number duplicated alignments removed:\s+\d+\s+\(([\d\.]+)%\)",
                'dedup_reads': r"^Total count of deduplicated leftover sequences:\s+(\d+)",
                'dedup_reads_percent': r"^Total count of deduplicated leftover sequences:\s+\d+\s+\(([\d\.]+)% of total\)"
            },
            'methextract': {
                # These calls are typically done after deduplication
                'me_total_c': r"^Total number of C's analysed:\s+(\d+)$",
                'me_meth_cpg': r"^Total methylated C's in CpG context:\s+(\d+)",
                'me_meth_cph': r"^Total methylated C's in CHG context:\s+(\d+)",
                'me_meth_chh': r"^Total methylated C's in CHH context:\s+(\d+)",
                'me_unmeth_cpg': r"^Total C to T conversions in CpG context:\s+(\d+)",
                'me_unmeth_cph': r"^Total C to T conversions in CHG context:\s+(\d+)",
                'me_unmeth_chh': r"^Total C to T conversions in CHH context:\s+(\d+)",
                'me_percent_cpg_meth': r"^C methylated in CpG context:\s+([\d\.]+)%",
                'me_percent_chg_meth': r"^C methylated in CHG context:\s+([\d\.]+)%",
                'me_percent_chh_meth': r"^C methylated in CHH context:\s+([\d\.]+)%"
            }
        }
        for sn, data in self.bismark_raw_data.iteritems():
            for report_type in regexes.keys():
                for k, r in regexes[report_type].iteritems():
                    try:
                        r_search = re.search(r, data[report_type], re.MULTILINE)
                        if r_search:
                            self.bismark_raw_data[sn][k] = float(r_search.group(1))
                    except KeyError:
                        pass # Missing report type

    def bismark_stats_table(self, report):
        """ Take the parsed stats from the Bismark reports and add them to the
        basic stats table at the top of the report """

        # Use several try blocks in case one of the report types is missing
        # If exception is triggered, header rows won't be added
        try:
            for sn, data in self.bismark_raw_data.iteritems():
                report['general_stats']['rows'][sn]['percent_cpg_meth'] = '<td class="text-right">{:.1f}%</td>'.format(data['me_percent_cpg_meth'])
                report['general_stats']['rows'][sn]['total_c'] = '<td class="text-right">{:.1f}</td>'.format(data['me_total_c']/1000000)
            report['general_stats']['headers']['percent_cpg_meth'] = '<th class="chroma-col" data-chroma-scale="BrBG" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: % Cytosines methylated in CpG context (meth&nbsp;extraction)">%&nbsp;Meth</span></th>'
            report['general_stats']['headers']['total_c'] = '<th class="chroma-col" data-chroma-scale="Purples" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Total number of C\'s analysed, in millions (meth&nbsp;extraction)">M&nbsp;C\'s</span></th>'
        except KeyError:
            # Use numbers from alignment instead
            try:
                for sn, data in self.bismark_raw_data.iteritems():
                    report['general_stats']['rows'][sn]['percent_cpg_meth'] = '<td class="text-right">{:.1f}%</td>'.format(data['aln_percent_cpg_meth'])
                    report['general_stats']['rows'][sn]['total_c'] = '<td class="text-right">{:.1f}</td>'.format(data['aln_total_c']/1000000)
                report['general_stats']['headers']['percent_cpg_meth'] = '<th class="chroma-col" data-chroma-scale="Greens" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: % Cytosines methylated in CpG context (alignment)">%&nbsp;Meth</span></th>'
                report['general_stats']['headers']['total_c'] = '<th class="chroma-col" data-chroma-scale="Purples" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Total number of C\'s analysed, in millions (alignment)">M&nbsp;C\'s</span></th>'
            except KeyError:
                pass

        try:
            for sn, data in self.bismark_raw_data.iteritems():
                report['general_stats']['rows'][sn]['bismark_dedup_reads_percent'] = '<td class="text-right">{:.1f}%</td>'.format(data['dup_reads_percent'])
                report['general_stats']['rows'][sn]['bismark_dedup_reads'] = '<td class="text-right">{:.1f}</td>'.format(data['dedup_reads']/1000000)
                report['general_stats']['rows'][sn]['bismark_aligned'] = '<td class="text-right">{:.1f}</td>'.format(data['aligned_reads']/1000000)
            report['general_stats']['headers']['bismark_dedup_reads_percent'] = '<th class="chroma-col" data-chroma-scale="RdYlGn-rev" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Percent Duplicated Alignments">%&nbsp;Dups</span></th>'
            report['general_stats']['headers']['bismark_dedup_reads'] = '<th class="chroma-col" data-chroma-scale="Greens" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Deduplicated Alignments (millions)">M&nbsp;Unique</span></th>'
        except KeyError:
            pass

        try:
            for sn, data in self.bismark_raw_data.iteritems():
                report['general_stats']['rows'][sn]['bismark_percent_aligned'] = '<td class="text-right">{:.1f}%</td>'.format((data['aligned_reads']/data['total_reads'])*100)
                report['general_stats']['rows'][sn]['bismark_aligned'] = '<td class="text-right">{:.1f}</td>'.format(data['aligned_reads']/1000000)
            report['general_stats']['headers']['bismark_percent_aligned'] = '<th class="chroma-col" data-chroma-scale="YlGn" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Percent Aligned Sequences">%&nbsp;Aligned</span></th>'
            report['general_stats']['headers']['bismark_aligned'] = '<th class="chroma-col" data-chroma-scale="PuRd" data-chroma-min="0"><span data-toggle="tooltip" title="Bismark: Total Aligned Sequences (millions)">M&nbsp;Aligned</span></th>'
        except KeyError:
            pass


    def parse_alignment_chart_data (self):
        """ Make a data structure suitable for HighCharts for the alignment plot """
        self.bismark_sn_categories = list()
        series = OrderedDict()
        series['No Genomic Sequence'] = list()
        series['Did Not Align'] = list()
        series['Aligned Ambiguously'] = list()
        series['Aligned Uniquely'] = list()
        series['Duplicated Unique Alignments'] = list()
        series['Deduplicated Unique Alignments'] = list()
        for sn in sorted(self.bismark_raw_data.keys()):
            self.bismark_sn_categories.append(sn)
            series['No Genomic Sequence'].append(int(self.bismark_raw_data[sn].get('discarded_reads', 0)))
            series['Did Not Align'].append(int(self.bismark_raw_data[sn].get('no_alignments', 0)))
            series['Aligned Ambiguously'].append(int(self.bismark_raw_data[sn].get('ambig_reads', 0)))
            try:
                series['Duplicated Unique Alignments'].append(int(self.bismark_raw_data[sn]['dup_reads']))
                series['Deduplicated Unique Alignments'].append(int(self.bismark_raw_data[sn]['dedup_reads']))
                series['Aligned Uniquely'].append(0)
            except KeyError:
                series['Aligned Uniquely'].append(int(self.bismark_raw_data[sn].get('aligned_reads', 0)))

        self.bismark_aln_plot_series = list()
        for cat in series:
            if(len(series[cat]) > 0 and max(series[cat]) > 0):
                self.bismark_aln_plot_series.append({
                    'name': cat,
                    'data': series[cat]
                })

    def bismark_alignment_chart (self):
        """ Make the HighCharts HTML to plot the alignment rates """

        return '<div class="btn-group switch_group"> \n\
			<button class="btn btn-default btn-sm active" data-action="set_numbers" data-target="#bismark_alignment_plot">Number of Reads</button> \n\
			<button class="btn btn-default btn-sm" data-action="set_percent" data-target="#bismark_alignment_plot">Percentages</button> \n\
		</div> \n\
        <div id="bismark_alignment_plot" class="fastqc-overlay-plot" style="height:500px;"></div> \n\
        <script type="text/javascript"> \n\
            bismark_alignment_cats = {};\n\
            bismark_alignment_data = {};\n\
            var bismark_alignment_pconfig = {{ \n\
                "colors": ["#f28f43", "#0d233a", "#492970", "#2f7ed8", "#8bbc21"], \n\
                "title": "Bismark Alignment Scores",\n\
                "ylab": "# Reads",\n\
                "ymin": 0,\n\
                "stacking": "normal" \n\
            }}; \n\
            $(function () {{ \
                plot_stacked_bar_graph("#bismark_alignment_plot", bismark_alignment_cats, bismark_alignment_data, bismark_alignment_pconfig); \
            }}); \
        </script>'.format(json.dumps(self.bismark_sn_categories), json.dumps(self.bismark_aln_plot_series));


    def parse_methylation_chart_data (self):
        """ Make a data structure suitable for HighCharts for the methylation plot """
        self.bismark_meth_helptext = "Numbers taken from methylation extraction report."
        self.bismark_meth_snames = list()
        series = OrderedDict()
        series['Unmethylated CpG'] = list()
        series['Methylated CpG'] = list()
        for sn in sorted(self.bismark_raw_data.keys()):
            self.bismark_meth_snames.append(sn)
            try:
                series['Unmethylated CpG'].append(int(self.bismark_raw_data[sn]['me_unmeth_cpg']))
                series['Methylated CpG'].append(int(self.bismark_raw_data[sn]['me_meth_cpg']))
            except KeyError:
                series['Unmethylated CpG'].append(int(self.bismark_raw_data[sn]['aln_unmeth_cpg']))
                series['Methylated CpG'].append(int(self.bismark_raw_data[sn]['aln_meth_cpg']))
                self.bismark_meth_helptext = "Numbers taken from Bismark alignment report"

        self.bismark_meth_plot_series = list()
        for cat in series:
            self.bismark_meth_plot_series.append({
                'name': cat,
                'data': series[cat]
            })

    def bismark_methlyation_chart (self):
        """ Make the HighCharts HTML to plot the methylation calls """

        return '<p class="text-muted">{}<p> \n\
        <div class="btn-group switch_group"> \n\
			<button class="btn btn-default btn-sm" data-action="set_numbers" data-target="#bismark_methylation_plot">Number of Calls</button> \n\
			<button class="btn btn-default btn-sm active" data-action="set_percent" data-target="#bismark_methylation_plot">Percentages</button> \n\
		</div> \n\
        <div id="bismark_methylation_plot" class="fastqc-overlay-plot" style="height:500px;"></div> \n\
        <script type="text/javascript"> \n\
            bismark_methylation_cats = {};\n\
            bismark_methylation_data = {};\n\
            var bismark_methylation_pconfig = {{ \n\
                "colors": ["#0d233a", "#2f7ed8", "#8bbc21", "#1aadce", "#910000", "#492970"], \n\
                "title": "Cytosine CpG Methylation",\n\
                "ylab": "% Calls",\n\
                "ymin": 0,\n\
                "stacking": "percent" \n\
            }}; \n\
            $(function () {{ \
                plot_stacked_bar_graph("#bismark_methylation_plot", bismark_methylation_cats, bismark_methylation_data, bismark_methylation_pconfig); \
            }}); \
        </script>'.format(self.bismark_meth_helptext, json.dumps(self.bismark_meth_snames), json.dumps(self.bismark_meth_plot_series));