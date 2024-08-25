#!/usr/bin/env python

__author__ = 'lenk'

import os
import sys
import shutil

rquast_dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

sys.path.insert(0, os.path.join(rquast_dirpath, 'quast_libs'))

import plugins.QuastRNA.rnaquast.quast_libs.fastaparser as fastaparser

from plugins.QuastRNA.rnaquast.general import rqconfig
from plugins.QuastRNA.rnaquast.general import log

from plugins.QuastRNA.rnaquast.general import UtilsGeneral
from plugins.QuastRNA.rnaquast.general import UtilsPipeline
from plugins.QuastRNA.rnaquast.general import UtilsTools
from plugins.QuastRNA.rnaquast.general import UtilsAlignment
from plugins.QuastRNA.rnaquast.general import UtilsAnnotations

from plugins.QuastRNA.rnaquast.objects import SortedExonsAttributes

from plugins.QuastRNA.rnaquast.metrics import TranscriptsMetrics
from plugins.QuastRNA.rnaquast.metrics import GeneDatabaseMetrics
from plugins.QuastRNA.rnaquast.metrics import ReadsCoverage

from plugins.QuastRNA.rnaquast.report import ShortReport
from plugins.QuastRNA.rnaquast.report import SeparatedReport
from plugins.QuastRNA.rnaquast.report import ComparisonReport

logger = log.get_logger(rqconfig.LOGGER_DEFAULT_NAME)

from plugins.QuastRNA.rnaquast.general.rqconfig import PRECISION

import PyIO
import PyPluMA
class QuastRNAPlugin:
  def input(self, inputfile):
        self.parameters = PyIO.readParameters(inputfile)
  def run(self):
        pass
  def output(self, outputfile):
    #program_name = sys.argv[0][:sys.argv[0].rfind('.')]
    program_name = "rnaQUAST"
    # parse running string of main program and get all arguments:
    #args = UtilsPipeline.get_arguments()
    #print(args.lower_threshold)
    #print(args.upper_threshold)
    if ("lower_threshold" not in self.parameters):
        self.parameters["lower_threshold"] = 0.5
    if ("upper_threshold" not in self.parameters):
        self.parameters["upper_threshold"] = 0.95
    WELL_FULLY_COVERAGE_THRESHOLDS = rqconfig.well_fully_coverage_thresholds(float(self.parameters["lower_threshold"]), float(self.parameters["upper_threshold"]))

    ALIGNMENT_THRESHOLDS = rqconfig.alignment_thresholds()

    # run rnaQUAST on test_data:
    
    #if ("test" not in self.parameters):
    #    self.parameters["test"] = True
    #if bool(self.parameters["test"]):
    #    UtilsPipeline.run_rnaQUAST_on_test_data(self.parameters, rquast_dirpath, program_name, logger)
    #    # UtilsPipeline.run_rnaQUAST_on_debug_data(args, rquast_dirpath, program_name, logger)
    #    sys.exit()

    #UtilsPipeline.get_abspath_input_data(self.parameters)

    # create output directory:
    output_dir = outputfile
    #args.output_dir = UtilsPipeline.create_output_folder(args.output_dir, program_name)
    # create temporary directory:
    tmp_dir = UtilsPipeline.create_empty_folder(os.path.join(output_dir, 'tmp'))
    # create directory for log files:
    log_dir = UtilsPipeline.create_empty_folder(os.path.join(output_dir, 'logs'))

    if ("debug" not in self.parameters):
        self.parameters["debug"] = True
    # SET LOGGER:
    if bool(self.parameters["debug"]):
        rqconfig.debug = True
        logger.set_up_console_handler(debug=True)
    else:
        logger.set_up_console_handler()
    logger.set_up_file_handler(log_dir)
    #logger.print_command_line([os.path.realpath(__file__)] + sys.argv[1:], wrap_after=None)
    if ("blat" not in self.parameters):
        self.parameters["blat"] = False
    if ("busco" not in self.parameters):
        self.parameters["busco"] = False
    if ("gene_mark" not in self.parameters):
        self.parameters["gene_mark"] = False
    logger.start(bool(self.parameters["blat"]), bool(self.parameters["busco"]), bool(self.parameters["gene_mark"]), tmp_dir)

    #UtilsPipeline.get_input_data_exist_error(args, logger)

    # THREADING:
    if ("threads" not in self.parameters):
        self.parameters["threads"] = 8
    #args.threads = UtilsPipeline.get_num_threads(args.threads, logger)

    if ("meta" not in self.parameters):
        self.parameters["meta"] = True
    if bool(self.parameters["meta"]):
        logger.info('\nYOU RUN QUALITY ASSESSMENT FOR METATRANSCRIPTOME ASSEMBLIES')
    # GET segregate FILES:
    #if args.reference and args.gtf and len(args.reference) != len(args.gtf):
    #    logger.error('Numbers of references and gene databases are different', exit_with_code=1)
    self.parameters["reference"] = [PyPluMA.prefix()+"/"+self.parameters["reference"]]
    self.parameters["GTF"] = [PyPluMA.prefix()+"/"+self.parameters["GTF"]]
    reference = \
        UtilsPipeline.get_single_file(self.parameters["reference"], tmp_dir, 'reference', rqconfig.list_ext_fa, self.parameters["meta"], logger)

    gtf = \
        UtilsPipeline.get_single_file(self.parameters["GTF"], tmp_dir, 'gene_database', rqconfig.list_ext_gtf, self.parameters["meta"], logger)

    # READ REFERENCE FROM MULTIFASTA:
    reference_dict = None
    ids_chrs = None
    if reference is not None:
        logger.print_timestamp()
        logger.info('Getting reference...')
        reference_dict = UtilsGeneral.list_to_dict(fastaparser.read_fasta(reference))
        logger.info('Done.')

        genome_len = UtilsGeneral.get_genome_len(reference_dict)

        ids_chrs = reference_dict.keys()

        # correction for fasta contained Y, W and etc:
        # for id_chr in ids_chrs:
        #     reference_dict[id_chr] = UtilsGeneral.correct_nucl_seq(reference_dict[id_chr])


    # for strand specific data we store + and - keys in dictionaries and only + for non strand specific data:
    if ("strand_specific" not in self.parameters):
        self.parameters["strand_specific"] = True
    strands = UtilsGeneral.get_strands(self.parameters, logger)


    if ("prokaryotes" not in self.parameters):
        self.parameters["prokaryote"] = True

    if bool(self.parameters["prokaryote"]):
        type_organism = 'prokaryotes'
    else:
        type_organism = 'eukaryotes'

    # USE ANNOTATION:
    sqlite3_db_genes = None
    sorted_exons_attr = None
    db_genes_metrics = None
    type_genes, type_isoforms, type_exons = \
        UtilsAnnotations.default_type_genes, \
        UtilsAnnotations.default_type_isoforms, \
        UtilsAnnotations.default_type_exons

    if gtf is not None or "gene_db" in self.parameters:
        if "gene_db" in self.parameters:
            gene_db_name = os.path.split(self.parameters["gene_db"])[1]
            label_db = gene_db_name[:gene_db_name.rfind('.db')]
        else:
            self.parameters["gene_db"] = None
            gtf_name = os.path.split(gtf)[1]
            label_db = gtf_name[:gtf_name.rfind('.g')]

            if ids_chrs is not None:
                gtf = UtilsAnnotations.clear_gtf_by_reference_chr(gtf, ids_chrs, tmp_dir, label_db, logger)

        if ("disable_infer_genes" not in self.parameters):
            self.parameters["disable_infer_genes"] = True
        if ("disable_infer_transcripts" not in self.parameters):
            self.parameters["disable_infer_transcripts`"] = True
        sqlite3_db_genes = \
            UtilsAnnotations.create_sqlite3_db(self.parameters["gene_db"], gtf, label_db,
                                               bool(self.parameters["disable_infer_genes"]), bool(self.parameters["disable_infer_transcripts"]),
                                               output_dir, tmp_dir, logger)

        type_genes, type_isoforms, type_exons = \
            UtilsAnnotations.get_type_features(sqlite3_db_genes, UtilsAnnotations.default_type_genes,
                                               UtilsAnnotations.default_type_isoforms,
                                               UtilsAnnotations.default_type_exons, self.parameters["prokaryote"], logger)

        # if UtilsAnnotations.default_type_exons == type_exons:
        #     type_organism = 'eukaryotes'
        # else:
        #     type_organism = 'prokaryotes'

        db_genes_metrics = GeneDatabaseMetrics.GeneDatabaseMetrics(sqlite3_db_genes, type_genes, type_isoforms, logger, self.parameters["prokaryote"])

        ALIGNMENT_THRESHOLDS.ERR_SPACE_TARGET_FAKE_BLAT = db_genes_metrics.max_intron_len + 100
        logger.info('\nSets maximum intron size equal {}. Default is 1500000 bp.\n'.format(ALIGNMENT_THRESHOLDS.ERR_SPACE_TARGET_FAKE_BLAT))

        # set exons starts / ends and ids for binning strategy:
        if ids_chrs is not None:
            sorted_exons_attr = \
                SortedExonsAttributes.SortedExonsAttributes(sqlite3_db_genes, type_exons, strands, ids_chrs, reference_dict, logger)

    reads_coverage = None
    if "reads_alignment" in self.parameters or \
            (("single_reads" in self.parameters or ("left_reads" in self.parameters and "right_reads" in self.parameters))
             and "reference" in self.parameters and "sqlite3_db_genes" in self.parameters):
        reads_coverage = \
            ReadsCoverage.ReadsCoverage(self.parameters["reads_alignment"], self.parameters["reference"], self.parameters["single_reads"],
                                        self.parameters["left_reads"], self.parameters["right_reads"], reference_dict, sqlite3_db_genes, type_isoforms,
                                        sorted_exons_attr, self.parameters["strand_specific"], db_genes_metrics.tot_isoforms_len,
                                        genome_len, tmp_dir, self.parameters["threads"], WELL_FULLY_COVERAGE_THRESHOLDS, logger, log_dir)


    if "transcripts" in self.parameters:
        mylist_t = []
        infile_t = open(PyPluMA.prefix()+"/"+self.parameters["transcripts"], 'r')
        for line in infile_t:
            mylist_t.append(PyPluMA.prefix()+"/"+line.strip())
        self.parameters["transcripts"] = mylist_t
        # GET TRANSCRIPTS:
        transcripts_dicts = []
        for i_transcripts in range(len(self.parameters["transcripts"])):
            logger.print_timestamp('  ')
            logger.info('  Getting transcripts from {}...'.format(self.parameters["transcripts"][i_transcripts]))
            transcripts_dicts.append(UtilsGeneral.list_to_dict(fastaparser.read_fasta(self.parameters["transcripts"][i_transcripts])))
            logger.info('  Done.')

        # get labels for folders names and names of transcripts in reports:
        all_labels_from_dirs = False
        if "labels" not in self.parameters:
            labels = UtilsPipeline.process_labels(self.parameters["transcripts"], None, all_labels_from_dirs)
    else:
        logger.warning('No transcripts. Use --transcripts option.')


    # GET PSL ALIGNMENT FILE:
    if "alignment" not in self.parameters and "reference" in self.parameters and "transcripts" in self.parameters:
        if self.parameters["blat"]:
            alignment = UtilsTools.run_blat(None, self.parameters["reference"], transcripts_dicts, labels,
                                                 self.parameters["threads"], tmp_dir, logger, log_dir)
        else:
            alignment = UtilsTools.run_gmap(self.parameters["reference"][0], genome_len, self.parameters["transcripts"], labels,
                                                 self.parameters["threads"], None, tmp_dir, logger, log_dir)

        #if args.fusion_misassemble_analyze:
        #    if not (args.left_reads is not None and args.right_reads is not None):
        #        logger.error('Usage: --left_reads LEFT_READS --right RIGHT_READS for analyse fusions and misassemblies',
        #                     exit_with_code=2, to_stderr=True)
        #        sys.exit(2)


    # FOR MISASSEMBLIES SEARCH:
    # GET DATABASE FOR FA ISOFORMS:
    blast = False
    if "reference" in self.parameters and sqlite3_db_genes is not None and "alignment" in self.parameters:
        blastn_run = os.path.join(rqconfig.rnaQUAST_LOCATION, '.', 'blastn')
        if not os.path.isfile(blastn_run):
            blastn_run = "blastn"

        if UtilsGeneral.which(blastn_run) is None:
            logger.warning('blastn not found! Please add blastn to PATH for better MISASSEMBLIES metrics.')
        else:
            blast = True

            isoforms_fa_path = os.path.join(tmp_dir, '{}.isoforms.fa'.format(label_db))
            isoforms_list = UtilsGeneral.dict_to_list(UtilsAnnotations.get_fa_isoforms(sqlite3_db_genes, type_isoforms, type_exons, reference_dict, logger))
            fastaparser.write_fasta(isoforms_fa_path, sorted(isoforms_list))

            isoforms_blast_db = UtilsTools.get_blast_db(isoforms_fa_path, label_db, tmp_dir, logger, log_dir)


    # LOGGING INPUT DATA:
    #logger.print_input_files(args)


    # INITIALIZATION TRANSCRIPTS METRICS AND REPORTS:
    transcripts_metrics = []
    separated_reports = []
    if "transcripts" in self.parameters:
        alignments_reports = []
        blast_alignments = []
        for i_transcripts in range(len(self.parameters["transcripts"])):
            # INITIALIZE TRANSCRIPTS METRICS:
            #if args.sam_file is not None:
            #    sam_file_tmp = args.sam_file[i_transcripts]
            #else:
            transcripts_metrics.append(
                TranscriptsMetrics.TranscriptsMetrics(self.parameters, labels[i_transcripts]))

            # INITIALIZE SEPARATED REPORTS:
            separated_reports.append(SeparatedReport.SeparatedReport(labels[i_transcripts], output_dir, transcripts_metrics[i_transcripts], WELL_FULLY_COVERAGE_THRESHOLDS))

            '''from joblib import Parallel, delayed

            n = len(self.parameters["transcripts"])
            run_n = n / self.parameters["threads"]
            for i_run in range(run_n):
                tmp = Parallel(n_jobs=self.parameters["threads"])(delayed(process_one_trascripts_file)(self.parameters, i_transcripts, reference_dict, annotation_dict,
                                                                                              annotated_exons, annotated_isoforms, strands, transcripts_metrics,
                                                                                              basic_isoforms_metrics, separated_reports)
                                                         for i_transcripts in range(i_run * self.parameters["threads"], self.parameters["threads"] * (i_run + 1), 1))
                for i in range(self.parameters["threads"]):
                    i_transcripts = i + i_run * self.parameters["threads"]
                    transcripts_metrics[i_transcripts] = tmp[i][0]
                    separated_reports[i_transcripts] = tmp[i][1]

            if n - run_n * self.parameters["threads"] != 0:
                tmp = Parallel(n_jobs=n - run_n * self.parameters["threads"])(delayed(process_one_trascripts_file)(self.parameters, i_transcripts, reference_dict, annotation_dict,
                                                                                                     annotated_exons, annotated_isoforms, strands, transcripts_metrics,
                                                                                                     basic_isoforms_metrics, separated_reports)
                                                                for i_transcripts in range(run_n * self.parameters["threads"], n, 1))
                for i in range(n - run_n * self.parameters["threads"]):
                    i_transcripts = i + run_n * self.parameters["threads"]
                    transcripts_metrics[i_transcripts] = tmp[i][0]
                    separated_reports[i_transcripts] = tmp[i][1]'''

            logger.info()
            logger.info('Processing transcripts from {}:'.format(self.parameters["transcripts"][i_transcripts]))

            if blast:
                blast_alignments.append\
                    (UtilsTools.align_transcripts_to_isoforms_by_blastn
                     (self.parameters["transcripts"][i_transcripts], isoforms_blast_db, tmp_dir, labels[i_transcripts], logger, log_dir))
            else:
                blast_alignments.append(None)

            # PROCESS TRANSCRIPTS ALIGNMENTS:
            if transcripts_metrics[i_transcripts].simple_metrics is not None:
                # GET FILES WITH ALIGNMENTS REPORTS:
                alignments_reports.append\
                    (UtilsAlignment.AlignmentsReport.get_alignments_report
                     (labels[i_transcripts], self.parameters["alignment"][i_transcripts], blast_alignments[i_transcripts],
                      transcripts_dicts[i_transcripts], tmp_dir, self.parameters["min_alignment"], logger, ALIGNMENT_THRESHOLDS))

                # UPDATE METRICS BY ASSEMBLED TRANSCRIPTS:
                transcripts_metrics[i_transcripts].processing_assembled_psl_file\
                    (alignments_reports[i_transcripts].blat_report.assembled_psl_file, sorted_exons_attr,
                     self.parameters["strand_specific"], logger, sqlite3_db_genes, type_isoforms, WELL_FULLY_COVERAGE_THRESHOLDS)

                # UPDATE METRICS BY MISASSEMBLED TRANSCRIPTS:
                # by blat:
                transcripts_metrics[i_transcripts].processing_misassembled_psl_file\
                    (alignments_reports[i_transcripts].blat_report.misassembled_psl_union_file, logger, True)
                # by blast:
                if blast:
                    transcripts_metrics[i_transcripts].processing_misassembled_psl_file\
                        (alignments_reports[i_transcripts].blast6_report.misassembled_blast6_union_file, logger, False)

            # GET METRICS:
            transcripts_metrics[i_transcripts].get_transcripts_metrics\
                (self.parameters, type_organism, reference_dict, self.parameters["transcripts"][i_transcripts], transcripts_dicts[i_transcripts],
                 labels[i_transcripts], self.parameters["threads"], sqlite3_db_genes, db_genes_metrics, reads_coverage, logger,
                 tmp_dir, log_dir, WELL_FULLY_COVERAGE_THRESHOLDS, rqconfig.TRANSCRIPT_LENS)

            # GET SEPARATED REPORT:
            separated_reports[i_transcripts].get_separated_report\
                (self.parameters, labels[i_transcripts], transcripts_dicts[i_transcripts], transcripts_metrics[i_transcripts],
                 db_genes_metrics, reads_coverage, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, rqconfig.TRANSCRIPT_LENS)

    # GET COMPARISON REPORT:
    comparison_report = None
    if len(separated_reports) != 1:
        comparison_report = ComparisonReport.ComparisonReport()
        comparison_report.get_comparison_report(self.parameters, output_dir, labels, transcripts_metrics,
                                                db_genes_metrics, reads_coverage, logger,
                                                WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION, rqconfig.TRANSCRIPT_LENS)

    # GET SHORT REPORT:
    short_report = \
        ShortReport.ShortReport(self.parameters, db_genes_metrics, transcripts_metrics, output_dir, separated_reports,
                                comparison_report, logger, WELL_FULLY_COVERAGE_THRESHOLDS, PRECISION,
                                rqconfig.TRANSCRIPT_LENS)

    # REMOVE TEMPORARY DIRECTORY FROM OUTPUT DIRECTORY:
    if os.path.exists(tmp_dir) and not self.parameters["debug"]:
        logger.debug('Remove temporary directory {}'.format(tmp_dir))
        shutil.rmtree(tmp_dir)
        logger.debug('Done.')

    # LOGGING RESULTS PATHES:
    logger.print_path_results(output_dir, separated_reports, comparison_report, short_report)

    #if self.parameters["debug"]:
    #    UtilsGeneral.profile_memory(self.parameters, reference_dict, db_genes_metrics, transcripts_metrics,
    #                                separated_reports, comparison_report, logger)

    # FINISH LOGGING:
    logger.finish_up()


#if __name__ == '__main__':
#    try:
#        return_code = main_utils()
#        exit(return_code)
#    except Exception:
#        _, exc_value, _ = sys.exc_info()
#        logger.exception(exc_value)
#        logger.error('Exception caught!', exit_with_code=1)
