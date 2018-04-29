from copy import copy
from ott_objects_defs.ClusterContext import ClusterContext, SeqType
from ott_objects_defs.ItsClusterContext import ItsClusterContext
from ott_objects_defs.ConcatClusterContext import ConcatClusterContext
from ott_objects_defs.ploidbCommon import *
import zlib
from ott_objects_defs.PathHelper import PathHelper

__author__ = 'moshe'


class PloiDbContext:
    PICKLED_FILE_EXTENSION = ".pkl"
    FIRST_NON_ITS_CLUSTER_IDX = 0
    # NUM_OF_CLUSTERS_FOR_CONCAT = 10
    MAX_NUMBER_OF_LOCUS_TREES = 7

    def __init__(self, id, working_dir, cluster_script_dir, cluster_method, taxa_list):
        self.p = PathHelper(working_dir, id)
        self.number_of_clusters = 0
        self.RerunID = 'None'
        self.OTTOptions = 'phylogenyTreeReconstruction' #default value
        self.its_support = False
        self.its_min_taxa = 'yes' #Check if ITS include the minimum number of species needed
        self.split_by_feature = False
        self.cluster_contexts = list()
        self.irrelevant_cluster_contexts = list()
        self.full_cluster_contexts = list()
        self.unresolved_species_names = set()
        self.rnr_concat_filtered_species = set()
        self.cluster_contexts_without_its = list()
        self.cluster_contexts_for_loci_trees = list()
        self.cluster_contexts_for_concat_tree_candidates = None
        self.cluster_contexts_for_concat_tree_full = list()  # Holds  all the clusters that are good for concat
        self.cluster_contexts_for_concat_tree = list()  # Holds  all the clusters that would actually be used for concat
        self.its1Cluster = None
        self.its2Cluster = None
        self.chloroplast_cluster = None
        self.Mt_cluster = None
        self.Nuc_cluster = None
        self.ClusterType_dict = dict()
        self.gb_seqs_dict = None
        self.gbDict = None
        self.standAlone_Flag = 'off' # default is running as webserver (optional to be Stand alone)
        self.top_outgroup_contexts = list()
        self.taxa_list = taxa_list
        self.species_list_names = list()
        self.species_list_ids = list()
        self.parent_species_tax_id_by_subsp_tax_id = dict()
        self.parent_species_names_by_parent_tax_id = dict()
        self.parent_species_names_by_subsp_tax_id = dict()
        self.seqs_number_by_taxId = dict()
        self.species_names_originals = list()
        self.list_final_species_names = list()
        self.names_vs_accepted_dict = dict()
        self.names_vs_CodedNames_dict = dict()
        # For cluster concat selection:
        self.Nuc_clusterMode = 1
        self.Mt_clusterMode = 1
        self.Chloro_clusterMode = 1

        self.TaxId_rank_dict=dict()
        self.Name_type_dict=dict()
        self.synTaxId_rank_dict=dict()
        self.TaxId_name_dict=dict()
        self.Name_TaxId_dict=dict()
        self.synTaxId_name_dict=dict()
        self.TaxId_Parent_dict=dict()
        self.synTaxId_Parent_dict=dict()
        self.original_vs_newName_dict=dict()
        self.spc_tax_afterMatch_dict=dict()

        #for user: some indication regarding species that were requested by the user but were not included in the final tree:
        self.Pipe_input_Spc_Names_vs_TaxId=dict()
        self.list_not_fount_in_ncbi=list()
        self.list_not_included_in_msa=list()

        self.large_taxid_seqs_number_dict = dict()
        self.taxIds_not_in_genbank_list = list()

        self.align_method_diff = 'false' # original alignment method is same as in rerun
        self.rerun_add_species_list = list()
        self.rerun_remove_species_list = list()
        self.rerun_chosen_clusters = list()
        self.rerun_new_its_files_list = list()
        self.rerun_new_non_its_files_list = list()


        self.taxon_Id_df = dict()
        self.Accession_df = dict()
        self.organism_df = dict()
        self.seq_length_df = dict()
        self.definition_df = dict()
        self.gene_name_location_df = dict()
        self.desc_df = dict()
        self.note_df = dict()
        self.mol_type_df = dict()
        self.product_df = dict()
        self.sequenceData_df = dict()
        self.outgroup_AccessionkeySeq_dict = dict()
        self.its_accession_ids = list()
        self.its_accession_vs_type = dict() #its1, its2,combined
        # In case blast results are missing-> rep seq will be chosen according to length
        self.its_accession_vs_length = dict()
        self.its_type_avg_len = dict()
        self.its_taxa_vs_counts = dict()


        self.CodedNames_TaxIds_dict = dict()
        self.MatchedNames_TaxIds_dict = dict()
        self.TaxId_Accepted_dict = dict()
        self.species_names_by_tax_id = dict()
        self.Species_names_vs_TaxId_AfterNR = dict()
        self.final_species_names_by_tax_id = dict()
        self.UserFlags_dict = dict()
        self.syn_acc_dict = dict()  #Keeps the Synonym vs Accepted names given by the user or the default TPL table
        self.syn_acc_TaxID_dict = dict() #Keeps the Synonym vs Accepted names given by the user or the default TPL table
        self.UserOutgroupInc = 'NO'
        self.User_Outgroup = 'None'
        self.UserOutgroupTaxId = '-'
        self.UserOutgroupName = '-'

        #default values for flags:
        self.UserFlag_mbUserModel = 'JmodelTest' # or Spesific Model, e.g. GTR+G
        self.UserFlag_runGuidance = 'None' #Yes/No
        self.UserFlag_OutgroupSel = 'None' #No Outgroup, Single or Multiple
        self.UserFlag_MSAmethod = 'mafft'   #mafft, Prank,ClustaOmega
        self.UserFlag_MSAfilterMethod = 'None' #Guidance, Gblocks,

        self.should_run_guidance = True
        self.should_ignore_subsp = False
        self.should_merge_subsp = False

        self.genus_id = zlib.crc32(bytes(id, 'UTF-8'))
        self.outgroupSelection = None

        self.id = id
        self.cluster_method = cluster_method
        if self.cluster_method == 'none':
            self.its_support = False
        if self.cluster_method == 'orthomclits':
            self.its_support = True
        if self.cluster_method == 'orthomclitsfeature':
            self.its_support = True
            self.split_by_feature = True

        self.cluster_script_dir = cluster_script_dir
        self.ott_scripts_dir = cluster_script_dir + '/ott_scripts/'

        #
        # General stats
        #
        self.start_time = None
        self.end_pre_tree_time = None
        self.end_time = None
        self.init_context_dir(working_dir)

    def __getstate__(self):
        odict = self.__dict__.copy()
        if 'gb_seqs_dict' in odict:
            del odict['gb_seqs_dict']
        return odict

    def init_context_dir(self, working_dir):
        working_dir = os.path.abspath(working_dir) + '/'
        self.pickled_name = working_dir + "/" + self.id + "-" + self.cluster_method + self.PICKLED_FILE_EXTENSION
        self.pickled_name_after_trees = working_dir + "/" + self.id + "-" + "-aftertrees" + self.PICKLED_FILE_EXTENSION
        self.pickled_json_name = working_dir + "/" + self.id + "-" + self.cluster_method + ".json"
        self.pickled_json_name_after_trees = working_dir + "/" + self.id + "-aftertrees" + ".json"

        self.working_dir = working_dir
        self.largeSpeciesDir = os.path.join(working_dir, "LargeSpecies")
        self.largeSpeciesITS = os.path.join(working_dir, "LargeSpecies/Large_ITS.fasta")
        self.largeSpeciesTaxIDsList = []
        self.largeSpeciesFileList = os.path.join(working_dir, "LargeSpecies/LargeFastaFiles_list.txt")
        self.largeSynToAdd_file = os.path.join(working_dir, "LargeSpecies/Syn_toAddto_LargeTaxon.fasta")
        self.pre_tree_done = os.path.join(working_dir, "pre_tree_done")
        self.create_tree_done = os.path.join(working_dir, "create_tree_done")
        self.db_file = os.path.join(working_dir, "ploidb.db")
        self.names_to_resolve = os.path.join(working_dir, "names_to_resolve.csv")
        self.name_resolve_results = os.path.join(working_dir, "resolved_names.csv")
        self.nr_analyze = os.path.join(working_dir, "NR_analyze.csv")
        self.subsp_variants_merge_names = os.path.join(working_dir, "subsp_var_merge_names.txt")
        self.ploidb_db_file = os.path.join(working_dir, "../../ploidb-full.db")
        self.fasta_seq_filename = os.path.join(working_dir, self.id + "-allseq" + ".fasta")
        self.fasta_cdhit_temp_input_f = os.path.join(working_dir, self.id + "_temp_input_cdhit" + ".fasta")
        self.fasta_cdhit_temp_output_f = os.path.join(working_dir, self.id + "_temp_output_cdhit" + ".fasta")
        self.fasta_seq_merged_subsp_filename = os.path.join(working_dir, self.id + "-merged-subsp-allseq" + ".fasta")
        self.fasta_seq_merged_syn_filename = os.path.join(working_dir, self.id + "-merged-syn-allseq" + ".fasta")
        self.fasta_seq_org_names_filename = os.path.join(working_dir, self.id + "-orgnames-allseq" + ".fasta")
        self.gb_seq_filename = os.path.join(working_dir, self.id + "-allseq" + ".gb")
        self.gb_new_seq_filename = os.path.join(working_dir, self.id + "-new-seq" + ".gb")   #For added species - Rerun options
        self.outgroup_workdir = os.path.join(working_dir, "outgroup-processing")
        self.mrbayes_input_dir = working_dir
        self.mrbayes_output_dir = working_dir
        self.clustering_dir = os.path.join(working_dir, "clustering")
        self.blast_results_filename = self.clustering_dir + "/blast_all-v-all_output.blastn"
        self.clustering_dir_inner = self.clustering_dir + "/inner"
        self.clustering_results_dir = os.path.join(working_dir, "cluster-results")
        self.cluter_its_dir = os.path.join(working_dir, "cluster_its")
        self.species_names_taxid_file = os.path.join(working_dir, "species_names_file.csv")
        self.outgroup_file = os.path.join(working_dir, "OutgroupSelection.txt")
        self.summary_dir = os.path.join(working_dir, "SummaryDir")
        self.debug_dir = os.path.join(working_dir, "DebugDir")
        self.summary_clusters_data_file = os.path.join(working_dir, "SummaryDir/All_Clusters_data.csv")
        self.summary_accession_matrix_file = os.path.join(working_dir, "SummaryDir/AccessionsMatrix.csv")
        self.accession_loci_matrix_dict = dict()
        self.loci_description_dict = dict()
        self.summary_Chosen_clusters_data_file = os.path.join(working_dir, "SummaryDir/Chosen_Clusters_data.txt")
        self.summary_file = os.path.join(working_dir, "SummaryDir/summary_file.txt")
        self.large_species_data = os.path.join(self.debug_dir, "large_species_data.txt")
        self.final_species_file = os.path.join(self.summary_dir,"FinalSpeciesList.txt")
        self.final_status = os.path.join(working_dir, "SummaryDir/FinalStatus.txt")
        self.TaxId_NoData_list_file = os.path.join(working_dir, "DebugDir/TaxIds_without_GenabnkSeqs.txt")

        # This part is for the concat part
        #
        self.mrbayes_concat_work_dir = os.path.join(self.mrbayes_input_dir, "concat")
        self.raxml_concat_work_dir = os.path.join(self.mrbayes_concat_work_dir, "raxml")
        self.clusters_to_concat_filename = os.path.join(self.mrbayes_concat_work_dir, "cluster_to_concat_list.txt")
        self.concat_workdir = self.mrbayes_concat_work_dir
        self.concat_seqs_fasta_filename = os.path.join(self.mrbayes_concat_work_dir,
                                                       self.id + "-concat-aligned.fasta")
        self.concat_seqs_report_filename = os.path.join(self.mrbayes_concat_work_dir,
                                                        self.id + "-concat-aligned.report")
        self.xml_dir = os.path.join(self.mrbayes_concat_work_dir,'XML_Run_dir')
        self.concat_xml_partition_filename = os.path.join(self.xml_dir,
                                                        self.id + "_xml_partition")
        self.fasta_files_list_to_concat = os.path.join(self.mrbayes_concat_work_dir, "fasta-files-to-concat.txt")
        self.concat_log_filename = os.path.join(self.mrbayes_concat_work_dir, "concat-alignment.out")
        self.concat_log_out_filename = os.path.join(self.mrbayes_concat_work_dir, "concat-alignment.out")
        self.concat_log_err_filename = os.path.join(self.mrbayes_concat_work_dir, "concat-alignment.err")

        self.rewrite_names_dict_filename = os.path.join(self.working_dir, "rewritten-names-concat.json")
        self.names_to_resolved_names = dict()

        self.concat_by_type_workdir = dict()
        self.concat_by_type_workdir[SeqType.Chloroplast] = os.path.join(self.mrbayes_concat_work_dir)
        self.concat_by_type_workdir[SeqType.Mt] = os.path.join(self.mrbayes_concat_work_dir)        # Check if needed
        self.concat_by_type_workdir[SeqType.Nuc] = os.path.join(self.mrbayes_concat_work_dir)
        self.concat_by_type_fastas = dict()
        self.concat_by_type_fastas[SeqType.Chloroplast] = os.path.join(self.concat_by_type_workdir[SeqType.Chloroplast],
                                                                       "chloroplast-all.fasta")
        self.concat_by_type_fastas[SeqType.Mt] = os.path.join(self.concat_by_type_workdir[SeqType.Mt],
                                                                       "Mt-all.fasta")
        self.concat_by_type_fastas[SeqType.Nuc] = os.path.join(self.concat_by_type_workdir[SeqType.Nuc],
                                                                       "Nuc-all.fasta")
        self.concat_by_type_outgroup_fasta = dict()
        self.concat_by_type_outgroup_fasta[SeqType.Chloroplast] = os.path.join(self.concat_by_type_workdir[SeqType.Chloroplast], "chloroplast-outgroup.fasta")
        self.concat_by_type_outgroup_fasta[SeqType.Mt] = os.path.join(self.concat_by_type_workdir[SeqType.Mt], "Mt-outgroup.fasta")
        self.concat_by_type_outgroup_fasta[SeqType.Nuc] = os.path.join(self.concat_by_type_workdir[SeqType.Nuc], "Nuc-outgroup.fasta")
        self.concat_by_type_seq_and_outgroup_fasta = dict()
        self.concat_by_type_seq_and_outgroup_fasta[SeqType.Chloroplast] = os.path.join(self.concat_by_type_workdir[SeqType.Chloroplast], "chloroplast-seq-and-out.fasta")
        self.concat_by_type_seq_and_outgroup_fasta[SeqType.Mt] = os.path.join(self.concat_by_type_workdir[SeqType.Mt], "Mt-seq-and-out.fasta")
        self.concat_by_type_seq_and_outgroup_fasta[SeqType.Nuc] = os.path.join(self.concat_by_type_workdir[SeqType.Nuc], "Nuc-seq-and-out.fasta")

        self.rerun_new_taxIds_dir = os.path.join(working_dir, "rerun_new_taxIds")

        #
        # This part is for the MB part
        #
        self.mb_config_file = os.path.join(self.mrbayes_concat_work_dir,
                                           "mb_config.nex")  # just for the configuration creation
        self.mb_config_file1 = os.path.join(self.mrbayes_concat_work_dir,
                                            "mb_config.nex1")  # run until 100,000 iterations
        self.mb_config_file2 = os.path.join(self.mrbayes_concat_work_dir,
                                            "mb_config.nex2")  # run from checkpoint until convergence/2,000,000 iterations
        self.mb_seq_file = os.path.join(self.mrbayes_concat_work_dir, "mb_final_seq.nex")
        self.mb_out_file = os.path.join(self.mrbayes_concat_work_dir, "mb.out")
        self.mb_out_consensus_file = os.path.join(self.mrbayes_concat_work_dir, "mb.out.con.tre")
        self.mb_out_consensus_file_newick = os.path.join(self.mrbayes_concat_work_dir, "mb.out.con.newick")

        self.mb_concat_log_out_filename1 = os.path.join(self.mrbayes_concat_work_dir, "mb-log1.out")
        self.mb_concat_log_out_filename2 = os.path.join(self.mrbayes_concat_work_dir, "mb-log2.out")
        self.mb_concat_log_err_filename1 = os.path.join(self.mrbayes_concat_work_dir, "mb-log1.err")
        self.mb_concat_log_err_filename2 = os.path.join(self.mrbayes_concat_work_dir, "mb-log2.err")


        self.tree_xml_namedate_config = os.path.join(self.summary_dir, "node_date_config.txt")
        self.NodeDatingFile = os.path.join(working_dir, "NodeDate.txt")

        # Just for cache:
        self.mraic_ouput_file = "seqs-organism-concat_phy.phy.MrAIC.txt"

        #
        # This part is for the chromosome numbers part
        #
        self.chr_counts_filename = os.path.join(working_dir, "/chromosome_counts.fa")
        # This part is for the chromevol part
        self.out_trees_filename = os.path.join(working_dir, "parsemb_trees.tre")
        self.out_consensus_tree_filename = os.path.join(working_dir, "parsemb_con_tree.tre")
        self.out_map_tree_filename = os.path.join(working_dir, "parsemb_map_tree.tre")
        self.mr_bayes_post_process_path = os.path.join(working_dir, "mrbayes")

    def reinit_context_dir(self, working_dir):
        self.init_context_dir(os.path.abspath(working_dir))

        for cc in self.cluster_contexts:
            cc.init_dirs(self)


    def get_concat_clusters_nuc(self):
        return self.get_concat_clusters_by_type(SeqType.Nuc)

    def get_cluster_by_id(self, cluster_id, clusters=None):
        if clusters is None:
            clusters = self.cluster_contexts
        m = dict((cc.cluster_id, cc) for cc in clusters)
        return m.get(cluster_id)

    def get_clusters_by_ids(self, cluster_ids, clusters=None):
        if clusters is None:
            clusters = self.cluster_contexts
        clusters_to_ret = list()
        for cc in clusters:
            if cc.cluster_id in cluster_ids:
                clusters_to_ret.append(cc)

        return clusters_to_ret

    def get_all_clusters_ids(self, order_func=None):
        logger.debug("cluster_contexts before sort %s" % ",".join(str(cc) for cc in self.cluster_contexts))
        if order_func is None:
            order_func = lambda cluster_context: cluster_context.get_data_matrix_size(True)
        # Sorting the cluster contexts according to the number of species in each cluster
        cc_list_copy = sorted(self.cluster_contexts, key=order_func, reverse=True)
        logger.debug("cluster_contexts after sort %s" % ",".join(str(cc) for cc in cc_list_copy))

        cluster_ids = list()
        for cc in cc_list_copy:
            cluster_ids.append(cc.cluster_id)

        logger.debug("Returning %s" % ",".join(cluster_ids))
        return cluster_ids

    def get_cluster_contexts_for_concat_tree_candidate_ids(self):
        ids = list()
        for cc in self.cluster_contexts_for_concat_tree_candidates:
            ids.append(cc.cluster_id)
        return ids

    def get_clusters_by_type(self, cluster_type, clusters):
        if not isinstance(cluster_type, list):
            cluster_type = [cluster_type]
        cluster_contexts_by_type_list = list()
        logger.debug("Looking for clusters of type %s out of a total of %d clusters" % (
            ','.join(str(v) for v in cluster_type), len(clusters)))
        for cc in clusters:
            if cc.cluster_type in cluster_type:
                cluster_contexts_by_type_list.append(cc)
        logger.debug("Found a total of %d clusters of type %s" % (
            len(cluster_contexts_by_type_list), ','.join(str(v) for v in cluster_type)))
        return cluster_contexts_by_type_list

    def get_clusters_contexts(self, clusters):
        logger.debug("----------   Get all clusters contexts for concat   ----------")
        ids = list()
        for cc in self.full_cluster_contexts:
            logger.debug(cc.cluster_id)
            ids.append(cc.cluster_id)
        return ids

    def get_clusters_for_concat_with_outgroup(self):
        cluster_contexts_with_outgroup_list = list()
        logger.debug("Getting clusters with outgroup")
        for cc in self.cluster_contexts_for_concat_tree:
            logger.debug(
                "Checking %s cc.no_of_outgroup_sequence_for_concat=%d" % (cc, cc.no_of_outgroup_sequence_for_concat))
            if cc.no_of_outgroup_sequence_for_concat > 0:
                cluster_contexts_with_outgroup_list.append(cc)
        logger.debug("Found a total for %d clusters with outgroup" % (len(cluster_contexts_with_outgroup_list)))
        return cluster_contexts_with_outgroup_list


    def get_locus_clusters_by_type(self, cluster_type):
        return self.get_clusters_by_type(cluster_type, self.cluster_contexts_for_loci_trees)

    def get_concat_clusters_by_type(self, cluster_type):
        return self.get_clusters_by_type(cluster_type, self.cluster_contexts_for_concat_tree)

    def get_concat_clusters(self):
        return self.get_clusters_contexts(self.cluster_contexts_for_concat_tree)

    def add_cluster(self, cluster_fasta_file):
        cluster_context = self.__add_cluster(cluster_fasta_file)
        if cluster_context is not None:
            self.cluster_contexts_without_its.append(cluster_context)
        return cluster_context

    def init_chloroplast_cluster(self):
        logger.info("Initializing chloroplast cluster")
        chloroplast_list = self.get_concat_clusters_by_type(SeqType.Chloroplast)
        if len(chloroplast_list) > 0:
            self.chloroplast_cluster = ConcatClusterContext(ConcatClusterContext.CHLOROPLAST_CLUSTER_INDEX,
                                                            self.concat_by_type_fastas[SeqType.Chloroplast], self,
                                                            chloroplast_list, SeqType.Chloroplast)
            self.chloroplast_cluster.init_cluster(self)
            logger.info("nexus file names: %s %s", self.chloroplast_cluster.mb_config_file1,
                        self.chloroplast_cluster.mb_config_file2)
            logger.info("DONE Initializing chloroplast cluster")
        else:
            logger.info("No clusters relevant for chloroplast")

    def init_mitochondrial_cluster(self):
        logger.info("Initializing mitochondrial cluster")
        mitochondrial_list = self.get_concat_clusters_by_type(SeqType.Mt)
        if len(mitochondrial_list) > 0:
            self.Mt_cluster = ConcatClusterContext(ConcatClusterContext.MITOCHONDRIAL_CLUSTER_INDEX,
                                                            self.concat_by_type_fastas[SeqType.Mt], self,
                                                            mitochondrial_list, SeqType.Mt)
            self.Mt_cluster.init_cluster(self)
            logger.info("nexus file names: %s %s", self.Mt_cluster.mb_config_file1,
                        self.Mt_cluster.mb_config_file2)
            logger.info("DONE Initializing mitochondrial cluster")
        else:
            logger.info("No clusters relevant for mitochondrial")

    def init_nuc_cluster(self):
        logger.info("Initializing nucleus cluster")
        nucleus_list = self.get_concat_clusters_by_type(SeqType.Nuc)
        if len(nucleus_list) > 0:
            self.Nuc_cluster = ConcatClusterContext(ConcatClusterContext.NUCLEUS_CLUSTER_INDEX,
                                                            self.concat_by_type_fastas[SeqType.Nuc], self,
                                                            nucleus_list, SeqType.Nuc)
            self.Nuc_cluster.init_cluster(self)
            logger.info("nexus file names: %s %s", self.Nuc_cluster.mb_config_file1,
                        self.Nuc_cluster.mb_config_file2)
            logger.info("DONE Initializing nucleus cluster")
        else:
            logger.info("No clusters relevant for nucleus")

    def add_its_cluster(self, cluster_fasta_file, its1_fasta_file, its2_fasta_file, combined_all_seqs_fasta_filename):
        allRecords = list(SeqIO.parse(cluster_fasta_file, "fasta"))
        num_of_cluster_seqs = len(allRecords)

        if num_of_cluster_seqs <= 0:
            return None
        else:
            self.number_of_clusters += 1
            its_context = ItsClusterContext(cluster_fasta_file, self, its1_fasta_file, its2_fasta_file,
                                            combined_all_seqs_fasta_filename)
            self.itsClusterContext = its_context
            return its_context


    def add_its_final_cluster(self):
        cluster_fasta_file = os.path.join(self.cluter_its_dir, "oneSeqPerSpecies.msa")
        if not os.path.exists(cluster_fasta_file):
            return None
        restandarize_description_after_mafft(cluster_fasta_file)

        its_context = self.__add_cluster(cluster_fasta_file)
        self.its_final_cluster = its_context
        if its_context is not None:
            self.its_final_cluster.is_its_cluster = True
            #self.init_its_context_seq_to_blast(its_context, self.cluter_its_dir)
            #fasta_forAlign_its = os.path.join(self.cluter_its_dir, "seqs-no-multiple-accessions.fasta")
            #shutil.copyfile(cluster_fasta_file, fasta_forAlign_its)
        return its_context

    # TODO: don't take the first - take according to blast results
    def init_its_context_seq_to_blast(self, its_context, its_clustering_results_dir):
        combined_filename = os.path.join(its_clustering_results_dir, "combined.fasta")
        its1_filename = os.path.join(its_clustering_results_dir, "ITS1_only.fasta")
        its2_filename = os.path.join(its_clustering_results_dir, "ITS2_only.fasta")
        cluster_seq_to_blast = get_first_seq_from_file(combined_filename)

        if cluster_seq_to_blast is None:
            cluster_seq_to_blast = get_first_seq_from_file(its1_filename)
        if cluster_seq_to_blast is None:
            cluster_seq_to_blast = get_first_seq_from_file(its2_filename)
        if cluster_seq_to_blast is None:
            raise Exception("Couldn't find seq to blast for ITS cluster at %s" % its_clustering_results_dir)

        its_context.cluster_seq_to_blast = cluster_seq_to_blast

    def add_its1_cluster(self, cluster_fasta_file):
        its_context = self.__add_cluster(cluster_fasta_file)
        self.its1Cluster = its_context
        return its_context

    def add_its2_cluster(self, cluster_fasta_file):
        its_context = self.__add_cluster(cluster_fasta_file)
        self.its2Cluster = its_context
        return its_context

    def __add_cluster(self, cluster_fasta_file):
        allRecords = list(SeqIO.parse(cluster_fasta_file, "fasta"))
        num_of_cluster_seqs = len(allRecords)

        if num_of_cluster_seqs <= 0:
            return None
        else:
            self.number_of_clusters += 1
            its_context = ClusterContext(self.number_of_clusters, cluster_fasta_file, self)
            its_context.init_cluster(self)
            self.cluster_contexts.append(its_context)
            return its_context

    def get_newick_concat_consensus_tree(self):
        return get_newick_tree_str(self.mb_out_consensus_file_newick)

    def init_cluster_contexts_for_loci_trees_list(self):
        ordered_clusters = sorted(self.cluster_contexts,
                                  key=lambda cluster: cluster.get_data_matrix_size(estimated=True), reverse=True)
        for cluster in ordered_clusters:
            if len(
                    self.cluster_contexts_for_loci_trees) <= self.MAX_NUMBER_OF_LOCUS_TREES and cluster.no_of_outgroup_sequence_for_concat > 0:
                cluster.is_used_for_locus_tree = True
                self.cluster_contexts_for_loci_trees.append(cluster)


    # Updating the relevant clusters, so that very small clusters won't be handled
    def optimize_active_clusters(self):
        relevant_clusters = list()
        irrelevant_clusters = list()
        data_matrix_size_threshold = None
        data_matrix_size_sum = 0
        NUM_OF_CLUSTERS_TO_AVERAGE_DATA_MATRIX_SIZE = 5
        MIN_NUMBER_OF_CLUSTERS = 10

        for ordinal, cluster_context in enumerate(self.cluster_contexts):
            # Calculating the data_matrix_size_threshold
            if ordinal <= NUM_OF_CLUSTERS_TO_AVERAGE_DATA_MATRIX_SIZE:
                data_matrix_size_sum += cluster_context.get_data_matrix_size(True)
            if ordinal > NUM_OF_CLUSTERS_TO_AVERAGE_DATA_MATRIX_SIZE and data_matrix_size_threshold is None:
                data_matrix_size_threshold = data_matrix_size_sum / NUM_OF_CLUSTERS_TO_AVERAGE_DATA_MATRIX_SIZE
                data_matrix_size_threshold *= 0.1

            if ordinal <= MIN_NUMBER_OF_CLUSTERS:
                relevant_clusters.append(cluster_context)
            else:
                if cluster_context.get_data_matrix_size(True) < data_matrix_size_threshold:
                    logger.info("Irrelevant cluster %s - data matrix size is only %d while the minimum is %d" %
                                (cluster_context, cluster_context.get_data_matrix_size(True),
                                 data_matrix_size_threshold))
                    irrelevant_clusters.append(cluster_context)
                else:
                    relevant_clusters.append(cluster_context)
        if len(relevant_clusters) + len(irrelevant_clusters) != len(self.cluster_contexts):
            raise Exception("Internal error - cluster division went wront len(relevant_clusters)=%d "
                            "len(irrelevant_clusters)=%d  len(cluster_context)=%d" % (
                                len(relevant_clusters), len(irrelevant_clusters), len(self.cluster_contexts)))
        self.full_cluster_contexts = self.cluster_contexts
        self.cluster_contexts = relevant_clusters
        self.irrelevant_cluster_contexts = irrelevant_clusters
        if len(irrelevant_clusters) == 0:
            logger.info("No irrelevant clusters were found - all will be used")
