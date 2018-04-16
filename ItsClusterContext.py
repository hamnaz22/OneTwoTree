from ClusterContext import ClusterContext

__author__ = 'moshe'

class ItsClusterContext(ClusterContext):
    ITS_CLUSTER_INDEX = 99


    def __init__(self,all_seqs_fasta_filename,genus_cluster,its1_all_seqs_fasta_filename,its2_all_seqs_fasta_filename,combined_all_seqs_fasta_filename):
        super(ClusterContext,self).__init__(self.ITS_CLUSTER_INDEX,all_seqs_fasta_filename,genus_cluster)
        self.init_cluster()
        self.its1_all_seqs_fasta_filename = its1_all_seqs_fasta_filename
        self.its2_all_seqs_fasta_filename = its2_all_seqs_fasta_filename
        self.combined_all_seqs_fasta_filename = combined_all_seqs_fasta_filename

        self.joint_its_all_seqs_fasta_filename = all_seqs_fasta_filename

