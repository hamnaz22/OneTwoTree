CREATE TABLE blast_results (outgroup_seq_id TEXT, taxonomy_last TEXT, organism TEXT,gi TEXT,gb TEXT, sacc TEXT,saccver TEXT,evalue REAL, bitscore REAL,score REAL,pident REAL,ppos REAL, qlen REAL, slen REAL, length REAL, is_cluster_gi INTEGER default 0, is_seq_in_species_list INTEGER default 0,is_min_score_for_acc INTEGER default 0);
CREATE VIEW blast_results_filtered AS
  SELECT *
  FROM blast_results
  WHERE evalue <= 10e-10 AND is_min_score_for_acc = 1
  AND (((qlen*1.0 / slen > 0.5) AND (length*1.0 / qlen > 0.5))
    OR ((slen*1.0 / qlen > 0.5) AND (length*1.0 / slen > 0.5)));

CREATE VIEW blast_results_max AS select gi, max(bitscore) as bitscore from blast_results group by gi;

CREATE VIEW blast_results_for_outgroup AS
  SELECT *
  FROM  blast_results_filtered
  where is_cluster_gi = 0 AND is_seq_in_species_list = 0;

-- Michal: changed id from INTEGER to TEXT since we cancelled the gi number and replaced it with accession id
create table genus_stats (id TEXT,
                          genus_name TEXT,
                          num_clusters_total INTEGER,
                          num_clusters_concat INTEGER,
                          num_clusters_locus_tree INTEGER,
                          num_filtered_by_name INTEGER,
                          PRIMARY KEY(id));

create table clusters (id TEXT,
                       genus_id INTEGER,
                       cluster_index INTEGER,
                       cluster_type TEXT,
                       is_used_for_locus_tree default 0,
                       is_used_for_concat_tree default 0,
                       is_its INTEGER default 0,
                       is_chloroplast INTEGER default 0,
                       chloroplast_loci INTEGER default 0,
                       num_of_taxons INTEGER default 0,
                       align_length INTEGER default 0,
                       data_matrix_size INTEGER default 0,
                       PRIMARY KEY(id));

CREATE VIEW outgroup_clusters AS
  select organism as outgroup_name,'organism' as outgroup_type, c.id as cluster_id, min(num_of_taxons) as num_of_taxons, min(align_length) as align_length,min(data_matrix_size) as data_matrix_size
  from blast_results_for_outgroup b, clusters c
  where 	b.outgroup_seq_id = c.id
  group by organism, cluster_id
  union
  select taxonomy_last,'genus', c.id, min(num_of_taxons), min(align_length),min(data_matrix_size)
  from blast_results_for_outgroup b, clusters c
  where 	b.outgroup_seq_id = c.id
  group by taxonomy_last, c.id
  ;

CREATE VIEW outgroup_stats AS
  select outgroup_name,outgroup_type, count(distinct cluster_id) as num_of_clusters,sum(data_matrix_size) as data_matrix_size
  from outgroup_clusters
  group by outgroup_name
  order by 4 desc,2 desc;

-- select * from outgroup_stats where data_matrix_size >= 0.5 * (select max(data_matrix_size) from outgroup_stats);

-- CREATE VIEW outgroup_stats AS
--   select organism as taxon, count(distinct c.id)
--   from blast_results_filtered b, clusters c
--   where 	b.outgroup_seq_id = c.id
--   group by organism
--   union
--   select taxonomy_last, count(distinct c.id)
--   from blast_results_filtered b, clusters c
--   where 	b.outgroup_seq_id = c.id
--   group by taxonomy_last
--   ;

