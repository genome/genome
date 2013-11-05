-- Verify model.model.index_genome_model_id_subject_class_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'm_m_id_subject_class_name_index';

ROLLBACK;
