-- Revert model.model.index_genome_model_id_subject_class_name

BEGIN;

DROP INDEX model.m_m_id_subject_class_name_index;

COMMIT;
