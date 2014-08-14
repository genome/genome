-- Verify config.subject_mapping_subject_index

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'subject_mapping_subject_subject_mapping_id_idx';

ROLLBACK;
