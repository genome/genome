-- Verify model.model.index_subject_class_name_subject_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'model_subject_index';

ROLLBACK;
