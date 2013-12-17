-- Verify model.model.index_subject_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'model_subject_id_index';

ROLLBACK;
