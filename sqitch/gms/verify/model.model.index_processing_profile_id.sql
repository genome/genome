-- Verify model.model.index_processing_profile_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'model_processing_profile_index';

ROLLBACK;
