-- Verify model.processing_profile_param.index_value_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'processing_profile_param_value_id_index';

ROLLBACK;
