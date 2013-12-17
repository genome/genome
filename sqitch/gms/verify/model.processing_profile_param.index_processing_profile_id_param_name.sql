-- Verify model.processing_profile_param.index_processing_profile_id_param_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'processing_profile_param_id_param_name_index';

ROLLBACK;
