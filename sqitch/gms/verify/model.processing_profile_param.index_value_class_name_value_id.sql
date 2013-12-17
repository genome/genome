-- Verify model.processing_profile_param.index_value_class_name_value_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'processing_profile_param_name_value_class_id_index';

ROLLBACK;
