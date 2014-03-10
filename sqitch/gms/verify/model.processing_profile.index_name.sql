-- Verify model.processing_profile.index_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'processing_profile_name_index';

ROLLBACK;
