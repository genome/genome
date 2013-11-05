-- Verify model.processing_profile.index_subclass_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'processing_profile_subclass_index';

ROLLBACK;
