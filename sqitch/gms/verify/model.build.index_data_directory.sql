-- Verify model.build.index_data_directory

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'build_directory_index';

ROLLBACK;
