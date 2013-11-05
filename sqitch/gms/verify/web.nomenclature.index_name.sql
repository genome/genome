-- Verify web.nomenclature.index_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'nomenclature_name_index';

ROLLBACK;
