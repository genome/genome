-- Verify web.nomenclature_field.index_nomenclature_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'nomenclature_field_nomenclature_id_index';

ROLLBACK;
