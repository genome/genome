-- Verify web.nomenclature_enum_value.index_nomenclature_field_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'nomenclature_enum_field_index';

ROLLBACK;
