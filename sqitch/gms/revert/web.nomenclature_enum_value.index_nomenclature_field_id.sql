-- Revert web.nomenclature_enum_value.index_nomenclature_field_id

BEGIN;

DROP INDEX web.nomenclature_enum_field_index;

COMMIT;
