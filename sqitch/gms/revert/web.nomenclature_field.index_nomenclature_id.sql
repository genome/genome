-- Revert web.nomenclature_field.index_nomenclature_id

BEGIN;

DROP INDEX web.nomenclature_field_nomenclature_id_index;

COMMIT;
