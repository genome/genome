-- Revert web.nomenclature.index_name

BEGIN;

DROP INDEX web.nomenclature_name_index;

COMMIT;
