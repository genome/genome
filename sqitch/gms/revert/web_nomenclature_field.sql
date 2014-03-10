-- Revert web_nomenclature_field

BEGIN;

DROP TABLE IF EXISTS web.nomenclature_field;

COMMIT;
